import vtk
import numpy as np
import scipy.sparse as sp
from .boundaryCondition import BoundaryCondition as bc
from .discretization import Discretization as disc
from .mesh import StructuredMesh, StructuredMesh1D
from .property import MaterialProperty as prop
from .solver import Solver as sol
from .visualization import MeshWriter, MeshWriter1D
from ..utils.utility import timing_decorator

class FVM:
    def __new__(cls, config):
        domain = config['simulation']['domain']
        divisions = domain['divisions']
        if len(divisions) == 1:
            instance = super(FVM, FVM1D).__new__(FVM1D)
        else:
            instance = super(FVM, FVM3D).__new__(FVM3D)
        instance.__init__(config)
        return instance

    def __init__(self, config):
        self.config = config
        self.mesh = None
        self.boundaryConditions = None
        self.materialProperties = None
        self.discretization = None
        self.solver = None
        self.visualization = None
        self.output = None
        self.nodalSolution = None

    def meshGeneration(self):
        raise NotImplementedError("meshGeneration must be implemented by subclass.")

    @timing_decorator
    def applyBoundaryConditions(self):
        boundary_config = self.config['simulation'].get('boundaryConditions', {}).get('parameters', {})
        
        # Initialize BoundaryCondition with default values if parameters are missing
        self.boundaryConditions = bc(
            self.mesh,
            **{key: boundary_config.get('temperature', {}).get(key, 0)
            for key in ['variableType', 'convectionCoefficient', 'emmissivity', 'dependentSource', 'independentSource', 'volumetricSource', 'ambientTemperature']}
        )

        conditions = self.config['simulation'].get('boundaryConditions', {})
        
        # Safely iterate over conditions if they exist
        apply_methods = {
            'x': lambda coord, bcValue: self.boundaryConditions.applyBoundaryCondition(coord, None, None, bcValue),
            'y': lambda coord, bcValue: self.boundaryConditions.applyBoundaryCondition(None, coord, None, bcValue),
            'z': lambda coord, bcValue: self.boundaryConditions.applyBoundaryCondition(None, None, coord, bcValue)
        }

        for axis, axis_conditions in conditions.items():
            if axis == 'parameters':
                continue  # Skip the parameters key itself during condition application
            for coord, bc_list in axis_conditions.items():
                for bcItem in bc_list if isinstance(bc_list, list) else [bc_list]:
                    apply_methods.get(axis, lambda c, b: None)(coord, bcItem['value'])
                    print(f"Applied {bcItem['type']} condition at {axis} = {coord} with value {bcItem['value']}")

        print("Boundary conditions applied.")

    @timing_decorator
    def loadMaterialProperty(self):
        material_config = self.config.get('simulation', {}).get('material', {})
        self.materialProperties = {}

        material_name = material_config.get('name', 'Unknown Material')
        material_property = prop(material_name)

        # Loop over properties inside 'properties' block
        for property_name, property_details in material_config.get('properties', {}).items():
            base_value = property_details.get('baseValue', 0)
            method = property_details.get('method', 'constant')
            coefficients = property_details.get('coefficients', [])
            # Ensure coefficients are converted to floats
            coefficients = [float(c) for c in coefficients]
            reference_temperature = property_details.get('referenceTemperature', 298.15)

            # Add the property
            material_property.add_property(
                propertyName=property_name,
                baseValue=base_value,
                method=method,
                referenceTemperature=reference_temperature,
                coefficients=coefficients
            )
            print(f"Added {property_name} to {material_name} with method {method}.")

        # Store the populated material property
        self.materialProperties[material_name] = material_property
        print(f"Material properties successfully initialized for {material_name}.")

    @timing_decorator
    def discretize(self):
        if not self.mesh:
            raise ValueError("Mesh must be generated before discretization.")
        material_name = self.config['simulation']['material']['name']  # Get material name dynamically
        self.discretization = disc(self.mesh, self.solver, self.materialProperties[material_name], self.boundaryConditions)
        self.discretization.discretizeHeatDiffusion()
        print("Discretization applied.")
    
    @timing_decorator
    def interpolateNodeFromCell(self):
        """
        Interpolate the cell‐centered solution onto the mesh nodes
        by averaging all adjacent cell values at each mesh point,
        then let subclasses overwrite Dirichlet BCs in-place.
        """
        # --- prepare in-place storage ---
        num_pts = self.mesh.GetNumberOfPoints()
        # reuse or (re)allocate the nodalSolution array
        self.nodalSolution = np.zeros(num_pts, dtype=float)
        counts = np.zeros(num_pts, dtype=int)

        # --- accumulate cell contributions ---
        cell_vals = self.solution[0]
        for cid in range(self.mesh.GetNumberOfCells()):
            v = cell_vals[cid]
            pts_ids = self.mesh.GetCell(cid).GetPointIds()
            for i in range(pts_ids.GetNumberOfIds()):
                pid = pts_ids.GetId(i)
                self.nodalSolution[pid] += v
                counts[pid] += 1

        # --- normalize to get averages ---
        counts[counts == 0] = 1
        # in-place division
        self.nodalSolution /= counts

        # --- apply any Dirichlet BCs in-place ---
        self._apply_nodal_bc(self.nodalSolution)

        return self.nodalSolution
    
    @timing_decorator
    def solveEquations(self):
        if not self.mesh:
            raise ValueError("Mesh must be generated before solving.")
        
        backend = self.config['simulation'].get('solver', {}).get('module')
        self.solver = sol(self.mesh.A, self.mesh.b, backend=backend)
        solver_type = self.config['simulation'].get('solver', {}).get('method')
        tolerance = self.config['simulation'].get('solver', {}).get('tolerance')
        maxIterations = self.config['simulation'].get('solver', {}).get('maxIterations')
        self.solution = self.solver.solve(method=solver_type, preconditioner="none")
        print(f"Solver {solver_type} completed with tolerance {tolerance} and max iterations {maxIterations}.")        
    
    @timing_decorator
    def visualizeResults(self):
        if not self.solver or self.solver.solution is None:
            raise ValueError("Solution must exist before visualization.")
            
        # Initialize MeshWriter with the mesh
        self.visualization = MeshWriter(self.mesh)

        # Read variable name from YAML or default to 'temperature_cell'
        variable_name = self.config['simulation'].get('visualization', {}).get('variableName', 'temperature') + "_cell"
        nodal_variable_name = self.config['simulation'].get('visualization', {}).get('nodalVariableName', 'temperature') + "_node"

        # Prepare the solution as a cell variable dictionary
        solution, err, info = self.solution
        variables = {
            variable_name: solution,
            nodal_variable_name: self.interpolateNodeFromCell()
        }
        
        # Read the output path from YAML or default to current directory
        output_path = self.config['simulation'].get('visualization', {}).get('path', './')
        
        # Write the VTS file
        self.visualization.writeVTS(output_path, variables)
        print(f"Visualization generated and saved at {output_path} with variable '{variable_name}'.")

    def simulate(self):
        self.meshGeneration()
        self.applyBoundaryConditions()
        self.loadMaterialProperty()
        self.discretize()
        self.solveEquations()
        self.visualizeResults()
        print("Simulation complete.")


class FVM3D(FVM):
    @timing_decorator
    def meshGeneration(self):
        domain = self.config['simulation']['domain']
        bounds = (
            tuple(domain['size']['x']),
            tuple(domain['size']['y']),
            tuple(domain['size']['z'])
        )
        divisions = (domain['divisions']['x'], domain['divisions']['y'], domain['divisions']['z'])
        self.mesh = StructuredMesh(bounds, divisions)
        print("3D Mesh initialized.")

    def _apply_nodal_bc(self, nodalSolution: np.ndarray):
        """
        Overwrite in-place any Dirichlet “temperature” BCs on the mesh nodes
        whose physical x,y or z coordinate matches the BC coordinate.
        """
        bc_conf = self.config['simulation'].get('boundaryConditions', {})
        tol = self.config['simulation'].get('boundaryConditions', {})\
                    .get('parameters', {})\
                    .get('tolerance', 1e-6)

        n_pts = self.mesh.GetNumberOfPoints()
        pts   = self.mesh.GetPoints()

        for axis, axis_cfg in bc_conf.items():
            if axis == 'parameters':
                continue

            for coord_key, bc_list in axis_cfg.items():
                # coord_key is the physical location (e.g. 0.0 or 1.0)
                coord = float(coord_key)

                # gather all pids on this plane
                if axis == 'x':
                    match_pids = [
                        pid for pid in range(n_pts)
                        if abs(pts.GetPoint(pid)[0] - coord) <= tol
                    ]
                elif axis == 'y':
                    match_pids = [
                        pid for pid in range(n_pts)
                        if abs(pts.GetPoint(pid)[1] - coord) <= tol
                    ]
                elif axis == 'z':
                    match_pids = [
                        pid for pid in range(n_pts)
                        if abs(pts.GetPoint(pid)[2] - coord) <= tol
                    ]
                else:
                    continue

                # apply all BCs on that plane
                for bc in (bc_list if isinstance(bc_list, list) else [bc_list]):
                    if bc.get('type') != 'temperature':
                        continue
                    T = float(bc['value'])
                    for pid in match_pids:
                        nodalSolution[pid] = T

class FVM1D(FVM):
    @timing_decorator
    def meshGeneration(self):
        domain = self.config['simulation']['domain']
        bounds = tuple(domain['size']['x'])
        divisions = (domain['divisions']['x'])
        self.mesh = StructuredMesh1D(bounds, divisions, faceArea=domain.get('area', 1.0))
        print("1D Mesh initialized.")

    @timing_decorator
    def visualizeResults(self):
        if not self.solver or self.solver.solution is None:
            raise ValueError("Solution must exist before visualization.")
            
        # Initialize MeshWriter with the mesh
        self.visualization = MeshWriter1D(self.mesh)

        # Read variable name from YAML or default to 'temperature_cell'
        variable_name = self.config['simulation'].get('visualization', {}).get('variableName', 'temperature_cell')

        # Prepare the solution as a cell variable dictionary
        solution, err, info = self.solution
        variables = {
            variable_name: solution
        }
        
        # Read the output path from YAML or default to current directory
        output_path = self.config['simulation'].get('visualization', {}).get('path', './')
        
        # Write the VTS file
        self.visualization.writeVTS(output_path, variables)
        print(f"Visualization generated and saved at {output_path} with variable '{variable_name}'.")

    def _apply_nodal_bc(self, nodalSolution: np.ndarray):
        """
        Overwrite in-place any Dirichlet “temperature” BCs on the 1D mesh
        nodes whose x-coordinate matches the BC coordinate.
        """
        bc_conf = self.config['simulation'].get('boundaryConditions', {})
        tol = self.config['simulation'].get('boundaryConditions', {})\
                    .get('parameters', {})\
                    .get('tolerance', 1e-6)

        n_pts = self.mesh.GetDimensions()  # returns int
        pts   = self.mesh.GetPoints()

        # pick any axis key (only 'x' makes sense in 1D)
        axis_cfg = bc_conf.get('x', {})
        for coord_key, bc_list in axis_cfg.items():
            coord = float(coord_key)
            match_pids = [
                pid for pid in range(n_pts)
                if abs(pts.GetPoint(pid)[0] - coord) <= tol
            ]

            for bc in (bc_list if isinstance(bc_list, list) else [bc_list]):
                if bc.get('type') != 'temperature':
                    continue
                T = float(bc['value'])
                for pid in match_pids:
                    nodalSolution[pid] = T