import vtk
import numpy as np
import scipy.sparse as sp
from .boundaryCondition import BoundaryCondition as bc
from .discretization import Discretization as disc
from .mesh import StructuredMesh, StructuredMesh1D
from .property import MaterialProperty as prop
from .solver import Solver as sol
from .visualization import MeshWriter, MeshWriter1D

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

    def meshGeneration(self):
        raise NotImplementedError("meshGeneration must be implemented by subclass.")

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

    def discretize(self):
        if not self.mesh:
            raise ValueError("Mesh must be generated before discretization.")
        material_name = self.config['simulation']['material']['name']  # Get material name dynamically
        self.discretization = disc(self.mesh, self.solver, self.materialProperties[material_name], self.boundaryConditions)
        self.discretization.discretizeHeatDiffusion()
        print("Discretization applied.")
    
    def solveEquations(self):
        if not self.mesh:
            raise ValueError("Mesh must be generated before solving.")
        
        self.solver = sol(self.mesh.A, self.mesh.b)
        solver_type = self.config['simulation'].get('solver', {}).get('method')
        tolerance = self.config['simulation'].get('solver', {}).get('tolerance')
        maxIterations = self.config['simulation'].get('solver', {}).get('maxIterations')
        self.solution = self.solver.solve(method=solver_type, preconditioner="none")
        print(f"Solver {solver_type} completed with tolerance {tolerance} and max iterations {maxIterations}.")
    
    def visualizeResults(self):
        if not self.solver or self.solver.solution is None:
            raise ValueError("Solution must exist before visualization.")
            
        # Initialize MeshWriter with the mesh
        self.visualization = MeshWriter(self.mesh)

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

    def simulate(self):
        self.meshGeneration()
        self.applyBoundaryConditions()
        self.loadMaterialProperty()
        self.discretize()
        self.solveEquations()
        self.visualizeResults()
        print("Simulation complete.")


class FVM3D(FVM):
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

class FVM1D(FVM):
    def meshGeneration(self):
        domain = self.config['simulation']['domain']
        bounds = tuple(domain['size']['x'])
        divisions = (domain['divisions']['x'])
        self.mesh = StructuredMesh1D(bounds, divisions)
        print("1D Mesh initialized.")

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