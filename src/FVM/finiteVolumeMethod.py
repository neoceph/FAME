from .boundaryCondition import BoundaryCondition as bc
from .discretization import Discretization as disc
from .mesh import StructuredMesh as mesh
from .property import MaterialProperty as prop
from .solver import Solver as sol
from .visualization import MeshWriter as vis


class FVM:
    def __init__(self, config):
        self.config = config
        self.mesh = None
        self.boundaryConditions = None
        self.material_properties = None
        self.discretization = None
        self.solver = None
        self.visualization = None

    def meshGeneration(self):
        domain = self.config['simulation']['domain']
        bounds = (
            tuple(domain['size']['x']),
            tuple(domain['size']['y']),
            tuple(domain['size']['z'])
        )
        div_x = domain['divisions']['x']
        div_y = domain['divisions']['y']
        div_z = domain['divisions']['z']
        self.mesh = mesh(bounds, (div_x, div_y, div_z))
        print("Mesh initialized.")

    def applyBoundaryConditions(self):
        boundary_config = self.config['simulation'].get('boundaryConditions', {}).get('parameters', {})
        
        # Initialize BoundaryCondition with default values if parameters are missing
        self.boundaryConditions = bc(
            self.mesh,
            **{key: boundary_config.get('temperature', {}).get(key, 0)
            for key in ['variableType', 'convectionCoefficient', 'emmissivity', 'ambientTemperature']}
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



    def _applyMultipleBoundaryConditions(self, x, y, z, bc):
        bcType = bc['type']
        bcValue = bc['value']
        self.boundaryConditions.applyBoundaryCondition(x, y, z, bcValue)
        print(f"Applied {bcType} of {bcValue} at {'x' if x else 'y' if y else 'z'} = {x or y or z}")

    def loadMaterialProperty(self):
        material_config = self.config.get('simulation', {}).get('material', {})
        self.material_properties = {}

        material_name = material_config.get('name', 'Unknown Material')

        for property_name, base_value in material_config.items():
            if property_name != 'name':
                material_property = prop(
                    materialName=material_name,
                    propertyName=property_name,
                    baseValue=base_value
                )
                self.material_properties[property_name] = material_property
                print(f"Loaded {property_name} for {material_name} with value: {base_value}")

        print("Material properties successfully initialized.")

    def discretization(self):
        if not self.mesh:
            raise ValueError("Mesh must be generated before discretization.")
        self.discretization = disc.Discretization(self.mesh)
        scheme = self.config['simulation'].get('discretization', {}).get('scheme', 'central')
        self.discretization.applyDiscretization(scheme)
        print(f"Discretization applied using {scheme} scheme.")

    def solveEquations(self):
        if not self.discretization:
            raise ValueError("Discretization must be performed before solving.")
        self.solver = sol.Solver(self.discretization)
        solver_type = self.config['simulation'].get('solver', {}).get('type', 'GaussSeidel')
        tolerance = self.config['simulation'].get('solver', {}).get('tolerance', 1e-8)
        maxIterations = self.config['simulation'].get('solver', {}).get('maxIterations', 1000)
        self.solver.solve(solver_type, tolerance, maxIterations)
        print(f"Solver {solver_type} completed with tolerance {tolerance} and max iterations {maxIterations}.")

    def visualizeResults(self):
        if not self.solver or not self.solver.solution:
            raise ValueError("Solution must exist before visualization.")
        self.visualization = vis.Visualization(self.mesh, self.solver.solution)
        output_type = self.config['simulation'].get('visualization', {}).get('type', 'contour')
        output_path = self.config['simulation'].get('visualization', {}).get('path', './testOutput')
        self.visualization.generateOutput(output_type, output_path)
        print(f"Visualization generated and saved at {output_path} as {output_type} plot.")

    def simulate(self):
        self.meshGeneration()
        self.applyBoundaryConditions()
        self.loadMaterialProperty()
        self.discretization()
        self.solveEquations()
        self.visualizeResults()
        print("Simulation complete.")
