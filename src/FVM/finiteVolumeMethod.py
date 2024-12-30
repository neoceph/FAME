import boundaryCondition as bc
import discretization as disc
import mesh as mesh
import property as prop
import solver as sol
import visualization as vis


class FVM:
    def __init__(self, config):
        self.config = config
        self.mesh = None
        self.boundary_conditions = None
        self.material_properties = None
        self.discretization = None
        self.solver = None
        self.visualization = None

    def meshGeneration(self):

        domain = self.config['simulation']['domain']
        
        # Extract and format bounds as tuples
        bounds = (
            tuple(domain['size']['x']),
            tuple(domain['size']['y']),
            tuple(domain['size']['z'])
        )
        
        # Extract divisions directly
        div_x = domain['divisions']['x']
        div_y = domain['divisions']['y']
        div_z = domain['divisions']['z']
        
        # Initialize mesh with tuple bounds and unpacked divisions
        self.mesh = mesh.StructuredMesh(bounds, (div_x, div_y, div_z))
        print("Mesh initialized.")

    def applyBoundaryConditions(self):
        # Initialize the BoundaryCondition object
        self.boundary_conditions = bc.BoundaryCondition(self.mesh, valueType='scalar')
        
        # Read boundary conditions from the config
        boundary_conditions = self.config['simulation']['boundary_conditions']

        # Apply boundary conditions at each specified axis and coordinate
        for axis, conditions in boundary_conditions.items():
            for coord, bc in conditions.items():
                bc_type = bc['type']
                bc_value = bc['value']

                # Apply boundary condition by directly passing x, y, z coordinates
                if axis == 'x':
                    self.boundary_conditions.applyBoundaryCondition(coord, None, None, bc_value)
                elif axis == 'y':
                    self.boundary_conditions.applyBoundaryCondition(None, coord, None, bc_value)
                elif axis == 'z':
                    self.boundary_conditions.applyBoundaryCondition(None, None, coord, bc_value)
                
                print(f"Applied {bc_type} of {bc_value} at {axis} = {coord}")
        
        print("Boundary conditions applied.")

    def loadMaterialProperty(self):
        """
        Load material properties from the configuration and initialize the MaterialProperty object.
        """
        # Extract material properties from the configuration
        material_config = self.config.get('simulation', {}).get('material', {})

        # Initialize the MaterialProperty object directly from YAML fields
        self.material_properties = {}

        for property_name, base_value in material_config.items():
            if property_name != 'name':  # Skip the material name, only use properties
                material_property = prop.MaterialProperty(propertyName=property_name, baseValue=base_value)
                self.material_properties[property_name] = material_property
                print(f"Loaded {property_name} with value: {base_value}")

        print("Material properties successfully initialized.")

