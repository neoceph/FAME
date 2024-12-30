import boundaryCondition as bc
import discretization as disc
import mesh as mesh
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
        bounds = {
            'x': domain['size']['x'],
            'y': domain['size']['y'],
            'z': domain['size']['z']
        }
        divisions = {
            'x': domain['divisions']['x'],
            'y': domain['divisions']['y'],
            'z': domain['divisions']['z']
        }

        self.mesh = mesh.StructuredMesh(bounds, divisions)
        print("Mesh initialized.")

    def applyBoundaryConditions(self):
        self.boundary_conditions = bc.BoundaryCondition(self.mesh)
        self.boundary_conditions.apply()
        print("Boundary conditions applied.")