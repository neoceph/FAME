import vtk
import numpy as np


class Vertex:
    def __init__(self, coordinate):
        self.coordinate = coordinate

    def getCoordinate(self):
        return self.coordinate


class CellFace:
    def __init__(self, points):
        self.points = points  # 4 points defining the face

    def getArea(self):
        """
        Calculate the area of the face.
        Assumes the face is a rectangle for simplicity.
        """
        if len(self.points) != 4:
            raise ValueError("Face must be defined by 4 points.")
        
        # Vector math to calculate area (simple planar quad area)
        vec1 = [self.points[1][i] - self.points[0][i] for i in range(3)]
        vec2 = [self.points[2][i] - self.points[1][i] for i in range(3)]
        
        # Cross product for area calculation
        cross_product = np.cross(vec1, vec2)
        
        # Magnitude of cross product vector gives the area of the rectangle or parrellopipe
        area = np.linalg.norm(cross_product)
        return abs(area)


class Cell:
    def __init__(self, unique_index, center, vertices):
        self.unique_index = unique_index
        self.center = center  # Cell center coordinate
        self.vertices = vertices  # 8 vertices of the cell
        self.faces = {}  # Dictionary to store faces by direction

    def getArea(self, direction):
        """
        Get the area of a face in a specific direction.
        :param direction: One of ['east', 'west', 'north', 'south', 'top', 'bottom']
        :return: Area of the specified face.
        """
        if direction not in self.faces:
            raise ValueError(f"Invalid direction: {direction}")
        return self.faces[direction].getArea()


class Mesh:
    def __init__(self, x_range, y_range, z_range, divisions):
        """
        Initialize the mesh with given domain ranges and divisions.
        :param x_range: Tuple (x_min, x_max) defining the range in the x direction.
        :param y_range: Tuple (y_min, y_max) defining the range in the y direction.
        :param z_range: Tuple (z_min, z_max) defining the range in the z direction.
        :param divisions: Tuple (nx, ny, nz) for the number of divisions in x, y, z.
        """
        self.x_min, self.x_max = x_range
        self.y_min, self.y_max = y_range
        self.z_min, self.z_max = z_range
        self.nx, self.ny, self.nz = divisions
        self.cells = []
        self.points = vtk.vtkPoints()
        self.structured_grid = vtk.vtkStructuredGrid()

    def generate_mesh(self):
        """
        Generate the structured grid and populate points.
        """
        dx = (self.x_max - self.x_min) / self.nx
        dy = (self.y_max - self.y_min) / self.ny
        dz = (self.z_max - self.z_min) / self.nz

        index = 0
        for z in range(self.nz + 1):
            for y in range(self.ny + 1):
                for x in range(self.nx + 1):
                    # Add point to vtkPoints
                    self.points.InsertNextPoint(
                        self.x_min + x * dx,
                        self.y_min + y * dy,
                        self.z_min + z * dz
                    )

        self.structured_grid.SetDimensions(self.nx + 1, self.ny + 1, self.nz + 1)
        self.structured_grid.SetPoints(self.points)

        # Create cells and calculate centers
        for z in range(self.nz):
            for y in range(self.ny):
                for x in range(self.nx):
                    # Calculate cell center
                    center = (
                        self.x_min + x * dx + dx / 2,
                        self.y_min + y * dy + dy / 2,
                        self.z_min + z * dz + dz / 2
                    )

                    # Calculate vertices
                    vertices = [
                        Vertex((self.x_min + x * dx, self.y_min + y * dy, self.z_min + z * dz)),
                        Vertex((self.x_min + (x + 1) * dx, self.y_min + y * dy, self.z_min + z * dz)),
                        Vertex((self.x_min + (x + 1) * dx, self.y_min + (y + 1) * dy, self.z_min + z * dz)),
                        Vertex((self.x_min + x * dx, self.y_min + (y + 1) * dy, self.z_min + z * dz)),
                        Vertex((self.x_min + x * dx, self.y_min + y * dy, self.z_min + (z + 1) * dz)),
                        Vertex((self.x_min + (x + 1) * dx, self.y_min + y * dy, self.z_min + (z + 1) * dz)),
                        Vertex((self.x_min + (x + 1) * dx, self.y_min + (y + 1) * dy, self.z_min + (z + 1) * dz)),
                        Vertex((self.x_min + x * dx, self.y_min + (y + 1) * dy, self.z_min + (z + 1) * dz))
                    ]

                    # Define faces
                    faces = {
                        'east': CellFace([v.getCoordinate() for v in [vertices[1], vertices[5], vertices[6], vertices[2]]]),
                        'west': CellFace([v.getCoordinate() for v in [vertices[0], vertices[3], vertices[7], vertices[4]]]),
                        'north': CellFace([v.getCoordinate() for v in [vertices[2], vertices[6], vertices[7], vertices[3]]]),
                        'south': CellFace([v.getCoordinate() for v in [vertices[0], vertices[1], vertices[5], vertices[4]]]),
                        'top': CellFace([v.getCoordinate() for v in [vertices[4], vertices[5], vertices[6], vertices[7]]]),
                        'bottom': CellFace([v.getCoordinate() for v in [vertices[0], vertices[1], vertices[2], vertices[3]]]),
                    }

                    # Create the cell and assign faces
                    cell = Cell(unique_index=index, center=center, vertices=vertices)
                    cell.faces = faces

                    self.cells.append(cell)
                    index += 1

    def write_to_vtk(self, filename="structured_mesh.vtk"):
        """
        Write the structured grid to a VTK file for visualization.
        """
        writer = vtk.vtkStructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self.structured_grid)
        writer.SetFileTypeToBinary()
        writer.Write()


# Example Usage
# mesh = Mesh(x_range=(0, 1), y_range=(0, 2), z_range=(0, 1), divisions=(2, 4, 2))
# mesh.generate_mesh()

# # Access a specific cell and face
# cell = mesh.cells[0]
# print("Cell Center:", cell.center)
# print("East Face Area:", cell.getArea('east'))

# # Write mesh to VTK file
# mesh.write_to_vtk()
# print("VTK file 'structured_mesh.vtk' written successfully.")
