import os
import unittest
from src.FVM.mesh import Vertex, CellFace, Cell, Mesh


class TestVertex(unittest.TestCase):
    def test_vertex_initialization(self):
        coordinate = (1.0, 2.0, 3.0)
        vertex = Vertex(coordinate)
        self.assertEqual(vertex.getCoordinate(), coordinate)


class TestCellFace(unittest.TestCase):
    def test_cellface_initialization(self):
        points = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]
        face = CellFace(points)
        self.assertEqual(face.points, points)

    def test_cellface_get_area_valid(self):
        points = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]
        face = CellFace(points)
        self.assertAlmostEqual(face.getArea(), 1.0)

    def test_cellface_get_area_invalid(self):
        points = [(0, 0, 0), (1, 0, 0), (1, 1, 0)]  # Only 3 points
        face = CellFace(points)
        with self.assertRaises(ValueError):
            face.getArea()


class TestCell(unittest.TestCase):
    def test_cell_initialization(self):
        center = (0.5, 0.5, 0.5)
        vertices = [Vertex((x, y, z)) for x, y, z in [
            (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
            (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)
        ]]
        cell = Cell(unique_index=0, center=center, vertices=vertices)
        self.assertEqual(cell.unique_index, 0)
        self.assertEqual(cell.center, center)
        self.assertEqual(len(cell.vertices), 8)

    def test_cell_get_area(self):
        center = (0.5, 0.5, 0.5)
        vertices = [Vertex((x, y, z)) for x, y, z in [
            (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
            (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)
        ]]
        faces = {
            'east': CellFace([(1, 0, 0), (1, 0, 1), (1, 1, 1), (1, 1, 0)]),
        }
        cell = Cell(unique_index=0, center=center, vertices=vertices)
        cell.faces = faces
        self.assertAlmostEqual(cell.getArea('east'), 1.0)

        with self.assertRaises(ValueError):
            cell.getArea('invalid_direction')


class TestMesh(unittest.TestCase):
    def test_mesh_initialization(self):
        mesh = Mesh(x_range=(0, 1), y_range=(0, 2), z_range=(0, 1), divisions=(2, 4, 2))
        self.assertEqual(mesh.x_min, 0)
        self.assertEqual(mesh.x_max, 1)
        self.assertEqual(mesh.y_min, 0)
        self.assertEqual(mesh.y_max, 2)
        self.assertEqual(mesh.z_min, 0)
        self.assertEqual(mesh.z_max, 1)
        self.assertEqual(mesh.nx, 2)
        self.assertEqual(mesh.ny, 4)
        self.assertEqual(mesh.nz, 2)

    def test_mesh_generate_mesh(self):
        mesh = Mesh(x_range=(0, 1), y_range=(0, 1), z_range=(0, 1), divisions=(1, 1, 1))
        mesh.generate_mesh()
        self.assertEqual(len(mesh.cells), 1)  # Only one cell in this simple case
        self.assertEqual(mesh.structured_grid.GetNumberOfPoints(), 8)  # 2x2x2 points

    def test_mesh_write_to_vtk(self):
        mesh = Mesh(x_range=(0, 1), y_range=(0, 1), z_range=(0, 1), divisions=(1, 1, 1))
        mesh.generate_mesh()
        mesh.write_to_vtk(filename="test_mesh.vtk")
        # Verify file is written
        self.assertTrue(os.path.exists("test_mesh.vtk"))
        # os.remove("test_mesh.vtk")


if __name__ == "__main__":
    unittest.main()
