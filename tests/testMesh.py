import os
import unittest
import vtk
from src.FVM.mesh import StructuredMesh  # Assuming your StructuredMesh is in the same directory as mesh.py


class TestStructuredMesh(unittest.TestCase):
    def setUp(self):
        """
        Set up a StructuredMesh object for each test case.
        """
        self.bounds = ((0, 10), (0, 5), (0, 2))
        self.divisions = (10, 5, 2)
        self.mesh = StructuredMesh(self.bounds, self.divisions)

    def testMeshInitialization(self):
        """
        Test that the StructuredMesh is initialized correctly.
        """
        # Verify dimensions
        nx, ny, nz = self.mesh.GetDimensions()
        self.assertEqual((nx, ny, nz), self.divisions)

        # Verify points are generated
        self.assertEqual(self.mesh.GetNumberOfPoints(), 11 * 6 * 3)

    def test_vertex_initialization(self):
        """
        Test that vertices (points in the mesh) are initialized correctly.
        """
        # Check the first point in the mesh
        point_id = 0
        point = self.mesh.GetPoint(point_id)
        expected_point = (0.0, 0.0, 0.0)  # First point should correspond to (0, 0, 0)
        self.assertEqual(point, expected_point)

        # Check the last point in the mesh
        last_point_id = self.mesh.GetNumberOfPoints() - 1
        last_point = self.mesh.GetPoint(last_point_id)
        expected_last_point = (10.0, 5.0, 2.0)  # Last point in the grid
        self.assertEqual(last_point, expected_last_point)

        # Check the number of points
        expected_num_points = 11 * 6 * 3
        self.assertEqual(self.mesh.GetNumberOfPoints(), expected_num_points)

    def testComputeCellCenters(self):
        """
        Test that cell centers are computed correctly.
        """
        self.assertEqual(len(self.mesh.cell_centers), self.mesh.GetNumberOfCells())

        # Check the center of the first cell
        first_cell_center = self.mesh.getCellCenter(0)
        self.assertIsInstance(first_cell_center, tuple)
        self.assertEqual(len(first_cell_center), 3)  # Should be an (x, y, z) coordinate

    def testComputeNeighbors(self):
        """
        Test that neighboring cells and shared vertices are computed correctly.
        """
        shared_info = self.mesh.getSharedCellsInfo(0)
        self.assertIn("shared_cells", shared_info)
        self.assertIn("shared_vertices", shared_info)

        # Verify the first cell has neighbors
        shared_cells = shared_info["shared_cells"]
        shared_vertices = shared_info["shared_vertices"]
        self.assertIsInstance(shared_cells, list)
        self.assertIsInstance(shared_vertices, list)

        # Ensure there is at least one neighbor
        self.assertGreater(len(shared_cells), 0)
        self.assertGreater(len(shared_vertices), 0)

    def testMeshSharedCellsIteration(self):
        """
        Test iteration over shared cells and vertices.
        """
        cell_id = 0
        shared_info = self.mesh.getSharedCellsInfo(cell_id)
        shared_cells = shared_info["shared_cells"]
        shared_vertices = shared_info["shared_vertices"]

        # Verify that shared cells and vertices can be iterated
        for shared_cell, vertices in zip(shared_cells, shared_vertices):
            self.assertIsInstance(shared_cell, int)
            self.assertIsInstance(vertices, list)
            self.assertGreater(len(vertices), 2)  # More than two vertices are shared

    def testTriangleArea(self):
        """
        Test that the area of a triangle is computed correctly using the calculateArea method.
        """
        points = vtk.vtkPoints()
        points.InsertNextPoint(0, 0, 0)
        points.InsertNextPoint(1, 0, 0)
        points.InsertNextPoint(0, 1, 0)

        area = self.mesh.calculateArea(points)
        expected_area = 0.5
        self.assertAlmostEqual(area, expected_area, places=5)
    
    def testRectangleAreaXYPlane(self):
        """
        Test that the area of a rectangle is computed correctly using the calculateArea method.
        """
        points = vtk.vtkPoints()
        points.InsertNextPoint(0, 0, 0)
        points.InsertNextPoint(1, 0, 0)
        points.InsertNextPoint(1, 1, 0)
        points.InsertNextPoint(0, 1, 0)

        area = self.mesh.calculateArea(points)
        expected_area = 1.0
        self.assertAlmostEqual(area, expected_area, places=5)

        points = vtk.vtkPoints()
        points.InsertNextPoint(0, 0, 0)
        points.InsertNextPoint(1, 1, 0)
        points.InsertNextPoint(1, 0, 0)
        points.InsertNextPoint(0, 1, 0)

        area = self.mesh.calculateArea(points)
        expected_area = 1.0
        self.assertAlmostEqual(area, expected_area, places=5)

        points = vtk.vtkPoints()
        points.InsertNextPoint(0, 1, 0)
        points.InsertNextPoint(0, 0, 0)
        points.InsertNextPoint(1, 1, 0)
        points.InsertNextPoint(1, 0, 0)

        area = self.mesh.calculateArea(points)
        expected_area = 1.0
        self.assertAlmostEqual(area, expected_area, places=5)

        points = vtk.vtkPoints()
        points.InsertNextPoint(0, 0, 0)
        points.InsertNextPoint(1, 0, 0)
        points.InsertNextPoint(0, 1, 0)
        points.InsertNextPoint(1, 1, 0)

        area = self.mesh.calculateArea(points)
        expected_area = 1.0
        self.assertAlmostEqual(area, expected_area, places=5)
    
    def testRectangleAreaXZPlane(self):
        """
        Test that the area of a rectangle is computed correctly using the calculateArea method.
        """
        points = vtk.vtkPoints()
        points.InsertNextPoint(0, 0, 0)
        points.InsertNextPoint(1, 0, 0)
        points.InsertNextPoint(1, 0, 1)
        points.InsertNextPoint(0, 0, 1)

        area = self.mesh.calculateArea(points)
        expected_area = 1.0
        self.assertAlmostEqual(area, expected_area, places=5)

        points = vtk.vtkPoints()
        points.InsertNextPoint(0, 0, 0)
        points.InsertNextPoint(1, 0, 1)
        points.InsertNextPoint(1, 0, 0)
        points.InsertNextPoint(0, 0, 1)

        area = self.mesh.calculateArea(points)
        expected_area = 1.0
        self.assertAlmostEqual(area, expected_area, places=5)

        points = vtk.vtkPoints()
        points.InsertNextPoint(1, 0, 0)
        points.InsertNextPoint(0, 0, 0)
        points.InsertNextPoint(1, 0, 1)
        points.InsertNextPoint(0, 0, 1)

        area = self.mesh.calculateArea(points)
        expected_area = 1.0
        self.assertAlmostEqual(area, expected_area, places=5)

        points = vtk.vtkPoints()
        points.InsertNextPoint(0, 0, 0)
        points.InsertNextPoint(1, 0, 0)
        points.InsertNextPoint(0, 0, 1)
        points.InsertNextPoint(1, 0, 1)

        area = self.mesh.calculateArea(points)
        expected_area = 1.0
        self.assertAlmostEqual(area, expected_area, places=5)

    def testRectangleAreaYZPlane(self):
        """
        Test that the area of a rectangle is computed correctly using the calculateArea method.
        """
        points = vtk.vtkPoints()
        points.InsertNextPoint(0, 0, 0)
        points.InsertNextPoint(0, 1, 0)
        points.InsertNextPoint(0, 1, 1)
        points.InsertNextPoint(0, 0, 1)

        area = self.mesh.calculateArea(points)
        expected_area = 1.0
        self.assertAlmostEqual(area, expected_area, places=5)

        points = vtk.vtkPoints()
        points.InsertNextPoint(0, 0, 0)
        points.InsertNextPoint(0, 1, 1)
        points.InsertNextPoint(0, 1, 0)
        points.InsertNextPoint(0, 0, 1)

        area = self.mesh.calculateArea(points)
        expected_area = 1.0
        self.assertAlmostEqual(area, expected_area, places=5)

        points = vtk.vtkPoints()
        points.InsertNextPoint(0, 1, 0)
        points.InsertNextPoint(0, 0, 0)
        points.InsertNextPoint(0, 1, 1)
        points.InsertNextPoint(0, 0, 1)

        area = self.mesh.calculateArea(points)
        expected_area = 1.0
        self.assertAlmostEqual(area, expected_area, places=5)

        points = vtk.vtkPoints()
        points.InsertNextPoint(0, 0, 0)
        points.InsertNextPoint(0, 1, 0)
        points.InsertNextPoint(0, 0, 1)
        points.InsertNextPoint(0, 1, 1)

        area = self.mesh.calculateArea(points)
        expected_area = 1.0
        self.assertAlmostEqual(area, expected_area, places=5)

    def testSharedCellArea(self):
        """
        Test that the area of a shared face between two cells is computed correctly using the calculateArea method.
        """
        cell_id = 0  # Choose a cell ID
        shared_cells_info = self.mesh.getSharedCellsInfo(cell_id)

        # Iterate over shared cells and compute area of the shared face
        for shared_cell, shared_vertices in zip(shared_cells_info["shared_cells"], shared_cells_info["shared_vertices"]):
            points = vtk.vtkPoints()

            # Add shared vertices to the vtkPoints object
            for vertex in shared_vertices:
                points.InsertNextPoint(vertex)

            # Calculate the area using the calculateArea method
            area = self.mesh.calculateArea(points)

            # For a structured grid with uniform spacing, validate expected area
            expected_area = 1.0  # Assume uniform cell face spacing for this test
            self.assertAlmostEqual(area, expected_area, places=5)
    


if __name__ == "__main__":
    unittest.main()
