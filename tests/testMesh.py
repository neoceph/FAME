import os
import unittest
import vtk
import numpy as np

from src.FVM.mesh import StructuredMesh  # Assuming your StructuredMesh is in the same directory as mesh.py


class TestStructuredMesh(unittest.TestCase):
    def setUp(self):
        """
        Set up a StructuredMesh object for each test case.
        """
        self.bounds = ((0, 10), (0, 5), (0, 3))
        self.divisions = (10, 5, 3)
        self.mesh = StructuredMesh(self.bounds, self.divisions)

    def testMeshInitialization(self):
        """
        Test that the StructuredMesh is initialized correctly.
        """
        # Verify dimensions
        nx, ny, nz = self.mesh.GetDimensions()
        div_x = nx - 1
        div_y = ny - 1
        div_z = nz - 1
        
        self.assertEqual((div_x, div_y, div_z), self.divisions)

        # Verify points are generated
        self.assertEqual(self.mesh.GetNumberOfPoints(), 11 * 6 * 4)

    def testVertexInitialization(self):
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
        expected_last_point = (10.0, 5.0, 3.0)  # Last point in the grid
        self.assertEqual(last_point, expected_last_point)

        # Check the number of points
        expected_num_points = 11 * 6 * 4
        self.assertEqual(self.mesh.GetNumberOfPoints(), expected_num_points)

    def testFaceInitialization(self):
        """
        Test that cell faces are computed correctly.
        """
        
        expectedFaces = 3*self.divisions[0]*self.divisions[1]*self.divisions[2] + self.divisions[0]*self.divisions[1] + self.divisions[1]*self.divisions[2] + self.divisions[0]*self.divisions[2]
        totalFaces = len(self.mesh.faces)
        self.assertEqual(totalFaces, expectedFaces)

    def testComputeCellCenters(self):
        """
        Test that cell centers are computed correctly.
        """
        self.assertEqual(len(self.mesh.cellCenters), self.mesh.GetNumberOfCells())

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

        # Verify the first cell has neighbors
        shared_cells = shared_info["shared_cells"]
        self.assertIsInstance(shared_cells, list)

        # Ensure there is at least one neighbor
        self.assertGreater(len(shared_cells), 0)

    def testMeshSharedCellsIteration(self):
        """
        Test iteration over shared cells and vertices.
        """
        cell_id = 0
        shared_info = self.mesh.getSharedCellsInfo(cell_id)
        shared_cells = shared_info["shared_cells"]

        # Verify that shared cells and vertices can be iterated
        for shared_cell in shared_cells:
            self.assertIsInstance(shared_cell, int)

    def test_shared_cells(self):
        """
        Test if sharedCells are correctly computed and populated.
        """
        # Ensure sharedCells list is not empty
        self.assertGreater(len(self.mesh.sharedCells), 0, "sharedCells should not be empty.")

        # Iterate through each cell and verify structure
        for cell_info in self.mesh.sharedCells:
            # Check required keys exist
            self.assertIn('cell_id', cell_info)
            self.assertIn('shared_cells', cell_info)
            self.assertIn('shared_faces', cell_info)
            self.assertIn('boundary_faces', cell_info)

            # Validate types of each key
            self.assertIsInstance(cell_info['cell_id'], int)
            self.assertIsInstance(cell_info['shared_cells'], list)
            self.assertIsInstance(cell_info['shared_faces'], list)
            self.assertIsInstance(cell_info['boundary_faces'], list)

            # Ensure the cell_id is within expected bounds
            self.assertGreaterEqual(cell_info['cell_id'], 0)
            self.assertLess(cell_info['cell_id'], self.mesh.GetNumberOfCells())

            # Check that shared cells and faces are integers
            for shared_cell in cell_info['shared_cells']:
                self.assertIsInstance(shared_cell, int)

            for face in cell_info['shared_faces']:
                self.assertIsInstance(face, int)

            for boundary_face in cell_info['boundary_faces']:
                self.assertIsInstance(boundary_face, int)            

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

        # Iterate over shared faces and compute area
        for face_id in shared_cells_info["shared_faces"]:
            face_points = vtk.vtkPoints()

            # Retrieve the points of the shared face
            for vertex_id in self.mesh.faces[face_id]:
                face_points.InsertNextPoint(self.mesh.GetPoint(vertex_id))

            # Calculate the area using the calculateArea method
            area = self.mesh.calculateArea(face_points)

            # For a structured grid with uniform spacing, validate expected area
            expected_area = 1.0  # Assume uniform cell face spacing for this test
            self.assertAlmostEqual(area, expected_area, places=5)

    def testCornerCellSharedFaces(self):
        """
        Test that corner cells have exactly 3 shared faces.
        """
        cell_id = 0  # Corner cell at (0, 0, 0)
        shared_cells_info = self.mesh.getSharedCellsInfo(cell_id)
        self.assertEqual(len(shared_cells_info["shared_faces"]), 3)

    def testEdgeIntersectionCellSharedFaces(self):
        """
        Test that cells at the intersection of two boundaries have exactly 4 shared faces.
        """
        cell_id = (self.divisions[0])  # Cell at the edge of two boundaries
        shared_cells_info = self.mesh.getSharedCellsInfo(cell_id)
        self.assertEqual(len(shared_cells_info["shared_faces"]), 4)

    def testBoundaryCellSharedFaces(self):
        """
        Test that boundary cells (not at corners) have exactly 5 shared faces.
        """
        cell_id = self.divisions[0] * self.divisions[1] + self.divisions[0]  # Mid-cell locating after 1st X-Y layer plus 1st row in X and advancing 2 more cells.
        shared_cells_info = self.mesh.getSharedCellsInfo(cell_id)
        self.assertEqual(len(shared_cells_info["shared_faces"]), 5)

    def testCenterCellSharedFaces(self):
        """
        Test that center cells have exactly 6 shared faces.
        """
        cell_id = self.divisions[0] * self.divisions[1] + self.divisions[0] + 4  # Mid-cell locating after 1st X-Y layer plus 1st row in X and advancing 2 more cells.
        shared_cells_info = self.mesh.getSharedCellsInfo(cell_id)
        self.assertEqual(len(shared_cells_info["shared_faces"]), 6)

    def testGetFaceById(self):
        """
        Test that getFaceById correctly retrieves face points.
        """
        face_id = 0  # Test first face
        face_points = self.mesh.getFaceById(face_id)
        
        self.assertIsNotNone(face_points)
        self.assertEqual(len(face_points), 4)  # A hexahedral face should have 4 points
        
        # Validate that the points correspond to valid mesh points
        for pointId in face_points:
            point = self.mesh.GetPoint(pointId)
            self.assertIsInstance(point, tuple)
            self.assertEqual(len(point), 3)

    def testGetFaceByCenter(self):
        """
        Test that getFaceByCenter correctly retrieves the nearest face ID(s).
        """
        firstFaceCenter = list(self.mesh.faceCenters.values())[0]
        faceIds = self.mesh.getFaceByCenter(firstFaceCenter)
        
        # Handle list of face IDs or single face ID
        if isinstance(faceIds, list):
            self.assertGreater(len(faceIds), 0)
            for faceId in faceIds:
                self.assertIn(faceId, self.mesh.faces)
                retrieved_center = np.array(self.mesh.faceCenters[faceId])
                np.testing.assert_array_almost_equal(firstFaceCenter, retrieved_center, decimal=5)
        else:
            self.assertIn(faceIds, self.mesh.faces)
            retrieved_center = np.array(self.mesh.faceCenters[faceIds])
            np.testing.assert_array_almost_equal(firstFaceCenter, retrieved_center, decimal=5)


    def testGetFaceByCenterWithTolerance(self):
        """
        Test that getFaceByCenter retrieves all faces within a given tolerance.
        """
        test_point = (0.5, 0.5, 0.5)
        tolerance = 0.6
        faceIDs = self.mesh.getFaceByCenter(test_point, tolerance=tolerance)
        
        # Ensure the result is always treated as a list
        if not isinstance(faceIDs, list):
            faceIDs = [faceIDs]
        
        self.assertGreater(len(faceIDs), 0)

        # Validate that all returned faces are within the tolerance distance
        for face_id in faceIDs:
            self.assertIn(face_id, self.mesh.faceCenters)  # Ensure the face exists
            face_center = np.array(self.mesh.faceCenters[face_id])
            distance = np.linalg.norm(face_center - np.array(test_point))
            self.assertLessEqual(distance, tolerance)


    def testCellVolumeCalculation(self):
        """
        Test that the volume of a cell is calculated correctly using the Gauss Divergence Theorem.
        """
        num_cells = self.mesh.GetNumberOfCells()
        
        # Iterate over all cells and calculate their volume
        for cell_id in range(num_cells):
            volume = self.mesh.getCellVolume(cell_id)
            
            # Expected volume for structured hexahedral grid cells
            dx = (self.bounds[0][1] - self.bounds[0][0]) / self.divisions[0]
            dy = (self.bounds[1][1] - self.bounds[1][0]) / self.divisions[1]
            dz = (self.bounds[2][1] - self.bounds[2][0]) / self.divisions[2]
            
            expected_volume = dx * dy * dz
            
            # Assert that the calculated volume matches the expected value
            self.assertAlmostEqual(volume, expected_volume, places=5, msg=f"Mismatch in volume for cell {cell_id}")

    def testGetFacesByX(self):
        
        expectedCellFaces = self.divisions[1] * self.divisions[2]
        result = self.mesh.getFacesByCoordinates(x=self.bounds[0][0], tolerance=0.1)
        self.assertEqual(len(result), expectedCellFaces)

        result = self.mesh.getFacesByCoordinates(x=self.bounds[0][1], tolerance=0.1)
        self.assertEqual(len(result), expectedCellFaces)

    def testGetFacesByY(self):
        
        expectedCellFaces = self.divisions[0] * self.divisions[2]
        result = self.mesh.getFacesByCoordinates(y=self.bounds[1][0], tolerance=0.1)
        self.assertEqual(len(result), expectedCellFaces)

        result = self.mesh.getFacesByCoordinates(y=self.bounds[1][1], tolerance=0.1)
        self.assertEqual(len(result), expectedCellFaces)

    def testGetFacesByZ(self):
        
        expectedCellFaces = self.divisions[0] * self.divisions[1]
        result = self.mesh.getFacesByCoordinates(z=self.bounds[2][0], tolerance=0.1)
        self.assertEqual(len(result), expectedCellFaces)

        result = self.mesh.getFacesByCoordinates(z=self.bounds[2][1], tolerance=0.1)
        self.assertEqual(len(result), expectedCellFaces)        