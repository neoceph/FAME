import unittest
import numpy as np
from scipy.sparse import lil_matrix
from src.FVM.boundaryCondition import BoundaryCondition

class MockMesh:
    def get_cells_connected_to_face(self, face_id):
        return [face_id]  # Mock: Assume one cell per face for simplicity

    def get_face_area(self, face_id):
        return 1.0  # Mock: Assume unit area for simplicity

    def get_face_center(self, face_id):
        return [0, 0, 0] if face_id % 2 == 0 else [1, 1, 1]  # Mock centers for simplicity

    def GetNumberOfFaces(self):
        return 10  # Mock: Assume 10 faces

class TestBoundaryCondition(unittest.TestCase):

    def setUp(self):
        self.mesh = MockMesh()
        self.bc = BoundaryCondition(self.mesh)
        self.num_cells = 10
        self.A = lil_matrix((self.num_cells, self.num_cells))
        self.b = np.zeros(self.num_cells)

    def test_apply_dirichlet(self):
        dirichlet_coordinates = {"x": 0, "y": None, "z": None}
        dirichlet_value = 100
        self.bc.apply_dirichlet(self.A, self.b, dirichlet_coordinates, dirichlet_value)

        # Check that the correct rows in A are set to identity
        for face_id in range(self.mesh.GetNumberOfFaces()):
            face_center = self.mesh.get_face_center(face_id)
            if np.isclose(face_center[0], 0):
                for cell_id in self.mesh.get_cells_connected_to_face(face_id):
                    self.assertTrue(np.all(self.A[cell_id, :].toarray() == 0))
                    self.assertEqual(self.A[cell_id, cell_id], 1)
                    self.assertEqual(self.b[cell_id], dirichlet_value)

    def test_apply_neumann(self):
        neumann_coordinates = {"x": 1, "y": None, "z": None}
        neumann_flux = 5
        self.bc.apply_neumann(self.b, neumann_coordinates, neumann_flux)

        # Check that the correct cells in b are updated
        for face_id in range(self.mesh.GetNumberOfFaces()):
            face_center = self.mesh.get_face_center(face_id)
            if np.isclose(face_center[0], 1):
                for cell_id in self.mesh.get_cells_connected_to_face(face_id):
                    self.assertEqual(self.b[cell_id], neumann_flux * self.mesh.get_face_area(face_id))
