import unittest
import numpy as np
from scipy.sparse import lil_matrix

from src.FVM.boundaryCondition import BoundaryCondition
from src.FVM.mesh import StructuredMesh  # Assuming StructuredMesh is defined in mesh.py

class TestBoundaryCondition(unittest.TestCase):
    def setUp(self):
        """
        Set up a StructuredMesh object for each test case.
        """
        self.bounds = ((0, 10), (0, 5), (0, 3))
        self.divisions = (10, 5, 3)
        self.mesh = StructuredMesh(self.bounds, self.divisions)

    def test_apply_scalar_boundary_conditions(self):
        bc = BoundaryCondition(self.mesh, valueType='scalar')
        
        # Apply scalar value to faces at x=0 and x=10
        faceIds_x = bc.applyBoundaryCondition(x=0, value=100)
        faceIds_x.extend(bc.applyBoundaryCondition(x=10, value=100))
        self.assertEqual(len(faceIds_x), 30)
        
        # Apply scalar value to faces at y=0 and y=5
        faceIds_y = bc.applyBoundaryCondition(y=0, value=100)
        faceIds_y.extend(bc.applyBoundaryCondition(y=5, value=100))
        self.assertEqual(len(faceIds_y), 60)

        # Apply scalar value to faces at z=0 and z=3
        faceIds_z = bc.applyBoundaryCondition(z=0, value=100)
        faceIds_z.extend(bc.applyBoundaryCondition(z=3, value=100))
        self.assertEqual(len(faceIds_z), 100)
        
        # Validate boundary matrix by iterating over faceIds
        for faceId in faceIds_x + faceIds_y + faceIds_z:
            print(faceId)
            self.assertEqual(bc.bcValues[faceId, 0], 100)

    def test_apply_vector_boundary_conditions(self):
        bc = BoundaryCondition(self.mesh, valueType='vector')
        vector_value = np.array([100, 50, 25])
        
        # Apply vector to faces at x=0 and x=10
        faceIds_x = bc.applyBoundaryCondition(x=0, value=vector_value)
        faceIds_x.extend(bc.applyBoundaryCondition(x=10, value=vector_value))
        self.assertEqual(len(faceIds_x), 30)
        
        # Apply vector to faces at y=0 and y=5
        faceIds_y = bc.applyBoundaryCondition(y=0, value=vector_value)
        faceIds_y.extend(bc.applyBoundaryCondition(y=5, value=vector_value))
        self.assertEqual(len(faceIds_y), 60)
        
        # Apply vector to faces at z=0 and z=3
        faceIds_z = bc.applyBoundaryCondition(z=0, value=vector_value)
        faceIds_z.extend(bc.applyBoundaryCondition(z=3, value=vector_value))
        self.assertEqual(len(faceIds_z), 100)

        # Validate boundary matrix by iterating over faceIds
        for faceId in faceIds_x + faceIds_y + faceIds_z:
            np.testing.assert_array_equal(bc.bcValues[faceId, :3].toarray().flatten(), vector_value)

    def test_apply_tensor_boundary_conditions(self):
        bc = BoundaryCondition(self.mesh, valueType='tensor')
        tensor_value = np.eye(3).flatten() * 100
        
        # Apply tensor to faces at x=0 and x=10
        faceIds_x = bc.applyBoundaryCondition(x=0, value=tensor_value)
        faceIds_x.extend(bc.applyBoundaryCondition(x=10, value=tensor_value))
        self.assertEqual(len(faceIds_x), 30)
        
        # Apply tensor to faces at y=0 and y=5
        faceIds_y = bc.applyBoundaryCondition(y=0, value=tensor_value)
        faceIds_y.extend(bc.applyBoundaryCondition(y=5, value=tensor_value))
        self.assertEqual(len(faceIds_y), 60)
        
        # Apply tensor to faces at z=0 and z=3
        faceIds_z = bc.applyBoundaryCondition(z=0, value=tensor_value)
        faceIds_z.extend(bc.applyBoundaryCondition(z=3, value=tensor_value))
        self.assertEqual(len(faceIds_z), 100)

        # Validate boundary matrix by iterating over faceIds
        for faceId in faceIds_x + faceIds_y + faceIds_z:
            np.testing.assert_array_equal(bc.bcValues[faceId, :9].toarray().flatten(), tensor_value)