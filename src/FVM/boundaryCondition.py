import numpy as np
from scipy.sparse import lil_matrix

class BoundaryCondition:
    def __init__(self, mesh, valueType='scalar'):
        """
        Initializes the BoundaryCondition object.
        
        Args:
            mesh (StructuredMesh): The mesh object.
            valueType (str): 'scalar', 'vector', or 'tensor'.
        """
        self.mesh = mesh
        self.valueType = valueType
        self.dof = 1 if valueType == 'scalar' else 3 if valueType == 'vector' else 9
        # Initialize bcValues based on number of face centers
        num_faces = len(mesh.faceCenters)
        self.bcValues = lil_matrix((num_faces, self.dof))

    def applyBoundaryCondition(self, x=None, y=None, z=None, value=0.0, tolerance=1e-6):
        """
        Apply boundary condition to a face based on x, y, z coordinates.

        Args:
            x, y, z (float, optional): Coordinates of the target face.
            value (float or np.array): Boundary condition value (scalar, vector, or tensor).
            tolerance (float): Tolerance to identify nearby faces.
        
        Returns:
            Faces to which the boundary condition was applied.
        """
        faceIds = self.mesh.getFacesByCoordinates(x=x, y=y, z=z, tolerance=tolerance)

        if not faceIds:
            raise ValueError("No matching face found within the specified tolerance.")

        # Determine value type based on input shape
        value = np.atleast_1d(value)
        
        for faceId in faceIds:
            if value.size == 1:  # Scalar
                self.bcValues[faceId, 0] = value[0]
            elif value.size == 3:  # Vector
                self.bcValues[faceId, :3] = value
            elif value.size == 9:  # Tensor
                self.bcValues[faceId, :9] = value.flatten()
            else:
                raise ValueError("Value size does not match scalar (1), vector (3), or tensor (9) dimensions.")
        
        return faceIds

    def getBoundaryMatrix(self):
        """
        Return the sparse boundary matrix.
        """
        return self.bcValues
