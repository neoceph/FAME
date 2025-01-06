import numpy as np
import scipy.sparse as sp
from .mesh import StructuredMesh1D, StructuredMesh

class BoundaryCondition:
    def __new__(cls, mesh, *args, **kwargs):
        if isinstance(mesh, StructuredMesh1D):
            return super(BoundaryCondition, BoundaryCondition1D).__new__(BoundaryCondition1D)
        else:
            return super(BoundaryCondition, BoundaryCondition3D).__new__(BoundaryCondition3D)

    def __init__(self, mesh, variableType='scalar', convectionCoefficient=0, emmissivity=0, dependentSource=0, independentSource=0, volumetricSource=0, ambientTemperature=0):
        """
        Initializes the BoundaryCondition object.
        Args:
            mesh (StructuredMesh): The mesh object.
            valueType (str): 'scalar', 'vector', or 'tensor'.
        """
        self.mesh = mesh
        self.valueType = variableType
        self.dof = 1 if variableType == 'scalar' else 3 if variableType == 'vector' else 9
        
        self.convectionCoefficient = convectionCoefficient
        self.emissivity = emmissivity
        self.ambientTemperature = ambientTemperature
        self.num_faces = len(self.mesh.faceCenters)
        self._initializeSources(dependentSource, independentSource, volumetricSource)

    def _initializeSources(self, dependentSource, independentSource, volumetricSource):
        raise NotImplementedError("This method must be implemented by subclasses.")

    def applyBoundaryCondition(self, x=None, y=None, z=None, value=0.0, tolerance=1e-6):
        raise NotImplementedError("applyBoundaryCondition must be implemented by subclass.")

    def getBoundaryMatrix(self):
        return self.bcValues

class BoundaryCondition3D(BoundaryCondition):
    def _initializeSources(self, dependentSource, independentSource, volumetricSource):
        vectorLength = self.mesh.divisions[0] * self.mesh.divisions[1] * self.mesh.divisions[2]
        self.bcValues = sp.lil_matrix((self.num_faces, self.dof))
        self.dependentSource = sp.lil_matrix((vectorLength, 1))
        self.independentSource = sp.lil_matrix((vectorLength, 1))
        self.volumetricSource = sp.lil_matrix((vectorLength, 1))
        
        for i in range(vectorLength):
            self.dependentSource[i, 0] = dependentSource
            self.independentSource[i, 0] = independentSource
            self.volumetricSource[i, 0] = volumetricSource

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
            if value.size == 1:
                self.bcValues[faceId, 0] = value[0]
            elif value.size == 3:
                self.bcValues[faceId, :3] = value
            elif value.size == 9:
                self.bcValues[faceId, :9] = value.flatten()
            else:
                raise ValueError("Value size does not match scalar (1), vector (3), or tensor (9) dimensions.")
        return faceIds

class BoundaryCondition1D(BoundaryCondition):
    def _initializeSources(self, dependentSource, independentSource, volumetricSource):
        vectorLength = self.mesh.divisions[0]
        self.bcValues = sp.lil_matrix((self.num_faces, self.dof))
        self.dependentSource = sp.lil_matrix((vectorLength, 1))
        self.independentSource = sp.lil_matrix((vectorLength, 1))
        self.volumetricSource = sp.lil_matrix((vectorLength, 1))
        
        for i in range(vectorLength):
            self.dependentSource[i, 0] = dependentSource
            self.independentSource[i, 0] = independentSource
            self.volumetricSource[i, 0] = volumetricSource

    def applyBoundaryCondition(self, x=None, y=None, z=None, value=0.0, tolerance=1e-6):
        faceIds = self.mesh.getFacesByCoordinates(x=x, tolerance=tolerance)
        if not faceIds:
            raise ValueError("No matching face found within the specified tolerance.")

        value = np.atleast_1d(value)
        for faceId in faceIds:
            if value.size == 1:
                self.bcValues[faceId, 0] = value[0]
            elif value.size == 3:
                self.bcValues[faceId, :3] = value
            elif value.size == 9:
                self.bcValues[faceId, :9] = value.flatten()
            else:
                raise ValueError("Value size does not match scalar (1), vector (3), or tensor (9) dimensions.")
        return faceIds
