import numpy as np
import scipy.sparse as sp
import vtk
from .mesh import StructuredMesh
from .solver import Solver
from .property import MaterialProperty
from .boundaryCondition import BoundaryCondition

class Discretization:
    def __init__(self, mesh, solver, property, boundaryCondition):
        """
        Initializes the Discretization class.

        Args:
            mesh (StructuredMesh): The mesh object containing cells and shared face information.
        """
        self.mesh = mesh
        self.solver = solver
        self.boundaryCondition = BoundaryCondition
        if property.propertyName == 'thermal_conductivity' and property.method == 'constant':
            self.property = property.baseValue
        else:
            print("other property is not implemented")
            self.property = 0

    def discretizeHeatDiffusion(self, thermalConductivity):
        """
        Discretize the 3D heat diffusion equation and populate the sparse matrix A and vector b from the solver class.

        Args:
            thermalConductivity (float): Thermal conductivity of the material provided from property module.
        """
        numCells = self.mesh.GetNumberOfCells()
        self.solver.A = sp.lil_matrix((numCells, numCells))  # Use LIL format for construction
        self.solver.b = np.zeros(numCells)

        for cellID in range(numCells):
            
            # off-diagonal matrix element construction
            for sharedCellID, sharedFace in zip (self.mesh.sharedCells[cellID]['shared_cells'], self.mesh.sharedCells[cellID]['shared_faces']):
                points = vtk.vtkPoints()
                for point in self.mesh.faces[sharedFace]:
                    points.InsertNextPoint(self.mesh.GetPoint(point))
                
                faceArea = self.mesh.calculateArea(points)
                cellDistance = np.linalg.norm(np.array(self.mesh.cellCenters[cellID]) - np.array(self.mesh.cellCenters[sharedCellID]))

                self.solver.A[cellID, sharedCellID] = self.property * (faceArea)/(cellDistance)

                # diagonal marix element construction
                self.solver.A[cellID, sharedCellID] += - self.property * (faceArea)/(cellDistance)
            
            # diagonal marix element construction for boundary faces
            for sharedBoundaryFace in self.mesh.sharedCells['boundary_faces']:
                points = vtk.vtkPoints()
                for point in self.mesh.faces[sharedBoundaryFace]:
                    points.InsertNextPoint(self.mesh.GetPoint(point))
                boundaryFaceArea = self.mesh.calculateArea(points)
                cell_to_boundary_face_distance = np.linalg.norm(np.array(self.mesh.cellCenters[cellID]) - np.array(self.mesh.faceCenters[sharedBoundaryFace]))
                self.solver.A[cellID, sharedCellID] += - self.property * (boundaryFaceArea)/(cell_to_boundary_face_distance) - self.boundaryCondition.convectionCoefficient
                self.solver.b[cellID] += -self.property * self.boundaryCondition.bcValues[sharedBoundaryFace, 0] * (boundaryFaceArea) / (cell_to_boundary_face_distance)    \
                                        -self.boundaryCondition.convectionCoefficient * boundaryFaceArea * self.boundaryCondition.ambientTemperature
                
            self.solver.A[cellID, sharedCellID] += -self.boundaryCondition.dependentSource[cellID]
            self.solver.b[cellID] += -self.boundaryCondition.independentSource[cellID] - self.boundaryCondition.volumetricSource[cellID] * self.mesh.getCellVolume[cellID]
