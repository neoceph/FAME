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
            solver (Solver): The solver object containing matrix A and vector b.
            property (MaterialProperty): MaterialProperty object with material and property information.
            boundaryCondition (BoundaryCondition): BoundaryCondition object for applying boundary conditions.
        """
        self.mesh = mesh
        self.solver = solver
        self.boundaryCondition = boundaryCondition
        self.property = property

    def discretizeHeatDiffusion(self):
        """
        Discretize the 3D heat diffusion equation and populate the sparse matrix A and vector b from the solver class.
        Uses temperature-dependent thermal conductivity from the MaterialProperty class.
        """

        for cellID in range(self.mesh.numCells):
            
            if 'thermalConductivity' in self.property.properties:
                thermalConductivity = self.property.evaluate('thermalConductivity', 298.15)  # Default temp used for evaluation
            else:
                raise ValueError("Material property must include 'thermalConductivity'")
            
            # off-diagonal matrix element construction
            for sharedCellID, sharedFace in zip (self.mesh.sharedCells[cellID]['shared_cells'], self.mesh.sharedCells[cellID]['shared_faces']):
                points = vtk.vtkPoints()
                for point in self.mesh.faces[sharedFace]:
                    points.InsertNextPoint(self.mesh.GetPoint(point))
                
                faceArea = self.mesh.calculateArea(points)
                cellDistance = np.linalg.norm(np.array(self.mesh.cellCenters[cellID]) - np.array(self.mesh.cellCenters[sharedCellID]))

                self.mesh.A[cellID, sharedCellID] = thermalConductivity * (faceArea)/(cellDistance)

                # diagonal marix element construction
                self.mesh.A[cellID, cellID] += - thermalConductivity * (faceArea)/(cellDistance)
            
            # diagonal marix element construction for boundary faces
            for sharedBoundaryFace in self.mesh.sharedCells[cellID]['boundary_faces']:
                points = vtk.vtkPoints()
                for point in self.mesh.faces[sharedBoundaryFace]:
                    points.InsertNextPoint(self.mesh.GetPoint(point))
                boundaryFaceArea = self.mesh.calculateArea(points)
                distance_cell_to_boundary_face = np.linalg.norm(np.array(self.mesh.cellCenters[cellID]) - np.array(self.mesh.faceCenters[sharedBoundaryFace]))
                self.mesh.A[cellID, cellID] += - thermalConductivity * (boundaryFaceArea)/(distance_cell_to_boundary_face) - self.boundaryCondition.convectionCoefficient
                self.mesh.b[cellID] += -thermalConductivity * self.boundaryCondition.bcValues[sharedBoundaryFace, 0] * (boundaryFaceArea) / (distance_cell_to_boundary_face)    \
                                        -self.boundaryCondition.convectionCoefficient * boundaryFaceArea * self.boundaryCondition.ambientTemperature
                
            self.mesh.A[cellID, cellID] += -self.boundaryCondition.dependentSource[cellID, 0]
            self.mesh.b[cellID] += -self.boundaryCondition.independentSource[cellID, 0] - self.boundaryCondition.volumetricSource[cellID, 0] * self.mesh.getCellVolume(cellID)
