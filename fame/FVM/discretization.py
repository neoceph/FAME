import numpy as np
import scipy.sparse as sp
import vtk
import warnings
from .mesh import StructuredMesh
from .solver import Solver
from .property import MaterialProperty
from .boundaryCondition import BoundaryCondition

class Discretization:
    def __init__(self, mesh, solver, property, boundaryCondition, solution=None):
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
        self.solution = solution

    def discretizeHeatDiffusion(self, timeStep: float = np.inf, timeIntegrationMethod: float = 0.0):
        """
        Discretize the 3D heat diffusion equation and populate the sparse matrix A and vector b from the solver class.
        Uses temperature-dependent thermal conductivity from the MaterialProperty class.
        Args:
            timeStep (float): Time step for the discretization, default is infinity for steady state case.
            timeIntegrationMethod (float): Method for time integration, default is 0.0 for explicit. For implicit, set to 1.0, and for Crank-Nicolson, set to 0.5.
        """

        for cellID in range(self.mesh.numCells):
            
            oldSolution = self.solution[cellID]
            
            thermalConductivity = self.property.evaluate('thermalConductivity', oldSolution)  # Default temp used for evaluation
            specificHeat = self.property.evaluate('specificHeat', oldSolution)
            density = self.property.evaluate('density', oldSolution)
            cellVolume = self.mesh.getCellVolume(cellID)
            temporalFlux = density * specificHeat * cellVolume / timeStep
            
            # off-diagonal matrix element construction
            for sharedCellID, sharedFace in zip (self.mesh.sharedCells[cellID]['shared_cells'], self.mesh.sharedCells[cellID]['shared_faces']): # the assumption is that sharedCells and sharedFaces are lists of the same length as only one face is shared between two cells.
                points = vtk.vtkPoints()
                for point in self.mesh.faces[sharedFace]:
                    points.InsertNextPoint(self.mesh.GetPoint(point))
                
                faceArea = self.mesh.calculateArea(points)
                cellDistance = np.linalg.norm(np.array(self.mesh.cellCenters[cellID]) - np.array(self.mesh.cellCenters[sharedCellID]))

                cellFlux = thermalConductivity * faceArea / cellDistance
                
                self.mesh.A[cellID, sharedCellID] = -timeIntegrationMethod * cellFlux
                self.mesh.A[cellID, cellID] += timeIntegrationMethod * cellFlux

                self.mesh.b[cellID] += (1-timeIntegrationMethod) * (cellFlux * oldSolution[sharedCellID] - cellFlux * oldSolution[cellID])
                           
            self.mesh.A[cellID, cellID] += temporalFlux - timeIntegrationMethod * self.boundaryCondition.dependentSource[cellID, 0]
            self.mesh.b[cellID] += (temporalFlux + (1-timeIntegrationMethod)*self.boundaryCondition.dependentSource[cellID, 0])*oldSolution[cellID] +  self.boundaryCondition.volumetricSource[cellID, 0] * cellVolume

            # diagonal marix element construction for boundary faces
            for sharedBoundaryFace in self.mesh.sharedCells[cellID]['boundary_faces']:
                points = vtk.vtkPoints()
                for point in self.mesh.faces[sharedBoundaryFace]:
                    points.InsertNextPoint(self.mesh.GetPoint(point))
                
                boundaryFaceArea = self.mesh.calculateArea(points)
                distance_cell_to_boundary_face = np.linalg.norm(np.array(self.mesh.cellCenters[cellID]) - np.array(self.mesh.faceCenters[sharedBoundaryFace]))
                
                self.mesh.A[cellID, cellID] += thermalConductivity * (boundaryFaceArea)/(distance_cell_to_boundary_face) + self.boundaryCondition.convectionCoefficient
                self.mesh.b[cellID] += thermalConductivity * self.boundaryCondition.bcValues[sharedBoundaryFace, 0] * (boundaryFaceArea) / (distance_cell_to_boundary_face)+self.boundaryCondition.convectionCoefficient * boundaryFaceArea * self.boundaryCondition.ambientTemperature
