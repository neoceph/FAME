import unittest
import scipy.sparse as sp
import numpy as np
import os
import shutil

from src.FVM.mesh import StructuredMesh
from src.FVM.property import MaterialProperty
from src.FVM.solver import Solver
from src.FVM.discretization import Discretization
from src.FVM.boundaryCondition import BoundaryCondition


class TestDiscretization(unittest.TestCase):

    def setUp(self):
        """
        Set up StructuredMesh, MaterialProperty, and Solver objects for each test case.
        """
        
        self.outputDir = "testOutput"
        
        # Initialize mesh with bounds and divisions
        self.mesh = StructuredMesh(((0, 10), (0, 5), (0, 3)), (10, 5, 3))
        
        # Initialize MaterialProperty for thermal conductivity
        self.prop = MaterialProperty('thermal_conductivity', baseValue=200, referenceTemperature=298.15, method='constant')
        
        # Initialize solver with default empty matrix and vector
        numCells = self.mesh.GetNumberOfCells()
        A = sp.lil_matrix((numCells, numCells))  # Empty sparse matrix
        b = np.zeros(numCells)  # Zero vector
        
        self.solver = Solver(A, b)

        # Initialize boundary condition
        self.bc = BoundaryCondition(self.mesh)
        
        # Initialize Discretization object
        self.discretization = Discretization(self.mesh, self.solver, self.prop, self.bc)
        os.makedirs(self.outputDir, exist_ok=True)

    def tearDown(self):
        """
        Clean up after each test by removing the output directory.
        """
        if os.path.exists(self.outputDir):
            shutil.rmtree(self.outputDir)  # Recursively remove the directory
            print(f"Removed directory: {self.outputDir}")

    def testInitialization(self):
        """
        Test if Discretization object initializes correctly with given material property.
        """
        self.assertEqual(self.discretization.property, 200, "Thermal conductivity should be set to 200.")
        
        # Test other property fallback
        invalid_prop = MaterialProperty('density', baseValue=100, referenceTemperature=300, method='variable')
        disc_invalid = Discretization(self.mesh, self.solver, invalid_prop)
        self.assertEqual(disc_invalid.property, 0, "Non-thermal conductivity property should default to 0.")

    def testDiscretizeHeatDiffusion(self):
        """
        Test if heat diffusion discretization correctly populates matrix A and vector b.
        """
        self.discretization.discretizeHeatDiffusion(self.prop.baseValue)

        num_cells = self.mesh.GetNumberOfCells()

        # Test matrix dimensions
        self.assertEqual(self.solver.A.shape, (num_cells, num_cells), "Matrix A dimensions should match the number of cells.")
        self.assertEqual(self.solver.b.shape[0], num_cells, "Vector b length should match the number of cells.")
        
        # Test matrix sparsity
        non_zero_elements = self.solver.A.count_nonzero()
        self.assertGreater(non_zero_elements, 0, "Matrix A should have non-zero elements after discretization.")

        # Test if the diagonal has non-zero values
        for i in range(num_cells):
            self.assertNotEqual(self.solver.A[i, i], 0, f"Diagonal element A[{i},{i}] should not be zero.")

        # Plot the sparse matrix after discretization
        outputFilename = os.path.join(self.outputDir, "test_matrix_plot.jpeg")
        self.solver.plotSparseMatrix(self.solver.A, filename=outputFilename)

        print("Sparse matrix visualization saved as 'test_matrix_plot.jpeg'.")