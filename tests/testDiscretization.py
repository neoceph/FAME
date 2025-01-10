import unittest
import scipy.sparse as sp
import numpy as np
import os
import shutil

from fame.FVM.mesh import StructuredMesh
from fame.FVM.property import MaterialProperty
from fame.FVM.solver import Solver
from fame.FVM.discretization import Discretization
from fame.FVM.boundaryCondition import BoundaryCondition


class TestDiscretization(unittest.TestCase):

    def setUp(self):
        """
        Set up StructuredMesh, MaterialProperty, and Solver objects for each test case.
        """
        
        self.outputDir = "testOutput"
        
        # Initialize mesh with bounds and divisions
        self.mesh = StructuredMesh(((0, 10), (0, 5), (0, 3)), (10, 5, 3))
        
        # Initialize MaterialProperty for Aluminum and add thermal conductivity
        self.prop = MaterialProperty('Aluminum')
        self.prop.add_property('thermalConductivity', baseValue=200, referenceTemperature=298.15, method='constant')
        
        self.solver = Solver(self.mesh.A, self.mesh.b)

        # Initialize boundary condition
        self.bc = BoundaryCondition(self.mesh)
               
        # Initialize Discretization object with evaluated property
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
        # Evaluate the thermal conductivity using the property object
        thermal_conductivity = self.prop.evaluate('thermalConductivity', 298.15)
        self.assertEqual(thermal_conductivity, 200, "Thermal conductivity should be set to 200.")
        
        # Test other property fallback
        invalid_prop = MaterialProperty('Steel')
        invalid_prop.add_property('density', baseValue=100, referenceTemperature=300, method='variable')
        
        disc_invalid = Discretization(self.mesh, self.solver, invalid_prop, self.bc)
        with self.assertRaises(ValueError):
            disc_invalid.discretizeHeatDiffusion()

    def testDiscretizeHeatDiffusion(self):
        """
        Test if heat diffusion discretization correctly populates matrix A and vector b.
        """
        # Access the property directly from the MaterialProperty object
        self.discretization.discretizeHeatDiffusion()

        # Test matrix dimensions
        self.assertEqual(self.mesh.A.shape, (self.mesh.numCells, self.mesh.numCells), "Matrix A dimensions should match the number of cells.")
        self.assertEqual(self.mesh.b.shape[0], self.mesh.numCells, "Vector b length should match the number of cells.")
        
        # Test matrix sparsity
        non_zero_elements = self.mesh.A.count_nonzero()
        self.assertGreater(non_zero_elements, 0, "Matrix A should have non-zero elements after discretization.")

        # Test if the diagonal has non-zero values
        for i in range(self.mesh.numCells):
            self.assertNotEqual(self.mesh.A[i, i], 0, f"Diagonal element A[{i},{i}] should not be zero.")

        # Plot the sparse matrix after discretization
        outputFilename = os.path.join(self.outputDir, "test_matrix_plot.jpeg")
        self.solver.plotSparseMatrix(self.mesh.A, filename=outputFilename)

        print("Sparse matrix visualization saved as 'test_matrix_plot.jpeg'.")
