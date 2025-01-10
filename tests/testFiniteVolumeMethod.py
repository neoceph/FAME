import unittest
import os
import yaml
import numpy as np
import scipy.sparse as sp

from fame.FVM.finiteVolumeMethod import FVM
from fame.FVM.solver import Solver


class TestDiscretizationBase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Placeholder to ensure child classes initialize FVM properly
        cls.fvm = None

    def test_materialProperties(self):
        # Ensure the FVM instance is initialized
        if self.fvm is None:
            self.skipTest("Skipping test because FVM instance is not initialized.")
        self.assertIn('Aluminum', self.fvm.materialProperties)
        aluminum = self.fvm.materialProperties['Aluminum']
        self.assertIn('thermalConductivity', aluminum.properties)
        thermal_conductivity = aluminum.properties['thermalConductivity']
        self.assertEqual(thermal_conductivity['referenceTemperature'], 298.15)


class TestDiscretization(TestDiscretizationBase):
    @classmethod
    def setUpClass(cls):
        yaml_path = os.path.join(os.path.dirname(__file__), '..', 'examples', 'FVM', 'HeatDiffusion', 'setup.yaml')
        with open(yaml_path, 'r') as file:
            cls.config = yaml.safe_load(file)
        cls.fvm = FVM(cls.config)
        cls.fvm.meshGeneration()
        cls.fvm.applyBoundaryConditions()
        cls.fvm.loadMaterialProperty()

        # Initialize solver after mesh generation
        cls.fvm.solver = Solver(cls.fvm.mesh.A, cls.fvm.mesh.b)
        
        cls.fvm.discretize()
        cls.fvm.solveEquations()
        
    def test_meshGeneration(self):
        
        self.assertIsNotNone(self.fvm.mesh)
        print("Mesh generation test passed.")

    def test_boundaryConditions(self):
        
        self.assertIsNotNone(self.fvm.boundaryConditions)
        print("Boundary conditions test passed.")

    # def test_materialProperties(self):
               
    #     # Check if Aluminum is loaded as the material
    #     self.assertIn('Aluminum', self.fvm.materialProperties)

    #     # Access the Aluminum material and its properties
    #     aluminum = self.fvm.materialProperties['Aluminum']
    #     self.assertIn('thermalConductivity', aluminum.properties)
        
    #     # Verify the base value, method, and reference temperature of thermal conductivity
    #     thermal_conductivity = aluminum.properties['thermalConductivity']
    #     self.assertEqual(thermal_conductivity['baseValue'], 237)
    #     self.assertEqual(thermal_conductivity['referenceTemperature'], 298.15)
    #     self.assertEqual(thermal_conductivity['method'], 'polynomial')
        
    def test_discretization_applied(self):
        
        self.assertIsNotNone(self.fvm.solver)
        self.assertIsNotNone(self.fvm.solver.A, "Discretization matrix A should not be None.")
        self.assertGreater(self.fvm.solver.A.count_nonzero(), 0, "Matrix A should have non-zero elements after discretization.")

    def test_solver(self):
        
        self.assertIsNotNone(self.fvm.solver)
        print("Solver test passed.")

    def test_visualization(self):
        output_path = self.fvm.config['simulation'].get('visualization', {}).get('path', './testOutput')
        self.fvm.visualizeResults()
        self.assertTrue(os.path.exists(output_path))
        print(f"Visualization output test passed. Output saved to: {output_path}")

    def test_fullSimulation(self):
        self.fvm.simulate()
        self.assertIsNotNone(self.fvm.mesh)
        self.assertIsNotNone(self.fvm.boundaryConditions)
        self.assertIn('thermalConductivity', self.fvm.materialProperties['Aluminum'].properties)
        self.assertIsNotNone(self.fvm.solver)
        output_path = self.fvm.config['simulation'].get('visualization', {}).get('path', './testOutput')
        self.assertTrue(os.path.exists(output_path))
        print("Full simulation test passed.")


class TestDiscretizationSmall(TestDiscretization):
    @classmethod
    def setUpClass(cls):
        yaml_path = os.path.join(os.path.dirname(__file__), '..', 'examples', 'FVM', 'HeatDiffusion', 'setup_small.yaml')
        with open(yaml_path, 'r') as file:
            cls.config = yaml.safe_load(file)
        cls.fvm = FVM(cls.config)
        cls.fvm.meshGeneration()
        cls.fvm.applyBoundaryConditions()
        cls.fvm.loadMaterialProperty()

        cls.fvm.solver = Solver(cls.fvm.mesh.A, cls.fvm.mesh.b)
        
        cls.fvm.discretize()
        cls.fvm.solveEquations()

    def test_fullSimulation(self):
        self.fvm.simulate()
        self.assertIsNotNone(self.fvm.mesh)
        self.assertIsNotNone(self.fvm.boundaryConditions)
        self.assertIn('thermalConductivity', self.fvm.materialProperties['Aluminum'].properties)
        self.assertIsNotNone(self.fvm.solver)
        output_path = self.fvm.config['simulation'].get('visualization', {}).get('path', './testOutput')
        self.assertTrue(os.path.exists(output_path))
        print("Full simulation test passed.")


class TestDiscretization1D(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        yaml_path = os.path.join(os.path.dirname(__file__), '..', 'examples', 'FVM', 'HeatDiffusion', 'setup_1D.yaml')
        with open(yaml_path, 'r') as file:
            cls.config = yaml.safe_load(file)

        # Dynamically create FVM instance based on the new static __new__ method
        cls.fvm = FVM(cls.config)
        
        # Initialize mesh, boundary conditions, and material properties
        cls.fvm.meshGeneration()
        cls.fvm.applyBoundaryConditions()
        cls.fvm.loadMaterialProperty()

        # Set up solver and discretization
        cls.fvm.solver = Solver(cls.fvm.mesh.A, cls.fvm.mesh.b)
        cls.fvm.discretize()
        cls.fvm.solveEquations()

    def test_fullSimulation(self):
        # Run full simulation and check that all parts are initialized correctly
        self.fvm.simulate()
        self.assertIsNotNone(self.fvm.mesh, "Mesh is not initialized.")
        self.assertIsNotNone(self.fvm.boundaryConditions, "Boundary conditions not applied.")
        self.assertIn('thermalConductivity', self.fvm.materialProperties['Aluminum'].properties)
        self.assertIsNotNone(self.fvm.solver, "Solver is not initialized.")
        
        # Verify output path
        output_path = self.fvm.config['simulation'].get('visualization', {}).get('path', './testOutput')
        self.assertTrue(os.path.exists(output_path), "Output path does not exist.")
        print("Full 1D simulation test passed.")

    def test_materialProperties(self):
        # Explicitly call the method from the base class
        TestDiscretizationBase.test_materialProperties(self)
        
        # Additional assertions specific to 1D
        thermal_conductivity = self.fvm.materialProperties['Aluminum'].properties['thermalConductivity']
        self.assertEqual(thermal_conductivity['baseValue'], 1000)
        self.assertEqual(thermal_conductivity['method'], 'constant')
