import unittest
import os
from src.FVM.finiteVolumeMethod import FVM
import yaml

class TestFVM(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        yaml_path = os.path.join(os.path.dirname(__file__), '..', 'examples', 'FVM', 'HeatDiffusion', 'setup.yaml')
        with open(yaml_path, 'r') as file:
            cls.config = yaml.safe_load(file)
        cls.fvm = FVM(cls.config)
        cls.fvm.meshGeneration()  # Generate mesh once for all tests

    def test_meshGeneration(self):
        self.fvm.meshGeneration()
        self.assertIsNotNone(self.fvm.mesh)
        print("Mesh generation test passed.")

    def test_boundaryConditions(self):
        self.fvm.applyBoundaryConditions()
        self.assertIsNotNone(self.fvm.boundaryConditions)
        print("Boundary conditions test passed.")

    def test_materialProperties(self):
        self.fvm.loadMaterialProperty()
        self.assertIn('thermalConductivity', self.fvm.material_properties)
        self.assertEqual(self.fvm.material_properties['thermalConductivity'].baseValue, 237)
        print("Material properties test passed.")

    def test_solver(self):
        self.fvm.solveEquations()
        self.assertIsNotNone(self.fvm.solver)
        print("Solver test passed.")

    def test_visualization(self):
        output_path = './testOutput'
        self.fvm.visualizeResults()
        self.assertTrue(os.path.exists(output_path))
        print("Visualization output test passed.")

    def test_fullSimulation(self):
        self.fvm.simulate()
        self.assertIsNotNone(self.fvm.mesh)
        self.assertIsNotNone(self.fvm.boundaryConditions)
        self.assertIn('thermalConductivity', self.fvm.material_properties)
        self.assertIsNotNone(self.fvm.solver)
        output_path = './testOutput'
        self.assertTrue(os.path.exists(output_path))
        print("Full simulation test passed.")
