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
        
        # Check if Aluminum is loaded as the material
        self.assertIn('Aluminum', self.fvm.material_properties)

        # Access the Aluminum material and its properties
        aluminum = self.fvm.material_properties['Aluminum']
        self.assertIn('thermalConductivity', aluminum.properties)
        
        # Verify the base value, method, and reference temperature of thermal conductivity
        thermal_conductivity = aluminum.properties['thermalConductivity']
        self.assertEqual(thermal_conductivity['baseValue'], 237)
        self.assertEqual(thermal_conductivity['referenceTemperature'], 298.15)
        self.assertEqual(thermal_conductivity['method'], 'polynomial')

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
