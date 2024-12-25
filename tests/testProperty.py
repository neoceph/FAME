import unittest
import numpy as np

from src.FVM.property import MaterialProperty

class TestMaterialProperty(unittest.TestCase):
    def test_constant_model(self):
        prop = MaterialProperty('thermal_conductivity', baseValue=200, referenceTemperature=298.15, method='constant')
        self.assertEqual(prop.evaluate(400), 200)
    
    def test_linear_model(self):
        thermalConductivity = MaterialProperty('thermal_conductivity', baseValue=200, referenceTemperature=298.15, method='linear', coefficients=[1e-3])
        self.assertAlmostEqual(thermalConductivity.evaluate(350), 210.37, places=1)

    def test_polynomial_model(self):
        thermalConductivity = MaterialProperty('thermal_conductivity', baseValue=200, referenceTemperature=298.15, method='polynomial', coefficients=[1e-3, -2e-6, 1e-9])
        self.assertAlmostEqual(thermalConductivity.evaluate(350), 537.66, places=1)

    def test_exponential_model(self):
        thermalConductivity = MaterialProperty('diffusivity', 1.0, method='exponential', coefficients=[1e-2])
        self.assertAlmostEqual(thermalConductivity.evaluate(350), 1.6795, places=3)

    def test_invalid_method(self):
        thermalConductivity = MaterialProperty('thermal_conductivity', baseValue=200, referenceTemperature=298.15, method='unknown')
        with self.assertRaises(ValueError):
            thermalConductivity.evaluate(350)