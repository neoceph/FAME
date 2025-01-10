import unittest
from fame.FVM.property import MaterialProperty

class TestMaterialProperty(unittest.TestCase):
    def setUp(self):
        # Initialize MaterialProperty instance with multiple properties
        self.aluminum = MaterialProperty('Aluminum')
        self.aluminum.add_property('thermal_conductivity', baseValue=200, method='constant')
        self.aluminum.add_property('thermal_conductivity_polynomial', baseValue=200, method='polynomial', coefficients=[1e-3, 1])
        self.aluminum.add_property('thermal_conductivity_complex', baseValue=200, method='polynomial', coefficients=[1e-3, -2e-6, 1e-9])
        self.aluminum.add_property('diffusivity', baseValue=1.0, method='exponential', coefficients=[1e-2])

    def test_constant_model(self):
        self.assertEqual(self.aluminum.evaluate('thermal_conductivity', 400), 200)
    
    def test_polynomial_model(self):
        result = self.aluminum.evaluate('thermal_conductivity_polynomial', 350)
        self.assertAlmostEqual(result, 210.37, places=1)

    def test_complex_polynomial_model(self):
        result = self.aluminum.evaluate('thermal_conductivity_complex', 350)
        self.assertAlmostEqual(result, 537.66, places=1)

    def test_exponential_model(self):
        result = self.aluminum.evaluate('diffusivity', 350)
        self.assertAlmostEqual(result, 1.6795, places=3)

    def test_invalid_method(self):
        self.aluminum.add_property('invalid_property', baseValue=200, method='unknown')
        with self.assertRaises(ValueError):
            self.aluminum.evaluate('invalid_property', 350)

    def test_missing_property(self):
        with self.assertRaises(ValueError):
            self.aluminum.evaluate('nonexistent_property', 300)

