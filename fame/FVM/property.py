import numpy as np

class MaterialProperty:
    def __init__(self, materialName, properties=None):
        """
        Initialize a material with multiple properties.

        :param materialName: Name of the material (e.g., 'Aluminum')
        :param properties: Dictionary of properties (e.g., {'thermal_conductivity': {...}, 'heat_capacity': {...}})
        """
        self.materialName = materialName
        self.properties = properties if properties else {}

    def add_property(self, propertyName, baseValue, referenceTemperature=298.15, method='constant', coefficients=None):
        """
        Add or update a property for the material.
        
        :param propertyName: Name of the property (e.g., 'thermal_conductivity')
        :param baseValue: Value at reference temperature (referenceTemperature)
        :param referenceTemperature: Reference temperature (default is 298.15 K)
        :param method: Method for temperature dependency ('constant', 'linear', 'polynomial', 'exponential')
            
            i. Constant: :math:`material property = baseValue`
            
            ii. Polynomial: :math:`material property = baseValue \cdot (1+c_0 \cdot {(\Delta T})^n + c_1 \cdot {(\Delta T})^{n-1} + ... + c_n)`

            iii. Exponential: :math:`material property = a_{0} \cdot e^{\\beta\Delta T}`
        
        :param coefficients: Coefficients for polynomial or exponential models
        """
        self.properties[propertyName] = {
            'baseValue': baseValue,
            'referenceTemperature': referenceTemperature,
            'method': method,
            'coefficients': coefficients if coefficients else []
        }

    def evaluate(self, propertyName, temperature):
        """
        Evaluate a specific property at a given temperature.
        
        :param propertyName: Property to evaluate (e.g., 'thermal_conductivity')
        :param temperature: Temperature at which to evaluate the property
        :return: Evaluated property value
        """
        if propertyName not in self.properties:
            raise ValueError(f"Property {propertyName} not found for material {self.materialName}.")

        prop = self.properties[propertyName]
        method = prop['method']
        baseValue = prop['baseValue']
        referenceTemperature = prop['referenceTemperature']
        coefficients = prop['coefficients']

        if method == 'polynomial':
            return self._polynomialModel(temperature, baseValue, referenceTemperature, coefficients)
        elif method == 'exponential':
            return self._exponentialModel(temperature, baseValue, referenceTemperature, coefficients)
        elif method == 'constant':
            return baseValue
        else:
            raise ValueError("Unknown method: choose 'linear', 'polynomial', 'exponential', or 'constant'")

    def _polynomialModel(self, temperature, baseValue, referenceTemperature, coefficients):
        delta_T = temperature - referenceTemperature
        return baseValue * np.polyval(coefficients, delta_T)

    def _exponentialModel(self, temperature, baseValue, referenceTemperature, coefficients):
        beta = coefficients[0] if coefficients else 0
        return baseValue * np.exp(beta * (temperature - referenceTemperature))

    def __repr__(self):
        return f"Material: {self.materialName}, Properties: {list(self.properties.keys())}"
