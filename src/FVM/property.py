import numpy as np

class MaterialProperty:
    def __init__(self, propertyName, baseValue, referenceTemperature=298.15, method='constant', coefficients=None):
        """
        Initialize a temperature-dependent property.

        :param propertyName: Name of the property (e.g., 'thermal_conductivity')
        :param baseValue: Value at reference temperature (referenceTemperature)
        :param referenceTemperature: Reference temperature (default is 298.15 K)
        :param method: Method for temperature dependency ('constant', 'linear', 'polynomial', 'exponential')
        :param coefficients: Coefficients for polynomial or exponential models
        """
        self.propertyName = propertyName
        self.baseValue = baseValue
        self.referenceTemperature = referenceTemperature
        self.method = method
        self.coefficients = coefficients if coefficients else []

    def evaluate(self, temperature):
        """
        Evaluate the property at a given temperature.
        
        :param temperature: Temperature at which to evaluate the property
        :return: Evaluated property value
        """
        if self.method == 'linear':
            return self._linearModel(temperature)
        elif self.method == 'polynomial':
            return self._polynomialModel(temperature)
        elif self.method == 'exponential':
            return self._exponentialModel(temperature)
        elif self.method == 'constant':
            return self.baseValue
        else:
            raise ValueError("Unknown method: choose 'linear', 'polynomial', 'exponential', or 'constant'")

    def _linearModel(self, temperature):
        """Linear temperature dependence model."""
        alpha = self.coefficients[0] if self.coefficients else 0
        return self.baseValue * (1 + alpha * (temperature - self.referenceTemperature))

    def _polynomialModel(self, temperature):
        """Polynomial temperature dependence model."""
        delta_T = temperature - self.referenceTemperature
        return self.baseValue * np.polyval(self.coefficients, delta_T)

    def _exponentialModel(self, temperature):
        """
        Exponential temperature dependence model in the form :math:`_{0} \cdot e^{\beta \delta T}`  
        """
        beta = self.coefficients[0] if self.coefficients else 0
        return self.baseValue * np.exp(beta * (temperature - self.referenceTemperature))

    def __repr__(self):
        return f"{self.propertyName} (Method: {self.method}, Base Value: {self.baseValue}, Reference Temp: {self.referenceTemperature})"

