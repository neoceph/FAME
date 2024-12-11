import os
import shutil
import unittest
import numpy as np
from src.FVM.mesh import StructuredMesh
from src.FVM.visualization import MeshWriter


class TestMeshWriter(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Set up a StructuredMesh and output directory for testing.
        """
        cls.outputDir = "testOutput"
        cls.mesh = StructuredMesh(bounds=((0, 1), (0, 1), (0, 1)), divisions=(10, 10, 10))
        cls.writer = MeshWriter(cls.mesh)
        if os.path.exists(cls.outputDir):
            shutil.rmtree(cls.outputDir)
        os.makedirs(cls.outputDir)

    @classmethod
    def tearDownClass(cls):
        """
        Clean up the test output directory.
        """
        if os.path.exists(cls.outputDir):
            shutil.rmtree(cls.outputDir)

    def testSteadyState(self):
        """
        Test writing steady-state scalar, vector, and tensor fields.
        """
        variables = {
            "temperature": np.random.rand(self.mesh.GetNumberOfPoints()),  # Scalar
            "velocity": np.random.rand(self.mesh.GetNumberOfPoints(), 3),  # Vector
            "stressTensor": np.random.rand(self.mesh.GetNumberOfPoints(), 9)  # Tensor
        }
        self.writer.writeVTS(self.outputDir, variables)

        pvdFile = os.path.join(self.outputDir, f"{os.path.basename(self.outputDir)}.pvd")
        vtsFile = os.path.join(self.outputDir, "output_0000.vts")

        self.assertTrue(os.path.exists(pvdFile), "PVD file not created for steady-state.")
        self.assertTrue(os.path.exists(vtsFile), "VTS file not created for steady-state.")

    def testTransientState(self):
        """
        Test writing transient-state scalar, vector, and tensor fields for multiple timesteps.
        """
        timeSteps = [0.0, 1.0, 2.0]
        for step, time in enumerate(timeSteps):
            variables = {
                "temperature": np.random.rand(self.mesh.GetNumberOfPoints()),  # Scalar
                "velocity": np.random.rand(self.mesh.GetNumberOfPoints(), 3),  # Vector
                "stressTensor": np.random.rand(self.mesh.GetNumberOfPoints(), 9)  # Tensor
            }
            self.writer.writeVTS(self.outputDir, variables, time=time, step=step)

        pvdFile = os.path.join(self.outputDir, f"{os.path.basename(self.outputDir)}.pvd")
        self.assertTrue(os.path.exists(pvdFile), "PVD file not created for transient-state.")

        for step in range(len(timeSteps)):
            vtsFile = os.path.join(self.outputDir, f"output_{step:04d}.vts")
            self.assertTrue(os.path.exists(vtsFile), f"VTS file not created for timestep {step}.")

        with open(pvdFile, 'r') as f:
            pvdContent = f.read()
        for step, time in enumerate(timeSteps):
            self.assertIn(f'timestep="{time}"', pvdContent, f"Time {time} not recorded in PVD file.")
            self.assertIn(f'output_{step:04d}.vts', pvdContent, f"File output_{step:04d}.vts not recorded in PVD file.")

    def testSteadyStateWithFixedValues(self):
        """
        Test writing steady-state data with fixed scalar, vector, and tensor fields.
        """
        scalar = np.ones(self.mesh.GetNumberOfPoints()) * 100  # Scalar: temperature = 100 everywhere
        vector = np.zeros((self.mesh.GetNumberOfPoints(), 3))  # Vector: velocity = (0, 0, 0)
        tensor = np.eye(3).flatten()  # Tensor: 3x3 identity matrix

        variables = {
            "temperature": scalar,
            "velocity": vector,
            "stressTensor": np.tile(tensor, (self.mesh.GetNumberOfPoints(), 1))  # Tensor repeated for all points
        }
        self.writer.writeVTS(self.outputDir, variables)

        pvdFile = os.path.join(self.outputDir, f"{os.path.basename(self.outputDir)}.pvd")
        vtsFile = os.path.join(self.outputDir, "output_0000.vts")

        self.assertTrue(os.path.exists(pvdFile), "PVD file not created for steady-state with fixed values.")
        self.assertTrue(os.path.exists(vtsFile), "VTS file not created for steady-state with fixed values.")

    def testTransientStateWithFixedValues(self):
        """
        Test writing transient-state data with fixed scalar, vector, and tensor fields for multiple timesteps.
        """
        scalar_values = [100, 200, 300]  # Different temperature for each timestep
        vector_values = [
            np.full((self.mesh.GetNumberOfPoints(), 3), i) for i in range(3)
        ]  # Velocity changes for each timestep
        tensor_values = [
            np.tile(np.eye(3).flatten() * i, (self.mesh.GetNumberOfPoints(), 1)) for i in range(1, 4)
        ]  # Tensor values scaled by timestep

        timeSteps = [0.0, 1.0, 2.0]
        for step, time in enumerate(timeSteps):
            variables = {
                "temperature": np.ones(self.mesh.GetNumberOfPoints()) * scalar_values[step],  # Scalar
                "velocity": vector_values[step],  # Vector
                "stressTensor": tensor_values[step]  # Tensor
            }
            self.writer.writeVTS(self.outputDir, variables, time=time, step=step)

        pvdFile = os.path.join(self.outputDir, f"{os.path.basename(self.outputDir)}.pvd")
        self.assertTrue(os.path.exists(pvdFile), "PVD file not created for transient-state with fixed values.")

        for step in range(len(timeSteps)):
            vtsFile = os.path.join(self.outputDir, f"output_{step:04d}.vts")
            self.assertTrue(os.path.exists(vtsFile), f"VTS file not created for timestep {step} with fixed values.")

        with open(pvdFile, 'r') as f:
            pvdContent = f.read()
        for step, time in enumerate(timeSteps):
            self.assertIn(f'timestep="{time}"', pvdContent, f"Time {time} not recorded in PVD file.")
            self.assertIn(f'output_{step:04d}.vts', pvdContent, f"File output_{step:04d}.vts not recorded in PVD file.")

    def testSteadyStateWithIncreasingValues(self):
        """
        Test writing steady-state data with monotonously increasing scalar, vector, and tensor fields.
        """
        num_points = self.mesh.GetNumberOfPoints()
        scalar = np.arange(1, num_points + 1)  # Scalar: monotonically increasing values
        vector = np.tile(np.arange(1, 4), (num_points, 1))  # Vector: same (1, 2, 3) for all points
        tensor = np.tile(np.arange(1, 10), (num_points, 1))  # Tensor: same [1, ..., 9] for all points

        variables = {
            "temperature": scalar,
            "velocity": vector,
            "stressTensor": tensor
        }
        self.writer.writeVTS(self.outputDir, variables)

        pvdFile = os.path.join(self.outputDir, f"{os.path.basename(self.outputDir)}.pvd")
        vtsFile = os.path.join(self.outputDir, "output_0000.vts")

        self.assertTrue(os.path.exists(pvdFile), "PVD file not created for steady-state with increasing values.")
        self.assertTrue(os.path.exists(vtsFile), "VTS file not created for steady-state with increasing values.")

    def testTransientStateWithIncreasingValues(self):
        """
        Test writing transient-state data with monotonously increasing scalar, vector, and tensor fields for multiple timesteps.
        """
        num_points = self.mesh.GetNumberOfPoints()
        timeSteps = [0.0, 1.0, 2.0]

        for step, time in enumerate(timeSteps):
            scalar = np.arange(step * num_points + 1, (step + 1) * num_points + 1)  # Scalar: different range for each step
            vector = np.tile(np.arange(step + 1, step + 4), (num_points, 1))  # Vector: step-dependent (step+1, step+2, step+3)
            tensor = np.tile(np.arange(step + 1, step + 10), (num_points, 1))  # Tensor: step-dependent [step+1, ..., step+9]

            variables = {
                "temperature": scalar,
                "velocity": vector,
                "stressTensor": tensor
            }
            self.writer.writeVTS(self.outputDir, variables, time=time, step=step)

        pvdFile = os.path.join(self.outputDir, f"{os.path.basename(self.outputDir)}.pvd")
        self.assertTrue(os.path.exists(pvdFile), "PVD file not created for transient-state with increasing values.")

        for step in range(len(timeSteps)):
            vtsFile = os.path.join(self.outputDir, f"output_{step:04d}.vts")
            self.assertTrue(os.path.exists(vtsFile), f"VTS file not created for timestep {step} with increasing values.")

        with open(pvdFile, 'r') as f:
            pvdContent = f.read()
        for step, time in enumerate(timeSteps):
            self.assertIn(f'timestep="{time}"', pvdContent, f"Time {time} not recorded in PVD file.")
            self.assertIn(f'output_{step:04d}.vts', pvdContent, f"File output_{step:04d}.vts not recorded in PVD file.")

    def testAppendPVD(self):
        """
        Test appending to an existing .pvd file and adding a new .vts file.
        """
        # Initial variables for timestep 0
        variables_step0 = {
            "temperature": np.ones(self.mesh.GetNumberOfPoints()) * 100,  # Scalar
            "velocity": np.zeros((self.mesh.GetNumberOfPoints(), 3)),  # Vector
            "stressTensor": np.tile(np.eye(3).flatten(), (self.mesh.GetNumberOfPoints(), 1))  # Tensor
        }

        # Write the initial .pvd and .vts files
        self.writer.writeVTS(self.outputDir, variables_step0)

        # Verify initial .pvd and .vts file creation
        pvdFile = os.path.join(self.outputDir, f"{os.path.basename(self.outputDir)}.pvd")
        vtsFileStep0 = os.path.join(self.outputDir, "output_0000.vts")
        self.assertTrue(os.path.exists(pvdFile), "PVD file not created for initial timestep.")
        self.assertTrue(os.path.exists(vtsFileStep0), "VTS file not created for initial timestep.")

        # New variables for timestep 1
        variables_step1 = {
            "temperature": np.ones(self.mesh.GetNumberOfPoints()) * 200,  # Scalar
            "velocity": np.ones((self.mesh.GetNumberOfPoints(), 3)) * 10,  # Vector
            "stressTensor": np.tile((np.eye(3) * 2).flatten(), (self.mesh.GetNumberOfPoints(), 1))  # Tensor
        }

        # Append to the existing .pvd file and add a new .vts file
        self.writer.writeVTS(self.outputDir, variables_step1, time=1.0, step=1)

        # Verify new .vts file creation
        vtsFileStep1 = os.path.join(self.outputDir, "output_0001.vts")
        self.assertTrue(os.path.exists(vtsFileStep1), "VTS file not created for appended timestep.")

        # Verify .pvd file content
        with open(pvdFile, 'r') as f:
            pvdContent = f.read()
        self.assertIn('output_0000.vts', pvdContent, "Initial timestep not recorded in PVD file.")
        self.assertIn('output_0001.vts', pvdContent, "Appended timestep not recorded in PVD file.")
        self.assertIn('timestep="0.0"', pvdContent, "Initial timestep not recorded in PVD file.")
        self.assertIn('timestep="1.0"', pvdContent, "Appended timestep not recorded in PVD file.")