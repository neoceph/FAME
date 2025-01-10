import os
import shutil
import unittest
import numpy as np
from fame.FVM.mesh import StructuredMesh
from fame.FVM.visualization import MeshWriter


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
        Test writing steady-state scalar, vector, and tensor fields for points and cells.
        """
        num_points = self.mesh.GetNumberOfPoints()
        num_cells = self.mesh.GetNumberOfCells()

        variables = {
            "temperature_point": np.random.rand(num_points),  # Scalar at points
            "velocity_point": np.random.rand(num_points, 3),  # Vector at points
            "stressTensor_point": np.random.rand(num_points, 9),  # Tensor at points
            "temperature_cell": np.random.rand(num_cells),  # Scalar at cells
            "velocity_cell": np.random.rand(num_cells, 3),  # Vector at cells
            "stressTensor_cell": np.random.rand(num_cells, 9)  # Tensor at cells
        }
        self.writer.writeVTS(self.outputDir, variables)

        pvdFile = os.path.join(self.outputDir, f"{os.path.basename(self.outputDir)}.pvd")
        vtsFile = os.path.join(self.outputDir, "output_0000.vts")

        self.assertTrue(os.path.exists(pvdFile), "PVD file not created for steady-state.")
        self.assertTrue(os.path.exists(vtsFile), "VTS file not created for steady-state.")

    def testTransientState(self):
        """
        Test writing transient-state scalar, vector, and tensor fields for points and cells.
        """
        timeSteps = [0.0, 1.0, 2.0]
        num_points = self.mesh.GetNumberOfPoints()
        num_cells = self.mesh.GetNumberOfCells()

        for step, time in enumerate(timeSteps):
            variables = {
                "temperature_point": np.random.rand(num_points),
                "velocity_point": np.random.rand(num_points, 3),
                "stressTensor_point": np.random.rand(num_points, 9),
                "temperature_cell": np.random.rand(num_cells),
                "velocity_cell": np.random.rand(num_cells, 3),
                "stressTensor_cell": np.random.rand(num_cells, 9)
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

    def testSteadyStateWithFixedValues(self):
        """
        Test writing steady-state data with fixed scalar, vector, and tensor fields for points and cells.
        """
        num_points = self.mesh.GetNumberOfPoints()
        num_cells = self.mesh.GetNumberOfCells()

        scalar_point = np.ones(num_points) * 100  # Scalar at points
        vector_point = np.zeros((num_points, 3))  # Vector at points
        tensor_point = np.tile(np.eye(3).flatten(), (num_points, 1))  # Tensor at points

        scalar_cell = np.ones(num_cells) * 200  # Scalar at cells
        vector_cell = np.zeros((num_cells, 3))  # Vector at cells
        tensor_cell = np.tile(np.eye(3).flatten(), (num_cells, 1))  # Tensor at cells

        variables = {
            "temperature_point": scalar_point,
            "velocity_point": vector_point,
            "stressTensor_point": tensor_point,
            "temperature_cell": scalar_cell,
            "velocity_cell": vector_cell,
            "stressTensor_cell": tensor_cell
        }
        self.writer.writeVTS(self.outputDir, variables)

        pvdFile = os.path.join(self.outputDir, f"{os.path.basename(self.outputDir)}.pvd")
        vtsFile = os.path.join(self.outputDir, "output_0000.vts")

        self.assertTrue(os.path.exists(pvdFile), "PVD file not created for steady-state with fixed values.")
        self.assertTrue(os.path.exists(vtsFile), "VTS file not created for steady-state with fixed values.")

    def testSteadyStateWithIncreasingValues(self):
        """
        Test writing steady-state data with monotonously increasing scalar, vector, and tensor fields for points and cells.
        """
        num_points = self.mesh.GetNumberOfPoints()
        num_cells = self.mesh.GetNumberOfCells()

        scalar_point = np.arange(1, num_points + 1)  # Increasing scalar at points
        vector_point = np.tile(np.arange(1, 4), (num_points, 1))  # Vector at points
        tensor_point = np.tile(np.arange(1, 10), (num_points, 1))  # Tensor at points

        scalar_cell = np.arange(1, num_cells + 1)  # Increasing scalar at cells
        vector_cell = np.tile(np.arange(1, 4), (num_cells, 1))  # Vector at cells
        tensor_cell = np.tile(np.arange(1, 10), (num_cells, 1))  # Tensor at cells

        variables = {
            "temperature_point": scalar_point,
            "velocity_point": vector_point,
            "stressTensor_point": tensor_point,
            "temperature_cell": scalar_cell,
            "velocity_cell": vector_cell,
            "stressTensor_cell": tensor_cell
        }
        self.writer.writeVTS(self.outputDir, variables)

        pvdFile = os.path.join(self.outputDir, f"{os.path.basename(self.outputDir)}.pvd")
        vtsFile = os.path.join(self.outputDir, "output_0000.vts")

        self.assertTrue(os.path.exists(pvdFile), "PVD file not created for steady-state with increasing values.")
        self.assertTrue(os.path.exists(vtsFile), "VTS file not created for steady-state with increasing values.")

    def testTransientStateWithIncreasingValues(self):
        """
        Test writing transient-state data with monotonously increasing scalar, vector, and tensor fields for points and cells.
        """
        num_points = self.mesh.GetNumberOfPoints()
        num_cells = self.mesh.GetNumberOfCells()
        timeSteps = [0.0, 1.0, 2.0]

        for step, time in enumerate(timeSteps):
            # Scalar field with increasing values per point/cell
            scalar_point = np.arange(1, num_points + 1) * (step + 1)
            scalar_cell = np.arange(1, num_cells + 1) * (step + 1)

            # Vector field with increasing values per point/cell
            vector_point = np.arange(1, num_points * 3 + 1).reshape(num_points, 3) * (step + 1)
            vector_cell = np.arange(1, num_cells * 3 + 1).reshape(num_cells, 3) * (step + 1)

            # Tensor field with increasing values per point/cell
            tensor_point = np.arange(1, num_points * 9 + 1).reshape(num_points, 9) * (step + 1)
            tensor_cell = np.arange(1, num_cells * 9 + 1).reshape(num_cells, 9) * (step + 1)


            variables = {
                "temperature_point": scalar_point,
                "velocity_point": vector_point,
                "stressTensor_point": tensor_point,
                "temperature_cell": scalar_cell,
                "velocity_cell": vector_cell,
                "stressTensor_cell": tensor_cell
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
        Test appending to an existing .pvd file and adding a new .vts file for points and cells.
        """
        num_points = self.mesh.GetNumberOfPoints()
        num_cells = self.mesh.GetNumberOfCells()

        # Initial timestep 0
        variables_step0 = {
            "temperature_point": np.ones(num_points) * 100,
            "velocity_point": np.zeros((num_points, 3)),
            "stressTensor_point": np.tile(np.eye(3).flatten(), (num_points, 1)),
            "temperature_cell": np.ones(num_cells) * 150,
            "velocity_cell": np.zeros((num_cells, 3)),
            "stressTensor_cell": np.tile(np.eye(3).flatten(), (num_cells, 1))
        }
        self.writer.writeVTS(self.outputDir, variables_step0)

        pvdFile = os.path.join(self.outputDir, f"{os.path.basename(self.outputDir)}.pvd")
        vtsFileStep0 = os.path.join(self.outputDir, "output_0000.vts")
        self.assertTrue(os.path.exists(pvdFile), "PVD file not created for initial timestep.")
        self.assertTrue(os.path.exists(vtsFileStep0), "VTS file not created for initial timestep.")

        # Timestep 1 - Append new data
        variables_step1 = {
            "temperature_point": np.ones(num_points) * 200,
            "velocity_point": np.ones((num_points, 3)) * 10,
            "stressTensor_point": np.tile((np.eye(3) * 2).flatten(), (num_points, 1)),
            "temperature_cell": np.ones(num_cells) * 250,
            "velocity_cell": np.ones((num_cells, 3)) * 5,
            "stressTensor_cell": np.tile((np.eye(3) * 2).flatten(), (num_cells, 1))
        }
        self.writer.writeVTS(self.outputDir, variables_step1, time=1.0, step=1)

        vtsFileStep1 = os.path.join(self.outputDir, "output_0001.vts")
        self.assertTrue(os.path.exists(vtsFileStep1), "VTS file not created for appended timestep.")

        # Verify PVD contents
        with open(pvdFile, 'r') as f:
            pvdContent = f.read()
        self.assertIn('output_0000.vts', pvdContent, "Initial timestep not recorded in PVD file.")
        self.assertIn('output_0001.vts', pvdContent, "Appended timestep not recorded in PVD file.")
        self.assertIn('timestep="0.0"', pvdContent, "Initial timestep not recorded in PVD file.")
        self.assertIn('timestep="1.0"', pvdContent, "Appended timestep not recorded in PVD file.")


class TestMeshWriter1D(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """
        Set up a 1D StructuredMesh and output directory for testing.
        """
        cls.outputDir = "testOutput1D"
        cls.mesh = StructuredMesh(bounds=(0, 5), divisions=(10,))  # 1D mesh
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
        Test writing steady-state scalar and vector fields for 1D mesh.
        """
        num_points = self.mesh.GetNumberOfPoints()
        num_cells = self.mesh.GetNumberOfCells()

        variables = {
            "temperature_point": np.random.rand(num_points),
            "velocity_point": np.random.rand(num_points, 3),  # Vector in 1D
            "tensor_point": np.random.rand(num_points, 9),  # tensor in 1D
            "temperature_cell": np.random.rand(num_cells),
            "velocity_cell": np.random.rand(num_cells, 3),
            "tensor_cell": np.random.rand(num_cells, 9),
        }
        self.writer.writeVTS(self.outputDir, variables)

        pvdFile = os.path.join(self.outputDir, f"{os.path.basename(self.outputDir)}.pvd")
        vtpFile = os.path.join(self.outputDir, "output_0000.vtp")

        self.assertTrue(os.path.exists(pvdFile), "PVD file not created for steady-state.")
        self.assertTrue(os.path.exists(vtpFile), "VTP file not created for steady-state.")

    def testTransientState(self):
        """
        Test writing transient-state scalar, vector (3D), and tensor (9D) fields for 1D mesh.
        """
        timeSteps = [0.0, 1.0, 2.0]
        num_points = self.mesh.GetNumberOfPoints()
        num_cells = self.mesh.GetNumberOfCells()

        for step, time in enumerate(timeSteps):
            variables = {
                "temperature_point": np.random.rand(num_points),  # Scalar at points
                "velocity_point": np.random.rand(num_points, 3),  # 3D Vector at points
                "stressTensor_point": np.random.rand(num_points, 9),  # 9D Tensor at points
                "temperature_cell": np.random.rand(num_cells),  # Scalar at cells
                "velocity_cell": np.random.rand(num_cells, 3),  # 3D Vector at cells
                "stressTensor_cell": np.random.rand(num_cells, 9)  # 9D Tensor at cells
            }
            self.writer.writeVTS(self.outputDir, variables, time=time, step=step)

        pvdFile = os.path.join(self.outputDir, f"{os.path.basename(self.outputDir)}.pvd")
        self.assertTrue(os.path.exists(pvdFile), "PVD file not created for transient-state.")

        # Check if VTP files are created for each timestep
        for step in range(len(timeSteps)):
            vtpFile = os.path.join(self.outputDir, f"output_{step:04d}.vtp")
            self.assertTrue(os.path.exists(vtpFile), f"VTP file not created for timestep {step}.")

        # Validate that the PVD file records each timestep and VTP file
        with open(pvdFile, 'r') as f:
            pvdContent = f.read()
        for step, time in enumerate(timeSteps):
            self.assertIn(f'timestep="{time}"', pvdContent, f"Time {time} not recorded in PVD file.")
            self.assertIn(f'output_{step:04d}.vtp', pvdContent, f"File output_{step:04d}.vtp not recorded in PVD file.")


    def testSteadyStateWithFixedValues(self):
        """
        Test writing steady-state scalar, vector (3D), and tensor (9D) fields for 1D mesh.
        """
        num_points = self.mesh.GetNumberOfPoints()
        num_cells = self.mesh.GetNumberOfCells()

        variables = {
            "temperature_point": np.ones(num_points) * 50,  # Scalar at points
            "velocity_point": np.zeros((num_points, 3)),  # 3D Vector at points (zero vector)
            "stressTensor_point": np.tile(np.eye(3).flatten(), (num_points, 1)),  # 9D Tensor at points
            "temperature_cell": np.ones(num_cells) * 75,  # Scalar at cells
            "velocity_cell": np.zeros((num_cells, 3)),  # 3D Vector at cells (zero vector)
            "stressTensor_cell": np.tile(np.eye(3).flatten(), (num_cells, 1))  # 9D Tensor at cells
        }
        
        self.writer.writeVTS(self.outputDir, variables)

        pvdFile = os.path.join(self.outputDir, f"{os.path.basename(self.outputDir)}.pvd")
        vtpFile = os.path.join(self.outputDir, "output_0000.vtp")

        # Check if the PVD and VTP files are created
        self.assertTrue(os.path.exists(pvdFile), "PVD file not created for steady-state with fixed values.")
        self.assertTrue(os.path.exists(vtpFile), "VTP file not created for steady-state with fixed values.")

        # Validate the PVD contents for the recorded file
        with open(pvdFile, 'r') as f:
            pvdContent = f.read()
        self.assertIn(f'output_0000.vtp', pvdContent, "VTP file not recorded in PVD file for steady-state.")
        self.assertIn(f'timestep="0.0"', pvdContent, "Timestep 0.0 not recorded in PVD file for steady-state.")


    def testTransientStateWithIncreasingValues(self):
        """
        Test writing transient-state data with increasing scalar, vector (3D), and tensor (9D) fields for 1D mesh.
        """
        num_points = self.mesh.GetNumberOfPoints()
        num_cells = self.mesh.GetNumberOfCells()
        timeSteps = [0.0, 1.0, 2.0]

        for step, time in enumerate(timeSteps):
            # Increasing scalar for points and cells
            scalar_point = np.arange(step * num_points + 1, (step + 1) * num_points + 1)
            scalar_cell = np.arange(step * num_cells + 1, (step + 1) * num_cells + 1)
            
            # Vector with increasing magnitude (3D vectors for 1D mesh, varying per point)
            vector_point = np.arange(1, num_points * 3 + 1).reshape(num_points, 3) * (step + 1)
            vector_cell = np.arange(1, num_cells * 3 + 1).reshape(num_cells, 3) * (step + 1)

            # Tensor (9D) with increasing values
            tensor_point = np.arange(1, num_points * 9 + 1).reshape(num_points, 9) * (step + 1)
            tensor_cell = np.arange(1, num_cells * 9 + 1).reshape(num_cells, 9) * (step + 1)


            variables = {
                "temperature_point": scalar_point,  # Scalar at points
                "velocity_point": vector_point,  # 3D Vector at points
                "stressTensor_point": tensor_point,  # 9D Tensor at points
                "temperature_cell": scalar_cell,  # Scalar at cells
                "velocity_cell": vector_cell,  # 3D Vector at cells
                "stressTensor_cell": tensor_cell  # 9D Tensor at cells
            }
            self.writer.writeVTS(self.outputDir, variables, time=time, step=step)

        pvdFile = os.path.join(self.outputDir, f"{os.path.basename(self.outputDir)}.pvd")
        self.assertTrue(os.path.exists(pvdFile), "PVD file not created for transient-state with increasing values.")

        # Check if VTP files are created for each timestep
        for step in range(len(timeSteps)):
            vtpFile = os.path.join(self.outputDir, f"output_{step:04d}.vtp")
            self.assertTrue(os.path.exists(vtpFile), f"VTP file not created for timestep {step}.")

        # Validate that the PVD file records each timestep and VTP file
        with open(pvdFile, 'r') as f:
            pvdContent = f.read()
        for step, time in enumerate(timeSteps):
            self.assertIn(f'timestep="{time}"', pvdContent, f"Time {time} not recorded in PVD file.")
            self.assertIn(f'output_{step:04d}.vtp', pvdContent, f"File output_{step:04d}.vtp not recorded in PVD file.")
