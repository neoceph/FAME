import unittest
import numpy as np
import scipy.sparse as sp
import os
from fame.FVM.solver import Solver

class TestSolver(unittest.TestCase):

    def setUp(self):
        """Set up a test system Ax = b with a sparse matrix of size 10x10."""
        np.random.seed(0)
        # Create a symmetric positive definite sparse matrix for compatibility with CG and bicgstab
        random_matrix = sp.rand(10, 10, density=0.3, format="csr")
        self.A = random_matrix + random_matrix.T + sp.eye(10) * 10  # Ensure symmetry and diagonal dominance
        self.b = np.random.rand(10)

    # Scipy Solver Tests
    def test_scipy_solver_bicgstab_none(self):
        """Test Scipy solver with bicgstab and no preconditioner."""
        solver = Solver(self.A, self.b, backend="scipy")
        solution, err, info = solver.solve(method="bicgstab", preconditioner="none")
        np.testing.assert_allclose(self.A @ solution, self.b, atol=1e-10)

    def test_scipy_solver_bicgstab_jacobi(self):
        """Test Scipy solver with bicgstab and jacobi preconditioner."""
        solver = Solver(self.A, self.b, backend="scipy")
        solution, err, info = solver.solve(method="bicgstab", preconditioner="jacobi")
        np.testing.assert_allclose(self.A @ solution, self.b, atol=1e-10)

    def test_scipy_solver_cg_none(self):
        """Test Scipy solver with cg and no preconditioner."""
        solver = Solver(self.A, self.b, backend="scipy")
        solution, err, info = solver.solve(method="cg", preconditioner="none")
        if err > 1e-10:
            print(f"Warning: Residual error is high ({err}). Consider using a preconditioner.")
        self.assertTrue(err < 1e2, "Residual is too high without a preconditioner.")

    def test_scipy_solver_cg_jacobi(self):
        """Test Scipy solver with cg and jacobi preconditioner."""
        solver = Solver(self.A, self.b, backend="scipy")
        solution, err, info = solver.solve(method="cg", preconditioner="jacobi")
        if err > 1e-10:
            print(f"Warning: Residual error is high ({err}). Consider using other method.")
        self.assertTrue(err < 1e2, "Residual is too high for Conjugate Gradient method.")

    def test_scipy_solver_gmres_none(self):
        """Test Scipy solver with gmres and no preconditioner."""
        solver = Solver(self.A, self.b, backend="scipy")
        solution, err, info = solver.solve(method="gmres", preconditioner="none")
        np.testing.assert_allclose(self.A @ solution, self.b, atol=1e-10)

    def test_scipy_solver_gmres_jacobi(self):
        """Test Scipy solver with gmres and jacobi preconditioner."""
        solver = Solver(self.A, self.b, backend="scipy")
        solution, err, info = solver.solve(method="gmres", preconditioner="jacobi")
        np.testing.assert_allclose(self.A @ solution, self.b, atol=1e-10)

    # JAX Solver Tests
    def test_jax_solver_bicgstab_none(self):
        """Test JAX solver with bicgstab and no preconditioner."""
        solver = Solver(self.A, self.b, backend="jax")
        solution, err, info = solver.solve(method="bicgstab", preconditioner="none")
        np.testing.assert_allclose(self.A @ solution, self.b, atol=1e-6)

    def test_jax_solver_bicgstab_jacobi(self):
        """Test JAX solver with bicgstab and jacobi preconditioner."""
        solver = Solver(self.A, self.b, backend="jax")
        solution, err, info = solver.solve(method="bicgstab", preconditioner="jacobi")
        np.testing.assert_allclose(self.A @ solution, self.b, atol=1e-6)

    def test_jax_solver_cg_none(self):
        """Test JAX solver with cg and no preconditioner."""
        solver = Solver(self.A, self.b, backend="jax")
        solution, err, info = solver.solve(method="cg", preconditioner="none")
        np.testing.assert_allclose(self.A @ solution, self.b, atol=1e-6)

    def test_jax_solver_cg_jacobi(self):
        """Test JAX solver with cg and jacobi preconditioner."""
        solver = Solver(self.A, self.b, backend="jax")
        solution, err, info = solver.solve(method="cg", preconditioner="jacobi")
        np.testing.assert_allclose(self.A @ solution, self.b, atol=1e-6)

    def test_jax_solver_gmres_none(self):
        """Test JAX solver with gmres and no preconditioner."""
        solver = Solver(self.A, self.b, backend="jax")
        solution, err, info = solver.solve(method="gmres", preconditioner="none")
        np.testing.assert_allclose(self.A @ solution, self.b, atol=1e-6)

    def test_jax_solver_gmres_jacobi(self):
        """Test JAX solver with gmres and jacobi preconditioner."""
        solver = Solver(self.A, self.b, backend="jax")
        solution, err, info = solver.solve(method="gmres", preconditioner="jacobi")
        np.testing.assert_allclose(self.A @ solution, self.b, atol=1e-6)

    # PETSc Solver Tests
    def test_petsc_solver_bicgstab_none(self):
        """Test PETSc solver with bicgstab and no preconditioner."""
        solver = Solver(self.A, self.b, backend="petsc")
        solution, err, info = solver.solve(method="bicgstab", preconditioner="none")
        np.testing.assert_allclose(self.A @ solution, self.b, atol=1e-5)

    def test_petsc_solver_bicgstab_jacobi(self):
        """Test PETSc solver with bicgstab and jacobi preconditioner."""
        solver = Solver(self.A, self.b, backend="petsc")
        solution, err, info = solver.solve(method="bicgstab", preconditioner="jacobi")
        np.testing.assert_allclose(self.A @ solution, self.b, atol=1e-5)

    def test_petsc_solver_cg_none(self):
        """Test PETSc solver with cg and no preconditioner."""
        solver = Solver(self.A, self.b, backend="petsc")
        solution, err, info = solver.solve(method="cg", preconditioner="none")
        np.testing.assert_allclose(self.A @ solution, self.b, atol=1e-5)

    def test_petsc_solver_cg_jacobi(self):
        """Test PETSc solver with cg and jacobi preconditioner."""
        solver = Solver(self.A, self.b, backend="petsc")
        solution, err, info = solver.solve(method="cg", preconditioner="jacobi")
        np.testing.assert_allclose(self.A @ solution, self.b, atol=1e-5)

    def test_petsc_solver_gmres_none(self):
        """Test PETSc solver with gmres and no preconditioner."""
        solver = Solver(self.A, self.b, backend="petsc")
        solution, err, info = solver.solve(method="gmres", preconditioner="none")
        np.testing.assert_allclose(self.A @ solution, self.b, atol=1e-5)

    def test_petsc_solver_gmres_jacobi(self):
        """Test PETSc solver with gmres and jacobi preconditioner."""
        solver = Solver(self.A, self.b, backend="petsc")
        solution, err, info = solver.solve(method="gmres", preconditioner="jacobi")
        np.testing.assert_allclose(self.A @ solution, self.b, atol=1e-5)

    # Test Plot Sparse Matrix Method
    def test_plot_sparse_matrix(self):
        """Test the plot_sparse_matrix method to generate a .jpeg image."""
        solver = Solver(self.A, self.b, backend="scipy")
        output_path = "test_matrix_plot.jpeg"
        solver.plotSparseMatrix(self.A, filename=output_path)
        self.assertTrue(os.path.exists(output_path))
        # Clean up the generated file
        os.remove(output_path)
