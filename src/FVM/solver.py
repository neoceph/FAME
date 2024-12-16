import numpy as np
import jax
import jax.numpy as jnp
import scipy.sparse as sp
import matplotlib.pyplot as plt

from petsc4py import PETSc
from jax.experimental.sparse import BCOO

class Solver:
    def __init__(self, A, b, backend="scipy"):
        """
        Initialize the solver with the matrix A, vector b, and backend.

        Parameters:
        A: scipy.sparse matrix (A in Ax = b)
        b: numpy array (b in Ax = b)
        backend: str, one of ["scipy", "jax", "petsc"]
        """
        if not sp.isspmatrix(A):
            raise TypeError("A must be a scipy sparse matrix.")
        if not isinstance(b, np.ndarray):
            raise TypeError("b must be a numpy array.")

        self.A = A
        self.b = b
        self.solution = None
        self.backend = backend.lower()

        if self.backend not in ["scipy", "jax", "petsc"]:
            raise ValueError("Unsupported backend. Choose from 'scipy', 'jax', or 'petsc'.")

    def solve(self, method="bicgstab", preconditioner="none"):
        """
        Solve the system Ax = b using the selected backend and method.

        Parameters:
            method: str, optional (default="bicgstab")
                The solver method to use (e.g., "bicgstab", "cg", "gmres").
            preconditioner: str, optional (default="none")
                Preconditioner type (e.g., "jacobi") or "none" for no preconditioning.
                For PETSc, passed directly to pc.setType().

        Returns:
            solution: numpy array
                The solution vector.
        """

        if self.backend == "scipy":
            self.solution = self._solve_scipy(method, preconditioner)
        elif self.backend == "jax":
            self.solution = self._solve_jax(method, preconditioner)
        elif self.backend == "petsc":
            self.solution = self._solve_petsc(method, preconditioner)
        return self.solution

    def _solve_scipy(self, method, preconditioner):
        """
        Solve using Scipy's iterative solvers with optional Jacobi preconditioning.
        """
        solverMethods = {
            "bicgstab": sp.linalg.bicgstab,
            "cg": sp.linalg.cg,
            "gmres": sp.linalg.gmres
        }
        if method not in solverMethods:
            raise ValueError(f"Unsupported method '{method}' for scipy backend.")

        if preconditioner == "jacobi":
        # Construct the Jacobi preconditioner
            jacobi_diag = self.A.diagonal()
            if np.any(jacobi_diag == 0):
                raise ValueError("Jacobi preconditioner cannot be constructed: zero diagonal entries.")
            preconditioner_fn = sp.linalg.LinearOperator(
                dtype=self.A.dtype,
                shape=self.A.shape,
                matvec=lambda x: x / jacobi_diag,
            )
        elif preconditioner == "none":
            preconditioner_fn = None
        else:
            raise ValueError(f"Unsupported preconditioner '{preconditioner}' for scipy backend.")

        solution, info = solverMethods[method](self.A, self.b, rtol=1e-10, atol=1e-10, maxiter=None, M=preconditioner_fn)
        err = np.linalg.norm(self.A @ solution - self.b)
        print(f"Scipy {method} solver residual: {err}")
        return solution, err, info

    def _solve_jax(self, method, preconditioner):
        """
        Solve using JAX's iterative solvers with optional Jacobi preconditioning.
        """
        solver_methods = {
            "bicgstab": jax.scipy.sparse.linalg.bicgstab,
            "cg": jax.scipy.sparse.linalg.cg,
            "gmres": jax.scipy.sparse.linalg.gmres
        }
        if method not in solver_methods:
            raise ValueError(f"Unsupported method '{method}' for JAX backend. Supported methods: {list(solver_methods.keys())}.")

        # Convert scipy sparse matrix to JAX sparse matrix
        A_jax = BCOO.from_scipy_sparse(self.A).sort_indices()

        # Prepare the right-hand side vector
        b_jax = jnp.array(self.b)

        # Create preconditioner (Jacobi diagonal scaling)
        if preconditioner == "jacobi":
            jacobi_diag = jnp.array(self.A.diagonal())
            if jnp.any(jacobi_diag == 0):
                raise ValueError("Jacobi preconditioner cannot be constructed: zero diagonal entries.")
            preconditioner_fn = lambda x: x / jacobi_diag
        elif preconditioner == "none":
            preconditioner_fn = None
        else:
            raise ValueError(f"Unsupported preconditioner '{preconditioner}' for JAX backend.")

        # Solve using the selected JAX method
        x0_jax = jnp.zeros_like(b_jax)
        solution, info = solver_methods[method](A_jax, b_jax, tol=1e-10, atol=1e-10, maxiter=None, M=preconditioner_fn, x0=x0_jax, )

        # Flatten the solution if necessary and convert to NumPy
        if isinstance(solution, (tuple, list)):
            solution = solution[0]  # Use the first element of the tuple if applicable
        solution = np.array(solution)

        # Verify convergence
        residual = jnp.linalg.norm(A_jax @ solution - b_jax)
        print(f"JAX {method} solver residual: {residual}")
        if info is not None and info != 0:
            raise RuntimeError(f"JAX solver failed to converge: info={info}")
        
        return solution, residual, info

    def _solve_petsc(self, method, preconditioner):
        """
        Solve using PETSc solver.

        Available Preconditioners:
        - jacobi: Diagonal scaling preconditioner.
        - ilu: Incomplete LU factorization.
        - sor: Successive over-relaxation.
        - none: No preconditioning.
        - asm: Additive Schwarz method.
        - bjacobi: Block Jacobi preconditioner.
        """
        petscMethods = {
            "bicgstab": "bcgs",  # Correct PETSc name for BiCGSTAB
            "cg": "cg",          # PETSc name for Conjugate Gradient
            "gmres": "gmres"     # PETSc name for GMRES
        }
        if method not in petscMethods:
            raise ValueError(f"Unsupported method '{method}' for petsc backend.")

        mat = PETSc.Mat().createAIJ(size=self.A.shape, csr=(self.A.indptr, self.A.indices, self.A.data))
        vec_b = PETSc.Vec().createWithArray(self.b)
        vec_x = PETSc.Vec().createWithArray(np.zeros_like(self.b))

        ksp = PETSc.KSP().create()
        ksp.setOperators(mat)
        ksp.setType(petscMethods[method])  # Use corrected PETSc solver type
        pc = ksp.getPC()
        pc.setType(preconditioner)

        ksp.solve(vec_b, vec_x)
        solution = vec_x.getArray()
        err = np.linalg.norm(self.A @ solution - self.b)
        iteration_number = ksp.getIterationNumber()

        print(f"PETSc {method} solver residual: {err}, Iterations: {iteration_number}")

        return solution, err, iteration_number


    # Utility method to visualize the matrix
    def plotSparseMatrix(self, matrix, filename="matrix.jpeg"):
        """
        Saves a visualization of the sparse matrix as a .jpeg image.

        Parameters:
            matrix (sp.spmatrix): Sparse matrix to visualize.
            filename (str): Output file name for the image.
        """
        if not sp.isspmatrix(matrix):
            raise TypeError("Matrix A must be a scipy sparse matrix.")

        plt.figure(figsize=(8, 8))
        try:
            plt.spy(matrix, markersize=1)
            plt.title("Sparse Matrix Visualization")
            plt.xlabel("Columns")
            plt.ylabel("Rows")
            plt.savefig(filename, format="jpeg")
            plt.close()
            print(f"Sparse matrix plot saved to {filename}")
        except Exception as e:
            plt.close()
            raise RuntimeError(f"Failed to generate the sparse matrix plot: {e}")