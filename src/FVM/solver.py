import numpy as onp
import jax
import jax.numpy as np
import jax.experimental.sparse as jsp
from petsc4py import PETSc
import scipy
import scipy.sparse as sp
import time
import os
gpu_idx = "3"
os.environ["CUDA_VISIBLE_DEVICES"] = gpu_idx

'''
solver module: 3 options

scipy cpu solver;
jax gpu iterative solver with Jacobi preconditioner
petsc cpu solver (mpi-cpu and mpi gpu)

'''
###scipy solver with Jacobi preconditioner
def jacobi_preconditioner(A_sp_scipy):
    print(f"Compute and use jacobi preconditioner")
    jacobi = np.array(A_sp_scipy.diagonal())
    ## for dirichlet bc nodes, might need to use 1 at the diagonal
    return jacobi


def get_jacobi_precond(jacobi):

    def jacobi_precond(x):
        return x * (1. / jacobi)

    return jacobi_precond

def scipy_solve(A_fn, b, x0, precond: bool, pc_matrix=None):
    """Solves the equilibrium equation using a JAX solver.
    Is fully traceable and runs on GPU.

    Parameters
    ----------
    A_fn: scipy sparse matrix
    precond
        Whether to calculate the preconditioner or not
    pc_matrix
        The matrix to use as preconditioner
    """
    pc = get_jacobi_precond(
        jacobi_preconditioner(A_fn)) if precond else None
    x, info = scipy.sparse.linalg.bicgstab(A_fn,
                                           b,
                                           x0=x0,
                                           M=pc,
                                           tol=1e-10,
                                           atol=1e-10,
                                           maxiter=10000)

    # Verify convergence
    err = np.linalg.norm(A_fn@x - b)
    print(f"JAX scipy linear solve res = {err}")
    return x


###jax solver with Jacobi preconditioner
def jax_solve(A_sp_scipy, A_fn, b, x0, precond: bool, pc_matrix=None):
    """Solves the equilibrium equation using a JAX solver.
    Is fully traceable and runs on GPU.

    Parameters
    ----------
    A_fn
        2D array or function that calculates the linear map (matrix-vector product) Ax when called like A(x) or A @ x
        may need special treatment for dirichlet bc
    precond
        Whether to calculate the preconditioner or not
    pc_matrix
        The matrix to use as preconditioner
    """
    pc = get_jacobi_precond(
        jacobi_preconditioner(A_sp_scipy)) if precond else None
    x, info = jax.scipy.sparse.linalg.bicgstab(A_fn,
                                               b,
                                               x0=x0,
                                               M=pc,
                                               tol=1e-10,
                                               atol=1e-10,
                                               maxiter=10000)

    # Verify convergence
    err = np.linalg.norm(A_fn(x) - b)
    print(f"JAX scipy linear solve res = {err}")

    # assert err < 0.1, f"JAX linear solver failed to converge with err = {err}"
    # x = np.where(err < 0.1, x, np.nan) # For assert purpose, some how this also affects bicgstab.

    return x


###petsc solver
'''
input: csr matrix
output: solution vector
dealing with Dirichlet bc: use mat.zeroRows(List(bc_nodes))
'''
def mat_cuda_from_scipy_csr(csr):
    """
    Convert a scipy.sparse.csr_matrix to a PETSc.Mat.
    csr: scipy.sparse.csr_matrix
    return PETSc.Mat
    """
    mat = PETSc.Mat().create()
    mat.setType(PETSc.Mat.Type.AIJCUSPARSE)
    mat.setSizes(csr.shape)
    mat.setValuesCSR(csr.indptr, csr.indices, csr.data)
    mat.setUp()
    mat.assemble()
    return mat

def mat_cuda_from_jax_csr(bcsr):
    """
    Convert a scipy.sparse.csr_matrix to a PETSc.Mat.
    csr: scipy.sparse.csr_matrix
    return PETSc.Mat
    """
    mat = PETSc.Mat().create()
    mat.setType(PETSc.Mat.Type.AIJCUSPARSE)
    mat.setSizes(bcsr.shape)
    mat.setValuesCSR(bcsr.indptr, bcsr.indices, bcsr.data)
    mat.setUp()
    mat.assemble()
    return mat

def mat_cpu_from_scipy_csr(csr):
    """
    Convert a scipy.sparse.csr_matrix to a PETSc.Mat.
    csr: scipy.sparse.csr_matrix
    return PETSc.Mat
    """
    mat = PETSc.Mat().create()
    mat.setType(PETSc.Mat.Type.AIJ)
    mat.setSizes(csr.shape)
    mat.setValuesCSR(csr.indptr, csr.indices, csr.data)
    mat.setUp()
    mat.assemble()
    return mat

def array_np_to_petsc(numpy_array):
    """
    Return a PETSc.Vec from a numpy array.
    """
    vec = PETSc.Vec().createWithArray(numpy_array)
    vec.setUp()
    vec.assemble()
    return vec

def array_jnp_to_petsc(numpy_array):
    """
    Return a PETSc.Vec from a numpy array.
    """
    vec = PETSc.Vec().create()
    vec = PETSc.Vec().createCUDAWithArrays(cpuarray = numpy_array)
    # vec.setType(PETSc.Vec.Type.CUDA)
    # vec = PETSc.Vec().setValues(jnp.arange(jnp_array.shape[0]), jnp_array)
    vec.setUp()
    vec.assemble()
    return vec

def solve_ax_b(mat, vec, ksp_type, pc_type):
    """
    Pre-conditioned iterative solver
    ksp_type: iterative scheme
    pc_type: precontioner type
    https://petsc.org/main/petsc4py/reference/petsc4py.PETSc.html
    Solve Ax = b using PETSc.
    mat: PETSc.Mat, A in Ax = b
    vec: PETSc.Vec, b in Ax = b
    """
    ksp = PETSc.KSP().create()
    ksp.setOperators(mat)
    ksp.setFromOptions()
    ksp.setType(ksp_type)
    ksp.pc.setType(pc_type)
    x = mat.createVecRight()
    x.assemble()
    ksp.solve(vec, x)
    return x

def solve_ax_b_pclu(mat, vec):
    """
    Solve Ax = b using PETSc.
    mat: PETSc.Mat, A in Ax = b
    vec: PETSc.Vec, b in Ax = b
    """
    pc = PETSc.PC().create()
    pc.setOperators(mat)
    x = pc(vec)
    pc.apply(vec, x)
    return x






