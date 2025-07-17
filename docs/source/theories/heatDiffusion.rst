Heat Diffusion
===============

This section discretizes the heat diffusion equation in three dimensions for any arbitrary shape mesh cell using the finite volume method (FVM) for transient conditions. The discretization is based on the general form of the heat diffusion equation, which includes both dependent and independent source terms.

Generic Form of the Heat Diffusion Equation
-------------------------------------------

The heat diffusion equation in three dimensions can be written as:

.. math::

    \rho c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial x_i}(k \frac{\partial T}{\partial x_i}) + S_T

where:

- :math:`T` is the temperature (dependent variable).
- :math:`\alpha` is the thermal diffusivity (:math:`\alpha = \frac{k}{\rho c_p}`) where

  - :math:`k` is the thermal conductivity.
  - :math:`\rho` is the density of the material.
  - :math:`c_p` is the specific heat capacity.
- :math:`S_T` is a source term.

In steady-state conditions, the time-dependent term vanishes:

.. math::

    \frac{\partial}{\partial x_i}(k \frac{\partial T}{\partial x_i}) + S_T = 0

Discretization in 3D
--------------------

The objective is to follow a generic discretization scheme for the heat diffusion equation in three dimensions, which can be applied to any arbitrary shape mesh cell. Ultimately, the goal is to generalize the transient cases as well allowing the formation of a sparse matrix system of equations for solution using iterative solvers.

Considering a control volume :math:`\Delta V` with a surface area :math:`A` and a normal vector :math:`\mathbf{n}`. The integral form of the heat diffusion equation over the control volume can be expressed as:

.. image:: ../media/images/FVM/discretization_1D.svg
   :alt: Example SVG
   :scale: 50%
   :align: center

Integrating the heat diffusion equation over the control volume, we can express it as:

.. math::

    \int_{\Delta V} \left[\frac{\partial}{\partial x}\left(k \frac{\partial T}{\partial x}\right) + S_T\right] dV = 0

Following the Gauss-divergence theorem and hence finite volume method, we can re-write the integral form of the heat diffusion equation over the control volume as:

.. math::
    \Rightarrow \int_A \left[k \frac{\partial T}{\partial x} \right] \cdot \mathbf{n} \ dA + S_T \Delta V = 0

Considering the i-th cell in the mesh, we can express the above equation as:

.. math::
    \Rightarrow \left[kA \frac{\partial T}{\partial x} \right]_{right} + \left[kA \frac{\partial T}{\partial x} \right]_{left} + S_T \Delta V = 0    

Considering :math:`S_T \Delta V = S_u + S_i T_i` for a dependent source variable

.. math::
    \Rightarrow \left[kA \frac{\partial T}{\partial x} \right]_{R} + \left[kA \frac{\partial T}{\partial x} \right]_{L} + S_u + S_i T_i = 0    

Assuming the convention of **fluxes moving out the cell as positive**, we can express the above equation as:

.. math::
    \Rightarrow \left[kA \frac{T_{i+1} - T_i}{||x_{i+1} - x_{i}||^2} \right]_{R} + \left[kA \frac{T_{i-1} - T_i}{||x_{i} - x_{i-1}||^2} \right]_{L} + S_u + S_i T_i = 0    

or


.. math::
    :label: eq:heatDiscretization
    :nowrap:

    \begin{align}
    \Rightarrow \left[kA \frac{T_{i+1} - T_i}{\Delta x} \right]_{i+1} + \left[kA \frac{T_{i-1} - T_i}{\Delta x} \right]_{i-1} + S_u + S_i T_i = 0   
    \end{align}

Here if needed

- :math:`k|_{i+1}` or :math:`k|_{i-1}` can be approximated as the average of the thermal conductivities at the current and next cell, i.e., :math:`k_{i+1} = \frac{k_i + k_{i+1}}{2}` and :math:`k_{i-1} = \frac{k_i + k_{i-1}}{2}`.
- If needed :math:`A|_{i+1}` can be approximated by considering the area of the face between the current cell :math:`ith` and the next cell :math:`(i+1)th`, i.e., :math:`A_{i, i+1}`.

From the above equation, we can rearrange the terms to form a system of equations suitable to be solved iteratively.:


.. math::
   :nowrap:

   \begin{align*}
   \Rightarrow\ 
   & \left[\left\{ \frac{kA}{\Delta x} \right\}_{i,i-1} \right] T_{i-1} \notag \\
   & -\left[\left\{ \frac{kA}{\Delta x} \right\}_{i,i+1} + \left\{ \frac{kA}{\Delta x} \right\}_{i,i-1} - S_i \right] T_i \notag \\
   & +\left[\left\{ \frac{kA}{\Delta x} \right\}_{i,i+1} \right] T_{i+1} + S_u = 0 
   \end{align*}

Changing the sign of the equation, and for a general case where the source term is not zero (:math:`S_u = q \Delta V = S_i^u \Delta V`), we can express the equation as:


.. math::
    :label: eq:organizedHeatDiffusionDiscretization
    :nowrap:

    \begin{align}
    \Rightarrow\ 
    & -\left[\left\{ \frac{kA}{\Delta x} \right\}_{i, i-1} \right] T_{i-1} \notag \\
    & +\left[\left\{ \frac{kA}{\Delta x} \right\}_{i, i+1} + \left\{ \frac{kA}{\Delta x} \right\}_{i, i-1} - S_i \right] T_i \notag \\
    & -\left[\left\{ \frac{kA}{\Delta x} \right\}_{i, i+1} \right] T_{i+1} = S_i^u \Delta V
    \end{align}

The equation above is actually a linear system of equations and in simplified form can be written as follows which can be expressed in a matrix form.

.. math::
    -a_{i, i-1} T_{i-1} + a_{i, i} T_i - a_{i, i+1} T_{i+1} = b_i

and

.. math::
    Ax = b 

Thus, the coefficients of the matrix :math:`A` and the vector :math:`b` can be defined as follows:



.. math:: 
    :label: eq:matrixCoefficients
    :nowrap:

    \begin{align}
        a_{i, i-1} &= -\left\{ \frac{kA}{\Delta x} \right\}_{i, i-1} \notag \\
        a_{i, i} &= \left( \left\{ \frac{kA}{\Delta x} \right\}_{i, i+1} + \left\{ \frac{kA}{\Delta x} \right\}_{i, i-1} - S_i\right) \notag \\
        a_{i, i+1} &= -\left\{ \frac{kA}{\Delta x} \right\}_{i, i+1} \notag \\
        b_i &= q \Delta V = S_u
    \end{align}

Where :math:`\left\{\Delta x\right\}_{i+1} = ||x_{i+1} - x_i||^2`, :math:`\left\{\Delta x\right\}_{i-1} = ||x_{i-1} - x_i||^2` is the distance between the two cell centers.

The above equations can be generalized for any arbitrary shape mesh cell in three dimensions, where the coefficients of the matrix :math:`A` and the vector :math:`b` are defined based on the connectivity of the cells and the properties of the material. A general assumption can be made that there are a total of :math:`m_i` cells and for a given cell indexed with :math:`i`, it can be connected to :math:`n` other cells indexed with :math:`j`. In such case, the coefficients can be defined as follows:

.. math::
    \begin{aligned}
        a_{ij} &= -\left\{ \frac{kA}{\Delta x} \right\}_{i, j} \notag \\
        a_{i, i} &= \left( \sum_{j=1}^{n} \left\{ \frac{kA}{\Delta x} \right\}_{i, j} - S_i\right) \notag \\
        b_i &= q \Delta V
    \end{aligned}

Boundary Conditions
---------------------

In the case of boundary conditions, the discretization can be modified to account for the specific conditions at the boundaries. For example, if a Dirichlet boundary condition is applied at either the left/right boundary (i.e., fixed temperature), the equation will be modified.

1. Dirichlet Boundary Condition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For a list of boundary faces defined with :math:`\mathcal{B}` and indexed with :math:`1 \le j \le n` for a cell indexed with :math:`i` equation :eq:`eq:heatDiscretization` can be modified as follows:

.. math::
    :nowrap:

    \begin{align*}
        \Rightarrow\ 
        & -\left[ 0 \right] T_{i-1} \notag \\
        & +\left[\left\{ \frac{kA}{\Delta x} \right\}_{i, i+1} + \left\{ \frac{kA}{\Delta x} \right\}_{i, B} - S_i \right] T_i \notag \\
        & -\left[\left\{ \frac{kA}{\Delta x} \right\}_{i, i+1} \right] T_{i+1} = S_i^u \Delta V + \left\{ \frac{kA}{\Delta x} \right\}_{i, B} T_B
    \end{align*}

From the equation, the coefficients of the matrix :math:`A` and the vector :math:`b` can be defined as follows:

.. math::
    :nowrap:

    \begin{align*}
        a_{ij} &= -\left\{ \frac{kA}{\Delta x} \right\}_{i, i+1} \notag \\
        a_{ii} &= \left\{ \frac{kA}{\Delta x} \right\}_{i, i+1} + \left\{ \frac{kA}{\Delta x} \right\}_{i, B} - S_i \notag \\
        b_{i} &= S_i^u \Delta V + \left\{ \frac{kA}{\Delta x} \right\}_{i, B} T_B \notag
    \end{align*}


2. Neumann Boundary Condition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The flux at the boundary is specified for the Neumann boundary condition, which can be expressed as:

#. Thermal Insulation :math:`\frac {\partial T}{\partial x} = 0` 
#. a convective heat transfer with coefficient :math:`h` and ambient temperature :math:`T_\infty` as :math:`\frac {\partial T}{\partial x} = hA(T_\infty - T_B)`, 
#. a radiative heat loss with emissivity :math:`\epsilon`, Stefan Boltzmann constant :math:`\epsilon`, and ambient temperature :math:`T_\infty` as :math:`\frac {\partial T}{\partial x} = \sigma \epsilon A(T_\infty^4 - T_B^4)`. 

We can consider the temperature of the boundary node and the cell almost equal as a mean to simplify where :math:`T_B = T_i`. In an attempt to develop a generic boundary condition that covers all three of the cases, we can add all of them together and replace the flux at the boundary. Then the equation :eq:`eq:organizedHeatDiffusionDiscretization` becomes:

.. math::
    :nowrap:

    \begin{align*}
    \Rightarrow\ 
    & -\left[ 0 \right] T_{i-1} \notag \\
    & +\left[\left\{ \frac{kA}{\Delta x} \right\}_{i, i+1} - S_i \right] T_i \notag \\
    & -\left[\left\{ \frac{kA}{\Delta x} \right\}_{i, i+1} \right] T_{i+1} \notag \\
    & = S_i^u \Delta V + hA(T_\infty - T_B) + \sigma \epsilon A(T_\infty^4 - T_B^4)
    \end{align*}


Generalization
----------------

Generalizing the equation :eq:`eq:matrixCoefficients` for any arbitrary shape mesh cell in n-dimension, we can express the coefficients of the matrix :math:`A` and the vector :math:`b` as follows considering the i-th cell is connected to n other cells:

.. math:: 
    :nowrap:

    \begin{align*}
        a_{i, j} &= - \sum_{j=1}^n \left\{ \frac{kA}{\Delta x} \right\}_{i, j} \notag \\
        a_{i, i} &= \left( \sum_{j=1}^n \left\{ \frac{kA}{\Delta x} \right\}_{i, j} - S_i\right) \notag \\
        b_i &= S_i^u \Delta V
    \end{align*}

For the boundary faces considering Dirichlet condition, the coefficients can be modified as follows:

.. math:: 
    :nowrap:

    \begin{align*}
        a_{i, i} &+= \left( \sum_{j=1}^n \left\{ \frac{kA}{\Delta x} \right\}_{i, j} - S_i\right) \notag \\
        b_i &+= \left\{ \frac{kA}{\Delta x} \right\}_{i, B} T_B
    \end{align*}

If the condition is Neumann, the coefficients can be modified as follows:

.. math:: 
    :nowrap:

    \begin{align*}
        a_{i, i} &+=  - S_i \notag \\
        b_i &+= hA(T_\infty - T_B) + \sigma \epsilon A(T_\infty^4 - T_B^4)
    \end{align*}


Summary
-------

The finite volume discretization of the heat diffusion equation is formulated in a comprehensive fashion considering Dirichlet, Neumann, Heat generation.
