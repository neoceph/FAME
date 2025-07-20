Heat Diffusion Transient
========================

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

Discretization in 3D
--------------------

Continuing from the steady state case 

.. image:: ../media/images/FVM/discretization_1D.svg
   :alt: Example SVG
   :scale: 50%
   :align: center

The equation now have the time derivative term which can be integrated over a time step :math:`\Delta t` to yield the following integral form:


.. math::
    :label: eq:heatDiffusionDiscretization-transient
    :nowrap:

    \begin{align}
    \int_{t}^{t+\Delta t} \int_{\Delta V} \rho c \frac{\partial T}{\partial t} dV dt & = \int_{t}^{t+\Delta t} \int_{\Delta V} \left[\frac{\partial}{\partial x}\left(k \frac{\partial T}{\partial x}\right) + S_T\right] dV dt \notag \\
    \rho c (T_{i}-T_{i}^{0}) \Delta V & = \int_{t}^{t+\Delta t} \left( \left[kA \frac{T_{i+1} - T_i}{\Delta x} \right]_{i+1} \right. \notag \\
    & \left. + \left[kA \frac{T_{i-1} - T_i}{\Delta x} \right]_{i-1} + S_u + S_i T_i \right)\ dt
    \end{align}

On the temperature integral of the right side of the equation,we can generalize the equation by means of a weighting parameter :math:`\theta` which can be used to approximate the temperature at the next time step as follows:

.. math::
    :nowrap:

    \begin{align*}
        T_{i}^{n+1} = T_{i}^{n} + \theta \left( T_{i}^{n+1} - T_{i}^{n} \right)
    \end{align*}

where :math:`\theta` is a weighting factor that can be set to 0 for implicit, 0.5 for the Crank-Nicolson method or 1 for the explicit method.

Dropping the superscript for the future time step and using :math:`n = 0` for the current time step, we can express the equation as:

.. math::
    :nowrap:

    \begin{align*}
        T_{i} & = T_{i}^{0} + \theta \left( T_{i} - T_{i}^{0} \right) \\
        T_{i} & = \theta T_{i} + (1-\theta) T_{i}^{0} \\
        \int_{t}^{t+\Delta t} T_{i} dt & = \left[ \theta T_{i} + (1-\theta) T_{i}^{0} \right] \Delta t
    \end{align*}

Thus the equation :eq:`eq:heatDiffusionDiscretization-transient` can be re-written as:

.. math::
    :label: eq:heatDiffusionDiscretized-transient-initial
    :nowrap:

    \begin{align}
        \rho c (T_{i}-T_{i}^{0}) \Delta V & = \theta\left(\left[kA \frac{T_{i+1} - T_i}{\Delta x} \right]_{i+1} + \left[kA \frac{T_{i-1} - T_i}{\Delta x} \right]_{i-1} + S_i T_{i}\right) \Delta t \notag \\
        + (1-\theta) & \left(\left[kA \frac{T_{i+1}^{0} - T_i^{0}}{\Delta x} \right]_{i+1} + \left[kA \frac{T_{i-1}^{0} - T_i^{0}}{\Delta x} \right]_{i-1} + S_i T_{i}^{0} \right) \Delta t \notag \\
        & + S_u \Delta V \Delta t \notag \\
        \Rightarrow \rho c (T_{i}-T_{i}^{0}) \frac{\Delta V}{\Delta t} & = \theta\left(\left[kA \frac{T_{i+1} - T_i}{\Delta x} \right]_{i+1} + \left[kA \frac{T_{i-1} - T_i}{\Delta x} \right]_{i-1} + S_i T_{i}\right) \notag \\
        + (1-\theta) & \left(\left[kA \frac{T_{i+1}^{0} - T_i^{0}}{\Delta x} \right]_{i+1} + \left[kA \frac{T_{i-1}^{0} - T_i^{0}}{\Delta x} \right]_{i-1} + S_i T_{i}^{0} \right) \notag \\
        & + S_u \Delta V
    \end{align}

Rearrenging them to organize all the unknowns on the left side and knowns on the right side, we can express the equation as:

.. math::
    :label: eq:genericTransientDiffusionDiscretization
    :nowrap:

    \begin{align}
        & - \theta \left[ \frac{kA}{\Delta x} \right]_{i-1} T_{i-1} \notag \\
        & + \left[ \frac {\rho c \Delta V}{\Delta t} + \theta \left\{ \left( \frac{kA}{\Delta x} \right)_{i+1} + \left( \frac{kA}{\Delta x} \right)_{i-1} -S_i \right\} \right] T_{i} \notag \\
        & -\theta \left[ \frac{kA}{\Delta x} \right]_{i+1} T_{i+1} \notag \\
        & = (1-\theta) \left[ \frac{kA}{\Delta x} \right]_{i-1} T^{0}_{i-1}\notag \\
        & + \left[ \frac {\rho c \Delta V}{\Delta t} - (1-\theta) \left\{ \left(\frac{kA}{\Delta x} \right)_{i+1} + \left( \frac{kA}{\Delta x} \right)_{i-1} -S_i \right\}\right] T_{i}^{0} \notag \\
        & +(1-\theta) \left[ \frac{kA}{\Delta x} \right]_{i+1} T^{0}_{i+1} + S_u \Delta V
    \end{align}


This is the generic discretization equation for the transient heat diffusion equation using the finite volume method. The equation can be applied to any arbitrary shape mesh cell by appropriately defining the coefficients based on the geometry, cell connectivity, boundary faces, and material properties of the mesh cell.

If :math:`\theta = 0`, the equation becomes explicit, meaning that the temperature at the next time step is calculated directly from the current temperature and source terms. If :math:`\theta = 1`, it becomes implicit, requiring a system of equations to be solved at each time step. For :math:`\theta = 0.5`, it represents the Crank-Nicolson method, which is a time-centered scheme providing a balance between stability and accuracy. The equation is organized such that all the unknowns (temperatures at the next time step) are on the left side, while all known values (temperatures at the current time step and source terms) are on the right side.

For :math:`\theta = 0`, the equation simplifies to an explicit form as follows

.. math::
    :label: eq:genericTransientDiffusionDiscretization-explicit
    :nowrap:

    \begin{align}
        & \left[ \frac {\rho c \Delta V}{\Delta t} \right] T_{i} \notag \\
        & = \left[ \frac{kA}{\Delta x} \right]_{i-1} T^{0}_{i-1}\notag \\
        & + \left[ \frac {\rho c \Delta V}{\Delta t} - \left\{ \left(\frac{kA}{\Delta x} \right)_{i+1} + \left( \frac{kA}{\Delta x} \right)_{i-1} -S_i\right\} \right] T_{i}^{0} \notag \\
        & +\left[ \frac{kA}{\Delta x} \right]_{i+1} T^{0}_{i+1} + S_u \Delta V
    \end{align}

For :math:`\theta = 1`, the equation simplifies to an implicit form as follows

.. math::
    :label: eq:genericTransientDiffusionDiscretization-implicit
    :nowrap:

    \begin{align}
        & - \left[ \frac{kA}{\Delta x} \right] T_{i-1} \notag \\
        & + \left[ \frac {\rho c \Delta V}{\Delta t} + \left\{ \left( \frac{kA}{\Delta x} \right)_{i+1} + \left( \frac{kA}{\Delta x} \right)_{i-1} -S_i \right\} \right] T_{i} \notag \\
        & -\left[ \frac{kA}{\Delta x} \right] T_{i+1} \notag \\
        & = \left[ \frac {\rho c \Delta V}{\Delta t} \right] T_{i}^{0} \notag \\
        & + S_u \Delta V
    \end{align}

For :math:`\theta = \frac{1}{2}`, the equation translates into a Crank-Nicolson form.

Most interestingly, the equation :eq:`eq:genericTransientDiffusionDiscretization` can be used to derive the steady state heat diffusion equation by setting :math:`\Delta t \to \infty` and :math:`\theta = 1` , which leads to the steady state form of the heat diffusion equation.

.. math::
    :label: eq:genericSteadyDiffusionDiscretization
    :nowrap:

    \begin{align}
        & - \left[ \frac{kA}{\Delta x} \right]_{i-1} T_{i-1} \notag \\
        & + \left[ \left( \frac{kA}{\Delta x} \right)_{i+1} + \left( \frac{kA}{\Delta x} \right)_{i-1} -S_i \right] T_{i} \notag \\
        & -\left[ \frac{kA}{\Delta x} \right]_{i+1} T_{i+1} \notag \\
        & = S_u \Delta V
    \end{align}

Boundary conditions
-------------------

Considering the equation :eq:`eq:heatDiffusionDiscretized-transient-initial`, we can estimate the discretization for the different types of boundary conditions.

1. Dirichlet Boundary Condition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a Dirichlet boundary condition, the temperature at the boundary is specified. For example, if the temperature at the boundary is fixed at :math:`T_B`, then :eq:`eq:heatDiffusionDiscretized-transient-initial` becomes:

.. math:: 
    :nowrap:

    \begin{align*}
    \rho c (T_{i}-T_{i}^{0}) \frac{\Delta V}{\Delta t} & = \theta\left(\left[kA \frac{T_{i+1} - T_i}{\Delta x} \right]_{i+1} + \left[kA \frac{T_{B} - T_i}{\Delta x} \right]_{B} + S_i T_{i}\right) \notag \\
            + (1-\theta) & \left(\left[kA \frac{T_{i+1}^{0} - T_i^{0}}{\Delta x} \right]_{i+1} + \left[kA \frac{T_{B}^{0} - T_i^{0}}{\Delta x} \right]_{B} + S_i T_{i}^{0} \right) \notag \\
            & + S_u \Delta V
    \end{align*}

Rearrenging the equation to organize all the unknowns on the left side and knowns on the right side similar to :eq:`eq:genericTransientDiffusionDiscretization`, we can express the equation as:

.. math::
    :nowrap:

    \begin{align}
        & - \theta \left[ 0 \right]_B T_{i-1} \notag \\
        & + \left[ \frac {\rho c \Delta V}{\Delta t} + \theta \left\{ \left( \frac{kA}{\Delta x} \right)_{i+1} + \left( \frac{kA}{\Delta x} \right)_B -S_i \right\} \right] T_{i} \notag \\
        & -\theta \left[ \frac{kA}{\Delta x} \right]_{i+1} T_{i+1} \notag \\
        & = (1-\theta) \left[ 0 \right]_B T^{0}_{i-1}\notag \\
        & + \left[ \frac {\rho c \Delta V}{\Delta t} - (1-\theta) \left\{ \left(\frac{kA}{\Delta x} \right)_{i+1} + \left( \frac{kA}{\Delta x} \right)_{B} -S_i \right\}\right] T_{i}^{0} \notag \\
        & +(1-\theta) \left[ \frac{kA}{\Delta x} \right]_{i+1} T^{0}_{i+1} \notag \\
        & + S_u \Delta V + \theta \left[ \frac{k A}{\Delta x} \right]_B T_B + (1-\theta) \left[ \frac{k A }{\Delta x} \right]_B T_B^0 \notag
    \end{align}

Since the temperature at the boundary is fixed and known, meaning :math:`T_B = T_B^0`, we can express the equation as:

.. math::
    :nowrap:

    \begin{align}
        & - \theta \left[ 0 \right]_B T_{i-1} \notag \\
        & + \left[ \frac {\rho c \Delta V}{\Delta t} + \theta \left\{ \left( \frac{kA}{\Delta x} \right)_{i+1} + \left( \frac{kA}{\Delta x} \right)_B -S_i \right\} \right] T_{i} \notag \\
        & -\theta \left[ \frac{kA}{\Delta x} \right]_{i+1} T_{i+1} \notag \\
        & = (1-\theta) \left[ 0 \right]_B T^{0}_{i-1}\notag \\
        & + \left[ \frac {\rho c \Delta V}{\Delta t} - (1-\theta) \left\{ \left(\frac{kA}{\Delta x} \right)_{i+1} + \left( \frac{kA}{\Delta x} \right)_{B} -S_i \right\}\right] T_{i}^{0} \notag \\
        & +(1-\theta) \left[ \frac{kA}{\Delta x} \right]_{i+1} T^{0}_{i+1} \notag \\
        & + S_u \Delta V + \theta \left[ \frac{k A }{\Delta x} \right]_B T_B^0 + (1-\theta) \left[ \frac{k A }{\Delta x} \right]_B T_B^0 \notag
    \end{align}

Therefore the coefficients with the term :math:`\theta` cancels out, and the equation simplifies to:

.. math::
    :label: eq:genericTransientDiffusionDiscretization-dirichlet
    :nowrap:

    \begin{align}
        & - 0 \notag \\
        & + \left[ \frac {\rho c \Delta V}{\Delta t} + \theta \left\{ \left( \frac{kA}{\Delta x} \right)_{i+1} + \left( \frac{kA}{\Delta x} \right)_B -S_i \right\} \right] T_{i} \notag \\
        & -\theta \left[ \frac{kA}{\Delta x} \right]_{i+1} T_{i+1} \notag \\
        & = 0\notag \\
        & + \left[ \frac {\rho c \Delta V}{\Delta t} - (1-\theta) \left\{ \left(\frac{kA}{\Delta x} \right)_{i+1} + \left( \frac{kA}{\Delta x} \right)_{B} -S_i \right\}\right] T_{i}^{0} \notag \\
        & +(1-\theta) \left[ \frac{kA}{\Delta x} \right]_{i+1} T^{0}_{i+1} \notag \\
        & + S_i^u \Delta V + \left[ \frac{k A }{\Delta x} \right]_B T_B^0 
    \end{align}

2. Neumann Boundary Condition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The flux at the boundary is specified for the Neumann boundary condition, which can be expressed as:

#. Thermal Insulation :math:`\frac {\partial T}{\partial x} = 0` 
#. a convective heat transfer with coefficient :math:`h` and ambient temperature :math:`T_\infty` as :math:`\frac {\partial T}{\partial x} = hA(T_\infty - T_B)`, 
#. a radiative heat loss with emissivity :math:`\epsilon`, Stefan Boltzmann constant :math:`\epsilon`, and ambient temperature :math:`T_\infty` as :math:`\frac {\partial T}{\partial x} = \sigma \epsilon A(T_\infty^4 - T_B^4)`. 

We can consider the temperature of the boundary node and the cell almost equal as a mean to simplify where :math:`T_B = T_i`. In an attempt to develop a generic boundary condition that covers all three of the cases, we can add all of them together and replace the flux at the boundary. Then the equation :eq:`eq:heatDiffusionDiscretized-transient-initial` becomes:

.. math:: 
    :nowrap:

    \begin{align*}
    \rho c (T_{i}-T_{i}^{0}) \frac{\Delta V}{\Delta t} & = \theta\left(\left[kA \frac{T_{i+1} - T_i}{\Delta x} \right]_{i+1} \notag \right. \\
    & \left. + \left[0 + hA(T_B-T_\infty) + \sigma \epsilon A(T_B^4 - T_\infty^4) \right]_{B} + S_i T_{i} \vphantom{\frac{\Delta V}{\Delta t}}\right) \notag \\
    + (1-\theta) & \left(\left[kA \frac{T^{0}_{i+1} - T^{0}_i}{\Delta x} \right]_{i+1} \notag \right. \\
    & \left. + \left[0 + hA(T_\infty - T_B) + \sigma \epsilon A(T_\infty^4 - T_B^4) \right]_{B} + S_i T_{i}^0 \vphantom{\frac{\Delta V}{\Delta t}}\right) \notag \\
    & + S_u \Delta V
    \end{align*}

Rearrenging the equation to organize all the unknowns on the left side and knowns on the right side similar to :eq:`eq:genericTransientDiffusionDiscretization` where the known boundary values having the coefficient :math:`\theta` cancels out and thus, we can express the equation as:

.. math::
    :label: eq:genericTransientDiffusionDiscretization-neumann
    :nowrap:

    \begin{align}
        & -0 + \left[ \frac {\rho c \Delta V}{\Delta t} + \theta \left\{ \left( \frac{kA}{\Delta x} \right)_{i+1} -S_i \right\} \right] T_{i} \notag \\
        & -\theta \left[ \frac{kA}{\Delta x} \right]_{i+1} T_{i+1} \notag \\
        & = 0 + \left[ \frac {\rho c \Delta V}{\Delta t} - (1-\theta) \left\{ \left(\frac{kA}{\Delta x} \right)_{i+1} -S_i \right\}\right] T_{i}^{0} \notag \\
        & +(1-\theta) \left[ \frac{kA}{\Delta x} \right]_{i+1} T^{0}_{i+1} \notag \\
        & + S_i^u \Delta V + \left[ hA(T_\infty - T_B) + \sigma \epsilon A(T_\infty^4 - T_B^4) \right]_{B}
    \end{align}

Generic Discretization for any dimensions
------------------------------------------

In scenarios of higher dimensions, the cell connectivity is more than two thus the generic equation for cells without any boundary faces are a modified form of the equation :eq:`eq:genericSteadyDiffusionDiscretization` as follows:

.. math::
    :label: eq:genericTransientDiffusionDiscretization-nD
    :nowrap:

    \begin{align}
        & - \theta \sum_{j=1}^n \left[ \frac{kA}{\Delta x} \right]_{i, j} T_{j} \notag \\
        & + \left[ \frac {\rho c \Delta V}{\Delta t} + \theta \left\{ \sum_{j=1}^n \left( \frac{kA}{\Delta x} \right)_{i, j} - S_i \right\} \right] T_{i} \notag \\
        & = (1-\theta) \sum_{j=1}^n \left[ \frac{kA}{\Delta x} \right]_{i-1} T^{0}_{i-1}\notag \\
        & + \left[ \frac {\rho c \Delta V}{\Delta t} - (1-\theta) \left\{ \sum_{j=1}^n \left(\frac{kA}{\Delta x} \right)_{i, j} - S_i \right\}\right] T_{i}^{0} \notag \\
        & + S_i^u \Delta V
    \end{align}

The above equation works for any arbitrary shape mesh cell in any dimensions, where :math:`n` is the number of shared cells (connected cells) with the cell indexed with :math:`i`. The coefficients are defined as follows:

.. math::
    :label: eq:coefficientDiscretization
    :nowrap:

    \begin{align}
        a_{i,j} & = - \theta \sum_{j=1}^n \left[ \frac{kA}{\Delta x} \right]_{i, j} \notag \\
        a_{i,i} & = \left[ \frac {\rho c \Delta V}{\Delta t} + \theta \left\{ \sum_{j=1}^n \left( \frac{kA}{\Delta x} \right)_{i, j} - S_i \right\} \right] \notag \\
        b_i & = (1-\theta) \sum_{j=1}^n \left[ \frac{kA}{\Delta x} \right]_{i, j} T^{0}_{i, j}\notag \\
        & + \left[ \frac {\rho c \Delta V}{\Delta t} - (1-\theta) \left\{ \sum_{j=1}^n \left(\frac{kA}{\Delta x} \right)_{i, j} - S_i \right\}\right] \notag \\
        & + S_i^u \Delta V
    \end{align}

for a generic system of linear algebraic equation of the form:

.. math::
    :nowrap:

    \begin{align*}
        A_{i,j} X_{i} = b_{i}
    \end{align*}

where :math:`A` is the coefficient matrix, :math:`X` is the vector of unknown temperatures, and :math:`b` is the vector of known values.

For a given cell if there are boundary faces, only then special equation is to be used. If it is a Dirichlet boundary conditions

.. math::
    :nowrap:

    \begin{align*}
        & \left[ \frac {\rho c \Delta V}{\Delta t} + \theta \left\{ \sum_{j=1}^n \left( \frac{kA}{\Delta x} \right)_{i, B_j} - S_i \right\} \right] T_{i} \notag \\
        & = \left[ \frac {\rho c \Delta V}{\Delta t} - (1-\theta) \left\{ \sum_{j=1}^n \left(\frac{kA}{\Delta x} \right)_{i, B_j} - S_i \right\}\right] T_{i}^{0} \notag \\
        & + S_i^u \Delta V + \sum_{j=1}^n \left[ \frac{k A T_{j,B}^0 }{\Delta x} \right]_B
    \end{align*}

The :math:`\frac{\rho c \Delta V}{\Delta t}` term should be already added to the coefficient matrix :math:`A` for the cell indexed with :math:`i` and thus, it is unnecessary/redundent to include this term in the equation above. Therefore, the co-efficient update can be simplified to:

.. math::
    :label: eq:coefficientDiscretization-dirichlet
    :nowrap:

    \begin{align}
        a_{i,i} & += \theta \left\{ \sum_{j=1}^n \left( \frac{kA}{\Delta x} \right)_{i, B_j} - S_i \right\} \notag \\
        b_i & += - (1-\theta) \left\{ \sum_{j=1}^n \left(\frac{kA}{\Delta x} \right)_{i, B_j} - S_i \right\} T_i^0 \notag \\
        & + \sum_{j=1}^n \left[ \frac{k A T_{j,B}^0 }{\Delta x} \right]_B
    \end{align}

For Neumann boundary conditions

.. math::
    :nowrap:

    \begin{align*}
        & \left[ \frac {\rho c \Delta V}{\Delta t} + \theta \left\{ - S_i \right\} \right] T_{i} \notag \\
        & = \left[ \frac {\rho c \Delta V}{\Delta t} - (1-\theta) \left\{- S_i \right\}\right] T_{i}^{0} \notag \\
        & + S_i^u \Delta V + \sum_{j=1}^n \left[ hA(T_\infty - T_{j,B}) + \sigma \epsilon A(T_\infty^4 - T_{j, B}^4) \right]_{B}
    \end{align*}

where

.. math::
    :label: eq:coefficientDiscretization-neumann
    :nowrap:

    \begin{align}
        a_{i,i} & += -\theta S_i \notag \\
        b_i & += (1-\theta) S_i T_i^0 \notag \notag \\
        & + \sum_{j=1}^n \left[ hA(T_\infty - T_{j,B}) + \sigma \epsilon A(T_\infty^4 - T_{j, B}^4) \right]_{B}
    \end{align}

Equation :eq:`eq:coefficientDiscretization`, :eq:`eq:coefficientDiscretization-dirichlet`, and :eq:`eq:coefficientDiscretization-neumann` are all we need for the construction of matrix :math:`A` and the known vector :math:`b` for FVM solution for any n-dimensional heat diffusion transient/steady state problem for any boundary conditions. We simply set :math:`\theta = 0` for explicit, :math:`\theta = 1` for implicit, and :math:`\theta = 0.5` for Crank-Nicolson method and :math:`\Delta t = \infty` with ``np.inf`` for steady state cases.

During coding, when considering the connected cells (cells with shared faces), we can use the following logic to determine the coefficients:


1. When iterating over cell

    .. math::
        :nowrap:

        \begin{align*}
            a_{i, i} & += \frac{\rho c}{\Delta t} - \theta S_i \\
            b_i & += \left [ \frac{\rho c}{\Delta t} + (1 - \theta) S_i \right ] T_i^0 + S_i^u \Delta V
        \end{align*}

2. When the inner loop for iteration over connected cells are considered, :math:`i \neq j` with :math:`n` being the number of connected cells
    
    .. math::
        :nowrap:

        \begin{align*}
            a_{i,j} & = -\theta \sum_{j=1}^n \left(\frac{kA}{\Delta x} \right)_{i,j} \\
            a_{i, i} & += \theta \sum_{j=1}^n \left( \frac{kA}{\Delta x} \right)_{i,j} \\
            b_i & += (1-\theta) \sum_{j=1}^n \left [ \left( \frac{kA}{\Delta x} \right)_{i,j} T^{0}_{j} - \sum_{j=1}^n \left( \frac{kA}{\Delta x} \right) T_i^0 \right ]
        \end{align*}

3. When the boundary faces are considered

    i. If Dirichlet boundary condition does not exists and flux is specified (either insulated or convective or radiative)

    .. math::
        :nowrap:

        \begin{align*}
            b_i += &\sum_{j=1}^n \left [ hA_j(T_\infty - T_{B_j}) + \sigma \epsilon A_j(T_\infty^4 - T_{B_j}^4) \right ] \\
        \end{align*}
    
    ii. If Dirichlet boundary condition also exists, them simply 

    .. math::
        :nowrap:

        \begin{align*}
            a_{i, i} & += \theta \sum_{j=1}^n \left( \frac{kA}{\Delta x} \right)_{i,B_j} \\
            b_i & += \sum_{j=1}^n \left( \frac{kA}{\Delta x} \right)_{i,B} T^{0}_{B_j} - (1-\theta) \sum_{j=1}^n \left( \frac{kA}{\Delta x} \right)_{B_j} T_i^0
        \end{align*}
