Heat Diffusion
==============

This section focuses on the heat diffusion equation and its discretization in three dimensions, including the generic form of the governing equation and the linearized source term.

Generic Form of the Heat Diffusion Equation
-------------------------------------------

The heat diffusion equation in three dimensions can be written as:

.. math::

    \rho c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial x_i}(k \frac{\partial T}{\partial x_i}) + S_T

where:

- :math:`T` is the temperature (dependent variable).
- :math:`\alpha` is the thermal diffusivity (:math:`\alpha = \frac{k}{\rho c_p}`).
- :math:`S_T` is a source term.

In steady-state conditions, the time-dependent term vanishes:

.. math::

    \frac{\partial}{\partial x_i}(k \frac{\partial T}{\partial x_i}) + S_T = 0

Discretization in 3D
--------------------

The discretization in 1-D follows the scheme as shown in the image. Using the finite volume method (FVM), the diffusion equation is discretized over a control volume. The discretization involves integrating the equation over the control volume and applying Gauss's divergence theorem:

.. image:: ../media/images/FVM/discretization_1D.svg
   :alt: Example SVG
   :scale: 50%
   :align: center
   
.. math::

    \int_{\Delta V} \left[\frac{\partial}{\partial x}\left(k \frac{\partial T}{\partial x}\right) + S_T\right] dV = 0

.. math::
    \Rightarrow \int_A \left[k \frac{\partial T}{\partial x} \right] \cdot \mathbf{n} \ dA + S_T \Delta V = 0

.. math::
    \Rightarrow \left[kA \frac{\partial T}{\partial x} \right]_{A_{right}} - \left[kA \frac{\partial T}{\partial x} \right]_{A_{left}} + S_T \Delta V = 0    

Considering :math:`S_T \Delta V = S_u + S_i T_i` for a dependent source variable

.. math::
    \Rightarrow \left[kA \frac{\partial T}{\partial x} \right]_{A_{right}} - \left[kA \frac{\partial T}{\partial x} \right]_{A_{left}} + S_u + S_i T_i = 0    

Now

- :math:`\left[\frac{\partial T}{\partial x} \right]_{A_{right}} = \frac{T_{i+1}-T_i}{||x_{i+1} - x_i||^2}`
- :math:`\left[\frac{\partial T}{\partial x} \right]_{A_{left}} = \frac{T_{i}-T_{i-1}}{||x_i-x_{i-1}||^2}`
- :math:`k|_{A_{right}}=\frac{k_{i+1}+k_i}{2}`
- :math:`k|_{A_{left}}=\frac{k_{i}+k_{i-1}}{2}`


.. math::
    
    \begin{align*}
        & \left[kA \frac{T_{i+1}-T_i}{||x_{i+1} - x_i||^2} \right] - \left[kA \frac{T_{i}-T_{i-1}}{||x_i-x_{i-1}||^2} \right] + S_u + S_i T_i = 0    \\
        & \Rightarrow \left[ \frac{kA_{right}}{||x_{i+1} - x_i||^2} \right]T_{i+1} + \left[-\frac{kA_{left}}{||x_{i} - x_{i-1}||^2} -\frac{kA_{right}}{||x_{i+1} - x_i||^2} + S_i \right]T_{i} + \left[ \frac{kA_{left}}{||x_{i} - x_{i-1}||^2} \right]T_{i-1} = -S_u   \\
        & \Rightarrow \sum_{j=1}^n\left[ \frac{k_{i \leftrightarrow j}A_{i \leftrightarrow j}}{||x_{i} - x_{j}||^2} \right]T_{j} -\sum_{j=1}^n\left[\left(\frac{k_{i \leftrightarrow j}A_{i \leftrightarrow j}}{||x_{i} - x_{j}||^2}\right) + S_i  \right]T_{i} = -S_u
    \end{align*}

Here considering total number of shared cells are :math:`n` for a given cell :math:`i` 

.. math::
    \begin{align*}
        & a_{ij} = \sum_{j=1}^n\left[ \frac{k_{i \leftrightarrow j}A_{i \leftrightarrow j}}{||x_{i} - x_{j}||^2} \right]    \\
        & a_{ii} = -\sum_{j=1}^n\left[\frac{k_{i \leftrightarrow j}A_{i \leftrightarrow j}}{||x_{i} - x_{j}||^2}\right] - S_i  \\
        & b_{i} = -S_u
    \end{align*}


Discretized Equation
---------------------

For a structured grid, the discrete form in 3D becomes:

.. math::

    a_P T_P = a_E T_E + a_W T_W + a_N T_N + a_S T_S + a_T T_T + a_B T_B + b

where:

- :math:`a_P` is the coefficient for the central node (:math:`P`).
- :math:`a_E, a_W, a_N, a_S, a_T, a_B` are the coefficients for the neighboring nodes (East, West, North, South, Top, Bottom).
- :math:`T_E, T_W, T_N, T_S, T_T, T_B` are the temperatures at neighboring nodes.
- :math:`b` is the linearized source term.

Coefficients
------------

The coefficients are defined as:

.. math::

    a_P = \sum_{nb} a_{nb} - S_P

.. math::

    a_{nb} = \frac{k A}{\Delta x_{nb}}

.. math::

    b = S_U \Delta V

where:

- :math:`k` is the thermal conductivity.
- :math:`A` is the face area of the control volume.
- :math:`\Delta x_{nb}` is the distance between node P and its neighbor.
- :math:`S_P` and :math:`S_U` are the linearized source term coefficients.
- :math:`\Delta V` is the control volume size.

Linearized Source Term
----------------------

The source term :math:`S_T` can be linearized as:

.. math::

    S_T = S_U + S_P T_P

where:

- :math:`S_U` represents the constant part of the source term.
- :math:`S_P` represents the coefficient of the temperature at the central node.

Substituting the linearized source term into the discretized equation modifies the central coefficient :math:`a_P` and the constant term :math:`b` as shown above.

Summary
-------

The finite volume discretization of the heat diffusion equation provides a robust framework for solving heat transfer problems in three dimensions, with linearized source terms ensuring computational efficiency and stability.
