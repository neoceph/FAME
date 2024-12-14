Heat Diffusion
==============

This section focuses on the heat diffusion equation and its discretization in three dimensions, including the generic form of the governing equation and the linearized source term.

Generic Form of the Heat Diffusion Equation
-------------------------------------------

The heat diffusion equation in three dimensions can be written as:

.. math::

    \frac{\partial T}{\partial t} = \alpha \nabla^2 T + S_T

where:

- :math:`T` is the temperature (dependent variable).
- :math:`\alpha` is the thermal diffusivity (:math:`\alpha = \frac{k}{\rho c_p}`).
- :math:`\nabla^2 T` is the Laplacian of the temperature field.
- :math:`S_T` is a source term.

In steady-state conditions, the time-dependent term vanishes:

.. math::

    \nabla^2 T + \frac{S_T}{\alpha} = 0

Discretization in 3D
--------------------

Using the finite volume method (FVM), the diffusion equation is discretized over a control volume. The discretization involves integrating the equation over the control volume and applying Gauss's divergence theorem:

.. math::

    \int_V \nabla^2 T \ dV = \int_A \nabla T \cdot \mathbf{n} \ dA

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
