==================
Beginners tutorial
==================

Quick start
===========

The first example is meant to quickly get into touch with MueLu.

Example problem
---------------

We generate a test matrix corresponding to the stencil of a 2D Laplacian operator on a structured Cartesian grid.
The matrix stencil is

.. math::
   \frac{1}{h^2}\begin{pmatrix} & -1 & \\ -1 & 4 & -1 \\ & -1 & \end{pmatrix}

where :math:`h` denotes the mesh size parameter.
The resulting matrix is symmetric positive definite.
We choose the right hand side to be the constant vector one
and use a random initial guess for the iterative solution process.
The problem domain is the unit square with a Cartesian (uniform) mesh.

Running this example
--------------------

The `s1_easy.xml` file has the following content:

.. literalinclude:: ../../src/packages/muelu/test/s1_easy.xml
  :language: xml

As one can easily find from the xml parameters,
a multigrid method with not more than 3 levels and a damped Jacobi method for level smoothing shall be used.

