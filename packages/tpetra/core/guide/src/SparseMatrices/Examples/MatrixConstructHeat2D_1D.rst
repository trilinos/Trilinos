:orphan:

..
  The code examples for this file are incomplete

.. _construct_ex_1:

Example: 2D Heat Transfer, One-Dimensional Constructor
######################################################

.. rubric:: Keywords

``CrsMatrix``, ``Map``, ``fillComplete``

Overview
========

The matrix associated with a 2D linear heat transfer operator is constucted using the ``Tpetra::CsrMatrix`` one-dimensional constructor.  The 2D linear heat transfer operator, discretized over a unit square by 8 linear three-node finite elements takes the form

.. math::

   \boldsymbol{A} = \begin{bmatrix}
     x & x &   & x \\
     x & x & x &   & x \\
       & x & x &   &   & x \\
     x &   &   & x & x &   & x \\
       & x &   & x & x & x &   & x \\
       &   & x &   & x & x &   &   & x \\
       &   &   & x &   &   & x & x &   \\
       &   &   &   & x &   & x & x & x \\
       &   &   &   &   & x &   & x & x \\
   \end{bmatrix}

where :math:`x` represents non-zero values and blank entries are :math:`0`.

In the following, suppose the operator :math:`\boldsymbol{A}` and the associated system of equations :math:`\boldsymbol{y}=\boldsymbol{A}\boldsymbol{x}` are distributed as follows amongst two processors:

.. math::

   \begin{Bmatrix}
     y_0 \\ y_1 \\ y_2 \\ y_3 \\ y_4 \\ \hline y_5 \\ y_6 \\ y_7 \\ y_8
   \end{Bmatrix}
   =
   \begin{bmatrix}
     x & x &   & x &   &   &   &   & P_0\\
     x & x & x &   & x \\
       & x & x &   &   & x \\
     x &   &   & x & x &   & x \\
       & x &   & x & x & x &   & x \\
       \hline
       &   & x &   & x & x &   &   & x \\
       &   &   & x &   &   & x & x &   \\
       &   &   &   & x &   & x & x & x \\
    P_1&   &   &   &   & x &   & x & x \\
   \end{bmatrix}
   \begin{Bmatrix}
     x_0 \\ x_1 \\ x_2 \\ x_3 \\ x_4 \\ \hline x_5 \\ x_6 \\ x_7 \\ x_8
   \end{Bmatrix}

For the distribution shown,  the row maps of :math:`\boldsymbol{A}` are

.. math::

   \begin{align}
   {\tt row\_map\_P0} &= \{0, 1, 2, 3, 4\} \\
   {\tt row\_map\_P1} &= \{5, 6, 7, 8\}
   \end{align}

The Sparse Matrix Construction Program
======================================

The following source code listing demonstrates the construction of the 2D heat
transfer operator derived.

.. only:: builder_html

   The source code can be downloaded from :download:`here </Examples/SourceCode/matrix_construct_heat2d_1.cpp>`.

.. literalinclude:: /Examples/SourceCode/matrix_construct_heat2d_1.cpp
   :language: c++
   :linenos:
   :lines: 41-

Notes
-----
