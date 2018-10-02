.. _matrix_fill:

Example: Matrix Fill and Complete
#################################

.. rubric:: Keywords

``CrsMatrix``, ``insertGlobalValues``, ``fillComplete``

Overview
========

In this example, a Tpetra sparse matrix object is created, filled, and
completed for a matrix of the form

.. math::

   \begin{bmatrix}
   2 & -1 \\
   -1 & 2 & -1 \\
   & -1 & 2 & -1 \\
   &&& \ddots \\
   &&& -1 & 2 & -1 \\
   &&&&& -1 & 2
   \end{bmatrix}

Full Source Code Listing
========================

The following source code listing demonstrates the creation of the contiguous
and uniform map from the preceding example.

.. only:: builder_html

   The source code can be downloaded from :download:`here </Examples/SourceCode/matrix_fill_1.cpp>`.

.. literalinclude:: /Examples/SourceCode/matrix_fill_1.cpp
   :language: c++
   :linenos:
   :lines: 41-

