.. _sparse_mat_redist:

Example: Sparse Matrix Redistribution
#####################################

.. rubric:: Keywords

Import, Export, data redistribution

Overview
========

This example demonstrates how to use ``Map``\s and Tpetra's ``Export`` class to redistribute data. In this case, a sparse matrix is built on one MPI process, and is redistributed to a sparse matrix stored in block row fashion, with an equal number of rows per process.

The matrix in this example is tridiagonal of the form:

.. math::

   \begin{bmatrix}
   2 & -1 \\
   & -1 & 2 & -1 \\
   &&& \ddots \\
   &&& -1 & 2 & -1 \\
   &&&& -1 & 2 \\
   \end{bmatrix}

The global number of rows in the matrix is scaled relative to the number of (MPI) processes, so that no matter how many MPI processes are used to run the example, every process will have 10 rows.

Two maps are created in the example.  The first places all rows of the matrix on
process 0.  The second distributes the rows of the matrix evenly amongst
available processes.  The matrix is created and filled on the first map and then
exported to the second.

The Sparse Matrix Redistribution Program
========================================

The following source code listing demonstrates redistribution of a sparse matrix
from one process to all processes.

.. only:: builder_html

   The source code can be downloaded from :download:`here </Examples/SourceCode/data_redist_1.cpp>`.

.. literalinclude:: /Examples/SourceCode/data_redist_1.cpp
   :language: c++
   :linenos:
   :lines: 41-
   :emphasize-lines: 112-137

Notes
-----

* Teuchos_ timers are used throughout the program to time different code
  segments.  Starting and stopping timer objects is local, but computing timer statistics (e.g., with ``TimeMonitor::summarize()``) is global.  There are ways to restrict the latter to any given MPI communicator; the default is ``MPI_COMM_WORLD``.

* Since both the source and target ``Map``\s are one-to-one, an ``Import`` or an ``Export`` could be used.  If only the source ``Map`` were one-to-one, an ``Import`` would have to be used; if only the target ``Map`` were one-to-one, then an ``Export`` would have to be used.  Redistribution using ``Import`` or ``Export`` is not allowed if neither the source nor target ``Map`` is one-to-one.

* The ``Export`` type has the same template parameters as a ``Map``.  Note that ``Export`` does not depend on the ``Scalar`` template parameter of the objects it redistributes.  You can reuse the same ``Export`` for different Tpetra object types, or for Tpetra objects of the same type but different Scalar template parameters (e.g., ``Scalar=float`` or ``Scalar=double``).

.. include:: /links.txt
