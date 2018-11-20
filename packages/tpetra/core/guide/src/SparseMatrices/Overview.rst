.. _sparse_matrices_overview:

Overview
########

The |tpetra_crsmatrix|_ class implements a distributed-memory parallel sparse
matrix, and provides sparse matrix-vector multiply (including transpose) and
sparse triangular solve operations. It provides access by rows to the elements
of the matrix, as if the local data were stored in compressed sparse row format.
(Implementations are not required to store the data in this way internally.)
This class has an interface like that of ``Epetra_CrsMatrix``, but also allows
insertion of data into nonowned rows, much like ``Epetra_FECrsMatrix``.

Rowwise Access, Not Rowwise Distribution
========================================

``CrsMatrix`` inherits from the ``RowMatrix`` interface; it "is a"
``RowMatrix``. Some people have the misconception that ``RowMatrix`` (and
therefore ``CrsMatrix``) refers to a rowwise distribution of the sparse matrix
over parallel processes. It does not. ``CrsMatrix`` supports an arbitrary
"two-dimensional" distribution of rows and columns over processes. The term
instead refers to rowwise access to the data. That is, the methods in this
interface and in ``CrsMatrix`` let you add or access entries on each process by
row.

This distinction matters because two or more processes might share entries in
a row. Asking for "all the entries in a row" on a particular process only
accesses the entries owned by that process, which is not necessarily all the
entries in a row.

Whether adding entries or modifying existing ones, you may always do so for any
number of entries in the row, such that their columns are owned by the calling
process. You should always modify as many entries with one method call as
possible, in order to amortize function call and data access overhead.

.. include:: /links.txt
