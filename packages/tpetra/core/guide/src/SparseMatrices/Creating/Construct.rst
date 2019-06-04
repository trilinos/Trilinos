.. _crsmatrix_construct:

Constructing a CrsMatrix
########################

A |tpetra_crsmatrix|_ is distributed over one or more parallel processes,
just like  ``Vector`` or other distributed objects. Also like a
``Vectors``, the ``CrsMatrix`` must be given its distribution on construction.
Unlike ``Vector``\s, though, a ``CrsMatrix`` has two dimensions over which to
be distributed: rows and columns.

Template Parameters
===================

The ``Tpetra::CrsMatrix`` takes the following template parameters:

* ``Scalar``
* ``LocalOrdinal``
* ``GlobalOrdinal``
* ``Node``

"One-Dimensional" Construction
==============================

For many users, distributing the matrix in "one-dimensional" fashion over rows
and ignoring the column distribution is sufficient. In that case, only a "row
``Map``" needs to be supplied to the constructor. This implies that for any row
which a process owns, that process may insert entries in any column of that row.

"Two-Dimensional" Construction
==============================

Other users may need the full flexibility of distributing both the rows and
columns of the matrix over processes. This "two-dimensional" distribution, if
chosen optimally, can significantly reduce the amount of communication needed
for distributed-memory parallel sparse matrix-vector multiply. Trilinos packages
like Zoltan_ and Zoltan2_ can help you compute this distribution. In that case,
both the "row ``Map``" and the "column ``Map``" may be given to the constructor.
This implies that for any row which a process owns, that process may insert
entries in any column in that row which that process owns in its column ``Map``.

Constructing a CrsMatrix from a CrsGraph
========================================

Still other users may already know the structure of the sparse matrix, and just
need to fill in values. These users should first create the graph (a
|tpetra_crsgraph|_) and then give the graph to the constructor of ``CrsMatrix``.
The graph may have either a "1-D" or "2-D" distribution, as mentioned above.

.. note::

   In future versions of Tpetra, the last method listed (constructing
   a ``CrsMatrix`` from an existing ``CrsGraph``) may become the preferred, and
   only, method of creating a ``CrsMatrix``.

.. |tpetra_crsmatrix| replace:: :code:`Tpetra::CrsMatrix`
.. _tpetra_crsmatrix: https://trilinos.org/docs/dev/packages/tpetra/doc/html/classTpetra_1_1CrsMatrix.html
.. |tpetra_crsgraph| replace:: :code:`Tpetra::CrsGraph`
.. _tpetra_crsgraph: https://trilinos.org/docs/dev/packages/tpetra/doc/html/classTpetra_1_1CrsGraph.html

.. _Zoltan: https://www.trilinos.org/packages/zoltan
.. _Zoltan2: https://www.trilinos.org/packages/zoltan2
