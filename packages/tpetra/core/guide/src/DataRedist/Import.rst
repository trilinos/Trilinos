.. _import:

Import
######

Overview
========

``Import`` sets up a communication plan for data redistribution from a uniquely-owned source distribution to a (possibly) multiply-owned target distribution.

Target Map Categorization
=========================

An ``Import`` object categorizes the elements of the target map into three sets
as follows:

#. All elements in the target map that have the same global ID as the corresponding element of the source map, starting with the first element in the target map, going up to the first element that is different from the source map.  The number of these IDs is returned by ``getNumSameIDs()``.

#. All elements that are local to the processor, but are not part of the first set of elements. These elements have global IDs that are owned by the calling processor, but at least the first element of this list is permuted. Even if subsequent elements are not permuted, they are included in this list.  The list of elements (local IDs) in the source map that are permuted can be found in the list ``getPermuteFromLIDs()``.  The list of elements (local IDs) in the target map that are the new locations of the source elements can be found in the list ``getPermuteToLIDs()``.

#. All remaining elements of the target map correspond to global IDs that are owned by remote processors.  The number of these elements is returned by ``getNumRemoteLIDs()`` and the list of these is returned by ``getRemoteLIDs()``.

Given the above information, the ``Import`` constructor builds a list of elements that must be communicated to other processors as a result of import requests. The number of exported elements (where multiple sends of the same element to different processors is counted) is returned by ``getNumExportIDs()``. The local IDs to be sent are returned by the list ``getExportLIDs()``. The processors to which each of the elements will be sent is returned in a list of the same length by ``getExportPIDs()``.

The total number of elements that will be sent by the calling processor is returned by ``numSend()``. The total number of elements that will be received is returned by ``numRecv()``.

Illustrative Example
--------------------

Assume we have 3 processors and 9 global elements with each processor owning 3 elements as follows

+-----------+---------+---------+---------+
| Processor | 0       |  1      |  2      |
+-----------+---------+---------+---------+
| Elements  | 0, 1, 2 | 3, 4, 5 | 6, 7, 8 |
+-----------+---------+---------+---------+

The above layout essentially defines the source map argument of the import object.

This could correspond to a 9 by 9 matrix with the first three rows on PE 0, and so on. Suppose that this matrix is periodic tridiagonal having the following sparsity pattern:

.. math::

   \begin{matrix}
     \\[-.2in] \text{PE 0} \\ \\
     \\ \text{PE 1} \\ \\
     \\ \text{PE 2} \\
   \end{matrix}
   \begin{bmatrix}
     X &  X &  0 &  0 &  0 &  0 &  0 &  0 &  X \\
     X &  X &  X &  0 &  0 &  0 &  0 &  0 &  0 \\
     0 &  X &  X &  X &  0 &  0 &  0 &  0 &  0 \\
     \hline
     0 &  0 &  X &  X &  X &  0 &  0 &  0 &  0 \\
     0 &  0 &  0 &  X &  X &  X &  0 &  0 &  0 \\
     0 &  0 &  0 &  0 &  X &  X &  X &  0 &  0 \\
     \hline
     0 &  0 &  0 &  0 &  0 &  X &  X &  X &  0 \\
     0 &  0 &  0 &  0 &  0 &  0 &  X &  X &  X \\
     X &  0 &  0 &  0 &  0 &  0 &  0 &  X &  X
   \end{bmatrix}

To perform a matrix vector multiplication operation :math:`y = Ax` (assuming that :math:`x` has the same distribution as the rows of the matrix :math:`A`) each processor will need to import elements of :math:`x` that are not local. To do this, we build a target map on each processor as follows:

+-----------+---------------+---------------+---------------+
| Processor | 0             |  1            |  2            |
+-----------+---------------+---------------+---------------+
| Elements  | 0, 1, 2, 3, 8 | 2, 3, 4, 5, 6 | 0, 5, 6, 7, 8 |
+-----------+---------------+---------------+---------------+

The above list is the elements that will be needed to perform the matrix vector multiplication locally on each processor. Note that the ordering of the elements on each processor is not unique, but has been chosen for illustration.

With these two maps passed into the ``Import`` constructor, we get the following attribute definitions:

+------------------------+-----------------+-----------------+-----------------+
| Processor              | 0               | 1               | 2               |
+========================+=================+=================+=================+
| Number of same IDs     | 3               | 0               | 0               |
+------------------------+-----------------+-----------------+-----------------+
| Number of permuted IDs | 0               | 3               | 3               |
+------------------------+-----------------+-----------------+-----------------+
| Permute to local IDs   | 0               | :math:`[0,1,2]` | :math:`[0,1,2]` |
+------------------------+-----------------+-----------------+-----------------+
| Permute from local IDs | 0               | :math:`[1,2,3]` | :math:`[2,3,4]` |
+------------------------+-----------------+-----------------+-----------------+
| Number of remote IDs   | 2               | 2               | 2               |
+------------------------+-----------------+-----------------+-----------------+
| Remote local IDs       | :math:`[3, 4]`  | :math:`[0,4]`   | :math:`[0,1]`   |
+------------------------+-----------------+-----------------+-----------------+
| Number of export IDs   | 2               | 2               | 2               |
+------------------------+-----------------+-----------------+-----------------+
| Export local IDs       | :math:`[0, 2]`  | :math:`[0,2]`   | :math:`[0,2]`   |
+------------------------+-----------------+-----------------+-----------------+
| Export processor IDs   | :math:`[1, 2]`  | :math:`[0,2]`   | :math:`[0,1]`   |
+------------------------+-----------------+-----------------+-----------------+
| Number of sends        | 2               | 2               | 2               |
+------------------------+-----------------+-----------------+-----------------+
| Number of receives     | 2               | 2               | 2               |
+------------------------+-----------------+-----------------+-----------------+

Using Import Objects
====================

Once a ``Import`` object has been constructed, it can be used by any of the Tpetra classes that support distributed global objects, namely ``Tpetra::Vector``, ``Tpetra::MultiVector``, ``Tpetra::CrsGraph``, and ``Tpetra::CrsMatrix``. All of these classes have ``Import`` and ``Export`` methods that will fill new objects whose distribution is described by the target map, taking elements from the source object whose distribution is described by the source map.

One use case of ``Import`` is bringing in remote source vector data for a distributed sparse matrix-vector multiply. The source vector itself is uniquely owned, but must be brought in into an overlapping distribution so that each process can compute its part of the target vector without further communication.

Note that the reverse operation, an export, using this importer is also possible and appropriate in some instances. For example, if we compute :math:`y = A^Tx`, the transpose matrix-multiplication operation, then we can use the importer we constructed in the above example to do an export operation to y, adding the contributions that come from multiple processors.
