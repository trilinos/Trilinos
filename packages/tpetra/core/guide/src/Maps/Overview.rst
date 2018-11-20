.. _data_dist:

Maps: Overview
##############

Tpetra, like Epetra, uses "``Map``" objects to encapsulate the details of the
distribution of data over MPI processes.  ``Map``\s make data distribution into
a first-class citizen. Each |tpetra_map|_ instance represents a particular
data distribution.  A ``Map``:

* has a ``Comm``\unicator;
* is like a vector space; and
* assigns entries of a data structure to (MPI) processes.

Data Encapsulation and Distribution
===================================

A ``Map`` assigns global and local indices to parallel processes.  An index is
global if it represents an entry of a distributed object uniquely over the
entire object and an index is local if it represents an entry of a distributed
object uniquely on the process that owns it.  Nearly any data structure
containing entries that can be assigned an integer index can be distributed
using a ``Map``.  For most Tpetra users, this means entries of a vector, rows of
a ``Tpetra::MultiVector``, or rows or columns of a sparse matrix. However, it is
not limited to these kinds of objects. ``Map`` can even be used to distribute
custom user objects.

Global indices
--------------

Global indices represent the entries of a distributed object (like columns of
a sparse matrix, or entries of a vector) uniquely over the entire object.  If
a ``Map`` assigns a global index ``G`` to a process ``P``, we say that process
``P`` owns global index ``G``. It is legal for multiple processes to own the
same global index ``G``. In fact, this is how many useful communication
patterns, including those in sparse matrix-vector multiply, are implemented.

Local indices
-------------

Local indices represent entries of a distributed object uniquely on the that
owns them. For efficiency, within a process, a global index is referred to by
using its "local index" on that process. Local indices are local to the process
that owns them. If process ``P`` owns global index ``G``, then there is a unique
local index ``L`` on process ``P`` corresponding to ``G``. If the local index
``L`` is valid on process ``P``, then there is a unique global index ``G`` owned
by ``P`` corresponding to the pair ``(L, P)``. However, multiple processes might
own the same global index, so a global index ``G`` might correspond to multiple
``(L, P)`` pairs.

Local indices matter to users of Tpetra because it may be more efficient to use
them to access or modify local data than it is to use global indices. This is
because distributed data structures must convert from global to local indices
every time a user asks for an element by its global index. This requires a table
lookup in general, since a process may own an arbitrary subset of all the global
indices, in an arbitrary order. Even though local indices are an implementation
detail, they are exposed because avoiding that table lookup on each access can
improve performance a lot.

Maps are themselves distributed data
====================================

If a ``Map`` has ``N`` global entries over ``P`` processes, and if no one process
owns all the global entries, we never store all ``N`` global indices on a single
process. Some kinds of ``Map``\s require storing all the global indices, but in this
case, the indices are themselves distributed over processes. This ensures memory
scalability (no one process has to store all the data).

.. include:: /links.txt
