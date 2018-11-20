Map compatibility
#################

A |tpetra_map|_ instance can be thought of abstractly as representing a vector
space. If two vectors have the `same` map, it's as if they come from the same
vector space and is, therefore, both legal and meaningful to add them together.
If, instead, they come from different vector spaces, then more information is
needed to know whether the two vectors are `compatible` and if it is legal to
add the vectors together.

Determining Map Compatibility
=============================

Two ``Map``\s are compatible when:

* they have the same global number of entries; and
* MPI processes in the Map's communicator that have the same MPI rank, own the
  same number of entries.

Two ``Map``\s are the same when:

* their minimum and maximum global indices are the same;
* they have the same global number of entries;
* the ``Map``\s are both distributed over multiple processes, or both not
  distributed over multiple processes;
* the ``Map``\s have the same index base (this means the smallest legal global
  index value, more or less);
* processes that have the same rank, own the same number of entries; and
* processes that have the same rank, own the same entries. That is, their
  entries have the same indices, in the same order.

Comparing Two Maps
------------------

|tpetra_map|_ member functions ``isCompatible`` and ``isSameAs`` provide methods
for determining compatibility and sameness of two maps, respectively.  For
example, ``Map``\s ``map1`` and ``map2`` are compatible if
``map1.isCompatible(map2)`` or the same if ``map1.isSameAs(map2)``.  Both
sameness and compatibility are commutative Boolean relations:
``map1.isCompatible(map2)`` means ``map2.isCompatible(map1)``.

Compatible Maps
===============

Compatibility of two ``Map``\s corresponds to isomorphism of two vector spaces.
Two ``Map``\s that are the same are always compatible. The ``isCompatible()``
criterion is less restrictive, and also less expensive to check (although
checking for compatibility requires a reduction on a Boolean over all processes
in the ``Map``'s communicator).

Adding together two vectors with compatible but not the same ``Map``\s is legal. It
might not make mathematical sense, depending on your application. This is
because entries of the vectors are ordered differently. (Also, just because two
vector spaces are isomorphic, doesn't necessarily mean that adding entries of
one to entries of another makes sense.) Adding together two vectors with the
same ``Map``\s is both legal and mathematically sensible.

.. include:: /links.txt
