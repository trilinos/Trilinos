.. _guide_overview:

Overview of this guide
######################

This guide is a reference to using Tpetra.  This guide should not be viewed as
exhaustive.  Indeed, Tpetra is a large library and this guide covers only a very
small portion of capabilities needed by most users.  The guide is divided into
two main parts:

.. raw:: html

   <h3> Part I: Overview </h3>

* :ref:`intro_and_overview`
* :ref:`obtaining_and_building`
* :ref:`initializing_tpetra`

.. raw:: html

   <h3> Part II: Key User Facing Classes
   </h3>

These are the classes that user's of Tpetra will want to become familiar with:

* :ref:`Parallel distributions <maps_index>`: ``Tpetra::Map``.  Maps contain information used to
  distribute vectors, matrices and other objects. This class is analogous to
  Epetra's ``Epetra_Map`` class.

* :ref:`Distributed dense vectors <vectors>`: ``Tpetra::MultiVector``, ``Tpetra::Vector``.
  Vectors rovides vector services such as scaling, norms, and dot products.

* :ref:`Distributed sparse matrices <sparse_matrices>`: ``Tpetra::RowMatrix``, ``Tpetra::CrsMatrix``.
  ``Tpetra::RowMatrix`` is a abstract interface for row-distributed sparse
  matrices. ``Tpetra::CrsMatrix`` is a specific implementation of
  ``Tpetra::RowMatrix``, utilizing compressed row storage format. Both of
  these classes derive from ``Tpetra::Operator``, the base class for linear
  operators.

* :ref:`Import/Export classes <data_redist>`: ``Tpetra::Import`` and ``Tpetra::Export``.  These
  classes allow efficient transfer of objects built using one mapping to a new
  object with a new mapping. Supports local and global permutations, overlapping
  Schwarz operations and many other data movement operations.

Examples
========

Most of the topics covered include commented (and functioning) examples.  See
:ref:`building_examples` for instructions on how to build and run the examples.

.. include:: /links.txt
