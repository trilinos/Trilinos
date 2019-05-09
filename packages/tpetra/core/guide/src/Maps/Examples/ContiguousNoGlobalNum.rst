.. _maps_example_contiguous_no_global_num:

Example: Contiguous Map, Unknown Number of Global Entries
#########################################################

.. rubric:: Keywords

``Map``, Contiguous map

Overview
========

In this example, a contiguous ``Map`` is created by scaling the global number of
entries in the ``Map`` with the number of MPI processes, like in
:ref:`maps_example_contiguous_and_uniform`.  However, in this example,
the number of global entries is not specified.  This is helpful for the case
that only the number of entries on each process is known, but not the global
number.

The Contiguous Map Program
==========================

The following source code listing demonstrates the creation of the contiguous
map with (initially) unknown number of global entries.

.. only:: builder_html

   The source code can be downloaded from :download:`here
   </Examples/SourceCode/map_contiguous_no_global_num.cpp>`.

.. literalinclude:: /Examples/SourceCode/map_contiguous_no_global_num.cpp
   :language: c++
   :linenos:
   :lines: 41-

Notes
-----

* Instead of ``numGlobalEntries`` for ``Tpetra::global_size_t``, an "invalid value" flag is used (``INVALID``) that signals to Tpetra that the global number of entries is unknown and must be computed.

* This ``Map`` constructor is helpful for the case that only the number of entries on each process is known, but not the global number.

* Since ``contigMap`` from :ref:`maps_example_contiguous_and_uniform` and ``contigMap3`` have the same number of entries on all MPI processes in their communicators, if they had the same communicators they would be "the same."  The construction of ``contigMap`` is repeated in lines 28-36 and the assertion tested.
