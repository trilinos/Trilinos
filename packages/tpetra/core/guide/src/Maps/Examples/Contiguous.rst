.. _maps_example_contiguous:

Example: Contiguous Map, Not Necessarily Uniform
################################################

.. rubric:: Keywords

``Map``, Contiguous map

Overview
========

In this example, a contiguous ``Map`` is created by scaling the global number of
entries in the ``Map`` with the number of MPI processes, like in
:ref:`maps_example_contiguous_and_uniform`.  However, in this example,
a different ``Map`` constructor is used that takes the number of entries on each
MPI process.  The resulting ``Map`` is "contiguous" but not necessarily uniform,
since the numbers of entries on different MPI processes may differ.

The Contiguous Map Program
==========================

The following source code listing demonstrates the creation of the contiguous
and uniform map from the preceding example.

.. only:: builder_html

   The source code can be downloaded from :download:`here </Examples/SourceCode/map_contiguous.cpp>`.

.. literalinclude:: /Examples/SourceCode/map_contiguous.cpp
   :language: c++
   :linenos:
   :lines: 41-

Notes
-----

* Since ``numLocalEntries`` was the same on all processors, this map is also uniform.  However, ``numLocalEntries`` could have been different on each process (the sum of ``numLocalEntries`` over all processes must still equal ``numGlobalEntries``) in which case the ``Map`` would not have been uniform.

* Since ``contigMap`` from :ref:`maps_example_contiguous_and_uniform` and ``contigMap2`` have the same number of entries on all MPI processes in their communicators, if they had the same communicators they would be "the same."  The construction of ``contigMap`` is repeated in lines 27-28 and the assertion tested.
