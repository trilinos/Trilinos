.. _maps_example_cyclic:

Example: Cyclic Map
###################

.. rubric:: Keywords

``Map``, Cyclic map

Overview
========

In this example, a ``Map`` is created which has the same number of global
entries per process and distributes them in a round-robin (1-D cyclic) fashion
instead of contiguously.

The Cyclic Map Program
======================

The following source code listing demonstrates the creation of the cyclic
map.

.. only:: builder_html

   The source code can be downloaded from :download:`here </Examples/SourceCode/map_cyclic.cpp>`.

.. literalinclude:: /Examples/SourceCode/map_cyclic.cpp
   :language: c++
   :linenos:
   :lines: 41-110

Notes
-----

* A version of the ``Map`` constructor is used that takes, on each MPI process, a list of the global entries in the ``Map`` belonging to that process (lines 28-38).  This constructor can be used to construct an overlapping (also called "not 1-to-1") ``Map``, in which one or more entries are owned by multiple processes.  This is not done here; in this example, a nonoverlapping (also called "1-to-1") ``Map`` is constructed.

* If there is more than one MPI process in the communicator, then cyclicMap is definitely NOT contiguous.  This condition is tested on lines 42-45.

* Since ``contigMap`` from :ref:`maps_example_contiguous_and_uniform` and ``cyclicMap`` have the same number of entries on all MPI processes in their communicators, if they have the same communicator they would be compatible.  However, if the communicator contains more than 1 process, then ``contigMap`` and ``cyclicMap`` are NOT the same.  The construction of ``contigMap`` is repeated on lines 55-56 and the assertion tested on lines 54-60.
