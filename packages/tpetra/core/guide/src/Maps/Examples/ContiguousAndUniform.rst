.. _maps_example_contiguous_and_uniform:

Example: Contiguous and Uniform Map
###################################

.. rubric:: Keywords

Contiguous map, Uniform map, ``Map``

Overview
========

In this example, a contiguous and uniform ``Map`` is created by scaling the global number of entries in the ``Map`` with the number of MPI processes.  That way, this example can be run with any number of MPI processes and every process will still have a positive number of entries.

The Contiguous and Uniform Map Program
======================================

The following source code listing demonstrates the creation of the contiguous
and uniform map from the preceding example.

.. only:: builder_html

   The source code can be downloaded from :download:`here </Examples/SourceCode/map_contiguous_and_uniform.cpp>`.

.. literalinclude:: /Examples/SourceCode/map_contiguous_and_uniform.cpp
   :language: c++
   :linenos:
   :lines: 41-

Notes
-----

* ``numLocalEntries`` is the local (on the calling MPI process) number of entries (indices) in the ``Map``.  Tpetra expects a ``size_t`` for this value.

* ``numGlobalEntries`` is the total (global, i.e., over all MPI processes) number of entries (indices) in the ``Map``.  Tpetra expects ``Tpetra::global_size_t`` for this value.  This type is at least 64 bits long on 64-bit machines.

* The ``Map`` constructor puts the same number of equations on each processor.  The resulting ``Map`` is "contiguous and uniform."

.. _indexbase:

* The ``indexBase`` argument to ``Map``\s constructor tells Tpetra the starting index of the entries of a ``Map``.  Common values of ``indexBase`` are 0 (C style), 1 (Fortran style), though any base can be used.  1-based indexing is handy when interfacing with Fortran.

* All ``Map`` constructors must be called as a collective over the input communicator.  Not all ``Map`` constructors necessarily require communication, but some do, so it's best to treat them all as collectives.

* ``Map``\s should be considered immutable objects.  This is why it is created as a "``const map_type``".  If a new data distribution is needed, create a new ``Map`` should be created for it.

* ``contigMap`` is contiguous by construction.  This assertion can be tested at run time with the ``TEUCHOS_TEST_FOR_EXCEPTION`` macro, which throws an exception of the given type (second argument) with the given message (third argument), if the first argument is true.
