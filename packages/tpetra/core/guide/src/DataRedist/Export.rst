.. _export:

Export
######

Overview
========

``Export`` sets up a communication plan for data redistribution from a (possibly) multiply-owned source distribution to a uniquely-owned target distribution.

Target Map Categorization
=========================

An ``Export`` object categorizes the elements of the target map into three sets
as follows:

#. All elements in the target map that have the same global ID as the corresponding element of the source map, starting with the first element in the target map, going up to the first element that is different from the source map.  The number of these IDs is returned by ``getNumSameIDs()``.

#. All elements that are local to the processor, but are not part of the first set of elements. These elements have global IDs that are owned by the calling processor, but at least the first element of this list is permuted. Even if subsequent elements are not permuted, they are included in this list.  The number of elements is returned by ``getNumPermuteIDs()``.  The list of elements (local IDs) in the source map that are permuted can be found in the list ``getPermuteFromLIDs()``.  The list of elements (local IDs) in the target map that are the new locations of the source elements can be found in the list ``getPermuteToLIDs()``.

#. All remaining elements of the target map correspond to global IDs that are owned by remote processors.  The number of these elements is returned by ``getNumRemoteLIDs()`` and the list of these is returned by ``getRemoteLIDs()``.

Given the above information, the ``Export`` constructor builds a list of elements that must be communicated to other processors as a result of export requests. The number of exported elements (where multiple sends of the same element to different processors is counted) is returned by ``getNumExportIDs()``. The local IDs to be sent are returned by the list ``getExportLIDs()``. The processors to which each of the elements will be sent is returned in a list of the same length by ``getExportPIDs()``.

The total number of elements that will be sent by the calling processor is returned by ``numSend()``. The total number of elements that will be received is returned by ``numRecv()``.

Illustrative Example
--------------------

Assume we have 3 processors and 9 global elements with each processor owning 3 elements as follows

+-----------+---------+---------+---------+
| Processor | 0       |  1      |  2      |
+-----------+---------+---------+---------+
| Elements  | 0, 1, 2 | 3, 4, 5 | 6, 7, 8 |
+-----------+---------+---------+---------+

The above layout essentially defines the target map argument of the export object.

This could correspond to a 9 entry forcing vector with the first three entries
on PE 0, and so on.  Suppose that the entries of this forcing vector are
computed by integrating over linear first-order finite element basis functions
("hat functions"):

.. image:: /Images/export_hatfun.pdf

In this case, PE 0 will make contributions to entries 0 through 3, PE 1 will make contributions to entries 2 through 6 and PE 2 will make contributions to entries 5 through 8. A convenient way to compute these contributions is to create a forcing vector with replicated entries for the shared contributions. Specifically the following source ``Map`` works for this scenario:

+-----------+---------------+---------------+---------------+
| Processor | 0             |  1            |  2            |
+-----------+---------------+---------------+---------------+
| Elements  | 0, 1, 2, 3    | 2, 3, 4, 5, 6 |    5, 6, 7, 8 |
+-----------+---------------+---------------+---------------+

A vector constructed using this source ``Map`` can be used to collect each processor's contributions to the forcing vector. Note that the ordering of the elements on each processor is not unique, but has been chosen for illustration.

With these two maps passed into the ``Export`` constructor, we get the following attribute definitions:

+------------------------+-----------------+-----------------+-----------------+
| Processor              | 0               | 1               | 2               |
+========================+=================+=================+=================+
| Number of same IDs     | 3               | 0               | 0               |
+------------------------+-----------------+-----------------+-----------------+
| Number of permuted IDs | 0               | 3               | 3               |
+------------------------+-----------------+-----------------+-----------------+
| Permute to local IDs   | 0               | :math:`[0,1,2]` | :math:`[0,1,2]` |
+------------------------+-----------------+-----------------+-----------------+
| Permute from local IDs | 0               | :math:`[1,2,3]` | :math:`[1,2,3]` |
+------------------------+-----------------+-----------------+-----------------+
| Number of remote IDs   | 1               | 2               | 1               |
+------------------------+-----------------+-----------------+-----------------+
| Remote local IDs       | :math:`[2]`     | :math:`[0,2]`   | :math:`[0]`     |
+------------------------+-----------------+-----------------+-----------------+
| Number of export IDs   | 1               | 2               | 1               |
+------------------------+-----------------+-----------------+-----------------+
| Export local IDs       | :math:`[3]`     | :math:`[0,4]`   | :math:`[0]`     |
+------------------------+-----------------+-----------------+-----------------+
| Export processor IDs   | :math:`[1]`     | :math:`[0,2]`   | :math:`[1]`     |
+------------------------+-----------------+-----------------+-----------------+
| Number of sends        | 1               | 2               | 1               |
+------------------------+-----------------+-----------------+-----------------+
| Number of receives     | 1               | 2               | 1               |
+------------------------+-----------------+-----------------+-----------------+

Using Export Objects
====================

Once a ``Export`` object has been constructed, it can be used by any of the Tpetra classes that support distributed global objects, namely ``Tpetra::Vector``, ``Tpetra::MultiVector``, ``Tpetra::CrsGraph``, and ``Tpetra::CrsMatrix``. All of these classes have ``Import`` and ``Export`` methods that will fill new objects whose distribution is described by the target map, taking elements from the source object whose distribution is described by the source map.

In the above example, if ``x_integrate`` is constructed using the source ``Map`` and then filled with local contributions, and ``x_force`` is constructed using the target map, the following operation will fill ``x_force`` with the combined results of ``x_integrate``, using the default template parameters for the ``Export`` type (the ``Export`` type has the same template parameters as a ``Map``):

.. code:: c++

   typedef Tpetra::Export<> export_type;
   export_type exporter(source_map, target_map);
   x_force.doExport(x_integrate, exporter, Tpetra::SUM);

The third argument above tells the export operation to sum results that come from multiple processors for the same global ID.

``Export`` objects can also be used by ``Import`` operations to perform the reverse operation. For example, if ``x_force`` in the above example had boundary conditions that should be sent to processors that share a boundary element, the following operation would send replicated values to ``x_integrate``:

.. code:: c++

   typedef Tpetra::Export<> export_type;
   export_type exporter(source_map, target_map);
   x_integrate.doImport(x_force, exporter, Tpetra::INSERT);

At the end of this operation, ``x_integrate`` would have replicated values from ``x_force`` of entries 2 and 3 on PEs 0 and 1, and entries 5 and 6 on PEs 1 and 2.
