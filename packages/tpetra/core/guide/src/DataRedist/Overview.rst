Overview
########

Tpetra's ``Map`` class describes a data distribution over one or more distributed-memory parallel processes. It "maps" global indices (unique labels for the elements of a data structure) to parallel processes. This ability to describe a data distribution calls for a redistribution capability, that is, to reorganize or remap data from one distribution to another. Tpetra provides this capability through the ``Import`` and ``Export`` classes.

``Import`` redistributes from a uniquely owned (one-to-one) ``Map`` to a possibly not uniquely owned ``Map``. ``Export`` redistributes from a possibly not uniquely owned to a uniquely owned ``Map``. We distinguish between these cases both for historical reasons and for performance reasons.

``Import`` and ``Export`` objects encapsulate and remember a communication pattern for reuse. Computing the computation pattern requires nontrivial work, but keeping around the ``Import`` or ``Export`` object lets you reuse that work. This is very important for operations that are performed frequently, such as the ``Import`` and ``Export`` operations in Tpetra's sparse matrix-vector multiply.

In both cases, ``Import`` and ``Export`` let the user specify how to combine incoming new data with existing data that has the same global index. For example, one may replace old data with new data or sum them together.

.. note::

   The names "``Import``" and "``Export``" have nothing to do with the direction in which data moves relative to the calling process; any process may do both receives and sends in an ``Import`` or ``Export``. Rather, the names suggest what happens in their most common use case, the communication pattern for sparse matrix-vector multiply. ``Import`` "brings in" remote source vector data (from the domain ``Map`` to the column ``Map``) for local computation, and ``Export`` "pushes" the result back (from the row ``Map`` to the range ``Map``). ``Import`` and ``Export`` have other uses as well.

.. note::

   Epetra separated Import and Export for performance reasons. The implementation is different, depending on which direction is the uniquely-owned Map. Tpetra retains this convention.

