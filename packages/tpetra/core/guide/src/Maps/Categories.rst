Categories of Maps
##################

Overview
========

A |tpetra_map|_ can be classified as belonging to one or more of the following
categories:

* one-to-one (or not);
* contiguous (or not); and/or
* uniform (or not).

One-to-one
==========

A ``Map`` is one-to-one if each global index in the ``Map`` is owned by only one
process. This means that the function from global index ``G`` to its local index
and process rank ``(L,P)`` is one to one in a mathematical sense ("injective").
In this case, the function is only onto ("surjective") if there is only one
process. Knowing whether a ``Map`` is one-to-one is important for data
redistribution, which Tpetra exposes as the :ref:`Import` and :ref:`Export`
operations.

An example of a one-to-one ``Map`` is a ``Map`` containing :math:`11` global indices
:math:`0 \ldots 10` and distributed over four processes, where

* Process :math:`0` owns :math:`0 \ldots 2`
* Process :math:`1` owns :math:`3 \ldots 5`
* Process :math:`2` owns :math:`6 \ldots 8`
* Process :math:`3` owns :math:`9 \ldots 10`

An example of a `not` one-to-one ``Map`` is a ``Map`` containing :math:`11`
global indices :math:`0 \ldots 10` and distributed over four processes, where

* Process :math:`0` owns :math:`0 \ldots 2`
* Process :math:`1` owns :math:`2 \ldots 5`
* Process :math:`2` owns :math:`5 \ldots 8`
* Process :math:`3` owns :math:`8 \ldots 10`

Note the overlap of one global index between each "adjacent" process. An example
of a mathematical problem with an overlapping distribution like this would be
a 1-D linear finite element or finite difference discretization, where entries
are distributed with unique ownership among the processes, but the boundary node
between two adjacent entries on different processes is shared among those two
processes.

Map Contiguity
==============

A ``Map`` is contiguous when each process' list of global indices forms an
interval and is strictly increasing, and the globally minimum global index
equals the :ref:`index base <indexbase>`. ``Map`` optimizes for the contiguous
case. In particular, noncontiguous ``Map``\s require communication in order to
figure out which process owns a particular global index.

Note that in Tpetra, "contiguous" is an optimization, not a predicate. Tpetra
may not necessarily work hard to check contiguity. The best way to ensure that
your Map is contiguous is to use one of the two constructors that always make
a contiguous Map.

An example of a contiguous Map is one containing :math:`11` global indices
:math:`0 \ldots 10` and distributed over four processes, where

* Process :math:`0` owns :math:`0 \ldots 2`
* Process :math:`1` owns :math:`3 \ldots 5`
* Process :math:`2` owns :math:`6 \ldots 8`
* Process :math:`3` owns :math:`9 \ldots 10`

Map Uniformity
--------------

Note that Process 3 in the previous example owns 2 global indices, whereas the
other processes each own 3. We say that a ``Map`` is uniform if each process
owns the same number of global indices. The above ``Map`` is not uniform.
``Map`` includes both a constructor for uniform contiguous ``Map``\s, where you
specify the total number of global indices, and a constructor for possibly
nonuniform contiguous Maps, where you specify the number of global indices owned
by each process.

Globally Distributed or Locally Replicated
==========================================

Globally distributed means that both of the following are true:

* The ``Map``'s communicator has more than one process.
* There is at least one process in the ``Map``'s communicator, whose local
  number of entries does not equal the number of global entries. (That is, not
  all the entries are replicated over all the processes.)

If at least one of the above are not true, then we call the ``Map`` locally
replicated. The two terms are mutually exclusive.

.. include:: /links.txt
