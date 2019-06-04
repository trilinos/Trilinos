.. _init_existing_mpi:

Example: Initialization for an existing MPI code
################################################

.. rubric:: Keywords

Initialization, MPI

Overview
========

For a code that initializes and finalizes MPI on its own by calling ``MPI_Init``
and ``MPI_Finalize``, an ``MPI_Comm`` (MPI communicator) must be made available
to Tpetra (either by using a predefined communicator such as ``MPI_COMM_WORLD``,
or by creating a new one).

The following example demonstrates how to initialize MPI and wrap an MPI communicator
in such a way that Trilinos understands that the program is responsible for
calling ``MPI_Comm_free`` on the ``MPI_Comm`` after use, if necessary.  (It's
not necessary for MPI_COMM_WORLD.)  There is a way to tell Trilinos to call
``MPI_Comm_free`` itself; though it is not shown here.  (It involves passing the
result of ``Teuchos::opaqueWrapper`` to ``MpiComm``'s constructor.)

The Initialization Program
==========================

The following source code listing demonstrates initializing Tpetra for an
existing MPI code.

.. only:: builder_html

   The source code can be downloaded from :download:`here </Examples/SourceCode/initializing_tpetra_with_standalone_mpi.cpp>`.

.. literalinclude:: /Examples/SourceCode/initializing_tpetra_with_standalone_mpi.cpp
   :language: c++
   :linenos:
   :lines: 43-83

If your code uses MPI on its own, as well as through Trilinos, consider giving
Trilinos a copy of your ``MPI_Comm`` (created via ``MPI_Comm_dup``) rather than
your ``MPI_Comm`` directly.  Trilinos may in the future duplicate the
``MPI_Comm`` automatically, but it does not currently do this.  Duplicating the
``MPI_Comm`` is not necessary, but may make it easier for you to overlap
asynchronous communication operations performed by Trilinos with those performed
by your code.
