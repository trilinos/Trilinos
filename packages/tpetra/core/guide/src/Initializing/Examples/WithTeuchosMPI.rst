.. _init_with_teuchos:

Example: Initialization for a code that only uses MPI through Trilinos
######################################################################

.. rubric:: Keywords

Initialization, MPI, ``Teuchos``

Overview
========

Trilinos provides an MPI interface that allows initializing, finalizing, and querying the global MPI session: ``Teuchos::GlobalMPISession``.  ``Teuchos::GlobalMPISession`` calls ``MPI_Init`` and ``MPI_Finalize`` for you in an MPI build, and does not call them if you did not build Trilinos with MPI support.

Though Tpetra was written for distributed-memory parallel programming using MPI, it will work correctly whether or not Trilinos is built with MPI support. It does so by interacting with MPI through the ``Teuchos::Comm`` interface. (If you are familiar with Epetra, this interface is analogous to ``Epetra_Comm``.) If MPI is enabled, then ``Teuchos::Comm`` wraps an ``MPI_Comm``. Otherwise, it is a "serial communicator" with one process, analogous to ``MPI_COMM_SELF``.

The Initialization Program
==========================

The following source code listing demonstrates how to initialize MPI (if available) through ``Teuchos::GlobalMPISession`` and get a ``Teuchos::Comm`` communicator corresponding to ``MPI_COMM_WORLD``. The example works whether or not Trilinos was built with MPI support.

.. only:: builder_html

   The source code can be downloaded from :download:`here </Examples/SourceCode/initializing_tpetra_with_teuchos_mpi.cpp>`.

.. literalinclude:: /Examples/SourceCode/initializing_tpetra_with_teuchos_mpi.cpp
   :language: c++
   :linenos:
   :lines: 41-74
