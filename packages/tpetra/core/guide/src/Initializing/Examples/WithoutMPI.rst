.. _init_without_mpi:

Example: Initialization for an existing non-MPI code
####################################################

.. rubric:: Keywords

Initialization, Serial, Non-MPI

Overview
========

If Trilinos was built with MPI disabled, ``Teuchos::GlobalMPISession`` can safely by invoked in the ``main()`` function, in which case it will do nothing. However, if Trilinos was built with MPI enabled, but you don't want to use MPI within an application, ``Teuchos::GlobalMPISession`` should not be used. Instead, a ``Teuchos::SerialComm`` should be created directly as the "communicator." The ``Teuchos::SerialComm`` doesn't actually "communicate," because it only has one process.  With a "serial" communicator, the rank is always 0, and the number of processes is always 1.

The Initialization Program
==========================

The following source code listing demonstrates how to create and use a ``Teuchos::SerialComm``.

.. only:: builder_html

   The source code can be downloaded from :download:`here </Examples/SourceCode/initializing_tpetra_with_teuchos_serial.cpp>`.

.. literalinclude:: /Examples/SourceCode/initializing_tpetra_with_teuchos_serial.cpp
   :language: c++
   :linenos:
   :lines: 41-74
