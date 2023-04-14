================================
Using MueLu in User Applications
================================

This tutorial demonstrates how to use MueLu from within user applications in **C++**.
We will use the **Tpetra** linear algebra stack.
We will read the configuration for MueLu from an xml file and then create a **MueLu::Hierarchy**.
This will then be used as a preconditioner within Belos as well as a standalone solver.

.. note::
   There is also support for Stratimikos.
   Please refer to the **examples** in the MueLu folder for more details.

In the next sections, we give some code snippets.
Most of them are borrowed form the **laplace2d.cpp** file in the tutorial.

Preparations
============
First of all, we extract the template parameters from Tpetra:

.. literalinclude:: ../../../test/tutorial/laplace2d.cpp
  :language: cpp
  :start-after: TpetraTemplateParameters begin
  :end-before: TpetraTemplateParameters end

Then, we define some abbreviations for better readibility of the code:

.. literalinclude:: ../../../test/tutorial/laplace2d.cpp
  :language: cpp
  :start-after: UsingStatements begin
  :end-before: UsingStatements end

Now, we need to grab a communicator object.
Therefore, it is easy to use some utilities from the Teuchos package:

.. literalinclude:: ../../../test/tutorial/laplace2d.cpp
  :language: cpp
  :start-after: CommunicatorObject begin
  :end-before: CommunicatorObject end

For the multigrid method, we need a linear operator :math:`A`.
The linear operator is usually assembled by the application code.
Hence, we will not covered it in detail in this tutorial.
For the sake of simplicity,
let us just assume that we have a Laplace operator in two dimensions readily available
and stored in a **Teuchos::RCP<Tpetra::CrsMatrix<SC,LO,GO,NO>> matrix**.

For aggregation-based algebraic multigrid methods,
one has to provide a valid set of near null space vectors to produce transfer operators.
In case of a Laplace problem, we just use a constant vector.

.. literalinclude:: ../../../test/tutorial/laplace2d.cpp
  :language: cpp
  :start-after: BuildNullSpaceVector begin
  :end-before: BuildNullSpaceVector end

Setup phase
===========
With a fine level operator :math:`A` available as **Xpetra::Matrix** object and a set of near null space vectors (available as **Xpetra::MultiVector**),
all minimum requirements are fulfilled for generating an algebraic multigrid hierarchy.
There are two different ways to setup a multigrid hierarchy in MueLu.
One can either use a parameter list driven setup process which accepts either **Teuchos::ParameterList** objects
or XML files in two different XML file formats.
Alternatively, one can use the MueLu C++ API to define the multigrid setup at compile time.
In the next sections we show both variants.

XML Interphase using CreateTpetraPreconditioner
-----------------------------------------------
The most comfortable way to declare the multigrid parameters for MueLu is using the XML interface.
In fact, MueLu provides two different XML interfaces.
There is a simplified XML interface for multigrid users and a more advanced XML interface for expert which allows to make use of all features of MueLu as a multigrid framework.
Both XML file formats are introduced in the previous sections of this hands on tutorial.
However, for the C++ code it makes no difference which type of XML interface is used.

We first read the MueLu configuration from the xml file:

.. literalinclude:: ../../../test/tutorial/laplace2d.xpp
  :language: cpp
  :start-after: ReadMueLuParamsFromXmlFile begin
  :end-before: ReadMueLuParamsFromXmlFile end

Then, we store the near null space vectors in the **"user data"** sublist of this parameter list:

.. literalinclude:: ../../../test/tutorial/laplace2d.xpp
  :language: cpp
  :start-after: InsertNullspaceInUserData begin
  :end-before: InsertNullspaceInUserData end

We then can create a MueLu object ready to be used as a preconditioner:

.. literalinclude:: ../../../test/tutorial/laplace2d.xpp
  :language: cpp
  :start-after: CreateTpetraPreconditioner begin
  :end-before: CreateTpetraPreconditioner end

.. _user_api/iteration phase:

Iteration Phase
===============
Once the setup phase is completed, the MueLu multigrid hierarchy is ready for being used.

There are several ways how to use the multigrid method.
One can apply the multigrid method as standalone solver for linear systems.
Multigrid methods are also known to be efficient preconditioners within iterative (Krylov) solvers such as CG or GMRES methods.

In the next subsections it is demonstrated how to use MueLu as standalone solver and as preconditioner for iterative solvers from the Belos and AztecOO package in Trilinos.

For solving a linear system :math:`Ax=b`, we need a right hand side vector :math:`b`.
When using iterative solvers we also need an initial guess for the solution vector.

.. literalinclude:: ../../../test/tutorial/laplace2d.cpp
  :language: cpp
  :start-after: RhsAndSolutionVector begin
  :end-before: RhsAndSolutionVector end

In this example we just **Tpetra::MultiVectors**.
The right hand side vector is initialized with one and the solution vector is filled with random values.

MueLu as multigrid solver
-------------------------
To use MueLu as standalone solver, we assume to have a **MueLu::Hierarchy** object.
Depending on the way of setting up the AMG preconditioner,
we might need to extract the hierarchy from the Tpetra preconditioner object:

.. literalinclude:: ../../../test/tutorial/laplace2d.xpp
  :language: cpp
  :start-after: ExtractHierarchyFromTpetraPrec begin
  :end-before: ExtractHierarchyFromTpetraPrec end

Then, we run the following code:

.. literalinclude:: ../../../test/tutorial/laplace2d.xpp
  :language: cpp
  :start-after: UseMultigridHierarchyAsSolver begin
  :end-before: UseMultigridHierarchyAsSolver end

The **MueLu::Hierarchy** object is set to the non-preconditioner mode and the **Iterate** routine is called
to perform **mgridSweeps** sweeps with the chosen multigrid cycle.
If successful, the **multigridSolVec** vector contains the solution.

.. _user_api/muelu as preconditioner for belos:

MueLu as preconditioner for Belos
---------------------------------
Belos is the Krylov package in Trilinos and works both for Epetra and Tpetra.
Here, we demonstrate how to use MueLu as preconditioner for Belos solvers using Tpetra.
First, we create a **Belos::LinearProblem**, equip it with the MueLu preconditioner object,
configure the iterative solver via a parameter list, and finally solver the linear system:

.. literalinclude:: ../../../test/tutorial/laplace2d.xpp
  :language: cpp
  :start-after: MueLuHierarchyAsPreconditionerWithinBelos begin
  :end-before: MueLuHierarchyAsPreconditionerWithinBelos end

Full example
============
The reader may refer to **laplace2d.cpp** for a working example to study the source code.
This demonstration program has some more features that are not discussed in this tutorial.

.. admonition:: Exercise 1

	Compile the example in **laplace2d.cpp** and then run the program in parallel using two processors

        ``mpirun -np 2 ./MueLu_tutorial_laplace2d.exe --help``

        Study the screen output and try to run the example with an XML file as input for the multigrid setup.

.. admonition:: Exercise 2

	Create large scale examples using the nx and ny parameters for a finer mesh.
  Choose reasonable numbers for nx and ny for your machine and make use of your knowledge about MueLu for generating efficient preconditioners.

Footnotes
=========
.. [1] L. Berger-Vergiat, C. A. Glusa, G. Harper, J. J. Hu, M. Mayr, P. Ohm, A. Prokopenko, C. M. Siefert, R. S. Tuminaro, and T. A. Wiesner. MueLu User's Guide. Technical Report SAND2023-12265, Sandia National Laboratories, Albuquerque, NM (USA) 87185, 2023.
