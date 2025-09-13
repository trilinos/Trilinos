================================
Using MueLu in User Applications
================================

This tutorial demonstrates how to use MueLu from within user applications in ``C++``.
We will use the ``Tpetra`` linear algebra stack.
We will read the configuration for MueLu from an XML file and then create a ``MueLu::Hierarchy``.
This will then be used as a preconditioner within Belos as well as a standalone solver.

.. note::
   There is also support for Stratimikos.
   Please refer to the ``examples`` in the MueLu folder for more details.

XML Interface using CreateTpetraPreconditioner
==============================================
The most comfortable way to declare the multigrid parameters for MueLu is using the XML interface.
In fact, MueLu provides two different XML interfaces.
There is a simplified XML interface for multigrid users and a more advanced XML interface for experts,
which allows to make use of all features of MueLu as a multigrid framework.
Both XML file formats are introduced in the previous sections of this hands on tutorial.
However, for the C++ code it makes no difference which type of XML interface is used.

.. note::
  In the next sections, we give some code snippets.
  THey are borrowed form the ``TutorialDriver.cpp`` file in the tutorial tests.

Preparations
------------
First of all, we need to grab a communicator object.
Therefore, it is easy to use some utilities from the Teuchos package:

.. literalinclude:: ../../../test/tutorial/TutorialDriver.cpp
  :language: cpp
  :start-after: CommunicatorObject begin
  :end-before: CommunicatorObject end

For the multigrid method, we need a linear operator :math:`A`.
We can generate a problem-dependent matrix operator using the Galeri package (see :ref:`quick_start/example problem`), with Tpetra as our underlying linear algebra framework (through Xpetra wrappers). Below is an example for the Laplace2D Problem:

.. literalinclude:: ../../../test/tutorial/Tutorial_cppInterface.cpp
  :language: cpp
  :start-after: 2DLaplacianOperator begin
  :end-before: 2DLaplacianOperator end

For aggregation-based algebraic multigrid methods,
one has to provide a valid set of near null space vectors to produce transfer operators.
In case of a Laplace problem, we just use a constant vector.

.. literalinclude:: ../../../test/tutorial/TutorialDriver.cpp
  :language: cpp
  :start-after: BuildNullSpaceVector begin
  :end-before: BuildNullSpaceVector end

Setup phase
-----------
With a fine level operator :math:`A` available as ``Tpetra::CrsMatrix`` object
and a set of near null space vectors (available as ``Tpetra::MultiVector``),
all minimum requirements are fulfilled for generating an algebraic multigrid hierarchy.
There are two different ways to setup a multigrid hierarchy in MueLu.
One can either use a parameter list driven setup process which accepts either ``Teuchos::ParameterList`` objects
or XML files in two different XML file formats.
Alternatively, one can use the MueLu ``C++`` API to define the multigrid setup at compile time.
In the next sections we show both variants.

The most comfortable way to declare the multigrid parameters for MueLu is using the XML interface.
In fact, MueLu provides two different XML interfaces.
There is a simplified XML interface for multigrid users and a more advanced XML interface for expert which allows to make use of all features of MueLu as a multigrid framework.
Both XML file formats are introduced in the previous sections of this hands on tutorial.
However, for the ``C++`` code it makes no difference which type of XML interface is used.

We first read the MueLu configuration from the XML file:

.. literalinclude:: ../../../test/tutorial/TutorialDriver.cpp
  :language: cpp
  :start-after: ReadMueLuParamsFromXmlFile begin
  :end-before: ReadMueLuParamsFromXmlFile end

Then, we store the near null space vectors in the ``"user data"`` sublist of this parameter list:

.. literalinclude:: ../../../test/tutorial/TutorialDriver.cpp
  :language: cpp
  :start-after: InsertNullspaceInUserData begin
  :end-before: InsertNullspaceInUserData end

We then can create a MueLu object ready to be used as a preconditioner:

.. literalinclude:: ../../../test/tutorial/TutorialDriver.cpp
  :language: cpp
  :start-after: CreateTpetraPreconditioner begin
  :end-before: CreateTpetraPreconditioner end

.. _user_api/iteration phase:

Iteration Phase
---------------
Once the setup phase is completed, the MueLu multigrid hierarchy is ready for being used.

There are several ways how to use the multigrid method.
One can apply the multigrid method as standalone solver for linear systems.
Multigrid methods are also known to be efficient preconditioners within iterative (Krylov) solvers such as CG or GMRES methods.

In the next subsections it is demonstrated how to use MueLu as standalone solver and as preconditioner for iterative solvers from the Belos and AztecOO package in Trilinos.

For solving a linear system :math:`Ax=b`, we need a right hand side vector :math:`b`.
When using iterative solvers we also need an initial guess for the solution vector.

.. literalinclude:: ../../../test/tutorial/TutorialDriver.cpp
  :language: cpp
  :start-after: SetRhsAndSolutionVector begin
  :end-before: SetRhsAndSolutionVector end

In this example we just create ``Tpetra::MultiVectors`` (with just a single vector each).
The right-hand side vector :math:`b` is initialized with ones and the solution vector :math:`x` is filled with random values.

.. _user_api/muelu as solver:

MueLu as multigrid solver
-------------------------
MueLu can be used as standalone solver.
First, we create and initialize a solution vector:

.. literalinclude:: ../../../test/tutorial/TutorialDriver.cpp
  :language: cpp
  :start-after: MueLuAsSolverCreateSolutionVector begin
  :end-before: MueLuAsSolverCreateSolutionVector end

If necessary, extract the ``MueLu::Hiearchy`` from the Tpetra preconditioner object:

.. literalinclude:: ../../../test/tutorial/TutorialDriver.cpp
  :language: cpp
  :start-after: ExtractHierarchyFromTpetraPrec begin
  :end-before: ExtractHierarchyFromTpetraPrec end

Then, the ``MueLu::Hierarchy`` object is set to the non-preconditioner mode:

.. literalinclude:: ../../../test/tutorial/TutorialDriver.cpp
  :language: cpp
  :start-after: MueLuAsSolverSetSolverMode begin
  :end-before: MueLuAsSolverSetSolverMode end

Finally, we solve the system by calling the ``Iterate()`` routine
to perform ``mgridSweeps`` sweeps with the chosen multigrid cycle.

.. literalinclude:: ../../../test/tutorial/TutorialDriver.cpp
  :language: cpp
  :start-after: MueLuAsSolverIterate begin
  :end-before: MueLuAsSolverIterate end

If successful, the ``multigridSolVec`` vector contains the solution.

.. _user_api/muelu as preconditioner for belos:

MueLu as preconditioner for Belos
---------------------------------
Belos is the Krylov package in Trilinos and works both for Epetra and Tpetra.
Here, we demonstrate how to use MueLu as preconditioner for Belos solvers using Tpetra.

First, we create and initialize a solution vector:

.. literalinclude:: ../../../test/tutorial/TutorialDriver.cpp
  :language: cpp
  :start-after: MueLuAsPrecCreateSolutionVector begin
  :end-before: MueLuAsPrecCreateSolutionVector end

Then, we create a ``Belos::LinearProblem`` and hand in the MueLu preconditioner object:

.. literalinclude:: ../../../test/tutorial/TutorialDriver.cpp
  :language: cpp
  :start-after: MueLuAsPrecSetupLinearSystem begin
  :end-before: MueLuAsPrecSetupLinearSystem end

Now, we define the linear solver configuration in a ``Teuchos::ParameterList``
and create the Belos solver with this configuration:

.. literalinclude:: ../../../test/tutorial/TutorialDriver.cpp
  :language: cpp
  :start-after: MueLuAsPrecConfigureAndCreateBelosSolver begin
  :end-before: MueLuAsPrecConfigureAndCreateBelosSolver end

Finally, we solve the linear system:

.. literalinclude:: ../../../test/tutorial/TutorialDriver.cpp
  :language: cpp
  :start-after: MueLuAsPrecSolve begin
  :end-before: MueLuAsPrecSolve end

Full example using the XML interface
------------------------------------
The reader may refer to ``TutorialDriver.cpp`` for a working example to study the source code.
This demonstration program has some more features that are not discussed in this tutorial.

.. admonition:: Exercise 1

	Compile the example in ``TutorialDriver.cpp`` and then run the program in parallel using two processors

        ``mpirun -np 2 ./MueLu_TutorialDriver.exe --help``

        Study the screen output and try to run the example with an XML file as input for the multigrid setup.

.. admonition:: Exercise 2

	Create large scale examples using the ``--nx`` and ``--ny`` parameters for a finer mesh.
  Choose reasonable numbers for ``--nx`` and ``--ny`` for your machine
  and make use of your knowledge about MueLu for generating efficient preconditioners.

C++ Interface
=============
As an alternative to the XML interfaces, the user can also define the multigrid hierarchy using the ``C++`` API directly.
In contrary to the XML interface, which allows to build the layout of the multigrid preconditioner at runtime,
the preconditioner is fully defined at compile time when using the ``C++`` interface.

First, a ``MueLu::Hierarchy`` object has to be defined, which manages the multigrid hierarchy including all multigrid levels.
It provides routines for the multigrid setup and the multigrid cycle algorithms (such as V-cycle and W-cycle).

.. literalinclude:: ../../../test/tutorial/Tutorial_cppInterface.cpp
  :language: cpp
  :start-after: CreateNewHierarchy begin
  :end-before: CreateNewHierarchy end

There are some member functions which can be used to describe the basic multigrid hierarchy.
The ``SetMaxCoarseSize`` member function is used to set the maximum size of the coarse level problem before the coarsening process can be stopped.

.. literalinclude:: ../../../test/tutorial/Tutorial_cppInterface.cpp
  :language: cpp
  :start-after: InstantiateNewHierarchyObject begin
  :end-before: InstantiateNewHierarchyObject end

Next, one defines an empty ``MueLu::Level`` object for the finest level.
The ``MueLu::Level`` objects represent a data container storing the internal variables on each multigrid level.
The user has to provide and fill the level container for the finest level only.
The ``MueLu::Hierarchy`` object then automatically generates the coarse levels using the multigrid parameters.
The absolute minimum requirements for the finest level that the user has to provide is the fine level operator :math:`A` which represents the fine level matrix.
MueLu is based on Xpetra. So, the matrix :math:`A` has to be of type ``Xpetra::Matrix``.
In addition, the user should also provide a valid set of near null space vectors.
For a Laplace problem we can just use the constant ``nullspace`` vector that has previously been defined.
Some routines need additional information.
For example, the user has to provide the node coordinates for repartitioning.

.. literalinclude:: ../../../test/tutorial/Tutorial_cppInterface.cpp
  :language: cpp
  :start-after: CreateFineLevelObject begin
  :end-before: CreateFineLevelObject end

.. note::
	When including the ``MueLu_UseShortNames.hpp`` header file,
  the template parameters usually can be dropped for compiling.
  The most important template parameters are ``SC`` for the scalar type,
  ``LO`` for the local ordinal type and ``GO`` for the global ordinal type.
  For a detailed description of the template parameters,
  the reader may refer to the Tpetra documentation.

A ``MueLu::FactoryManager`` object is used for the internal management of data dependencies and generating algorithms of the multigrid setup.
Even though not absolutely necessary,
we show the usage of the ``MueLu::FactoryManager`` object as it allows for user-specific enhancements of the multigrid code.

.. literalinclude:: ../../../test/tutorial/Tutorial_cppInterface.cpp
  :language: cpp
  :start-after: DefineFactoryManager begin
  :end-before: DefineFactoryManager end

The user can define its own factories for performing different tasks in the setup process.
The following code shows how to define a smoothed aggregation transfer operator and a restriction operator.
The ``MueLu::RAPFactory`` is used for the (standard) Galerkin product to generate the coarse level matrix :math:`A`.

.. literalinclude:: ../../../test/tutorial/Tutorial_cppInterface.cpp
  :language: cpp
  :start-after: DeclareSomeFactories begin
  :end-before: DeclareSomeFactories end

The user-defined factories have to be registered in the ``FactoryManager`` using the lines

.. literalinclude:: ../../../test/tutorial/Tutorial_cppInterface.cpp
  :language: cpp
  :start-after: ConfigureFactoryManager begin
  :end-before: ConfigureFactoryManager end

.. warning::
	If you forget to register the new factories, the ``FactoryManager`` will use some internal default factories for being responsible to create the corresponding variables.
  Then your user-specified factories are just ignored during the multigrid setup!

.. note::
	The ``FactoryManager`` is also responsible for resolving all dependencies between different factories.
  That is, after the user-defined factories have been registered,
  all factories that request variable :math:`P` are provided with the prolongation operator :math:`P` that has been generated by the registered factory ``PFact``.
  If there is some data requested for which no factory has been registered by the user,
  the ``FactoryManager`` manages an internal list for reasonable default choices and default factories.

Next, the user has to declare a level smoother.
The following code can be used to define a symmetric Gauss-Seidel smoother.
Other methods can be set up in a similar way.

.. literalinclude:: ../../../test/tutorial/Tutorial_cppInterface.cpp
  :language: cpp
  :start-after: DefineSmootherObject begin
  :end-before: DefineSmootherObject end

Before the level smoother can be used, a ``MueLu::SmootherFactory`` has to be defined for the smoother factory.
The ``SmootherFactory`` is used in the multigrid setup to generate level smoothers for the corresponding levels using the prototyping design pattern.
Note, that the ``SmootherFactory`` has also to be registered in the ``FactoryManager`` object.
If the user forgets this, the multigrid setup will use some kind of default smoother, i.e., the user-chosen smoother options are just ignored.

.. literalinclude:: ../../../test/tutorial/Tutorial_cppInterface.cpp
  :language: cpp
  :start-after: CreateSmootherFactory begin
  :end-before: CreateSmootherFactory end

Once the ``FactoryManager`` is set up, it can be used with the ``Hierarchy::Setup`` routine to initiate the coarsening process and set up the multigrid hierarchy.

.. literalinclude:: ../../../test/tutorial/Tutorial_cppInterface.cpp
  :language: cpp
  :start-after: SetupMultigridHierarchy begin
  :end-before: SetupMultigridHierarchy end

Footnotes
=========
.. [1] L. Berger-Vergiat, C. A. Glusa, G. Harper, J. J. Hu, M. Mayr, P. Ohm, A. Prokopenko, C. M. Siefert, R. S. Tuminaro, and T. A. Wiesner. MueLu User's Guide. Technical Report SAND2023-12265, Sandia National Laboratories, Albuquerque, NM (USA) 87185, 2023.
