================================
Using MueLu in User Applications
================================

This tutorial demonstrates how to use MueLu from within user applications in ``C++``.
In [[1]_, Section 2.6] it is explained how to use MueLu through the ``MueLu::CreateE/T/XpetraPreconditioner`` interface.
This interface is designed for beginners which want to try MueLu through standard Trilinos interfaces.

.. note::
   There is also support for Stratimikos.
   Please refer to the ``examples`` in the MueLu folder for more details.

This tutorial aims at more advanced methods to use MueLu
such as creating an explicit instance of some MueLu classes like ``MueLu::Hierarchy`` and ``MueLu::HierarchyManager``.
In the next sections, we give some code snippets.
Most of them are borrowed form the ``laplace2d.cpp`` file in the tutorial.

Preparations
============
First of all, we have to define a communicator object.

.. literalinclude:: ../../../test/tutorial/laplace2d.cpp
  :language: cpp
  :start-after: CommunicatorObject begin
  :end-before: CommunicatorObject end 

For the multigrid method we need a linear operator :math:`A`. For demonstration purposes, here we hust generate a 2D Laplacian operator using the Galeri package (see :ref:`quick_start/example problem`). In this example we use Epetra for the underlying linear algebra framework, but it shall be mentioned that it works for Tpetra in a similar way (refer to the code examples in the MueLu examples folder).

.. literalinclude:: ../../../test/tutorial/laplace2d.cpp
  :language: cpp
  :start-after: 2DLaplacianOperator begin
  :end-before: 2DLaplacianOperator end 
	
Muelu is based on Xpetra which provides a common interface both for Epetra and Tpetra. Therefore we have to enxapsulate our Epetra objects into Xpetra wrapper objects. This is done using the following code.

.. literalinclude:: ../../../test/tutorial/laplace2d.cpp
  :language: cpp
  :start-after: EpetraToXpetra begin
  :end-before: EpetraToXpetra end 
	
.. note::
	The MueLu setup routines require a ``Xpetra::Matrix`` object. 

The wrapper class ``Xpetra::CrsMatrixWrap`` is just a wrapper derived from ``Xpetra::Matrix`` which manages a ``Xpetra::CrsMatrix`` object which is the common base class for both Epetra and Tpetra CRS matrix classes. The details are not really important as long as one understands that one needs a ``Xpetra::Matrix`` object for MueLu in the end. With the ``SetFixedBlockSizeroutine`` we state that there is only one degree of freedom per node (pure Laplace problem). For aggregation-based algebraic multigrid methods one has to provide a valid set of near null space vectors to produce transfer operators. In case of a Laplace problem we just use a constant vector.

.. literalinclude:: ../../../test/tutorial/laplace2d.cpp
  :language: cpp
  :start-after: BuildNullSpaceVector begin
  :end-before: BuildNullSpaceVector end 
	
Setup phase
===========
With a fine level operator :math:`A` available as ``Xpetra::Matrix`` object and a set of near null space vectors (available as ``Xpetra::MultiVector``) all minimum requirements are fulfilled for generating an algebraic multigrid hierarchy. There are two different ways to setup a multigrid hierarchy in MueLu. One can either use a parameter list driven setup process which accepts either ``Teuchos::ParameterList`` objects or XML files in two different XML file formats. Alternatively, one can use the MueLu C++ API to define the multigrid setup at compile time. In the next sections we show both variants.

XML Interphase
--------------
The most comfortable way to declare the multigrid parameters for MueLu is using the XML interface. In fact, MueLu provides two different XML interfaces. There is a simplified XML interface for multigrid users and a more advanced XML interface for expert which allows to make use of all features of MueLu as a multigrid framework. Both XML file formats are introduced in the previous sections of this hands on tutorial.
However, for the C++ code it makes no difference which type of XML interface is used.

Assuming that we have a ``Teuchos::ParameterList`` object with valid MueLu parameters we can create a ``MueLu::HierarchyManager`` object

.. literalinclude:: ../../../test/tutorial/ScalingTestParamList.cpp
  :language: cpp
  :start-after: HierarchyManager begin
  :end-before: HierarchyManager end 
	
For an example how to fill the parameter list the reader may refer to [[1]_, Section 2.3]. Note that there are routines to fill the parameter lists with the information from XML files. You can also directly provide a file name of a XML file to the ``Muelu::ParameterListInterpreter``. For details you may refer to the doxygen documentation or the example in ``laplace2d.cpp``.

Next a new ``MueLu::Hierarchy`` object is generated

.. literalinclude:: ../../../test/tutorial/ScalingTestParamList.cpp
  :language: cpp
  :start-after: HierarchyObject begin
  :end-before: HierarchyObject end 
  
The ``CreateHierarchy`` creates a new empty multigrid hierarchy with a finest level only. The user has to feed in the linear operator :math:`A` and the near null space vector. If further information is available, such as the node coordinates, they can be also stored in the finest level. The coordinates are needed, e.g., for rebalancing the coarse levels.
Finally, the ``SetupHierarchy`` call initiates the coarsening process and the multigrid hierarchy is built according to the parameters from the ``mueluList`` parameters.

.. literalinclude:: ../../../test/tutorial/ScalingTestParamList.cpp
  :language: cpp
  :start-after: SetupHierarchyCall begin
  :end-before: SetupHierarchyCall end 
  
As XML parameter file any of the files shown in the previous tutorials can be used.
  
.. note::
 	As one can see from the last code snippet, the ``Hierarchy`` allows access to all important parts of the multigrid method
 	before setup. So, if you have to feed in some non-standard information, this is the way how it works. Using the 
 	``XreateE/TpetraPreconditioner`` interface may be easier but does not allow to access the finest level before setup.

Once the ``SetupHierarchy`` call is completed, the multigrid hierarchy is ready to use. The reader can skip the next section about the C++ interface and proceed with :ref:`user_api/MueLu as preconditioner for Belos` for an example how to use the multigrid method as preconditioner within a Krylov subspace method from the Belos package.

C++ Interface
-------------
As an alternative to the XML interfaces, the user can also define the multigrid hierarchy using the C++ API directly. In contrary to the XML interface which allows to build the layout of the multigrid preconditioner at runtime, the preconditioner is fully defined at compile time when using the C++ interface.

First, a ``MueLu::Hierarchy`` object has to be defined, which manages the multigrid hierarchy including all multigrid levels. It provides routines for the multigrid setup and the multigrid cycle algorithms (such as V-cycle and W-cycle).
 	
.. literalinclude:: ../../../test/tutorial/ScalingTest.cpp
  :language: cpp
  :start-after: CreateNewHierarchy begin
  :end-before: CreateNewHierarchy end 

There are some member functions which can be used to describe the basic multigrid hierarchy. The ``SetMaxCoarseSize`` member function is used to set the maximum size of the coarse level problem before the coarsening process can be stopped.

.. literalinclude:: ../../../test/tutorial/ScalingTest.cpp
  :language: cpp
  :start-after: InstantiateNewHierarchyObject begin
  :end-before: InstantiateNewHierarchyObject end 

Next, one defines an empty ``MueLu::Level`` object for the finest level. The ``MueLu::Level`` objects represent a data container storing the internal variables on each multigrid level. The user has to provide and fill the level container for the finest level only. The ``MueLu::Hierarchy`` object then automatically generates the coarse levels using the multigrid parameters. The absolute minimum requirements for the finest level that the user has to provide is the fine level operator :math:`A` which represents the fine level matrix. Muelu is based on Xpetra. So, the matrix :math:`A` has to be of type ``Xpetra::Matrix``. In addition, the user should also provide a valid set of near null space vectors. For a Laplace problem we can just use the constant ``nullspace`` vector that has previously been defined. Some routines need additional information. For example, the user has to provide the node coordinates for repartitioning. 

.. literalinclude:: ../../../test/tutorial/ScalingTest.cpp
  :language: cpp
  :start-after: CreateFineLevelObject begin
  :end-before: CreateFineLevelObject end 

.. note::
	When including the ``MueLu\_UseShortNames.hpp`` header file the template parameters usually can be dropped for compiling. The most important template parameters are ``SC`` for the scalar type, ``LO`` for the local ordinal type (usually ``int}``) and ``GO`` for the global ordinal type. For a detailed description of the template parameters the reader may refer to the Tpetra documentation.

A ``MueLu::FactoryManager`` object is used for the internal management of data dependencies and generating algorithms of the multigrid setup. Even though not absolutely necessary, we show the usage of the ``MueLu::FactoryManager`` object as it allows for  user-specific enhancements of the multigrid code.

.. literalinclude:: ../../../test/tutorial/ScalingTest.cpp
  :language: cpp
  :start-after: DefineFactoryManager begin
  :end-before: DefineFactoryManager end 

The user can define its own factories for performing different tasks in the setup process. The following code shows how to define a smoothed aggregation transfer operator and a restriction operator. The ``MueLu::RAPFactory`` is used for the (standard) Galerkin product to generate the coarse level matrix :math:`A`.

.. literalinclude:: ../../../test/tutorial/ScalingTest.cpp
  :language: cpp
  :start-after: DeclareSomeFactories begin
  :end-before: DeclareSomeFactories end
  
The user-defined factories have to be registered in the ``FactoryManager`` using the lines
 
.. literalinclude:: ../../../test/tutorial/ScalingTest.cpp
  :language: cpp
  :start-after: ConfigureFactoryManager begin
  :end-before: ConfigureFactoryManager end  
  
.. warning::
	If you forget to register the new factories, the ``FactoryManager`` will use some internal default factories for being responsible to create the corresponding variables. Then your user-specified factories are just ignored during the multigrid setup!
	
.. note::
	The ``FactoryManager`` is also responsible for resolving all dependencies between different factories. That is, after the user-defined factories have been registered, all factories that request variable :math:`P` are provided with the prolongation operator :math:`P` that has been generated by the registered factory ``PFact``. If there is some data requested for which no factory has been registered by the user, the ``FactoryManager`` manages an internal list for reasonable default choices and default factories.

Next, the user has to declare a level smoother. The following code can be used to define a symmetric Gauss-Seidel smoother. Other methods can be set up in a similar way.

.. literalinclude:: ../../../test/tutorial/ScalingTest.cpp
  :language: cpp
  :start-after: DefineSmootherObject begin
  :end-before: DefineSmootherObject end 
  
Before the level smoother can be used, a ``MueLu::SmootherFactory`` has to be defined for the smoother factory. The ``SmootherFactory`` is used in the multigrid setup to generate level smoothers for the corresponding levels using the prototyping design pattern. Note, that the ``SmootherFactory`` has also to be registered in the ``FactoryManager`` object. If the user forgets this, the multigrid setup will use some kind of default smoother, i.e., the user-chosen smoother options are just ignored. 

.. literalinclude:: ../../../test/tutorial/ScalingTest.cpp
  :language: cpp
  :start-after: CreateSmootherFactory begin
  :end-before: CreateSmootherFactory end 
  
Once the ``FactoryManager`` is set up, it can be used with the ``Hierarchy::Setup`` routine to initiate the coarsening process and set up the multigrid hierarchy.

.. literalinclude:: ../../../test/tutorial/ScalingTest.cpp
  :language: cpp
  :start-after: SetupMultigridHierarchy begin
  :end-before: SetupMultigridHierarchy end 
  
.. _user_api/iteration phase:

Iteration Phase
===============
Once the setup phase is completed, the MueLu multigrid hierarchy is ready for being used.

There are several ways how to use the multigrid method. One can apply the multigrid method as standalone solver for linear systems. Multigrid methods are also known to be efficient preconditioners within iterative (Krylov) solvers such as CG or GMRES methods.

In the next subsections it is demonstrated how to use MueLu as standalone solver and as preconditioner for iterative solvers from the Belos and AztecOO package in Trilinos.

For solving a linear system :math:`Ax=b` we need a right hand side vector :math:`b`. When using iterative solvers we also need an initial guess for the solution vector.

.. literalinclude:: ../../../test/tutorial/laplace2d.cpp
  :language: cpp
  :start-after: SetRhsAndSolutionVector begin
  :end-before: SetRhsAndSolutionVector end 
 
In this example we just create Epetra vectors and wrap them into Xpetra objects. The right hand side vector is initialized with one and the solution vector is filled with random values.

MueLu as multigrid solver
-------------------------
To use MueLu as standalone solver one can use the following code

.. literalinclude:: ../../../test/tutorial/laplace2d.cpp
  :language: cpp
  :start-after: UseMultigridHierarchyAsSolver begin
  :end-before: UseMultigridHierarchyAsSolver end 

In this code snippet a solution vector is created using the ``Xpetra::VectorFactory`` and initialized with the content from the solution vector :math:`xX` containing the initial guess. Then, the ``MueLu::Hierarchy`` object is set to the non-preconditioner mode and the ``Iterate`` routine is called
to perform ``mgridSweeps`` sweeps with the chosen multigrid cycle. If successful, the ``mgridLsgVec`` vector contains the solution.

MueLu as preconditioner for AztecOO
-----------------------------------
Commonly, multigrid methods are used as preconditioners for iterative linear solvers. Here, we show how to use the ``MueLu::Hierarchy`` as preconditioner within an AztecOO solver (using Epetra).
After an Epetra solution vector has been created by

.. literalinclude:: ../../../test/tutorial/laplace2d.cpp
  :language: cpp
  :start-after: EpetraSolutionVector begin
  :end-before: EpetraSolutionVector end 
  
the following code can be used to apply the MueLu hierarchy as preconditioner within the AztecOO CG solver

.. literalinclude:: ../../../test/tutorial/laplace2d.cpp
  :language: cpp
  :start-after: MueLuHierarchyAsPreconditionerWithinAztecOO begin
  :end-before: MueLuHierarchyAsPreconditionerWithinAztecOO end 
 
Basically, the MueLu hierarchy is put into an ``MueLu::EpetraOperator`` object, which implements the Epetra interface for preconditioners.
With the ``SetPrecOperator`` routine from the AztecOO solver the ``MueLu::EpetraOperator`` object then is defined as preconditioner.

.. _user_api/muelu as preconditioner for belos:

MueLu as preconditioner for Belos
---------------------------------
Belos is the successor package of AztecOO for linear solvers in Trilinos and works both for Epetra and Tpetra. Here we demonstrate how to use MueLu as preconditioner for Belos solvers using Xpetra.
First, we have to declare objects for the solution vector and the right hand side vector in Xpetra. The following code just uses a random vector for the initial guess and solution variable.

.. literalinclude:: ../../../test/tutorial/ScalingTest.cpp
  :language: cpp
  :start-after: DefineXB begin
  :end-before: DefineXB end 

In the following we demonstrate how to use the MueLu hierarchy as preconditioner within a Belos solver. There are special wrapper objects for wrapping the Xpetra matrix and the MueLu hierarchy to Belos compatible objects. These can be used to define a linear problem for use with Belos.

.. literalinclude:: ../../../test/tutorial/ScalingTest.cpp
  :language: cpp
  :start-after: OperatorAndMultivectorTypeBelos begin
  :end-before: OperatorAndMultivectorTypeBelos end 

Then, one can set up the Belos solver. For a Belos GMRES solver one uses

.. literalinclude:: ../../../test/tutorial/ScalingTest.cpp
  :language: cpp
  :start-after: BelosParameterList begin
  :end-before: BelosParameterList end 
  
Finally, we can solve the linear system using Belos with the MueLu multigrid preconditioner (left-preconditioning) by calling

.. literalinclude:: ../../../test/tutorial/ScalingTest.cpp
  :language: cpp
  :start-after: SolveLinearSystem begin
  :end-before: SolveLinearSystem end 
  
and perform some convergence checks

.. literalinclude:: ../../../test/tutorial/ScalingTest.cpp
  :language: cpp
  :start-after: CheckConvergence begin
  :end-before: CheckConvergence end 
  
Full example
============
The reader may refer to ``laplace2d.cpp`` for a working example to study the source code. This demonstration program has some more features that are not discussed in this tutorial.

.. admonition:: Exercise 1

	Compile the example in ``laplace2d.cpp`` and then run the program in parallel using two processors
	
        mpirun -np 2 ./MueLu_tutorial_laplace2d.exe --help
        \end{verbatim}
        Study the screen output and try to run the example with an XML file as input for the multigrid setup.
.. admonition:: Exercise 2

	Create large scale examples using the nx and ny parameters for a finer mesh. Choose reasonable numbers for nx and ny for your machine and make use of your knowledge about MueLu for generating efficient preconditioners.



Footnotes
=========
.. [1] A. Prokopenko, J.J Hu, T.A. Wiesner, C.M. Siefert and R.S. Tuminaro ``MueLu User's Guide 1.0 (Trilinos Version 11.12)``, SAND2014-18874, 2014

