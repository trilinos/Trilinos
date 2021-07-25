============================
ML ParameterList interpreter
============================

Backwards compatibility
=======================

ML [1]_ is the predecessor multigrid package of MueLu in Trilinos and widely used in the community for smoothed aggregation multigrid methods. ML is implemented in C and known for its good performance properties. However, the disadvantage is that ML is harder to adapt to new applications and non-standard problems. Furthermore, ML uses its own internal data structure and is somewhat limited to the use with Epetra objects only. In contrast, MueLu provides a fully flexible multigrid framework which is designed to be adapted to any kind of new application with non-standard requirements. Furthermore, it is based on Xpetra and therefore can be used both with Epetra or Tpetra. Nevertheless, it is an important point to provide some kind of backwards compatibility to allow ML users to easily migrate to MueLu (or make experiments with MueLu without having to write to much new code).

In this tutorial we present the ``MueLu::MLParameterListInterpreter`` which provides support for the most important ML parameters to be used with MueLu.

C++ part
========
Preparations
------------
In order to use MueLu (instead or aside of ML) you first have to add it to your application. Please refer to the MueLu user guide for information about compilation and linking (see [2]_). Basically, if your application is already working with ML you should only need to compile and install MueLu and make sure that the MueLu libraries are found by the linker.

C++ interface
-------------
In the following we assume that the linear operator :math:`A` is available as ``RCP<Xpetra::Matrix> A``.

Then we create a parameter list and fill it with ML parameters. Please refer to the ML guide [1]_ for a complete list of available parameters.

.. literalinclude:: ../../../test/tutorial/MLParameterList.cpp
  :language: cpp
  :start-after: ParameterList begin
  :end-before: ParameterList end

.. note::
	Be aware that the MLParameterListInterpreter does not support all ML parameters but only the most important ones (e.g., smoothers, transfer operators, rebalancing, ...). There is, e.g., no support for the Maxwell specific enhancements in ML.

Instead of defining the ML parameters by hand in the ParameterList you can also read in XML files with ML parameters using

.. literalinclude:: ../../../test/tutorial/MLParameterList.cpp
  :language: cpp
  :start-after: GetParametersFromXMLFile begin
  :end-before: GetParametersFromXMLFile end
  
Next, you create a MLParameterListInterpreter object using the parameters and create a new ``MueLu::Hierarchy`` from it.

.. literalinclude:: ../../../test/tutorial/MLParameterList.cpp
  :language: cpp
  :start-after: MultigridHierarchy begin
  :end-before: MultigridHierarchy end

Of course, we have to provide all necessary information for the multigrid setup routine. This does not only include the fine level operator but also the set of near null space vectors. Assuming that ``numPDEs`` stores the number of equations (and near null space vectors) the following code allows to produce piecewise constant standard near null space vectors (which should be valid for many PDE discretizations).

.. literalinclude:: ../../../test/tutorial/MLParameterList.cpp
  :language: cpp
  :start-after: BuildDefaultNullSpace begin
  :end-before: BuildDefaultNullSpace end

Then we just feed in the information to the finest level

.. literalinclude:: ../../../test/tutorial/MLParameterList.cpp
  :language: cpp
  :start-after: FeedInInformation begin
  :end-before: FeedInInformation end

Finally we call the ``Setup`` routine which actually builds the multigrid hierarchy.

.. literalinclude:: ../../../test/tutorial/MLParameterList.cpp
  :language: cpp
  :start-after: CallSetupRoutine begin
  :end-before: CallSetupRoutine end

Once we have the multigrid hierarchy set up we can use it the same way as described in :ref:`Iteration Phase <user_api/iteration phase>`.

.. admonition:: Exercise 1

	Study the source code of ``../../../test/tutorial/MLParameterList.cpp`` and compile it. Run the executable ``MueLu_tutorial_MLParameterList.exe`` with the ``--help`` command line parameter to get an overview of all available command line parameters. Run the example using

		``./MueLu_tutorial_MLParameterList.exe --ml=1 --muelu=0``
		``--xml=xml/ml_ParameterList.xml --linAlgebra=Epetra``

	and study the ML output. Compare the output and results when switching to MueLu using the same input file

		``./MueLu_tutorial_MLParameterList.exe --ml=0 --muelu=1``
		``--xml=xml/ml_ParameterList.xml --linAlgebra=Epetra``
	

.. admonition:: Exercise 2

	Play around with the parameters from ``MueLu_tutorial_MLParameterList.exe``. Change, e.g., the problem type to a 2D Laplace problem (``--matrixType=Laplace2D``) and adapt the ``--nx`` and ``--ny`` parameters accordingly.
	
	Try to run both ML and MueLu and compare the results. Do you find significant differences?


Footnotes
=========
.. [1] M.W. Gee, C.M. Siefert, J.J. Hu, R.S. Tuminaro and M.G. Sala, ML 5.0 Smoothed Aggre-gation User’s Guide, Sandia National Laboratories, 2006, SAND2006-2649
.. [2] A. Prokopenko, J.J. Hu, T.A. Wiesner, C.M. Siefert and R.S. Tuminaro ``MueLu User’sGuide 1.0 (Trilinos Version 11.12)``, SAND2014-18874, 2014
