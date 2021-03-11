============================
ML ParameterList interpreter
============================

Backwards compatibility
=======================

ML [1]_ is the predecessor multigrid package of MueLu in Trilinos and widely used in the community for smoothed aggregation multigrid methods.
ML is implemented in C and known for its good performance properties.
However, the disadvantage is that ML is harder to adapt to new applications and non-standard problems.
Furthermore, ML uses its own internal data structure and is somewhat limited to the use with Epetra objects only.
In contrast, MueLu provides a fully flexible multigrid framework which is designed to be adapted to any kind of new application with non-standard requirements.
Furthermore, it is based on Xpetra and therefore can be used both with Epetra or Tpetra.
Nevertheless, it is an important point to provide some kind of backwards compatibility to allow ML users to easily migrate to MueLu
(or make experiments with MueLu without having to write to much new code).

In this tutorial we present the **MueLu::MLParameterListInterpreter** which provides support for the most important ML parameters to be used with MueLu.

C++ part
========
Preparations
------------
In order to use MueLu (instead or aside of ML) you first have to add it to your application.
Please refer to the MueLu user guide for information about compilation and linking (see [2]_).
Basically, if your application is already working with ML,
you should only need to compile and install MueLu and make sure that the MueLu libraries are found by the linker.

C++ interface
-------------
In the following we assume that the linear operator :math:`A` is available as **RCP<Xpetra::Matrix> A**.

Then we create a parameter list and fill it with ML parameters.

.. literalinclude:: ../../../test/tutorial/MLParameterList.cpp
  :language: cpp
  :start-after: ParameterList begin
  :end-before: ParameterList end

Please refer to the ML guide [1]_ for a complete list of available parameters.

.. note::
	Be aware that the MLParameterListInterpreter does not support all ML parameters,
  but only the most important ones (e.g., smoothers, transfer operators, rebalancing, ...).

Instead of defining the ML parameters by hand in the ParameterList,
you can also read in XML files with ML parameters using

.. literalinclude:: ../../../test/tutorial/MLParameterList.cpp
  :language: cpp
  :start-after: ReadParametersFromXMLFile begin
  :end-before: ReadParametersFromXMLFile end

Next, you create a **MLParameterListInterpreter** object using the parameters
and create a new **MueLu::Hierarchy** from it.

.. literalinclude:: ../../../test/tutorial/MLParameterList.cpp
  :language: cpp
  :start-after: MultigridHierarchy begin
  :end-before: MultigridHierarchy end

Of course, you have to provide all necessary information for the multigrid setup routine.
This does not only include the fine level operator, but also the set of near null space vectors.
Assuming that **numPDEs** stores the number of equations (and near null space vectors),
the following code allows to produce piecewise constant standard near null space vectors
(which should be valid for many PDE discretizations).

.. literalinclude:: ../../../test/tutorial/MLParameterList.cpp
  :language: cpp
  :start-after: BuildDefaultNullSpace begin
  :end-before: BuildDefaultNullSpace end

Then, we just feed in the information to the finest level:

.. literalinclude:: ../../../test/tutorial/MLParameterList.cpp
  :language: cpp
  :start-after: FeedInInformation begin
  :end-before: FeedInInformation end

Finally, we call the **Setup** routine which actually builds the multigrid hierarchy:

.. literalinclude:: ../../../test/tutorial/MLParameterList.cpp
  :language: cpp
  :start-after: CallSetupRoutine begin
  :end-before: CallSetupRoutine end

Once we have the multigrid hierarchy set up,
we can use it the same way as described in :ref:`Iteration Phase <user_api/iteration phase>`.

.. admonition:: Exercise 1

	Study the source code of **../../../test/tutorial/MLParameterList.cpp** and compile it.
  Run the executable **MueLu_tutorial_MLParameterList.exe** with the **--help** command line parameter
  to get an overview of all available command line parameters.

  Run the example on a 1D mesh with 256 elements using

		**./MueLu_tutorial_MLParameterList.exe --xml=xml/ml_ParameterList.xml --nx=256**

	and study the MueLu output.

.. note::
  You will see a warning by the **MLParameterListInterpreter**,
  that the parameter list could not be validated.
  The reason is the follwing:
  Since this tutorial example runs with the Tpetra backend,
  ML, which is purely Epetra-based, cannot validate the parameter list.

.. admonition:: Exercise 2

	Play around with the parameters from **MueLu_tutorial_MLParameterList.exe**.
  Change, e.g., the problem type to a 2D Laplace problem (**--matrixType=Laplace2D**) and adapt the **--nx** and **--ny** parameters accordingly.

Footnotes
=========
.. [1] M. W. Gee, C. M. Siefert, J. J. Hu, R. S. Tuminaro, and M. G. Sala. ML 5.0 Smoothed Aggregation User's Guide, Sandia National Laboratories, 2006, SAND2006-2649
.. [2] L. Berger-Vergiat, C. A. Glusa, G. Harper, J. J. Hu, M. Mayr, P. Ohm, A. Prokopenko, C. M. Siefert, R. S. Tuminaro, and T. A. Wiesner. MueLu User's Guide. Technical Report SAND2023-12265, Sandia National Laboratories, Albuquerque, NM (USA) 87185, 2023.
