.. First_Steps documentation master file, created by
   sphinx-quickstart on Wed Mar 18 00:30:46 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

==================
The MueLu Tutorial
==================

.. epigraph::

   This is the MueLu Tutorial.
   Additional resources can be found at the `MueLu package web page <https://trilinos.github.io/muelu.html>`_,
   in the `MueLu User's Guide <https://github.com/trilinos/Trilinos/tree/master/packages/muelu/doc/UsersGuide>`_
   and the `Doxygen source code documentation <https://trilinos.org/docs/dev/packages/muelu/doc/html/index.html>`_.

Preface
=======

The MueLu tutorial is written as a hands-on tutorial for MueLu,
the next generation multigrid framework in Trilinos.
It covers the whole spectrum from absolute beginners' topics to expert level.
Since the focus of this tutorial is on practical and technical aspects
of multigrid methods in general and MueLu in particular,
the reader should already have a basic understanding of multigrid methods and their general underlying concepts.
Please refer to multigrid textbooks for the theoretical background.

Content
=======

The tutorial is split into three parts.
The first part contains four tutorials for beginners who are interested in using multigrid methods.
No knowledge about ``C++`` is required if the programs are used that come with the tutorial (in the Trilinos repository).
If one uses the virtual box image one can even avoid the Trilinos compilation process.
So, the tutorials in the first part can also be used for teaching purposes.
One can easily study the smoothing effect of multigrid smoothers and perform some very basic experiments,
which helps to gain a better understanding of multigrid methods.
In the quick start tutorial, all steps are documented step by step such that it should be very easy to follow the tutorial.
Different exercises may encourage the reader for performing some more experiments and tests.
The following tutorials give an overview of the existing level smoothers and transfer operators,
that can easily be used with the simple XML format, that MueLu uses for defining the multigrid hierarchies.
In addition, it is explained how to visualize the aggregates and export the multigrid levels for a more in-depth analysis.

The second part consists of tutorials for users which are interested in some more background on the underlying techniques that are used in MueLu.
The user still does not need explicit knowledge of ``C++`` or any other programming language,
but some interest in object-oriented design concepts may be helpful to understand the factory concept.
The focus of the second part is on the introduction of the advanced XML interface for MueLu,
which describes all internal building blocks of the multigrid setup procedures with its internal dependencies.
In context of transfer operator smoothing,
a brief introduction of the theory is given with some in-depth details on the algorithmic design in MueLu.
More advanced topics are handled as well such as rebalancing or aggregation strategies.
Additional exercises help the reader to perform some experiments in practice.

The third part is meant for expert users and for users, who want to use MueLu within their own software.
Many detailed ``C++`` examples show how to use MueLu from an user application as preconditioner for a Krylov subspace method or as a standalone multigrid solver.
We expect the reader to be familiar with Trilinos,
especially with the linear algebra packages Epetra and Tpetra as well as the linear solver packages AztecOO or Belos.
For users who are already using ML, the predecessor multigrid package of MueLu in Trilinos,
we provide a chapter describing the migration process from ML to MueLu.

.. toctree::
   :maxdepth: 2
   :caption: Beginners Tutorial

   pages/quick_start.rst
   pages/level_smoothers.rst
   pages/multigrid_for_non_symmetric_problems.rst
   pages/useful_tools_for_analysis.rst
   pages/cd_example.rst

.. toctree::
   :maxdepth: 2
   :caption: Advanced Topics

   pages/xml_interface_for_advanced_users.rst
   pages/muelu_factories_for_transfer_operators.rst
   pages/rebalancing.rst
   pages/advanced_concepts.rst
   pages/aggregation.rst
   pages/useful_commands_and_debugging.rst
   pages/challenge.rst

.. toctree::
   :maxdepth: 3
   :caption: Expert Tutorials

   pages/multigrid_for_multiphysics.rst
   pages/user_api.rst
   pages/ml_parameterlist_interpreter.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
