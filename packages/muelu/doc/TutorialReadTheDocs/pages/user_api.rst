================================
Using MueLu in User Applications
================================

This tutorial demonstrates how to use MueLu from within user applications in ``C++``.
In \cite[Section 2.6]{Mue} it is explained how to use MueLu through the ``MueLu::CreateE/T/XpetraPreconditioner`` interface.
This interface is designed for beginners which want to try MueLu through standard Trilinos interfaces.

.. note::
   There is also support for Stratimikos.
   Please refer to the ``examples`` in the MueLu folder for more details.

This tutorial aims at more advanced methods to use MueLu
such as creating an explicit instance of some MueLu classes like ``MueLu::Hierarchy`` and ``MueLu::HierarchyManager``.
In the next sections, we give some code snippets.
Most of them are borrowed form the ``laplace2d.cpp`` file in the tutorial.