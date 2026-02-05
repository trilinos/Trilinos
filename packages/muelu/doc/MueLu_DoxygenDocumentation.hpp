// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_DOXYGEN_DOCUMENTATION_HPP
#define MUELU_DOXYGEN_DOCUMENTATION_HPP

/*!

\mainpage

@image html muelu.png

\section muelu_index Table of Contents

- \ref muelu_overview
- \ref muelu_documentation
- \ref muelu_tutorial
- \ref muelu_questions
- \ref muelu_cite
- \ref muelu_developers
- \ref muelu_authors_and_contributors

\section muelu_overview MueLu Overview

MueLu is designed for the parallel solution of large sparse linear systems of equations arising from PDE discretizations. MueLu provides easy-to-use multigrid solvers and preconditioners based on smoothed aggregation algorithms. As a multigrid framework, MueLu supports the design of application-specific multigrid preconditioners.  MueLu is extensible and allows for the research and development of novel multigrid preconditioning methods.  Finally, MueLu is able to run on virtually all modern CPU and accelerator based architectures.

If you're looking for something ready-to-use, MueLu provides scalable multigrid preconditioners for these problem classes:

    Poisson
    Elasticity
    convection-diffusion
    Maxwell’s equations (eddy current formulation)


Features:
    - Easy-to-use interface: MueLu has a user-friendly parameter input deck which covers most important use cases, with reasonable defaults provided for common problem types.
    - Modern object-oriented software architecture: MueLu is written completely in C++ as a modular object-oriented multigrid framework, which provides flexibility to combine and reuse existing components to develop novel multigrid methods.
    - Extensibility: Due to its flexible design, MueLu is an excellent toolkit for research on novel multigrid concepts. Experienced multigrid users have full access to the underlying framework through an advanced XML based interface. Expert users may use and extend the C++ API directly.
    - Integration with Trilinos: As a package of Trilinos, Muelu is well integrated into the Trilinos environment and uses the Tpetra-based sparse linear algebra stack: Krylov solvers (Belos), algebraic solvers (Ifpack2), sparse direct solvers (Amesos2), and load-balancing (Zoltan2).
    - Broad range of supported platforms: MueLu runs on wide variety of architectures, from desktop workstations to parallel Linux clusters and accelerator-based supercomputers.
    - Open source: MueLu is freely available through a simplified BSD license.

If you are interested in using MueLu for your research or want to contribute, please contact the MueLu team through github.


\section muelu_documentation Documentation

Useful views of the Doxygen documentation are:
  - Browsing the entire class documentation. Choose the <a href="classes.html">Class Index subtab</a>.
  - Browsing the class documentation organized by logical groups.  Choose the <a href="modules.html">Modules tab</a>.

The MueLu User's Guide is located in muelu/doc/UsersGuide and at the
<a href=https://trilinos.org/packages/muelu/muelu-documentation>MueLu home page</a>.

\section muelu_tutorial Tutorial

The <a href=https://muelu.github.io>MueLu tutorial</a> describes the most important features of MueLu, from beginning information to advanced usage. It contains accompanying exercises where the interested user can do their own experiments with multigrid parameters. The MueLu tutorial is part of the Trilinos repository, and its examples are automatically tested to ensure that they compile and run.

\section muelu_questions For All Questions and Comments...

    Open an issue: https://github.com/trilinos/Trilinos/issues.
    Start a discussion: https://github.com/orgs/trilinos/discussions.

\section muelu_cite Citing

To cite MueLu, please use the following bibliography entry.

@techreport{MueLu,
title={Mue{L}u User’s Guide},
author={Luc Berger-Vergiat and Christian A. Glusa and Graham Harper and Jonathan J. Hu and Matthias Mayr and Peter Ohm and Andrey Prokopenko and Christopher M. Siefert and Raymond S. Tuminaro and Tobias A. Wiesner}, number={SAND2023-12265}, year={2023}, institution = {Sandia National Laboratories}, }


\section muelu_developers Current Developers

  - Luc Berger-Vergiat, Sandia National Labs
  - Max Firmbach, University of the Bundeswehr Munich
  - Christian Glusa, Sandia National Labs
  - Graham Harper, Sandia National Labs
  - Jonathan Hu, Sandia National Labs
  - Matthias Mayr, University of the Bundeswehr Munich
  - Malachi Phillips, Sandia National Labs
  - Chris Siefert, Sandia National Labs
  - Ray Tuminaro, Sandia National Labs

\section muelu_authors_and_contributors Authors and Contributors

  - Tom Benson, LLNL
  - Emily Furst, University of Washington (summer intern, 2015)
  - Jeremie Gaidamour, INRIA
  - Axel Gerstenberger, Rolls Royce
  - Brian Kelley, Sandia National Labs
  - Andrey Prokopenko, ORNL
  - Paul Tsuji, LLNL
  - Peter Ohm, RIKEN
  - Jerry Watkins, Sandia National Labs
  - Tobias Wiesner, Leica

*/

/* ************************************************************************ */
/* ************************************************************************ */

#endif  // ifndef MUELU_DOXYGEN_DOCUMENTATION_HPP
