//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

/*! 
  \mainpage %NOX: An Object-Oriented Nonlinear Solver Package

  \section intro Introduction
  
  %NOX is short for "Nonlinear Object-Oriented Solutions". Its 
  objective is to enable efficient solution of the equation:
  \f[
  F(x)=0.
  \f]
  Here \f$F:\Re^n \rightarrow \Re^n\f$.

  \section started Getting Started

  See \ref configuration.

  Developers, also see \ref devinfo

  \section classes Class Overview

  See \ref overview.

  \section params Parameter Specifications

  See \ref parameters.

  \section license License

  Our intention is to make NOX available for download under the terms
  of the <a href="http://www.gnu.org/copyleft/lesser.html">GNU Lesser
  General Public License</a>.

  \section download Downloads

  Downloadable versions of NOX will be posted shortly. Until then,
  please contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
  (rppawlo@sandia.gov) for a copy.

  \section bugs Reporting Bugs and Making Enhancement Requests

  To reports bugs or make enhancement reuqests, visit 
  <A HREF="http://software.sandia.gov/bugzilla/">NOX's Bugzilla (Bug Tracking) Database</A>,
  and use the following instructions.
      <UL>
      <LI>Click on "Enter a new bug report"
      <LI>Choose "%NOX"
      <LI>Either login or create a new account
      <LI>Submit your bug report
      </UL>
  </UL>

  \section email Mailing Lists

<p>We recommend that users and developers subscribe to the following mailing lists as appropriate.
<ul>
<li><a href="http://software.sandia.gov/mailman/listinfo/nox-announce">NOX-Announce Mailing List</a> - 
Low-volume, primarily for release announcements. (Highly Recommended)
<li><a href="http://software.sandia.gov/mailman/listinfo/nox-users">NOX-Users Mailing List</a> - 
Discussion forum for users of NOX. (Recommended)
<li><a href="http://software.sandia.gov/mailman/listinfo/nox-developers">NOX-Developers Mailing List</a> - 
Discussion forum for NOX developers. (Developers)
<li><a href="http://software.sandia.gov/mailman/listinfo/nox-checkins">NOX-Checkins Mailing List</a> - 
CVS Commit Message (Active  Developers)
</ul>

\section contributors APPSPACK Contributors

The following people have contributed code to APPSPACK:

<ul>
<li> Russ Hooper, Sandia National Labs, rhoope@sandia.gov
<li> <a href="http://csmr.ca.sandia.gov/~tgkolda/">Tammy Kolda</a>,  
     Sandia National Labs, 
     tgkolda@sandia.gov
<li> Roger Pawlowski, Sandia National Labs, rppawlo@sandia.gov
<li> Eric Phipps, Sandia National Labs, etphipp@sandia.gov
<li> Andy Salinger, Sandia National Labs, agsalin@sandia.gov
<li> Brett Bader (Visiting researcher, Summer 2002)
</ul>



\section questions For All Other Questions and Comments...

Please contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski (rppawlo@sandia.gov).

*/

/*! \page configuration NOX Configuration and Compliation Instructions

\section section1 Basic Setup

To just compile the basic NOX library, type

\verbatim
configure
make
\endverbatim

This only compiles the library \c libnox.a - it does not build any of the
examples. See below for more information on compiling the examples.

\section section2  Generating Documentation

After running configure, type 
\verbatim
make dox
\endverbatim
to generate the documentation.

\section section3 Compiling the Examples

NOX is distributed with a simple set of examples that have no other
external dependencies besides LAPACK and BLAS. These are primarily useful to
allow some basic testing of NOX without requiring an interface to a
more complicated linear algebra system. To enable these examples, type

\verbatim
configure --enable-examples
make
\endverbatim

The compiled examples are located in the \c examples
subdirectory.  To run each example, follow the instructions given below.

<ol>
<li>To run a simple 2D example based on the Rosenbrock function, type
\verbatim
rosenbrock
\endverbatim
</ol>

\section ss1 Building the NOX/Epetra Library

To compile NOX's Epetra library \c libnoxepetra.a, use one
of the following sets of configuration options.

The following will automatically use the environmental variable 
\c TRILINOS_HOME
to find the necessary headers and libraries.
\verbatim 
configure --enable-epetra
\endverbatim


The following will automatically search the directory
\c /home/jdoe/Trilinos to find the necessary headers and
libraries. 
\verbatim
configure --with-trilinos=/home/jdoe/Trilinos
\endverbatim

More specifically, it will look in the following locations
for \c Epetra_Vector.h:
<ul>
<li>directory optionally specified by \c --with-epetra-include
<li>\c /home/jdoe/Trilinos/include
<li>\c /home/jdoe/Trilinos/include/epetra/
<li>\c /home/jdoe/Trilinos/packages/epetra/src
</ul>

Similarily, it will look in the following locations
for \c AztecOO.h:
<ul>
<li>directory optionally specified by \c --with-aztecoo-include
<li>\c /home/jdoe/Trilinos/include
<li>\c /home/jdoe/Trilinos/include/aztecoo/
<li>\c /home/jdoe/Trilinos/packages/aztecoo/src
</ul>
These include files are used in building \c libnoxepetra.a.

Further, if we are also building the Epetra examples, we will look
for epetra, aztecoo, and ifpack libraries. For example, we would
search for \c libepetra.a in the following
directories:
<ul>
<li>directory optionally specified by \c --with-epetra-libdir
<li>\c /home/jdoe/Trilinos/lib
<li>\c /home/jdoe/Trilinos/lib/epetra
<li>\c /home/jdoe/Trilinos/lib/$TRILINOS_TARGET
</ul>
In the above, if \c TRILINOS_TARGET is not defined as an
environmental variable, configure will define it.
</ul>

\section ss2 Building the NOX/Epetra Examples

To compile the examples for NOX/Epetra, do the following:

\verbatim
configure --enable-epetra --enable-epetra-examples.
\endverbatim

\section section5 Using MPI with NOX/Epetra

NOX itself does not call any MPI commands, but some of the NOX/Epetra
examples do make use of MPI. So, we include some configuration options
to support MPI.

Most MPI
packages come with a wrapper to the compiler that specifies the
correct flags, include files, libraries, and so on to compile with
MPI. This is often the easiest way to compile with MPI.
<ul>
<li>  Use the default MPI C++ compiler (\c mpiCC). 
\verbatim
configure --with-mpi-cxx 
\endverbatim
<li> Specify the exact name and location of the MPI C++ compiler. This is
useful if you have multiple installations of MPI or if \c mpiCC
is not in your default path.
\verbatim
configure --with-mpi-cxx=/home/jdoe/mpich-1.2.3/bin/mpiCC
\endverbatim
</ul>

An alternative to using \c mpiCC is to specify the
include directory, library directory, and libraries specifically,
There are several flags that you can use to do this. These are
generally incompatible with \c --with-mpi-cxx
<ul>
<li> This specifies the general location of the MPI installation. We
then assume that the headers are located in
\c /home/jdoe/lam/include, that the libraries are
located in \c /home/jdoe/lam/lib, and that the
library is specifed by \c -lmpi.
\verbatim
configure --with-mpi=/home/jdoe/lam
\endverbatim

<li> This is the same as above except that the library to be linked with
is explicitly specified.
\verbatim
configure --with-mpi=/home/jdoe/mpich --with-mpi-libs="-lmpich"
\endverbatim

<li>
Specify the exact locations of the MPI include files and the MPI
libraries. Assumes that the
library is specifed by \c -lmpi.
\verbatim
configure --with-mpi-include="/home/jdoe/lam/include" \
          --with-mpi-libdir="/home/jdoe/lam/lib"
\endverbatim

<li>
This is the same as above except that the library to be linked with
is explicitly specified.
\verbatim
configure --with-mpi-include="/home/jdoe/mpich/include" \
          --with-mpi-libdir="/home/jdoe/mpich/lib" \
          --with-mpi-libs="-lmpich"
\endverbatim

<li> 
Explicitly specify the MPI library without specifying the library directory.

\verbatim
configure --with-mpi-include="/home/jdoe/mpich/include" \
          --with-mpi-libs="/home/jdoe/mpich/lib/libmpich.a"
\endverbatim

</ul>

\section petsclub Building the NOX/PETSc Library</h5>

To compile the NOX/PETSc library \c libnoxpetsc.a, use 
the following configuration option.

This will automatically use the environmental variables 
\c PETSC_DIR and \c PETSC_ARCH 
to find the necessary headers.
\verbatim
configure --enable-petsc
\endverbatim

For example, to build \c libnoxpetsc.a using a pre-built PETSc
version 2.1.3 installed in the user's home directory on a linux machine:
\verbatim
setenv PETSC_DIR ${HOME}/petsc-2.1.3
setenv PETSC_ARCH linux
configure --enable-petsc
\endverbatim

\section petscexamples Building the NOX/PETSc Examples</h5>

To compile the examples for NOX/PETSc, an additional environment variable,
\c PETSC_BOPT, must be defined and set to the value of BOPT used to build the 
user's PETSc installation.  For example, for a debug-enabled C++ PETSc build, 
BOPT=g_c++ so that PETSC_BOPT=g_c++.  The configure script then allows
the examples to be built using the following:

\verbatim
configure --enable-petsc --enable-petsc-examples.
\endverbatim


Depending on the user's PETSc build, dynamically linked libraries may be used.
If an error occurs when trying to run the example executable, make sure the
environment variable \c LD_LIBRARY_PATH points to the PETSc library directory.
From the example above, setting the path would likely require:
\verbatim
setenv LD_LIBRARY_PATH ${PETSC_DIR}/lib/lib${PETSC_BOPT}/${PETSC_ARCH}
\endverbatim


*/

/*! \page devinfo NOX Developer's Information

\section cvsrepos CVS Repository for NOX

The CVS repository is located on <a
href="http://software.sandia.gov/>software.sandia.gov</a>. See Paul
Sery (pgsery@sandia.gov) or Tammy Kolda (tgkolda@sandia.gov) for an
account.

You should not do any development of software.sandia.gov. Instead,
check out a copy to your local machine as follows.

First, be sure to correctly set your \c CVS_SSH environment
variable. For example, in (t)csh, you would type:

\verbatim
setenv CVS_RSH ssh
\endverbatim

Next, to check out nox, type
\verbatim
cvs -d :ext:username@hopper.ca.sandia.gov:/space/CVS checkout -P nox
\endverbatim

To update your copy of nox to the latest version, type
\verbatim
cvs update -d -P
\endverbatim

*/

/*! \page overview NOX Class Overview

  %NOX's flexibility is based on the fact that it does not depend on
  any particular linear algebra package. In order to interface to %NOX,
  the user needs to supply methods that derive from the following
  abstract classes:

  - NOX::Abstract::Vector
  - NOX::Abstract::Group

  The \b Vector supports basic vector operations such as dot products
  and so on. The \b Group supports the linear algebra functionality as
  well as the interface to evaluate the function and,
  optionally, the Jacobian. Complete details are provided in
  the class decription.

  Included with %NOX are three such interfaces:

  - NOX::Example::Vector and NOX::Example::Group
  - NOX::Epetra::Vector and NOX::Epetra::Group
  - NOX::Petsc::Vector and NOX::Petsc::Group

  The \b Example interface is a very simple interface that is not
  intended for actual computations, but to serve as an
  easy-to-understand example of how one might interface to %NOX.
  
  The \b Epetra interface is an interface to the
  Trilinos/Epetra library.

  The \b PETSc interface is an interface with the PETSc library. PETSc
  was developed at Argonne National Labs; for more information see the
  <a href="http://www-fp.mcs.anl.gov/petsc">PETSc Home Page</a>

  All solvers are in the NOX::Solver namespace. The solvers are
  accessed via the NOX::Solver::Manager. The recommended solver is
  NOX::Solver::LineSearchBased, which is a basic nonlinear solver
  based on a line search. Each solver has a number of options that can
  be specified, as documented in each class or on the \ref parameters.

  The search directions are in the NOX::Direction namespace and
  accessed via the NOX::Direction::Manager. The default search
  direction for a line-search based method in NOX::Direction::Newton.

  Several line searches are available, as defined in the
  NOX::LineSearch, and accessed via the NOX::LineSearch::Manager
  class. Examples include

  - NOX::LineSearch::FullStep
  - NOX::LineSearch::Backtrack
  - NOX::LineSearch::MoreThuente

  Convergence or failure of a given solver method is determined by 
  the status tests defined in the NOX::StatusTest namespace. Various
  status tests may be combined via the NOX::StatusTest::Combo
  object. Users are free to create additional status tests that derive
  from the NOX::StatusTest::Generic class.


 */
 
/*! \page parameters NOX Parameter Reference Page

  All parameters are specified via parameter lists as defined by the
  NOX::Parameter::List class.

  Here is the basic hierarchical structure of the top level list that
  is passed to a NOX::Solver::Manager object constructor reset call.

<ul>
<li> NOX::Solver::Manager - \e Nonlinear \e Solver
<li> NOX::Solver::LineSearchBased, NOX::Solver::TrustRegionBased - \b Direction
     <ul>
     <li> NOX::Direction::Manager - \e Method
     <li> NOX::Direction::Newton - \b Linear \b Solver
          <ul>
	  <li> NOX::Direction::Newton - Tolerance
	  </ul>
     <li> NOX::Direction::Newton - \e Forcing \e Term \e Method
     <li> NOX::Direction::Newton - \e Forcing \e Term \e Initial \e Tolerance
     <li> NOX::Direction::Newton - \e Forcing \e Term \e Minimum \e Tolerance
     <li> NOX::Direction::Newton - \e Forcing \e Term \e Maximum \e Tolerance
     <li> NOX::Direction::Newton - \e Forcing \e Term \e Alpha
     <li> NOX::Direction::Newton - \e Forcing \e Term \e Gamma
     <li> NOX::Direction::SteepestDescent - \e Scaling \e Type
     <li> NOX::Direction::NonlinearCG - \e Orthogonalize
     <li> NOX::Direction::NonlinearCG - \e Precondition
     <li> NOX::Direction::NonlinearCG - \e Restart Frequency
     </ul>
<li> NOX::Solver::LineSearchBased - \b Line \b Search
     <ul>
     <li> NOX::LineSearch::Manager - \e Method
     <li> NOX::LineSearch::FullStep - \e Full \e Step
     <li> NOX::LineSearch::Backtrack, NOX::LineSearch::Polynomial, NOX::LineSearch::Quadratic, NOX::LineSearch::MoreThuente, NOX::LineSearch::Secant - \e Default \e Step
     <li> NOX::LineSearch::Backtrack, NOX::LineSearch::Polynomial, NOX::LineSearch::Quadratic, NOX::LineSearch::MoreThuente, NOX::LineSearch::Secant   - \e Minimum \e Step
     <li> NOX::LineSearch::Backtrack, NOX::LineSearch::Polynomial, NOX::LineSearch::Quadratic, NOX::LineSearch::MoreThuente, NOX::LineSearch::Secant   - \e Recovery \e Step
     <li> NOX::LineSearch::Backtrack, NOX::LineSearch::Polynomial, NOX::LineSearch::Quadratic, NOX::LineSearch::MoreThuente, NOX::LineSearch::Secant  - \e Max \e Iters
     <li> NOX::LineSearch::Backtrack - \e Decrease \e Condition
     <li> NOX::LineSearch::Polynomial, NOX::LineSearch::Quadratic - \e Convergence \e Criteria
     <li> NOX::LineSearch::Polynomial, NOX::LineSearch::Quadratic - \e Alpha \e Factor
     <li> NOX::LineSearch::Polynomial, NOX::LineSearch::Quadratic - \e Min \e Bounds \e Factor
     <li> NOX::LineSearch::Polynomial, NOX::LineSearch::Quadratic - \e Max \e Bounds \e Factor
     <li> NOX::LineSearch::Polynomial, NOX::LineSearch::Quadratic - \b Output
          <ul>
	  <li> NOX::LineSearch::Polynomial, NOX::LineSearch::Quadratic - Total Number of Line Search Calls
	  <li> NOX::LineSearch::Polynomial, NOX::LineSearch::Quadratic - Total Number of Non-trivial Line Searches
	  <li> NOX::LineSearch::Polynomial, NOX::LineSearch::Quadratic - Total Number of Failed Line Searches
	  <li> NOX::LineSearch::Polynomial, NOX::LineSearch::Quadratic - Total Number of Line Search Inner Iterations
	  </ul>
     <li> NOX::LineSearch::MoreThuente - \e Sufficient \e Decrease
     <li> NOX::LineSearch::MoreThuente - \e Curvature \e Condition
     <li> NOX::LineSearch::MoreThuente - \e Interval \e Width
     <li> NOX::LineSearch::MoreThuente - \e Maximum \e Step
     </ul>
<li> NOX::Solver::LineSearchBased, NOX::Solver::TrustRegionBased  - \b Output
     <ul> 
     <li> NOX::Solver::LineSearchBased - Nonlinear Iterations
     <li> NOX::Solver::LineSearchBased - 2-Norm or Residual
     </ul>
<li> NOX::Solver::TrustRegionBased - \b Cauchy \b Direction 
     <ul>
     <li> (see \b Direction above)
     </ul>
<li> NOX::Solver::TrustRegionBased - \e Minimum \e Trust \e Region \e Radius
<li> NOX::Solver::TrustRegionBased - \e Maximum \e Trust \e Region \e Radius
<li> NOX::Solver::TrustRegionBased - \e Minimum \e Improvement \e Ratio
<li> NOX::Solver::TrustRegionBased - \e Contraction \e Trigger \e Ratio
<li> NOX::Solver::TrustRegionBased - \e Contraction \e Factor
<li> NOX::Solver::TrustRegionBased - \e Expansion \e Trigger \e Ratio
<li> NOX::Solver::TrustRegionBased - \e Expansion \e Factor
<li> NOX::Solver::TrustRegionBased - \e Recovery \e Step
<li> NOX::Utils - \e MyPID
<li> NOX::Utils - \e Output \e Information
<li> NOX::Utils - \e Output \e Processor
<li> NOX::Utils - \e Output \e Precision
     
</ul>



*/
