//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! 
  \mainpage %NOX & %LOCA: Object-Oriented Nonlinear Solver and Continuation Packages

  \section introduction Introduction
  
  %NOX is short for <em>Nonlinear Object-Oriented Solutions</em>, and its 
  objective is to enable the efficient solution of the equation:
  \f$
  F(x)=0
  \f$,
  where \f$F:\Re^n \rightarrow \Re^n\f$.  
  %NOX is designed to
  work with any linear algebra package and to be easily customized.
  %NOX is part of Sandia's <A
  HREF="http://software.sandia.gov/trilinos/">Trilinos project</a>.

  %LOCA, distributed as part of %NOX,
  is short for <em>Library of Continuation Algorithms</em>, and its
  objective is to compute families of solutions to
  \f$
  F(x,p)=0
  \f$
  and their bifurcations, where
  \f$F:\Re^n\times\Re^m\rightarrow\Re^n\f$.  For %LOCA specific information, 
  see \ref loca_overview.

  \section user_information User Information

  See \ref nox_user_information

  \section developer_information Developer Information

  See \ref nox_developer_information

  \section license License

  %NOX is available for download under the terms
  of the <a href="http://www.gnu.org/copyleft/lesser.html">GNU Lesser
  General Public License</a>.

  \section downloads Downloads

  You may download %NOX (which includes %LOCA) in one of two ways.

  %NOX only:

  <ul>
  <li> 
  <a href="http://software.sandia.gov/nox/downloads/nox-3.1a-dev.tar.gz">Latest Development Version</a>, updated nightly
  </ul>

  Trilinos (includes %NOX as well as other packages such as Epetra, Ifpack, and more):

  <ul>
  <li> Download site to be announced soon!
  </ul>


  \section bugs Reporting Bugs and Making Enhancement Requests

  To reports bugs or make enhancement requests, visit 
  <A HREF="http://software.sandia.gov/bugzilla/">NOX's Bugzilla (Bug Tracking) Database</A>,
  and use the following instructions.
      <UL>
      <LI>Click on "Enter a new bug report"
      <LI>Choose "NOX"
      <LI>Either login or create a new account
      <LI>Submit your bug report
      </UL>

  \section email Mailing Lists

  <p>We recommend that users and developers subscribe to the following mailing lists as appropriate.
  <ul>
  <li><a href="http://software.sandia.gov/mailman/listinfo/nox-announce">NOX-Announce Mailing List</a> - 
  Low-volume, primarily for release announcements. (Highly Recommended)
  <li><a href="http://software.sandia.gov/mailman/listinfo/nox-users">NOX-Users Mailing List</a> - 
  Discussion forum for users of %NOX. (Recommended)
  <li><a href="http://software.sandia.gov/mailman/listinfo/nox-developers">NOX-Developers Mailing List</a> - 
  Discussion forum for %NOX developers. (Developers)
  <li><a href="http://software.sandia.gov/mailman/listinfo/nox-checkins">NOX-Checkins Mailing List</a> - 
  CVS Commit Messages (Active  Developers)
  </ul>
  
  \section contributors NOX & LOCA Contributors
  
  The following people have contributed code to %NOX and/or %LOCA:
  
  <ul>
  <li> Brett Bader, Sandia National Labs, bwbader@sandia.gov
  <li> Russ Hooper, Sandia National Labs, rhoope@sandia.gov
  <li> <a href="http://csmr.ca.sandia.gov/~tgkolda/">Tammy Kolda</a>,  
  Sandia National Labs, 
  tgkolda@sandia.gov (%NOX Project Lead)
  <li> Roger Pawlowski, Sandia National Labs, rppawlo@sandia.gov (%NOX Project Lead)
  <li> Eric Phipps, Sandia National Labs, etphipp@sandia.gov
  <li> Andy Salinger, Sandia National Labs, agsalin@sandia.gov (LOCA Project Lead)
  </ul>
  
  \section questions For All Other Questions and Comments...
  
  Please contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski (rppawlo@sandia.gov).
  
*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page nox_user_information NOX User Information

\section overview Overview

To get started with %NOX, you'll want to follow these four steps.
<ul>
<li> \ref step1
<li> \ref step2
<li> \ref step3
<li> \ref step4
</ul>

\section documentation  Generating Documentation

The documentation is available online at 
<a href="http://software.sandia.gov/nox/">http://software.sandia.gov/nox/</a>.

To generate a local copy of the documentation, you'll need to have doxygen installed. 
Then follow the instructions in \c nox/doc/README.
If you've already successfully run configure, type 
\verbatim
make dox
\endverbatim
to generate the documentation. The main page is \c nox/doc/html/index.html.

\section Quicklinks

- \ref nox_configuration_options
- \ref parameters
- \ref downloads
- \ref bugs

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*!
\page step1 Step 1: Download, Configure, Compile, and Install NOX

\section step1_download Downloading and Unpacking NOX

If you have not already done so, you'll need to download %NOX; see
\ref downloads. You can download %NOX as a standalone package or as
part of Trilinos. 

\note If you download it as part of Trilinos, you should
follow the Trilinos instructions for configuring, compiling and
installing.

To unpack %NOX, type

\verbatim
gunzip -c nox-\<version\>.tar.gz | tar -xzvf -
\endverbatim

This will create a directory named nox-\<version\> containing the source
code for %NOX. Rename this direction to nox; i.e.,

\verbatim
mv nox-\<version\> nox
\endverbatim

\section step1_configure Configuring NOX

Change directories to the nox directory:

\verbatim
cd nox
\endverbatim

Run the configure script, which will create a Makefile customized to
your machine.

\verbatim
./configure
\endverbatim

There are many options that can be passed to the configure script, and
they will be discussed as appropriate. For a full listing, see \ref
nox_configuration_options.

\note If you have previously run the configure script and compiled the
code, be sure to type "make clean" before re-running the configure
script.

\section step1_compile Compiling NOX

To compile %NOX, simply type:
\verbatim
make
\endverbatim

\section step1_install Installing NOX

To install %NOX, type
\verbatim
make install
\endverbatim

By default, %NOX will be installed in /usr/local. If you prefer that
it be installed in a different directory, that can be specified using
the --prefix=PREFIX configuration option.

\section step1_problems If you have problems...

If you have problems configuring and compiling %NOX, please report
them using
<A HREF="http://software.sandia.gov/bugzilla/">Bugzilla</A>; see \ref bugs for more information.

\section step1_moving_on Moving on...

Go on to \ref step2.

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*!
\page step2 Step 2: Create concrete implementations of the NOX::Abstract classes

\section step2_abstract NOX abstract classes

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
  the class descriptions.

  In order to link your code to NOX, you have to write your own
  instantiation of the NOX::Abstract::Vector and NOX::Abstract::Group
  classes (see \ref step2_writing) or use one of the predefined
  implementations (see \ref step2_instantiations).

\section step2_writing Writing your own instantiation

  We recommend using the NOX::LAPACK framework as a guide.  In this
  case, the underlying vectors of C++ STL vector<double> objects, and
  the matrices are stored in our own NOX::LAPACK::Matrix class.  Both
  the vectors and matrices are manipulated using LAPACK.  The
  NOX::LAPACK::Vector is straightforward; see NOX_LAPACK_Vector.H and
  NOX_LAPACK_Vector.C.  The NOX::LAPACK::Group uses the
  NOX::LAPACK::Vector and NOX::LAPACK::Matrix objects. The interface
  with the application is handled in a separate class
  (NOX::LAPACK::Interface) that is passed to the NOX::LAPACK::Group
  when it is constructed.

  \note The NOX header files should be installed in /usr/local/include/nox.

\section step2_instantiations Instantiations provided with NOX

  %NOX includes three ready-made instantiations.

  <ul>
  <li>
  The NOX::LAPACK instantiation is an interface to the BLAS/LAPACK library.
  It is not intended for large-scale computations, but to serve as an
  easy-to-understand example of how one might interface to %NOX.  
  To compile the NOX::LAPACK library and examples, use the
  --enable-nox-lapack and --enable-nox-lapack-examples options to
  configure; see \ref nox_configuration_options.  For instructions on 
  interfacing your code to the lapack implementation see \ref lapack_interface.
  
  <li>
  The NOX::Epetra instantiation is an interface to the Trilinos/Epetra
  library developed by 
  <a href="http://www.sandia.gov">Sandia National Labs</a>; 
  for more information see the 
  <a href="http://software.sandia.gov/Trilinos">Trilinos Home Page</a>. 
  If you download %NOX with Trilinos (see \ref downloads),
  then Epetra is included in the distribution.
  To compile the NOX::Epetra library and examples, use the
  --enable-nox-epetra and --enable-nox-epetra-examples options to
  configure; see \ref nox_configuration_options.  For instructions on 
  interfacing your code to the epetra implementation see \ref epetra_interface.

  <li>
  The NOX::PETSc instantiation is an interface with the PETSc library. PETSc
  was developed at 
  <a href="http://www.anl.gov">Argonne National Labs</a>; 
  for more information see the
  <a href="http://www-fp.mcs.anl.gov/petsc">PETSc Home Page</a>
  To compile the NOX::Petsc library and examples, use the
  --enable-nox-petsc and --enable-nox-petsc-examples options to
  configure; see \ref nox_configuration_options.  For instructions on 
  interfacing your code to the petsc implementation see \ref petsc_interface.
  </ul>


\section step2_problems If you have problems...

If you have problems configuring and compiling %NOX, please report
them using
<A HREF="http://software.sandia.gov/bugzilla/">Bugzilla</A>; see \ref bugs for more information.

\section step2_moving_on Moving on...

Go on to \ref step3.

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*!
\page step3 Step 3: Call NOX from your code

A small example based on the NOX::LAPACK Vector and Group is included
with %NOX; see Rosenbrock.C.

There are four steps to calling the solver.

\section step3_stepa Step A: Set up a Group object containing the initial guess.

The first step is to somehow construct an object deriving from the
Abstract::Group which contains the initial guess as its x-vector. See
\ref step2 for more details.

\code

// Construct an object that derives from the abstract Group object. Here
// we call it the ExampleGroup. It should contain the initial guess.
ExampleGroup grp;

\endcode

\section step3_stepb Step B: Construct the status tests

The NOX::StatusTest objects are used to test for convergence or
failure. For example...

- The NOX::StatusTest::MaxIters test is used to control the maximum
  number of iterations that the nonlinear solver is allowed to take.

- The NOX::StatusTest::NormF test checks the size of \f$\| F(x)
  \|\f$. This check can be done based on relative or absolute norms,
  scaled by \f$ \sqrt{n} \f$, and more. See the description for more
  details.

These and other status tests can be combined using the
NOX::StatusTest::Combo object. It takes an arbitrary number of status
tests and either AND's or OR's them together. 

Finally, the user can create their own status tests, so long as they
derive from NOX::StatusTest::Generic.

Here is some sample code.

\code
  // Set up the status tests
  NOX::StatusTest::NormF statusTestA(grp, 1.0e-4);
  NOX::StatusTest::MaxIters statusTestB(20);
  NOX::StatusTest::Combo statusTestsCombo(NOX::StatusTest::Combo::OR, statusTestA, statusTestB);
\endcode

\section step3_stepc Step C: Set the solver parameters

Next, we set up a parameter list containing the information on what
types of solvers and so forth should be used. 

\code
// Create the list of solver parameters
NOX::Parameter::List solverParameters;

// Set the solver (this is the default)
solverParameters.setParameter("Nonlinear Solver", "Line Search Based");

// Create the line search parameters sublist
NOX::Parameter::List& lineSearchParameters = solverParameters.sublist("Line Search");

// Set the line search method
lineSearchParameters.setParameter("Method","More'-Thuente");
\endcode

For a full list of parameters; see \ref parameters.

\section step3_stepd Step D: Construct the solver and solve

The last step is to create the solver, passing in the group with the
initial guess, the status tests, and the solver parameters.

\code
// Create the solver
NOX::Solver::Manager solver(grp, statusTestsCombo, solverParameters);

// Solve the nonlinear system
NOX::StatusTest::StatusType status = solver.solve();

// Print the answer
cout << "\n" << "-- Parameter List From Solver --" << "\n";
solver.getParameterList().print(cout);

// Get the answer
grp = solver.getSolutionGroup();

// Print the answer
cout << "\n" << "-- Final Solution From Solver --" << "\n";
grp.print();
\endcode

\section step3_problems If you have problems...

If you have problems configuring and compiling %NOX, please report
them using
<A HREF="http://software.sandia.gov/bugzilla/">Bugzilla</A>; see \ref bugs for more information.

\section step3_moving_on Moving on...

Go on to \ref step4.

**/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*!
\page step4 Step 4: Link your code to NOX

\section step4_basics The Basics

To compile your code, your Makefile should look something like the
following, assuming %NOX was installed in /usr/local.

\verbatim
NOX_PREFIX = /usr/local
LDFLAGS = -L$(NOX_PREFIX)/lib 
LIBS = -lnox
CXXFLAGS = -I$(NOX_PREFIX)/include/
\endverbatim

\section step4_lapack Compiling with LAPACK

\verbatim
NOX_PREFIX = /usr/local
LDFLAGS = -L$(NOX_PREFIX)/lib 
LIBS = -lnox -lnoxlapack -llapack -lblas 
CXXFLAGS = -I$(NOX_PREFIX)/include/
\endverbatim

\section step4_epetra Compiling with Epetra

We assume that Epetra, Aztec, and Ifpack are installed in the same place as %NOX.

\verbatim
NOX_PREFIX = /usr/local
LDFLAGS = -L$(NOX_PREFIX)/lib 
LIBS = -lnox -lnoxepetra -laztecoo -lifpack -lepetra  -llapack -lblas 
CXXFLAGS = -I$(NOX_PREFIX)/include/
\endverbatim

\section step4_problems If you have problems...

If you have problems configuring and compiling %NOX, please report
them using
<A HREF="http://software.sandia.gov/bugzilla/">Bugzilla</A>; see \ref bugs for more information.

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page nox_developer_information NOX Developer Information

To become a %NOX developer, contact Tammy Kolda (tgkolda@sandia.gov) or
Roger Pawlowski (rppawlo@sandia.gov).

The following pages are relevant to developers.
<ul> 
<li>\ref cvsrepos
<li>\ref coding
<li>\ref prerelease
</ul>


*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page nox_configuration_options NOX Configuration Options 

Here is a detailed list and explanation of the %NOX configuration options.

\section basic_config Basic Options

<ul>
  <li> \c --enable-debug          
    <ul>
      <li> This turns on compiler debugger flags. It has not been fully
           tested. As an alternate, specify CXXFLAGS on the configure line.
    </ul>
  <li> \c --enable-opt            
    <ul>
      <li> This turns on compiler optimization flags. It has not been fully
           tested. As an alternate, specify CXXFLAGS on the configure line.
    </ul>
  <li> \c --with-cppflags
    <ul>
      <li> Specify additional preprocessor flags (e.g., "-Dflag -Idir")
    </ul>
  <li> \c --with-cxxflags
    <ul>
      <li> Specify additional C++ flags
    </ul>
  <li> \c --with-ldflags
    <ul>
      <li> Specify additional linker flags (e.g., "-Ldir")
    </ul>
  <li> \c --with-ar
    <ul>
      <li> Specify a special archiver command, the default is "ar cru".
    </ul>


</ul>

\section influential_env Influential Environmental Variables

<ul>
   <li> \c CXX         
      <ul>
         <li> C++ compiler command
      </ul>
   <li> \c CXXFLAGS    
      <ul> 
         <li> C++ compiler flags 
      </ul>
   <li> \c LDFLAGS
      <ul>
         <li> Specify linker flags
      </ul>
   <li> \c CPPFLAGS
      <ul>
         <li> C/C++ preprocessor flags
      </ul>
   <li> \c CXXCPP  
      <ul>
         <li> C++ preprocessor
      </ul>
   <li> \c F77     
      <ul>
         <li> Fortran 77 compiler command <em>(used only in configure tests)</em>
      </ul>
   <li> \c FFLAGS  
      <ul>
         <li> Fortran 77 compiler flags <em>(used only in configure tests)</em>
      </ul>
</ul>

\section mpi_config MPI-Related Options

The %NOX library itself has no MPI dependencies. However, the %NOX-Epetra
library and related examples do depend on MPI. 

<ul>
  <li> \c --enable-mpi
    <ul>
      <li> Enables MPI mode. Defines HAVE_MPI in the NOX_Config.h
           file. Will test for the ability to preprocess the MPI header file
           and may test ability to link with MPI. Note that this does not
	   automatically use the MPI compilers.
    </ul>
  <li> \c --with-mpi-compilers[=PATH]
    <ul>
      <li> Use the MPI compilers (mpicxx or mpiCC, mpicc, and
      mpif77). If the optional PATH is specified, then it uses those
      compilers that it finds in the specified directories. For the
      C++ compiler, it looks for mpicxx first. If it doesn't find
      that, it defaults to mpiCC. Alternatively, each compiler may be
      specified precisely using the CXX, CC, and F77 environmental
      variables.
    </ul>
  <li> \c --with-mpi=MPIROOT
    <ul>
      <li> Specify the MPI root directory. Then
           <tt>MPIROOT/include</tt> is assumed to be the include
           directory and and <tt>MPIROOT/lib</tt> is assumed to be the
           lib directory. Automatically enables MPI mode.
    </ul>
  <li> \c --with-mpi-libs="LIBS"
    <ul>
      <li> Specify the MPI libraries. 
           If this is not specified, defaults to \c "-lmpi" if either
           \c --with-mpi or \c --with-mpi-libdir is specified.
    </ul>
  <li> \c --with-mpi-incdir=DIR
    <ul>
      <li> Specify the MPI include files location. 
           Defaults to <tt>MPIROOT/include</tt> if <tt>--with-mpi</tt> is specified. 
	   If multiple directories must be specified, 
	   try <tt>--with-cppflags="-I\<dir1\> -I\<dir2\>"</tt> instead.
    </ul>
  <li> \c --with-mpi-libdir=DIR
    <ul>
      <li> Specify the MPI libraries location. 
           Defaults to <tt>MPIROOT/lib</tt> if <tt>--with-mpi</tt> is specified.
	   If multiple directories must be specified, 
	   try <tt>--with-ldflags="-L\<dir1\>  -L\<dir2\></tt>" instead.
    </ul>



</ul>

\section lapack_options LAPACK-Related Options

<ul>
  <li> \c --enable-nox-lapack         
    <ul>
      <li> Compile noxlapack library (libnoxlapack.a). Default is not to compile the library.
    </ul>
  <li> \c --enable-nox-lapack-examples
    <ul>
      <li> Compile lapack examples. 
    </ul>
</ul>

\section epetra_options Epetra-Related Options

<ul>
  <li> \c --enable-nox-epetra         
    <ul>
      <li> Compile noxepetra library (libnoxepetra.a). Default is not to compile the library.
    </ul>
  <li> \c --enable-nox-epetra-examples
    <ul>
      <li> Compile epetra examples. 
    </ul>
</ul>

\section petsc_options PETSc-Related Options

<ul>
  <li> \c --enable-nox-petsc      
    <ul>
      <li> Compile noxpetsc library (libnoxpetsc.a). Default is not to compile the library.
      <li> Requires existing PETSc source specified by environment variables, 
           \c PETSC_DIR and \c PETSC_ARCH
    </ul>
  <li> \c --enable-nox-petsc-examples
    <ul>
      <li> Compile petsc examples. 
           Requires \c --enable-nox-petsc.
      <li> Requires existing PETSc libraries specified by environment variables,
           \c PETSC_DIR,  \c PETSC_ARCH and \c PETSC_BOPT
    </ul>

  See \ref petsc_interface for more information regarding building 
       %NOX with PETSc
</ul>

\section loca_config LOCA-Related Options

<ul>
  <li> \c --enable-loca
    <ul>
      <li> Compile the loca library (libloca.a). Default is not to compile the library.
           Also compiles appropriate interfaces if they have been enabled (e.g., lapack and epetra)
    </ul>
</ul>

\section developer_config Developer-Related Options

<ul>
<li>  \c --enable-maintainer-mode 
  <ul>
    <li> Enable make rules and dependencies not useful
         (and sometimes confusing) to the casual installer
  </ul>
<li>  \c --enable-prerelease
  <ul>
    <li> Cannot be used with the distribution version of the code - only works withe CVS version.
         Compiles code that is not included in the releases but is part of the CVS repository. 
         See \ref prerelease.
  </ul>
</ul>

*/

/*! \page prerelease NOX Developer's Guide to Prerelease Code


We define <em>prerelease code</em> to be code that should not be
distributed, documented, or tested as part of the automatic nightly
scripts. 

Prerelease code will only be used if configure is executed with the \c
--enable-prerelease option. In that case, 
<ul>
<li>all prerelease code will be compiled (\c -DWITH_PRERELEASE is added to \c CPPFLAGS), 
<li>doxygen will document prerelease code (\c WITH_PRERELEASE is added to the \c PREDEFINED variable in Doxyfile),
<li>\c BUILD_PRERELEASE is true for automake/autoconf Makefile generation.
</ul>

There are two steps to adding prerelease code.

<ol>

<li> Add ifdef's to the source and header files. 

In other words, the code should be surrounded by a \c WITH_PRERELEASE
ifdef as follows.

\verbatim
#ifdef WITH_PRERELEASE
// Insert All Code Here
#endif
\endverbatim

<li> Modify Makefile.am.

For example, the following modifications were made to \c
nox/src/Makefile.am to add NOX_Direction_QuasiNewton as prerelease code
in the libnox.a library.

\verbatim
if BUILD_PRERELEASE

nodist_pkginclude_HEADERS = \
	$(srcdir)/NOX_Direction_QuasiNewton.H 

nodist_libnox_a_SOURCES = \
	$(srcdir)/NOX_Direction_QuasiNewton.C

endif
\endverbatim

</ol>

Remember to run configure with the \c --enable-prerelease option!!

\note
Prerelease code must still be able to compile.

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page cvsrepos NOX Developer's Guide to the CVS Repository for NOX

The CVS repository is located on 
<A HREF="http://software.sandia.gov/>software.sandia.gov</A>. See Paul
Sery (pgsery@sandia.gov) or Tammy Kolda (tgkolda@sandia.gov) for an
account.

You should not do any development of software.sandia.gov. Instead,
check out a copy to your local machine.  Be sure to correctly set your
\c CVS_RSH environment variable to \c ssh.

<ul>
<li>
To check out nox:
\verbatim
cvs -d :ext:username@software.sandia.gov:/space/CVS checkout -P nox
\endverbatim

<li>
To update your copy of nox to the latest version:
\verbatim
cvs update -d -P
\endverbatim
</ul>

We have installed a <a
href="http://software.sandia.gov/bonsai/cvsqueryform.cgi?cvsroot=/space/CVS&module=nox">web
front end for browsing and querying the CVS repository</a>. 

<ul>
<li> \b Browsing: From the
"CVS Query Form" select the "Trilinos Tree" under "Browse" in the
yellow box on the right hand side. Next, choose "Trilinos", then
"packages", and finally "nox".
<li> \b Querying: Select "nox" as the "Module"
</ul>

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page coding NOX Developer's Coding Guidelines

\htmlonly
Adapted largely from <a
href="http://www.doc.ic.ac.uk/lab/cplus/c++.rules/">Programming in
C++, Rules and Recommendations, by Mats Henricson and Erik
Nyquist</a>.


<H3>Structure of the Files</H3>

<H5>Naming Conventions</H5>
<ul>
<li>C++ header files end in <code>.H</code> and source files end in
<code>.C</code>
<br>&nbsp;
<li>The name of the files should correspond to the name of the class
they define, with double-colons replaced by underscores. For example,
the definition of the class <code>NOX::Abstract::Group</code> is in
the file <code>NOX_Abstract_Group.H</code>.
<br>&nbsp;
</ul>

<H5>General File Structure</H5>
<ul>
<li>Each file should begin as follows:

<pre>
// $Id$
// $Source$
//@HEADER
# ************************************************************************
# 
#            NOX: An Object-Oriented Nonlinear Solver Package
#                 Copyright (2002) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#   
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#   
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
# 
# Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
# (rppawlo@sandia.gov).
# 
# ************************************************************************
//@HEADER
</pre>

Once the file is committed to the CVS repository, the first two lines
will look something like the following:

<pre>
// $Id$ 
// $Source$ 
</pre>

The header information is automatically filled in between the two
<code>//@HEADER</code> keys when we run the <code>nox/maintenance/autoheader</code> command.

</ul>


<H5>Include File Structure</H5>
<ul>
<li>No include file should define more than one class.
<ul>
<li> The exception is for utility classes that are <em>only</em> used by the main class defined in the header file.
</ul>
<br>&nbsp;
<li>Every include file must contain a mechanism that prevents multiple
inclusions of the file. For example, the following should follow the
header information for the <code>NOX_Abstract_Vector.H</code> header
file.

<pre>
#ifndef NOX_ABSTRACT_VECTOR_H
#define NOX_ABSTRACT_VECTOR_H

<em>...body of include file goes here...</em>

#endif
</pre>

<li>Do not include system files (e.g., <CODE>iostream</CODE>) directly
in your files. Instead, include <CODE>NOX_Common.H</CODE>. The goal is
to better enable system portability since some machines have
<CODE>math.h</CODE> and others have <CODE>cmath</CODE> and so
on. Currently, we have the following system headers:
<br>&nbsp;
<ul>
<li><code>iostream</code>
<li><code>iomanip</code>
<li><code>string</code>
<li><code>cmath</code>
<li><code>vector</code>
<li><code>map</code>
<li><code>deque</code>
<li><code>algorithm</code>
<br>&nbsp;
</ul>
<li>Definitions of classes that are only accessed via pointers
(<CODE>*</CODE>) or references (<CODE>&amp;</CODE>) should be declared
using forward declarations, and <EM>not</EM> by including the header files.
<br>&nbsp;

<li>These are the cases when header files should be included in the
header file:
<br>&nbsp;
<UL>
<LI>classes that are used as <EM>base classes</EM>,
<LI>classes that are used as <EM>member variables</EM>,
<LI>classes that appear as <EM>return types</EM> or as <EM>argument
types</EM> in function/member function prototypes.
</UL>
</ul>

<H3>Naming Conventions</H3>
<ul>
<li>Everything must be declared within the <code>NOX</code>
namespace. No exceptions! 
<br>&nbsp;
<li>Furthermore, each class should be within the appropriate
sub-namespace. The choices are:
<br>&nbsp;
<ul>
<li><code>Abstract</code>
<li><code>Parameter</code>
<li><code>Solver</code>
<li><code>StatusTest</code>
<li><code>LineSearch</code>
<li><code>Direction</code>
<li><code>Epetra</code>
<br>&nbsp;
</ul>
<li>Class names should begin with an uppercase letter. Variable and
function names should begin with a lowercase letter.
<br>&nbsp;
<li>In names (function, class, variable, etc) 
    which consist of more than one word, the words are
    written together and each word that follows the first is begun
    with an uppercase letter. (e.g., <code>NOX::Linesearch::MoreThuente</code>).
<br>&nbsp;
<li>Do not use identifiers which begin
with one or two underscores (`<CODE>_</CODE>' or `<CODE>__</CODE>').
<br>&nbsp;
<li>Names should not include abbreviations that are not generally
    accepted.
<br>&nbsp;
<li>Choose variable names that suggest the usage. 
</ul>

<H3>Style</H3>

<h5>Classes</h5>
<ul>
<li>The public, protected, and private sections of a class are to be
    declared in that order (the public section is declared before the
    protected section which is declared before the private section).
<br>&nbsp;
<li>No inline functions, except the empty <code>{}</code>
function. The reason for this is that if the function does not inline
correctly, it can actually lead to slower code rather than faster
code. 
</ul>
<h5>Functions</h5>
<ul>
<li>Always provide the
<EM>return type</EM> of a function explicitly.
<br>&nbsp;


<li>When declaring functions, the leading parenthesis and the first
argument (if any) are to be written on the <EM>same line</EM> as the
function name. If space permits, other arguments and the closing
parenthesis may also be written on the same line as the function
name. Otherwise, each additional argument is to be written on a
separate line (with the closing parenthesis directly after the last
argument).
<br>&nbsp;

<li>Always write the left parenthesis
directly after a function name.
<PRE>
   void foo ();    // No!!
   void foo();     // Better
</PRE>
</ul>

<h5>Variable declarations</h5>

<ul>
<li>Only one variable per line.
<PRE>
   int i,j;   // No!!

   int i;     // Yes   
   int j;   
</PRE>

<li>
The characters `<CODE>*</CODE>' and `<CODE>&amp;</CODE>' should
be written together with the types of variables instead of with the
names of variables in order to emphasize that they are part of the type
definition. Instead of saying that <CODE>*i</CODE> is an <CODE>int</CODE>,
say that <CODE>i</CODE> is an <CODE>int*</CODE>.
<PRE>
   int *i;   // No!!
   int* i;   // Yes   
</PRE>
</ul>

<h5>Loops and conditionals: <code>if</code>, <code>for</code>,
<code>while</code>, etc.</h5>

<ul>

<li>User parens to make code readable.
<pre>
  if (a == b && c < d || e == f)  // No!
  {
    /* Stuff */
  }

  if (((a == b) && (c < d)) || (e == f))  // Yes
  {
    /* Stuff */
  } 
</pre>
<li>The block of any <code>if</code> statement should always follow on
a separate line.
<PRE>
   if ( /*Something*/ ) i++; // No!!

   if ( /*Something*/ )      // Yes!
     i++; 
</PRE>
<li>Braces ("{}") which enclose a block should be aligned as follows:
<PRE>
   if ( /*Something*/ ) // Yes!
   {
     i++; 
     j++;
   }

   if ( /*Something*/ ) { // Okay
     i++; 
     j++;
   }

   if ( /*Something*/ ) // No!
     {
       i++; 
       j++;
     }
</PRE>

Adding the following line to your \c .emacs file will help:

<PRE>
(c-set-offset 'substatement-open 0)
</PRE>

</ul>


<h5>Miscellaneous</h5>
<UL>
<li>Always provide a space on both sides of <code>=</code> signs and all logical
operators.
<br>&nbsp;
<li>Each statement shall always be in a separate line, however small it may appear. 
<br>&nbsp;
 
<Li>Do not use <EM>spaces</EM> around
`<CODE>.</CODE>' or `<CODE>-&gt;</CODE>', nor between unary operators
and operands.
<br>&nbsp;
<li>Use the c++ mode in GNU Emacs to
format code.
</UL>

<H3>Coding Rules</H3>


<ul>
<li>A public member function must never return a non-const reference
or pointer to member data. 
<br>&nbsp;
<li>Constants are to be defined using <code>const</code> or
<code>enum</code>; never using <code>#define</code>. 
<br>&nbsp;
<li>A switch statement must always contain a default branch which
handles unexpected cases.  
</ul>


<H3>Output</H3>

The <code>NOX::Utils</code> class has static utility functions related
to printing. To use it, include <code>NOX_Utils.H</code>.

<ul>
<li>For <b>any</b> non-error print statement, call the
<code>NOX::Utils::doPrint()</code> function with the appropriate
MsgType flag. The flags are:
<br>&nbsp;
<ul>
<li><code>NOX::Utils::Error  </code>
<li><code>NOX::Utils::Warning</code>
<li><code>NOX::Utils::OuterIteration</code>
<li><code>NOX::Utils::InnerIteration</code>
<li><code>NOX::Utils::Parameters</code>
<li><code>NOX::Utils::Details</code>
</ul>

</ul>


<H3>Error Handling</H3>

<ul>
<li>Always check return values of functions for errors.
<li>In general, try to recover from errors.
<li>If you must throw an exception, always print an explanation with
the function name to <code>cerr</code> and then throw an exception with the string
<code>"NOX Error"</code>. For example,

<pre>
  if (/* Error Condition */) 
  {
    cerr << "ERROR: NOX::Epetra::Group::getNewton() - invalid Newton vector" << endl;
    throw "NOX Error";
  }
</pre>

</ul>

<H3>Comments</H3>

We use Doxygen for the comments. To generate the documentation, do the
following:

<pre>
cd nox/doc
doxygen 
</pre>

<ul>
<li> Document each class, function, and enum in the header
files.
<br>&nbsp;
<ul>
<li> The one exception is that functions in derived objects do not
need to be documented <b>if</b> the documentation is inherited from
the base class. This should be tested in Doxygen to be sure that it
works correctly.
<br>&nbsp;
</ul>
<li> Here's an example of documented a class. Note the formatting of
the comments. It's a C-style comment. The open comment marker
(<code>/*</code>) is followed by an exclamation mark to indicate that
it's a Doxygen comment. The open and close comment markers are on
lines by themselves, and the text of the comment is indented two
spaces. Always include a <code>\brief</code> description. The long
description follows. Observe the use of the formatting tags
<code>\c</code> and <code>\e</code>. The <code>\note</code> tag is
used for any special notes. The <code>\author</code> tag is
recommended.

<pre>
/*!
  \brief Arbitrary combination of status tests.

  In the \c AND (see NOX::Status::Combo::ComboType) combination, the
  result is \c Unconverged (see NOX::Status::StatusType) if \e any of
  the tests is \c Unconverged. Otherwise, the result is equal to the
  result of the \e first test in the list that is either \c Converged
  or \c Failed. It is not recommended to mix \c Converged and \c
  Failed tests in an \c AND combination.

  In the \c OR combination, the result is \c Unconverged if \e all of
  the tests are \c Unconverged. Otherwise, it is the result of the \e
  first test in the list that is either \c Converged or \c
  Failed. Therefore, it will generally make sense to put the \c Failed
  -type tests at the end of the \c OR list.

  \note We always runs through all tests, even if we don't need
  to. This is useful so that the user knows which tests have and have
  not be satisfied.

  \author Tammy Kolda (SNL 8950)
*/
class Combo : public Test {
...
}; // class Combo
</pre>

<li>Any parameters that are used within the class <b>must</b> be
documented in the class description <b>and</b> in the file NOX_Description.H on the parameters "page".
Note that the name, a brief description, and the default value for
each parameter is listed.
<pre>
/*!
  \brief %Newton-like solver with a line search.
  
  The following parameters are valid for this solver:

  - "Line Search" - Sublist of the line search parameters, passed to
    the NOX::Linesearch::Manager constructor. Defaults to an empty list.

  - "Linear %Solver" - Sublist of the linear solver parameters, passed
    to Abstract::Group::computeNewton(). Furthermore, the "Tolerance"
    within this list may be modified by the
    resetForcingTerm(). Defaults to an empty list.

  - "Forcing Term Method" - Method to compute the forcing term, i.e.,
    the tolerance for the linear solver. Defaults to ""
    (nothing). Choices are "Type 1" and "Type 2".

  - "Forcing Term Minimum Tolerance" - Minimum acceptable linear
    solver tolerance. Defaults to 1.0e-6.

  - "Forcing Term Maximum Tolerance" = Maximum acceptable linear
    solver tolerance. Default to 0.01.

  - "Forcing Term Alpha" - Used for the "Type 2" forcing term
    calculation. Defaults to 1.5.

  - "Forcing Term Gamma" - Used for the "Type 2" forcing term
    calculation. Defaults to 0.9.

  \author Tammy Kolda (SNL 8950), Roger Pawlowski (SNL 9233)
*/
</pre>
Here's a more complicated example to produce a two-tiered list.
<pre>
/*! 
   The parameters must specify the type of line search as well as all
   the corresponding parameters for that line search.
 
   &lt;ul&gt;
   &lt;li&gt; "Method" - Name of the line search. Valid choices are
     &lt;ul&gt; 
     &lt;li&gt; "Full Step" (NOX::Linesearch::FullStep)
     &lt;li&gt; "Interval %Halving" (NOX::Linesearch::Halving)
     &lt;li&gt; "Polynomial" (NOX::Linesearch::Polynomial)
     &lt;li&gt; "More'-Thuente" (NOX::Linesearch::MoreThuente)
     &lt;/ul&gt;
   &lt;/ul&gt;
 */
</pre>


<li>Constants and enums can generally be described with simple
<code>\brief</code> comments. Those can be formatted in either of two
ways, as follows.

<pre>
  /*! 
    \brief The test can be either the AND of all the component
    tests, or the OR of all the component tests.
  */
  enum ComboType {AND, OR};

  //! Constructor
  Combo(ComboType t = OR);
</pre>

<li>Doxygen does automatically cross-linking, which is very
convenient. However, sometimes it cross-links when you don't intend
for it to. For example, the following line would automatically
generate a link from the word <code>Newton</code> to the
NOX::Solver::Newton class. 

<pre>
//! Newton-like solver with a line search.
</pre>
 
To prevent that automatic link, insert a percent sign (<code>%</code>)
immediately before the word that is causing the link. For example,

<pre>
//! %Newton-like solver with a line search.
</pre>
</ul>
\endhtmlonly

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/* \page nox_class_overview NOX Class Overview

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
  the class description.

  Included with %NOX are three such interfaces:

  - NOX::LAPACK
  - NOX::Epetra
  - NOX::Petsc

  The \b LAPACK interface is an interface to the BLAS/LAPACK library.
  It is not intended for large-scale computations, but to serve as an
  easy-to-understand example of how one might interface to %NOX.  
  
  The \b Epetra interface is an interface to the
  Trilinos/Epetra library developed by 
  <a href="http://www.sandia.gov">Sandia National Labs</a>; 
  for more information see the 
  <a href="http://software.sandia.gov/Trilinos">Trilinos Home Page</a>.

  The \b PETSc interface is an interface with the PETSc library. PETSc
  was developed at 
  <a href="http://www.anl.gov">Argonne National Labs</a>; 
  for more information see the
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
 
/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page parameters NOX Parameter Reference Page

  All parameters are specified via parameter lists as defined by the
  NOX::Parameter::List class.

  Here is the basic hierarchical structure of the top level list that
  is passed to a NOX::Solver::Manager object constructor reset call.

<ul>
<li> "Nonlinear Solver" - see NOX::Solver::Manager for options

<li> "Direction" sublist - used by NOX::Solver::LineSearchBased,
     NOX::Solver::TrustRegionBased

     <ul>

     <li> "Method" - see NOX::Direction::Manager for options

     <li> "User Defined Constructor" - see NOX::Direction::Manager for details

     <li> "Newton" sublist - used by NOX::Direction::Newton,
          see class documentation for all details.

	  <ul>

	  <li> "Forcing Term Method"
	  <li> "Forcing Term Initial Tolerance"
	  <li> "Forcing Term Minimum Tolerance"
	  <li> "Forcing Term Maximum Tolerance"
	  <li> "Forcing Term Alpha"
	  <li> "Forcing Term Gamma"
	  <li> "Forcing Term User Defined Norm"
	  <li> "Rescue Bad Newton Solver"

	  <li> "Linear Solver" sublist, passed to
	       NOX::Abstract::Group::computeNewton() by
	       NOX::Direction::Newton

	       <ul>

	       <li> "Tolerance"

	       </ul>

	  </ul>

     <li> "Steepest Descent" sublist, used by
          NOX::Direction::SteepestDescent, see class documentation for
          details

	  <ul>

	  <li> "Scaling Type"

	  </ul>

     </ul>

<li> "Cauchy Direction" sublist - used by
     NOX::Solver::TrustRegionBased, same options as "Direction" sublist

<li> "Line Search" sublist - used by NOX::Solver::LineSearchBased

     <ul>

     <li> "Method" - see NOX::LineSearch::Manager for options

     <li> "User Defined Constructor" - see NOX::LineSearch::Manager for details

     <li> "Backtrack" sublist - used by NOX::LineSearch::Backtrack,
          see class documentation for all details.

	  <ul>

	  <li> "Minimum Step"
	  <li> "Default Step"
	  <li> "Recovery Step"
	  <li> "Maximum Step"
	  <li> "Max Iters"
	  <li> "Decrease Condition"

	  </ul>

     <li>"Full Step" sublist - used by NOX::LineSearch::FullStep, 
         see class documentation for all details

	 <ul>

	 <li> "Full Step"

	 </ul>

     <li> "More'-Thuente" sublist - used by
          NOX::LineSearch::MoreThuente, see class documentation for all details

	  <ul>

	  <li> "Sufficient Decrease"
	  <li> "Curvature Condition"
	  <li> "Interval Width"
	  <li> "Minimum Step"
	  <li> "Maximum Step"
	  <li> "Max Iters"
	  <li> "Default Step"
	  <li> "Recovery Step"

	  </ul>


     <li> "Polynomial" sublist - used by NOX::LineSearch::MoreThuente,
          see class documentation for details

	  <ul>
	  <li> "Default Step"
	  <li> "Max Iters"
	  <li> "Minimum Step"
	  <li> "Recovery Step"
	  <li> "Interpolation Type"
	  <li> "Min Bounds Factor"
	  <li> "Max Bounds Factor"
	  <li> "Sufficient Decrease Condition"
	  <li> "Alpha Factor"
	  <li> "Force Interpolation"
	  <li> "Use Counters"
	  <li> "Maximum Iteration for Increase"
	  <li> "Allowed Relative Increase"
	  <li> "User Defined Merit Function"
	  <li> "User Defined Norm"
	  </ul>

     </ul>


<li> "Trust Region" sublist, used by NOX::Solver::TrustRegionBased,
     see class documentation for details

     <ul>

     <li> "Minimum Trust Region Radius"
     <li> "Maximum Trust Region Radius"
     <li> "Minimum Improvement Ratio"
     <li> "Contraction Trigger Ratio"
     <li> "Contraction Factor"
     <li> "Expansion Trigger Ratio"
     <li> "Expansion Factor"
     <li> "Recovery Step"

     </ul>

<li> "Solver Options" sublist used by all solvers for basic options that all solvers should support (see individual solvers to make sure the option is supported).  
     <ul>

     <li> "User Defined Pre/Post Operator" - a
     NOX::Parameter::PrePostOperator derived object that allows users
     to insert pre and post method calls into nox's solvers (before
     and after the NOX::Solver::Generic::iterate() and before and
     after the NOX::Solver::Generic::solve() methods).  This object
     can be set in the parameter list since it is derived from
     NOX::Parameter::Arbitrary.  See NOX::Parameter::PrePostOperator
     for details on how to implement this object.

     </ul>

<li> "Printing" sublist - used by NOX::Utils, see class documentation
     for details

     <ul>

     <li> "Output Information"
     <li> "MyPID"
     <li> "Output Processor"
     <li> "Output Precision"
     
     </ul>

</ul>

*/


/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page epetra_interface The Epetra interface

 \section epetra_interface_derived Step 2a: Implementing the interface

 To link your code to NOX using the epetra interface, you must
 write a concrete class derived from the NOX::Epetra::Interface
 object.  This object defines four methods:

 <ul>
 <li> computeF(const Epetra_Vector& x, Epetra_Vector& RHS, FillType flag) - 
      Computes the set of nonlinear equations, \f$ F \f$, to be solved 
      by NOX. This method must be supplied by the user.  
      - x - solution vector specified from NOX.
      - RHS - the Epetra_Vector to be filled with the \f$ F \f$ values that 
        correspond to the input solution vector x.
      - flag - enumerated type (see NOX::Epetra::FillType) that 
        tells a users interface why computeF() was called.  NOX 
        has the ability to generate Jacobians based on numerical 
        differencing using calls to computeF().  In this case, the 
	user may want to compute an inexact (and hopefully cheaper) 
        \f$ F \f$ since it is only used in the Jacobian (or preconditioner).  
      
 <li> computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac) - this
      is an optional method that the user can implement if they wish to 
      supply their own evaluation of the Jacobian. If the user does not 
      wish to supply their own Jacobian, they should implement this method 
      so that it throws an error if it is called.  This method should update 
      the Jac operator so that subsequent Epetra_Operator::Apply() calls 
      on that operator correspond to the Jacobian at the current solution 
      vector x.
      - x - solution vector specified from NOX.
      - Jac - a reference to the Jacobian operator that the user supplied 
        in the NOX::Epetra::Group constructor.  

 <li> computePrecMatrix(const Epetra_Vector& x, Epetra_RowMatrix& M) - This 
      method fills a Preconditioner matrix (Epetra_RowMatrix) object that the 
      user supplied in the NOX::Epetra::Group constructor.  The Matrix 
      should NOT be inverted, It should only contain an estiamte of the 
      Jacobian.  After the user fills the Row Matrix epetra/aztecOO 
      preconditioners will compute the inverse.  
      If the user does not wish to supply their own Preconditioning matrix, 
      they should implement this method such that if called, it will 
      throw an error.
      - x - solution vector specified from NOX.
      - M - a reference to the preconditioner operator that is to be 
        filled.  This operator should be updated to correspond to the 
	current solution vector x.

 <li> computePreconditioner(const Epetra_Vector& x, Epetra_Operator& M) - This
      method allows a user to supply their own preconditioner.  The method 
      should compute a preconditioner based upon the solution vector 
      x and store it in the Epetra_Operator M.  Subsequent calls to the 
      Epetra_Operator::Apply method will apply this user supplied 
      preconditioner to epetra vectors. 
      - x - solution vector specified from NOX.
      - M - a reference to the operator that is to be filled.  This operator 
	should be updated to correspond to the current solution vector x.

 </ul>

 \section epetra_interface_construction Step 2b: Construct an Epetra Group

  Once the interface is written, the user can construct a
  NOX::Epetra::Group for use with the NOX solvers.  All constructors
  require (1) a parameter list for the printing output, (2) a
  parameter list for the linear solver options, (3) an initial guess
  for the solution, and (4) an operator for the Jacobian, J.
  Optionally (5), the user may supply a second operator, M for the
  preconditioner.

  (1) The printing parameter list is a NOX::Parameter::List that
  descibes the type of output requested from NOX.  See \ref parameters
  for more details.  A simple example is:
  \verbatim
  NOX::Parameter::List& printParams = nlParams.sublist("Printing");
  printParams.setParameter("MyPID", MyPID); 
  printParams.setParameter("Output Precision", 3);
  printParams.setParameter("Output Processor", 0);
  printParams.setParameter("Output Information", 
			NOX::Utils::OuterIteration + 
			NOX::Utils::OuterIterationStatusTest + 
			NOX::Utils::InnerIteration +
			NOX::Utils::Parameters + 
			NOX::Utils::Details + 
			NOX::Utils::Warning);
  \endverbatim

  (2) The linear solver parameter list is a NOX::Parameter::List
  object.  A list of all parameters available can be found in
  NOX::Epetra::Group.  Below is a simple example of how to use the
  aztecOO iterative solver in epetra to use GMRES with ilut
  preconditioning (WARNING - not all preconditioning options are
  compatible with every Jacobian and preconditioner operator - see
  NOX::Epetra::Group for more details):
  \verbatim
  NOX::Parameter::List& lsParams = newtonParams.sublist("Linear Solver");
  lsParams.setParameter("Aztec Solver", "GMRES");  
  lsParams.setParameter("Max Iterations", 800);  
  lsParams.setParameter("Tolerance", 1e-4);
  lsParams.setParameter("Output Frequency", 50); 
  lsParams.setParameter("Preconditioning", "AztecOO: Jacobian Matrix");
  lsParams.setParameter("Aztec Preconditioner", "ilut"); 
  lsParams.setParameter("Overlap", 2);   
  lsParams.setParameter("Fill Factor", 2);   
  lsParams.setParameter("Drop Tolerance", 1.0e-12);   
  \endverbatim

  (3) The initial guess is stored in an Epetra_Vector.  Below is an
  example of creating an Epetra_Vector and setting the initial guess
  to one.  The epetra block map is an epetra object that sets the size
  of the problem and has information on how the unknowns are
  distributed across processors.  For examples on how to construct this
  see the examples-epetra directory in the nox source code.
  \verbatim 
  Epetra_Vector initialGuess(EpetraBlockMap); 
  initialGuess.PutScalar(1.0);
  \endverbatim

  (4) and (5)The user can write their own operators for the Jacobian or
  preconditioner or use one of the NOX supplied operators.  For
  example, here are three ways to create a Jacobian operator:
  - User defined (user returns their own operator from an object they created):
    \verbatim
    Epetra_Operator& J = UserProblem.getJacobian();
    \endverbatim
  - Create a NOX matrix-free operator (see NOX::Epetra::MatrixFree):
    \verbatim
    NOX::Epetra::MatrixFree J(userDerivedInterface, solutionVector);
    \endverbatim
  - Create a NOX finite difference operator 
    (see NOX::Epetra::FiniteDifference):
    \verbatim
    NOX::Epetra::FiniteDifference J(interface, soln);
    \endverbatim

  Create the NOX::Epetra::Interface derived object described in 
  \ref epetra_interface_derived:
  \verbatim
  UserInterface interface(.......);
  \endverbatim

  Finally, construct the group:
  \verbatim
  NOX::Epetra::Group group(printParams, lsParams, interface, initialGuess, J);
  \endverbatim

  You are ready to move on to \ref step3

 \section epetra_interface_notes Notes on NOX/Epetra:

  - The concrete Epetra implementations for the NOX Group and Vector
    are in a separate library from the nox solver algorithms.  To build
    the library for nox epetra support use the flag --enable-nox-epetra
    in NOX's configure script.  This will generate a separate library called
    libnoxepetra.a

  - In addition to the NOX headers, be sure to include the NOX Epetra 
    support specific headers in your interface:
    \verbatim
    // NOX headers for an epetra interface
    #include "NOX.H"
    #include "NOX_Epetra.H"
    \endverbatim

*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page petsc_interface The PETSc interface

At present, %NOX supports a concrete implementation of the
<b>NOX::Abstract::Group</b> and <b>NOX::Abstract::Vector</b> using
PETSc data structures in the namespace NOX::Petsc.  The current
implementation has been tested against PETSc 2.1.6.

\subsection build_nox_petsc_library Building the Library

To compile the %NOX/PETSc library \c libnoxpetsc.a, use the following 
configuration sequence:

\verbatim
setenv PETSC_DIR ${HOME}/petsc-2.1.6
setenv PETSC_ARCH linux
configure --enable-nox-petsc
\endverbatim
This will build \c libnoxpetsc.a using a pre-built PETSc
version 2.1.6 installed in the user's home directory on a linux machine:

If the environment variables \c PETSC_DIR and \c PETSC_ARCH have already been
defined, simply use:
\verbatim
configure --enable-nox-petsc
\endverbatim

A third option is to specify \c PETSC_DIR and \c PETSC_ARCH as arguments
to the configure script, eg
\verbatim
configure --enable-nox-petsc PETSC_DIR=${HOME}/petsc-2.1.6 PETSC_ARCH=linux
\endverbatim


\subsection build_nox_petsc_examples Building the Examples

To compile the examples for %NOX/PETSc, an additional environment variable,
\c PETSC_BOPT, must be defined and set to the value of BOPT used to build the 
user's PETSc installation.  For example, for PETSc previously built 
debug-enabled for C++ using BOPT=g_c++, simply set PETSC_BOPT=g_c++.  
The configure script then allows the examples to be built using the following:
\verbatim
configure --enable-nox-petsc --enable-nox-petsc-examples.
\endverbatim


\subsection additional_nox_petsc_build_information Additional Information

Depending on the user's PETSc build, dynamically linked libraries may be used.
If an error occurs when trying to run the example executable, make sure the
environment variable \c LD_LIBRARY_PATH points to the PETSc library directory.
From the example above, setting the path would likely require:
\verbatim
setenv LD_LIBRARY_PATH ${PETSC_DIR}/lib/lib_${PETSC_BOPT}/${PETSC_ARCH}
\endverbatim


*/

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*
\section ss1 Building the NOX/Epetra Library

To compile NOX's Epetra library \c libnoxepetra.a, use one
of the following sets of configuration options.

--- These instructions are out-of-date. They will be updated soon. ---

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

Similarly, it will look in the following locations
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


\section ss2 Building the NOX/Epetra Examples

To compile the examples for NOX/Epetra, do the following:

\verbatim
configure --enable-epetra --enable-epetra-examples.
\endverbatim

\section section5 Using MPI with NOX/Epetra

%NOX itself does not call any MPI commands, but some of the %NOX/Epetra
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
library is specified by \c -lmpi.
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
library is specified by \c -lmpi.
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


*/


