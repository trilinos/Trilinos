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

  Since our code has not yet been approved for general release, \b Do
  \b Not \b Distribute \b Outside \b of \b Sandia \b Without \b
  Express \b Permission \b from \b Tammy \b Kolda \b or \b Roger \b
  Pawlowski.


  <ul>
  <li> 
  <a href="http://software.sandia.gov/nox/downloads/nox-1.0a-dev.tar.gz">Latest Development Version</a>, updated nightly
  </ul>


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

\section contributors NOX Contributors

The following people have contributed code to NOX:

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

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page confoptions NOX Configuration Options 

Here is a detailed list and explanation of the NOX configuration options.

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

\section mpi_config MPI-Related Options

The NOX library itself has no MPI dependencies. However, the NOX-Epetra
library and related examples do depend on MPI. 

<ul>
  <li> \c --enable-mpi
    <ul>
      <li> Enables MPI mode. Defines HAVE_MPI in the NOX_Config.h
           file. Will test for the ability to preprocess the MPI header file
           and may test ability to link with MPI.
    </ul>
  <li> \c --with-mpi-cxx
    <ul>
      <li> If no argument is specified, uses "mpiCC". Otherwise, uses
           specified MPI C++ compiler. Automatically enables MPI mode.
    </ul>
  <li> \c --with-mpi=MPIROOT
    <ul>
      <li> Specify the MPI root directory. Automatically enables MPI mode.
    </ul>
  <li> \c --with-mpi-libs="LIBS"
    <ul>
      <li> Specify the MPI libraries. 
           Defaults to \c "-lmpi" if either \c --with-mpi or \c --with-mpi-libdir is specified.
    </ul>
  <li> \c --with-mpi-incdir=DIR
    <ul>
      <li> Specify the MPI include files location. 
           Defaults to \c MPIROOT/include if \c --with-mpi is specified. 
	   If multiple directories must be specified, try \c --with-cppflags="-I<dir1> \c -I<dir2>" instead.
    </ul>
  <li> \c --with-mpi-libdir=DIR
    <ul>
      <li> Specify the MPI libraries location. 
           Defaults to \c MPIROOT/lib if \c --with-mpi is specified.
	   If multiple directories must be specified, try \c --with-ldflags="-L<dir1> \c -L<dir2>" instead.
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
    </ul>
  <li> \c --enable-nox-petsc-examples
    <ul>
      <li> Compile petsc examples. 
           Requires \c --enable-nox-petsc.
    </ul>
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
    <li> Compile code that is not included in the nightly release but is part of the CVS repository. 
         See \ref prerelease.
  </ul>
</ul>

*/


/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page configuration NOX Configuration and Compliation Instructions

See \ref confoptions

\section section1 Basic Setup

To just compile the basic NOX library, type

\verbatim
configure
make
\endverbatim

This only compiles the library \c libnox.a - it does not build any of the
examples. See below for more information on compiling the examples.

\section section2  Generating Documentation

The documentation is available online at 
<a href="http://software.sandia.gov/nox/">http://software.sandia.gov/nox/</a>.

To generate a local copy of the documentation, you'll need to have Doxygen installed. 
After running configure, type 
\verbatim
make dox
\endverbatim
to generate the documentation. The main page is \c nox/doc/html/index.html.

\section section3 Compiling the Examples

NOX is distributed with a simple set of examples that have no other
external dependencies besides LAPACK and BLAS. These are primarily useful to
allow some basic testing of NOX without requiring an interface to a
more complicated linear algebra system. To enable these examples, type

\verbatim
configure --enable-lapack --enable-lapack-examples
make
\endverbatim

The compiled examples are located in the \c example-lapack
subdirectory.  To run each example, follow the instructions given below.

<ol>
<li>To run a simple 2 equation example based on the Rosenbrock function, type
\verbatim
rosenbrock
\endverbatim
</ol>

*/

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

/* ************************************************************************ */
/* ************************************************************************ */
/* ************************************************************************ */

/*! \page devinfo NOX Developer's Information

To become a NOX developer, contact Tammy Kolda (tgkolda@sandia.gov) or
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

/*! \page prerelease NOX Developer's Guide to Prerelease Code


We define \b prerelease \b code to be code that should not be
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
nox/src/Makefile.am to add NOX_Direction_QuasiNewton as prereease code
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

The CVS repository is located on <a
href="http://software.sandia.gov/>software.sandia.gov</a>. See Paul
Sery (pgsery@sandia.gov) or Tammy Kolda (tgkolda@sandia.gov) for an
account.

You should not do any development of software.sandia.gov. Instead,
check out a copy to your local machine.  Be sure to correctly set your
\c CVS_SSH environment variable to \c ssh.

<ul>
<li>
To check out nox:
\verbatim
cvs -d :ext:username@hopper.ca.sandia.gov:/space/CVS checkout -P nox
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
<li>Furthermore, each class should be within the approiate
sub-namespace. The choices are:
<br>&nbsp;
<ul>
<li><code>Abstract</code>
<li><code>Parameter</code>
<li><code>Solver</code>
<li><code>Status</code>
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
correctly, it can actaully lead to slower code rather than faster
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

<h5>Variable declations</h5>

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
  if (a == b && c < d || e == f) { // No!
    /* Stuff */
  }

  if (((a == b) && (c < d)) || (e == f)) { // Yes
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
<li>Braces ("{}") which enclose a block should be aligned in the
K&amp;R format. 
<PRE>
   if ( /*Something*/ ) // No!
   {
     i++; 
     j++;
   }

   if ( /*Something*/ ) { // Yes!
     i++; 
     j++;
   }
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

<br>&nbsp;

<li><em>We need to add a mechanism for controlling the
output streams. Currently, all output goes to cout and cerr.</em>
<br>&nbsp;

</ul>


<H3>Error Handling</H3>

<ul>
<li>Always check return values of functions for errors.
<li>In general, try to recover from errors.
<li>If you must throw an exception, always print an explanation with
the function name to <code>cerr</code> and then throw an exception with the string
<code>"NOX Error"</code>. For example,

<pre>
  if (/* Error Condition */) {
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

  - "Linear %Solver" - Sublist of the linear solver paramaters, passed
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
    calcuation. Defaults to 1.5.

  - "Forcing Term Gamma" - Used for the "Type 2" forcing term
    calcuation. Defaults to 0.9.

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
   &lt;li&gt; "%Polynomial" (NOX::Linesearch::Polynomial)
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
convienent. However, sometimes it cross-links when you don't intend
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

     <li> "Newton" sublist - used by NOX::Direction::Newton,
          see class documentation for all details.

	  <ul>

	  <li> "%Linear Solver" sublist, passed to
	       NOX::Abstract::Group::computeNewton() by
	       NOX::Direction::Newton

	       <ul>

	       <li> "Tolerance"

	       </ul>

	  <li> "Rescue Bad Newton Solver"
	  <li> "Forcing Term Method"
	  <li> "Forcing Term Initial Tolerance"
	  <li> "Forcing Term Minimum Tolerance"
	  <li> "Forcing Term Maximum Tolerance"
	  <li> "Forcing Term Alpha"
	  <li> "Forcing Term Gamma"

	  </ul>

     <li> "Quasi-Newton" sublist, used by
          NOX::Direction::QuasiNewton, see class documentation for
          details

	  <ul>

	  <li> "Memory"

	  </ul>

	  

     <li> "Nonlinear CG" sublist, used by
          NOX::Direction::NonlinearCG, see class documentation for
          details

	  <ul>

	  <li> "Orthogonalize"
	  <li> "Precondition"
	  <li> "Restart Frequency"

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

	  <li> "Convergence Criteria"
	  <li> "Interpolation Type"
	  <li> "Minimum Step"
	  <li> "Default Step"
	  <li> "Recovery Step"
	  <li> "Max Iters"
	  <li> "Alpha Factor"
	  <li> "Min Bounds Factor"
	  <li> "Max Bounds Factor"
	  <li> "Force Interpolation"

	  </ul>

     <li> "Secant" sublist - used by NOX::LineSearch::Secant, see
          class documentation for details

	  <ul>

	  <li> "Minimum Step"
	  <li> "Default Step"
	  <li> "Recovery Step"
	  <li> "Max Iters"

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

