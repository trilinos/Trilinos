=====================================
A HelloWorld TriBITS Project
=====================================

:Author: Joe Frye (jfrye@sandia.gov)
:Date: |date|

.. |date| date::

.. sectnum::
   :depth: 2

.. Sections in this document use the underlines:
..
.. Level-1 ==================
.. Level-2 ------------------
.. Level-3 ++++++++++++++++++
.. Level-4 ..................

.. contents::


Introduction
=============

This tutorial is designed to give you a basic overveiw of the
structure of a TriBITS project at a high level.  See the tutorial on
TriBITS Example Project for more detail and concrete examples. TriBITS
defines a few structural units hat are then put together to form a
tribits project.  In this tutorial we will discuss what each of these
unts are, and how to assemble them to create a tribits project.  One
key advantage of TriBITS is that once a project is using TriBITS, it
is east for another project to use the packages defined in another
TriBITS project.


TriBITS Packages
=================

The simplest way to think of a package is as a collection of source
files being built into one (or more, but ideally one) library or
executable and a set of tests.  If you have componetns of your project
that fit this description then consider making into a package.  A
package does not have to build just one library or executable and if
there are targets that are always built together then make them all
part of the same package.  

The goal of a tribits package is straight forward.  Take some source
files and build a target and the tests for that target.  This will be
the main function of the CMakeLists.txt file in top level package
directory.  In that CmakeLists.txt file you need to tell tribits which
source files to use to build each target and which targets to build.
You will also need to define any dependencies this package may have on
other packages in the project.

In order top define a package, you must have the following files
defined::



TriBITS Projects
================

A tribits project is just a collection of Tribits packages, there is
no limit to how many packages a project can contain but it must
contain at least one. A project can include packages that are defined
in the same repository as the project or it can use packages from
TriBITS repositories (see below for more detail). Source files and
build targets cannot be specified at the project level, that must be
done in a package. Once you have Tribits packages defined then you can
put them together in a Tribits project.  In order to do this you need
to define some things at the project level such as packages your
project wishes to use and where to find them, any TPL dependencies you
may have, Project name and other project options, and the current
version of the project.

In order top define a project, you must have the following files
defined::

CTestConfig.cmake
PackagesList.cmake
ProjectName.cmake
TPLsList.cmake
Version.cmake


TriBITS Sub-packages
=====================

Similar to the way a TriBITS project contains packages that can
optionally be turned on and off, a package can defein subpackages that
can be optionally enabled. Subpackages allow for a finer level of
control over dependencies as each subpackage can define it's own
dependencies.  A subpackage cannot have subpackages of it's own.

In order top define a subpackage, you must have the following files
defined::


TriBITS Repositories
=====================

A TriBITS repository is a collection of TriBITS packages.  Like a
project a repository must define which packages it contains, any TPL
dependenecies, and options.  Dependencies between packages are still
defined by the package.  A TriBITS project can also be used as a
TriBITS repository for another TriBITS project.

In order top define a tribits repository, you must have the following
files defined::

TriBITS TPLS
============

A TPL is software that your project depends on that is not built by
your project.  You will have to tell TriBITS which TPLs to find and
how/where to look for them.


Summary
========

The basic unit of a TriBITS build system is the TriBITS package where
source files are identified and targets are built.  Tribits has a
standard way to handle dependencies within packages and will
automatically enable/disable all appropriate packages for a specified
configuration.  A TriBITS project is made up of one or more TriBITS
packages that is either defined in the project's repo or in a TriBITS
repository.  This document is intended to be brief and give just a
high level summary of the structure of a TriBITS project.  Please see
the tutorial on converting a project to TriBITS or the TriBITS Example
Project tutorial for much more detail and examples.
