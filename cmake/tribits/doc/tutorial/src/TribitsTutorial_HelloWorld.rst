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


The Simplest HelloWorld project
================================

This short tutorial will walk you through setting up a basic
HelloWorld project built by TriBITS to intdroduce basic concepts in
TriBITS.  To begin you will need cmake and TriBITS installed on your
machine.  You will also need a working c and c++ compiler.

**Before you begin make sure you have:**
- CMake installed
- Tribits installed 
- C compiler
- C++ compiler

Initial Setup
-------------------

TriBITS projects have a specific structure ie a project is a
collection of packages.  Packages may or may not depend on other
packages, and may or may not be required.  For this example we will be
creating a project that has just one package, the "HelloPackage"
package. First lets create all the directories for our project.  We
will need a top level directory for the project which I will call
tribits_hello_world. We need a directory for the "HelloPackage"
package.  We will also need a directory for the build which I call
"build".  Under the hello_package_dir also create the directories
"cmake" and "src"

You should have the following directories::

  tribits_hello_world/
  tribits_hello_wolrd/build
  tribits_hello_world/hello_package_dir
  tribits_hello_world/hello_package_dir/cmake
  tribits_hello_world/hello_package_dir/src

  $ tree
  .
  |__ tribits_hello_world
      |__ build
      |__ hello_package_dir
          |__ cmake
          |__ src


Create a TriBITS package
------------------------

Any TriBITS package needs to have at least 3 files.

- a top level CMakeLists file
- a file that track package dependencies
- source files


First lets create a single source file which is the classic
HelloWorld.cpp.  Just copy this to HelloWorld.cpp in the src
directory::

  #include <iostream>

  int main()
  {
    std::cout << "Hello World!\n";
    return 0;

  }

Second lets create the package dependencies file which should be
placed in the cmake directory.  Copy the below text into a file called
Dependencies.cmake::


  tribits_package_define_dependencies()


In this case the package we are creating has no dependencies but we
still need this file.  The lack of arguments to the
tribits_package_define_dependencies() call reflects that this package
does not have dependencies.  The last and most interesting file we
will create in the package directory is the CMakeLists.txt file.  Copy
the following into CMakeLists.txt::

  tribits_package(HelloPackage)
  
  tribits_add_executable(Hello-Executable-Name NOEXEPREFIX SOURCES
    src/HelloWorld.cpp INSTALLABLE)
  
  tribits_package_postprocess()

**tribits_package(HelloPackage)** Sets up a TriBITS package with the
name "HelloPackage"

**tribits_add_executable(Hello-Executable-Name NOEXEPREFIX SOURCES src/HelloWorld.cpp INSTALLABLE)** 
  tells TriBITS that we want to build an executable named
  "Hello-Executable-Name" from the source file src/HelloWorld.cpp.
  NOEXEPREFIX and INSTALLABLE are options to tribits_add_executable()
  that I will not go into right now.

**tribits_package_postprocess()** Must be at the end of any
packages top level CMakeLists file


Create a Tribits Project
------------------------

Recall that a TriBITS project is made up of TriBITS packages.  We have
just defeined a package now we will create a project that consists of
just that one package.  In order to do this we are going to create 4
files in the top level directory and they are named:

- CMakeLists.txt
- PackageList.cmake
- ProjectName.cmake
- TPLsList.cmake


**TPLsList.cmake** this file tells Tribits about TPLs needed for the
project.  In this case, the package does not depend on any TPLs so
this file will be very simple.  It should contain just the following
single line::

  tribits_repository_define_tpls()

**ProjectName.cmake** this file sets the name of the project.  Some
 other options can be specified in this file but we will just set the
 project name. It should contain the following::
  
  set(PROJECT_NAME TribitsHelloWorld)

**PackageList.cmake** defines which packages are in the project.  We
 will just need to tell it the name and location of our one package::

  tribits_repository_define_packages(
    HelloPackage  hello_package_dir  PT
  )

**CMakeLists.txt** This is the most interesting file in this example.
 Here we will set a minimum cmake version, load some options, and tell
 cmake that this is a Tribits project.  The CMakeLists.txt file should
 have the following contents::

  # To be safe, define your minimum CMake version
  cmake_minimum_required(VERSION 3.23.0 FATAL_ERROR)
  
  # Make CMake set WIN32 with CYGWIN for older CMake versions
  set(CMAKE_LEGACY_CYGWIN_WIN32 1 CACHE BOOL "" FORCE)
  
  # Get PROJECT_NAME (must be in file for other parts of system)
  include(${CMAKE_CURRENT_SOURCE_DIR}/ProjectName.cmake)
  
  # CMake requires that you declare the CMake project in the top-level file 
  project(${PROJECT_NAME} NONE)

  # This needs to be set to the path to the installation of TriBITS on your machine 
  set(${PROJECT_NAME}_TRIBITS_DIR 
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/tribits CACHE PATH "TriBITS base
  directory (default assumes in TriBITS source tree).")

  # Include the TriBITS system
  include("${${PROJECT_NAME}_TRIBITS_DIR}/TriBITS.cmake")
  
  # MPI and Fortran are enabled by default, turn them off for this project
  set(TPL_ENABLE_MPI OFF CACHE BOOL "" FORCE)
  # Turn off Fortran support by default
  set(${PROJECT_NAME}_ENABLE_Fortran_DEFAULT OFF)
  
  # Only one package in this simple project so just enable it :-)
  set(${PROJECT_NAME}_ENABLE_HelloPackage ON CACHE BOOL "" FORCE)
  
  # Do all of the processing for this Tribits project
  tribits_project()

**${PROJECT_NAME}_TRIBITS_DIR** Make sure you set this to your Tribits
Installation path it may not be the same as this path.  Now you should
have a directory structure that looks like this::

  .
  |__ CMakeLists.txt
  |__ PackagesList.cmake
  |__ ProjectName.cmake
  |__ TPLsList.cmake
  |__ build
  |__ hello_package_dir
      |__ CMakeLists.txt
      |__ cmake
      |__ |__ Dependencies.cmake
      |__ src
          |__ HelloWorld.cpp


Build your TriBITS project
----------------------------

Go to the build directory and type the following to configure your
project::

  cmake ../

The configure step will have created several files inside your build
directory, most notably it will have created necessary make files to
actually build your project.  The other file I will mention here is
the CMakeCache.txt which stores information about how the project was
configured. To build your project just type::

  make

you should see::

  [ 50\%] Building CXX object
   hello_package_dir/CMakeFiles/Hello-Executable-Name.dir/src/HelloWorld.cpp.o
  [100\%] Linking CXX executable Hello-Executable-Name.exe
  [100\%] Built target Hello-Executable-Name

now in build/hello\_package\_dir you will see an executable named
"Hello-Executable-Name" and if you run that executable you will see::

  $ ./hello_package_dir/Hello-Executable-Name.exe 
  Hello World!


Adding other targets
======================

Types of targets
------------------

Previously we had just one source file and we compiled it into one
executable.  In addition to executables we may also want to create
other targets such as libraries abd tests.  In the
hello_package_dir/src directory create the following files:
 
**hello_world_main.cpp**::

  #include <iostream>
  #include "hello_world_lib.hpp"
  int main() {
    std::cout << HelloWorld::gethelloworld() << "\n";
    return 0;
  }

**hello_world_lib.hpp**::

  #include <string>
  
  namespace HelloWorld { std::string gethelloworld(); }

**hello_world\_lib.cpp**::

  #include "hello_world_lib.hpp"
  std::string HelloWorld::gethelloworld()
  { return "Hello World!"; }

**hello_world_unit_tests.cpp**::

  #include <iostream>
  #include "hello_world_lib.hpp"
  
  int main() {
  
    bool success = true;
  
    const std::string rtn = HelloWorld::gethelloworld();
    std::cout << "HelloWorld::gethelloworld() = '"<<rtn<<"' == 'Hello World'? ";
    if (rtn == "Hello World!") {
       std::cout << "passed\n";
    }
    else {
      std::cout << "FAILED\n";
      success = false;
    }
  
    if (success) {
      std::cout << "All unit tests passed :-)\n";
    }
    else {
      std::cout << "At least one unit test failed :-(\n";
    }
  
  }

We will use these files to build an executalbe, a library, and tests.
Remember in the CMakeLists.txt file for the HelloPackage
(hello_package_dir/CMakeList.txt) we have the line::

  tribits_add_executable(Hello-Executable-Name NOEXEPREFIX SOURCES
  src/HelloWorld.cpp INSTALLABLE)

lets now modify that line to build an executable of the same name but
using hello_world_main.cpp instead of HelloWorld.cpp::

  tribits_add_executable(Hello-Executable-Name NOEXEPREFIX SOURCES
  src/hello_world_main.cpp INSTALLABLE)

to create a library we need to call tribits_add_library() and give it
a name, headers and sources.  add this the CMakeLists.txt::

  tribits_add_library(hello_world_lib HEADERS src/hello_world_lib.hpp
  SOURCES src/hello_world_lib.cpp)

we can also add tests.  You can add a test based on an executable you
have already specified for example::

  tribits_add_test(Hello-Executable-Name NOEXEPREFIX
  PASS_REGULAR_EXPRESSION "Hello World")

will run "Hello-Executable-Name" and verify that the output is "Hello
World".  You can also add a test and an executable att he same
time. for example::

  tribits_add_executable_and_test(unit_tests SOURCES
  src/hello_world_unit_tests.cpp PASS_REGULAR_EXPRESSION "All unit
  tests passed")

will create an executable named "unit_tests" from the source file
hello_world_unit_tests.cpp.  This executable will be used in a test
that will be marked as passing if the output of that executable is
"All unit tests passed".  After making these changes and additions to
the CMakeLists.txt file it should read::

  tribits_package(HelloPackage)

  tribits_add_library(hello_world_lib HEADERS src/hello_world_lib.hpp
   SOURCES src/hello_world_lib.cpp)

  tribits_add_executable(Hello-Executable-Name NOEXEPREFIX SOURCES
   hello_world_main.cpp INSTALLABLE)

  tribits_add_test(Hello-Executable-Name NOEXEPREFIX
   PASS_REGULAR_EXPRESSION "Hello World")

  tribits_add_executable_and_test(unit_tests SOURCES
   hello_world_unit_tests.cpp PASS_REGULAR_EXPRESSION "All unit tests
   passed")

  tribits_package_postprocess()

now reconfigure and rebuild in the build directory with::

  cmake ../
  make


What did we build?
====================

In the build directory there are many new files created by
TriBITS/CMake lets look at a few that are important for understanding
how TriBITS is building your project.

Build Targets
----------------

In the last section we built a library, an executable, and two tests.
Where are they? look in::

  build/hello_package_dir

among other things you will see::
  

  Hello-Executable-Name.exe
  HelloPackage_unit_tests.exe
  libhello_world_lib.a

by default, TriBITS will place the targets inside a directory with the
same name as the package directory.  If you have more than one package
then the files will be in separate directories::

    build
    |__ package_one
        |__ build_target_A
        |__ build_target_B
    |__  package_two
        |__  build_target_C
        |__  build_target_D

You can install the built targets to the default location
(/usr/local/bin) with::


  make install

You may want to install somewhere other than the default.  In this
case you want to set a CMamke variable called CMAKE_INSTALL_PREFIX. If
this is set then the files will be installed to the directory
specified.  For example in the top level CMakeLists set this variable
to a diecroyr called "Install" in the current source tree::

  set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/Install)

now clear the contents of the build directory and reconfigure, build,
and install the project with::

  cmake ../
  make install

Now you should see a directory called "Install" in the top level of the
project with contents::

  tree
  .
  |__ bin
  |   |__ Hello-Executable-Name.exe
  |__ include
  |   |__ hello_world_lib.hpp
  |__ lib
      |__ cmake
      |   |__ TribitsGreetings
      |       |__ TribitsGreetingsConfigVersion.cmake
      |__ libhello_world_lib.a


Summary
========

This tutorial has covered the most basic concepts in a TriBITS
project. A TriBITS project a collection of TriBITS packages and each
package defines its build targets (executables, tests, and libraries)
and source files.  A package also must define its dependencies. See
the TriBITS example project tutarial for a more complicate example of
a project and more detail about Tribits packages, TPLs, and
dependencies
