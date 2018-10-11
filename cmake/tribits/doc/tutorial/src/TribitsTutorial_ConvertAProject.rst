=====================================
 Convert your project to use TriBITS
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
============

This document will help you convert your project that is currently
being built with cmake to build with Tribits instead.  I will be
introducing a subset of tribits features in the hopes of keeping this
document easy to follow.  There is a wealth of detailed documentation
for tribits.  Once you have some ecperience with Tribits, I recommend
that you use the "TriBITS Developers Guide" as a reference.  Tribits
is desigend to take advantage of modularity in projects and excels at
managing complex dependencies between software packages.  Tribits
projects should be structured a certain in order to work best. For the
purposes of this tutorial, the most basic unit of Tribits is a Tribits
package.  A Tribits project is a collection of Tribits packages.
(This is not the whole story but for now just consider packages and
projects)


TriBITS Packages
================

The first thing to do when converting a projec to tribits is to think
about the structure of the project and to figure out the best way to
break it into packages.  In this section, I hope to give you a better
idea of what exactly constitutes a tribits package and give you some
guidlines on how to split them up.  Once you have your project
building with Tribits, you will be able to build different
configurations by simply turning packages on or off.

Identify your packages
----------------------

The simplest way to think of a package is as a collection of source
files being built into one (or more, but ideally one) library or
executable and a set of tests.  If you have componetns of your project
that fit this description then consider making into a package.  A
package does not have to build just one library or executable and if
there are targets that are always built together then make them all
part of the same package.  For example consider the following
situations where you are building 3 different libraries.

examples
~~~~~~~~

if libA libB, and liC do not depend on each other and it makes sense
to build them separately, then make them each ther own package.

libA ==> libB ==> libC 
libA depends on libB which depends on libC.
Again consider making them all their own packages

libA <==> libB ==> libC

libA and libB depend on eachother and libB depends on libC.  Make two
packages, one that builds A and B, and one that builds C


Components of a Tribits Package
-------------------------------

The goal of a tribits package is straight forward.  Take some source
files and build a target and the tests for that target.  This will
probably be the main function of the CMakeLists.txt file in top level
package directory.  In that CmakeLists.txt file you need to have a few
commands if you follow this list, you will be well on your way to a
Tribits package.

1. TRIBITS_PACKAGE(<package_name>) - must be called before you add any
   targets you want to build in this package.  It is a good idea to
   make it the first line in the top level CMakLists.txt file
#. Add subdirectories 
#. Identify your source files 
#. Add targets to be built using TRIBITS_ADD_LIBRARY(),
   TRIBITS_ADD_EXECUTABLE(), and TRIBITS_ADD_TEST()
#. TRIBITS_PACKAGE_POSTPROCESS() - do not add any new targets or
   include new directories after this command.  This should be the
   last line in the file

You will also need to define any dependencies this package may have on
other packages in the project.  This is done in a File called
Dependencies.camke in a cmake directory.  All that is required in this
file is a call to TRIBITS_PACKAGE_DEFINE_DEPENDENCIES().  Even if the
package does not depand on another you still need to have this.


Example of simple package CMakeLists.txt files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a simple example of a CMakelists.txt file for the top level of
a package.  In this example the targets are being specified in this
top level file::

  #
  # A) Define the package
  #
  TRIBITS_PACKAGE( ExamplePackageName )

  #
  # B) Set up package-specific options
  #
  SET(example_srcs example.cpp)
  SET(example_headers example.h)

  #
  # C) Add the libraries, tests, and examples
  #
  TRIBITS_ADD_LIBRARY(library_name SOURCES ${example_srcs} HEADERS ${example_headers})

  #
  # D) Do standard post processing
  #
  TRIBITS_PACKAGE_POSTPROCESS()

You may have source files grouped together into files used to build
the target(s) and files used to build tests of the targets.  Suppose
you have a src/ directory and a test/ directory.  Then you may want
your top level file to just include those subdirectories and let the
CMakeLists.txt in each subdirectory do B) and C) from the above
example.  So you would have something like::

  #  
  # A) Define the package
  #
  TRIBITS_PACKAGE( SimpleCxx  ENABLE_SHADOWING_WARNINGS  CLEANED )

  #
  # B) Add the libraries, tests, and examples
  #
  ADD_SUBDIRECTORY(src)
  TRIBITS_ADD_TEST_DIRECTORIES(test)
  
  #
  # C) Do standard post processing
  #
  TRIBITS_PACKAGE_POSTPROCESS()

In this case you would also have a CMakeLists.txt file in the src/
directory that looks like.  Note there are no calls to
TRIBITS_PACKAGE() or TRIBITS_PACKAGE_POSTPROCESS() in this lower level
CMakeLists file.  These functions must be called in the top level
CMakLists file but not in any others ::

  #
  # A) Set up package-specific options
  #
  SET(example_srcs example.cpp)
  SET(example_headers example.h)

  #
  # B) Add the libraries, tests, and examples
  #
  TRIBITS_ADD_LIBRARY(library_name SOURCES ${example_srcs} HEADERS ${example_headers})


Examples of Dependencies.cmake files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In Addition to the CMakeLists files described above, you will also
need to tell Tribits about the dependencies this package has on other
packages.  This is done through a call to
TRIBITS_PACKAGE_DEFINE_DEPENDENCIES() in a package's
Dependencies.camke file. If there are no dependencies then this file
will contain a call to TRIBITS_PACKAGE_DEFINE_DEPENDENCIES() with no
arguments::

  # Dependencies.camke 
  TRIBITS_PACKAGE_DEFINE_DEPENDENCIES()

Suppose this package has an optional and a required dependency on
other packagages then the call would look something like::

  TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
    LIB_REQUIRED_TPLS  name_of_required_package
    LIB_OPTIONAL_TPLS  name_of_otional_package
    )


Linking libraries
-----------------

If you are building more than one target in a package you may be
accustomed to calling TARGET_LINK_LIBRARY() to tell camke to build the
target against some library.  If the library is being built in the
same package as your target, you do not need to do this because
Tribits will automatically link against any libraies built in the same
package.  Additionally you do not need to do this for targets built in
other packages because tribits will link against any libraries built
in packages that the current package depends on.  If you are calling
TARGET_LINK_LIBRARY() then it is either redundant, or it indicates
there is a dependancy that needs to be defined in the
Dependencies.camke file of your package.  


TriBITS Projects
================

A tribits project is a collection of Tribits packages.  If you have
Tribits packages defined then you can put them together in a Tribits
project.  In order to do this you need to define some things at the
project level.  The top level project will have a CMakeLists file as
well as a few .camke files to define packages that are in the project,
TPLs that he project may depend on, software version, and other
project infrormation.  In the tribits project directory you need to
have the following files:

1. *CMakeLists.txt* - top level CMakeLists for the project.  here you
   will initialize your Tribits project, define which packages will be
   built by defult, and define some setting for your project
#. *PackagesList.cmake* - Tells tribits which packages are part of
   this projec t and where to find them
#. *TPLsList.cmke* - Tells tribits which TPLs that packages my depend
   on and how to find them
#. *ProjectName.cmake* - defines the projet name and possibly some
   other project level settings
#. *Version.cmake* - set the version of the software being built
#. *Package Directories* - A directory for each package that contains
   everything nesseesary for a Tribits package described above. Often
   Projects will have a packages directory that contains all of the
   individual package directories in the project

An example direcory structure could look like this::

  my_tribits_project/
  |__CMakeLists.txt
  |__PackageLists.cmake
  |__TPLSsList.cmake
  |__ProjectName.cmake
  |__Version.cmkae
  |__packages/
     |__my_first_tribists_package/
     |  |__CMakeList.txt
     |  |__cmake/
     |  |  |__Dependencies.cmake
     |  |__src/
     |  |  |__ a bunch of source files
     |  |__test/ 
     |     |__ more source files to build tests
     |
     |__my_second_tribists_package/
        |__CMakeList.txt
        |__cmake/
        |  |__Dependencies.cmake
        |__src/
        |  |__ a bunch of source files
        |__test/ 
           |__ more source files to build tests

next we will go over what each of the project level cmake files need
to contain and some of the option you may want to use for your project

ProjectName.cmake
-----------------

This file simply needs to set the variable PROJECT_NAME.  For
example::

  SET(PROJECT_NAME Your_Project_Name)

you may also want to do other stuff in here (ask Ross What?/why?)


PackageList.cmake
-----------------

Here you will define all of the packages in the project with a name, a
location, and some options.  This is done with a call to
TRIBITS_REPOSITORY_DEFINE_PACKAGES().  For example::

  TRIBITS_REPOSITORY_DEFINE_PACKAGES(
    TriBitsPackageA      packages/package_a         PT
    TribitsPackageB      packages/package_b         PT
  )

Each line will contain information about one package.  The foirst
entry is the name of the package, the second is the path to the
package, and the third is for package classification.  Package
classifications you can specify are:

Testing category (Required)

- *PT (Primary Tested)* - The code is the highest priority to keep
  working.  This package package is essential to developer
  productivity and would adversly effect customers if broken.
- *ST (Seconday Tested)* - This code is very important to the project
  but will not nessesarily halt developement if it breaks.  Consider
  making a package as ST if it depends on difficult to install TPLs or
  TPLs that are not available no all deveopment platforms.
- *EX (Experimental)* - This code is unstable and difficult to
  maintain.  It is not portible or not important enough to be tested
  at teh same level as other code

Package Maturity (Defaults to UM) (ask Ross what this is? to include?)
- EP
- RS
- PG
- PM
- GRS
- GPG
- GPM
- UM

In this file you may also choose to call
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS() if you have packages that you
know will not work on certain platform and you want to disable them on
those platforms::

  TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(package_name
    platform_0 platform_1 ...
  )

This will cause package_name to be disabled by default on platform_0,
platform_1, ...


TPLsList.cmake
--------------

Here you will define all of the tpls in the project.  The function
call is very similar to defining packages above. You do this by
calling TRIBITS_REPOSITORY_DEFINE_TPLs() with a name, a path to a tpl
cmake find module, and a classification for each tpl.  For example::

 TRIBITS_REPOSITORY_DEFINE_TPLS(
   MPI  "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/FindTPLMPI.cmake"  PT
   SomeTplA   "cmake/tpls/"         PT
   SomeTplB   "cmkae/tpls/"         PT
 )

In this example you can see a path to the tribits findTPLMPI.camke
moduel which will find mpi on the sytem as well as exaples of reletive
paths to a directory where the project has defiend some cmkae find
modules for the required TPLs.  Each line will contain information
about one tpl.  The first entry is the name of the tpl, the second is
the path to the tpl find module, and the third is for tpl
classification.  Tpl classifications you can specify are:

Testing category (Required)

- *PT (Primary Tested)* - This tpl is essential to developer
  productivity and would adversly effect customers if broken.
- *ST (Seconday Tested)* - This tpl is important to the project but
  mat be difficult to install or the TPL is not available no all
  deveopment platforms.
- *EX (Experimental)* - TPL is experimental, unstable and/or difficult to
  maintain.

The recommendation is to list all TPLs as "PT"


Version.cmake
-------------

This file just contains version infromation for the code example::

  SET(${REPOSITORY_NAME}_VERSION 1.1)
  SET(${REPOSITORY_NAME}_MAJOR_VERSION 01)
  SET(${REPOSITORY_NAME}_MAJOR_MINOR_VERSION 010100)
  SET(${REPOSITORY_NAME}_VERSION_STRING "1.1 (Dev)")
  SET(${REPOSITORY_NAME}_ENABLE_DEVELOPMENT_MODE_DEFAULT ON) # Change to 'OFF' for a release

CMakeList.txt
-------------

Here you will tell tribits some basic information it need to build as
atribits project.  You need to specify where Tribits is located on the
system (many projects choose to snapshot tribits into their
repository) You will also be able to specify if packages are turned
on/off by default. Here is the order of commandsthat you should have
in this project level CMakeLists file:

1. CMAKE_MINIMUM_VERSION() - set the minimum version of cmake required
   for this project o build.  If you try and run with a lower version
   then there wil be an error. You cannot specify a version lower than
   3.10.0
#. Include ProjectNmae.cmake and call PROJECT() with argument PROJECT_NAME
#. specify the directory to tribits and include TriBITS.cmake
#. specify which packages are turned on/off by default
#. call TRIBITS_PROJECT()

Here is an examlpe of a project CMakeLists::

  # Deefine your minimum CMake version
  CMAKE_MINIMUM_REQUIRED(VERSION 3.10.0 FATAL_ERROR)


  # Define your project name and set up major project options
  INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/ProjectName.cmake")
  PROJECT(${PROJECT_NAME} NONE)


  # Pull in the TriBITS system and execute
  SET(${PROJECT_NAME}_TRIBITS_DIR
     "${CMAKE_CURRENT_LIST_DIR}/../.."  CACHE  STRING
    "TriBITS base directory (default assumes in TriBITS source tree)")
  INCLUDE("${${PROJECT_NAME}_TRIBITS_DIR}/TriBITS.cmake")


  # Do all of the processing for this Tribits project
  TRIBITS_PROJECT()



TriBITS Repositories
--------------------

In the simplest case your project will use packages that are in the
same repository as your project but his does not have to be the case.
Packages may be defined in a TriBITS repository that can then be used
by your project by adding the extra repositories in a file
"<projectDir>/cmake/ExtraRepositoriesList.cmake" which sets up the
repositories with a call to::

  TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES()

such as ::

  TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES(
   <REPO_NAME> <REPO_DIR> <REPO_VCTYPE> <REPO_URL> <REPO_PACKSTAT> <REPO_CLASSIFICATION>
    ...
  )

where each line is one repo and

- **REPO_NAME** is the name ofthe repo

- **REPO_DIR** is the relative path to the repo (asssumed to be
  ./REPO_NAME/ if it is blank)

- **REPO_VCTYPE** the type of version control used for this repo (must
  be: "GIT", "SVN", or "")

- **REPO_URL** the url to the repo (can be "". if REPO_VCTYPE is ""
  then this must be "")

- **REPO_PACKSTAT** indicates if this is a TriBITS repository with
  packages or not.  "NOPACKAGES" means this repo does not contain
  TriBITS packages.  "HASPACKAGES, PRE" means this repo does have
  packages and you would like them to be processed before the packages
  listed by your project because packages in your project depend on
  the packages in this repo.  "HASPACKAGES, POST" means this repo has
  packages and you would like for them to be processed after the
  packages listed in your project because the packages in this repo
  depend on the packages in your project

- **REPO_CLASSIFICATION** indicates when this repo should be included
  for testing must be: "Continuous", "Nightly", or "Experimental"

