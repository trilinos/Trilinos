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
guidelines on how to split them up.  Once you have your project
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
to build them separately, then make them each their own package.

libA ==> libB ==> libC 
libA depends on libB which depends on libC.
Again consider making them all their own packages

libA <==> libB ==> libC

libA and libB depend on each other and libB depends on libC.  Make two
packages, one that builds A and B, and one that builds C


Components of a Tribits Package
-------------------------------

The goal of a tribits package is straight forward.  Take some source
files and build a target and the tests for that target.  This will
probably be the main function of the CMakeLists.txt file in top level
package directory.  In that CmakeLists.txt file you need to have a few
commands if you follow this list, you will be well on your way to a
Tribits package.

1. tribits_package(<package_name>) - must be called before you add any
   targets you want to build in this package.  It is a good idea to
   make it the first line in the top level CMakLists.txt file
#. Add subdirectories 
#. Identify your source files 
#. Add targets to be built using tribits_add_library(),
   tribits_add_executable(), and tribits_add_test()
#. tribits_package_postprocess() - do not add any new targets or
   include new directories after this command.  This should be the
   last line in the file

You will also need to define any dependencies this package may have on
other packages in the project.  This is done in a File called
Dependencies.camke in a cmake directory.  All that is required in this
file is a call to tribits_package_define_dependencies().  Even if the
package does not depand on another you still need to have this.


Example of simple package CMakeLists.txt files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a simple example of a CMakelists.txt file for the top level of
a package.  In this example the targets are being specified in this
top level file::

  #
  # A) Define the package
  #
  tribits_package( ExamplePackageName )

  #
  # B) Set up package-specific options
  #
  set(example_srcs example.cpp)
  set(example_headers example.h)

  #
  # C) Add the libraries, tests, and examples
  #
  tribits_add_library(library_name SOURCES ${example_srcs} HEADERS ${example_headers})

  #
  # D) Do standard post processing
  #
  tribits_package_postprocess()

You may have source files grouped together into files used to build
the target(s) and files used to build tests of the targets.  Suppose
you have a src/ directory and a test/ directory.  Then you may want
your top level file to just include those subdirectories and let the
CMakeLists.txt in each subdirectory do B) and C) from the above
example.  So you would have something like::

  #  
  # A) Define the package
  #
  tribits_package( SimpleCxx  ENABLE_SHADOWING_WARNINGS  CLEANED )

  #
  # B) Add the libraries, tests, and examples
  #
  add_subdirectory(src)
  tribits_add_test_directories(test)
  
  #
  # C) Do standard post processing
  #
  tribits_package_postprocess()

In this case you would also have a CMakeLists.txt file in the src/
directory that looks like.  Note there are no calls to
tribits_package() or tribits_package_postprocess() in this lower level
CMakeLists file.  These functions must be called in the top level
CMakLists file but not in any others ::

  #
  # A) Set up package-specific options
  #
  set(example_srcs example.cpp)
  set(example_headers example.h)

  #
  # B) Add the libraries, tests, and examples
  #
  tribits_add_library(library_name SOURCES ${example_srcs} HEADERS ${example_headers})


Examples of Dependencies.cmake files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In Addition to the CMakeLists files described above, you will also
need to tell Tribits about the dependencies this package has on other
packages.  This is done through a call to
tribits_package_define_dependencies() in a package's
Dependencies.camke file. If there are no dependencies then this file
will contain a call to tribits_package_define_dependencies() with no
arguments::

  # Dependencies.camke 
  tribits_package_define_dependencies()

Suppose this package has an optional and a required dependency on
other packagages then the call would look something like::

  tribits_package_define_dependencies(
    LIB_REQUIRED_TPLS  name_of_required_package
    LIB_OPTIONAL_TPLS  name_of_otional_package
    )


Linking libraries
-----------------

If you are building more than one target in a package you may be
accustomed to calling target_link_library() to tell camke to build the
target against some library.  If the library is being built in the
same package as your target, you do not need to do this because
Tribits will automatically link against any libraries built in the same
package.  Additionally you do not need to do this for targets built in
other packages because tribits will link against any libraries built
in packages that the current package depends on.  If you are calling
target_link_library() then it is either redundant, or it indicates
there is a dependency that needs to be defined in the
Dependencies.camke file of your package.  


TriBITS Projects
================

A tribits project is a collection of Tribits packages.  If you have
Tribits packages defined then you can put them together in a Tribits
project.  In order to do this you need to define some things at the
project level.  The top level project will have a CMakeLists file as
well as a few .camke files to define packages that are in the project,
TPLs that he project may depend on, software version, and other
project information.  In the tribits project directory you need to
have the following files:

1. *CMakeLists.txt* - top level CMakeLists for the project.  here you
   will initialize your Tribits project, define which packages will be
   built by default, and define some setting for your project
#. *PackagesList.cmake* - Tells tribits which packages are part of
   this projec t and where to find them
#. *TPLsList.cmke* - Tells tribits which TPLs that packages my depend
   on and how to find them
#. *ProjectName.cmake* - defines the project name and possibly some
   other project level settings
#. *Version.cmake* - set the version of the software being built
#. *Package Directories* - A directory for each package that contains
   everything necessary for a Tribits package described above. Often
   Projects will have a packages directory that contains all of the
   individual package directories in the project

An example directory structure could look like this::

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

  set(PROJECT_NAME Your_Project_Name)

you may also want to do other stuff in here (ask Ross What?/why?)


PackageList.cmake
-----------------

Here you will define all of the packages in the project with a name, a
location, and some options.  This is done with a call to
tribits_repository_define_packages().  For example::

  tribits_repository_define_packages(
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
- *ST (Secondary Tested)* - This code is very important to the project
  but will not necessarily halt development if it breaks.  Consider
  making a package as ST if it depends on difficult to install TPLs or
  TPLs that are not available no all deveopment platforms.
- *EX (Experimental)* - This code is unstable and difficult to
  maintain.  It is not portible or not important enough to be tested
  at the same level as other code

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
tribits_disable_package_on_platforms() if you have packages that you
know will not work on certain platform and you want to disable them on
those platforms::

  tribits_disable_package_on_platforms(package_name
    platform_0 platform_1 ...
  )

This will cause package_name to be disabled by default on platform_0,
platform_1, ...


TPLsList.cmake
--------------

Here you will define all of the tpls in the project.  The function
call is very similar to defining packages above. You do this by
calling tribits_repository_define_tpls() with a name, a path to a tpl
cmake find module, and a classification for each tpl.  For example::

 tribits_repository_define_tpls(
   MPI  "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/FindTPLMPI.cmake"  PT
   SomeTplA   "cmake/tpls/"         PT
   SomeTplB   "cmkae/tpls/"         PT
 )

In this example you can see a path to the tribits findTPLMPI.camke
module which will find mpi on the system as well as examples of relative
paths to a directory where the project has defined some cmkae find
modules for the required TPLs.  Each line will contain information
about one tpl.  The first entry is the name of the tpl, the second is
the path to the tpl find module, and the third is for tpl
classification.  Tpl classifications you can specify are:

Testing category (Required)

- *PT (Primary Tested)* - This tpl is essential to developer
  productivity and would adversly effect customers if broken.
- *ST (Secondary Tested)* - This tpl is important to the project but
  mat be difficult to install or the TPL is not available no all
  deveopment platforms.
- *EX (Experimental)* - TPL is experimental, unstable and/or difficult to
  maintain.

The recommendation is to list all TPLs as "PT"


Version.cmake
-------------

This file just contains version information for the code example::

  set(${REPOSITORY_NAME}_VERSION 1.1)
  set(${REPOSITORY_NAME}_MAJOR_VERSION 01)
  set(${REPOSITORY_NAME}_MAJOR_MINOR_VERSION 010100)
  set(${REPOSITORY_NAME}_VERSION_STRING "1.1 (Dev)")
  set(${REPOSITORY_NAME}_ENABLE_DEVELOPMENT_MODE_DEFAULT ON) # Change to 'OFF' for a release

CMakeList.txt
-------------

Here you will tell tribits some basic information it need to build as
atribits project.  You need to specify where Tribits is located on the
system (many projects choose to snapshot tribits into their
repository) You will also be able to specify if packages are turned
on/off by default. Here is the order of commandsthat you should have
in this project level CMakeLists file:

1. cmake_minimum_version() - set the minimum version of cmake required
   for this project o build.  If you try and run with a lower version
   then there will be an error. You cannot specify a version lower than
   3.23.0
#. Include ProjectName.cmake and call project() with argument PROJECT_NAME
#. specify the directory to tribits and include TriBITS.cmake
#. specify which packages are turned on/off by default
#. call tribits_project()

Here is an example of a project CMakeLists::

  # Define your minimum CMake version
  cmake_minimum_required(VERSION 3.23.0 FATAL_ERROR)

  # Define your project name and set up major project options
  include("${CMAKE_CURRENT_SOURCE_DIR}/ProjectName.cmake")
  project(${PROJECT_NAME} NONE)

  # Pull in the TriBITS system and execute
  set(${PROJECT_NAME}_TRIBITS_DIR
     "${CMAKE_CURRENT_LIST_DIR}/../.."  CACHE  STRING
    "TriBITS base directory (default assumes in TriBITS source tree)")
  include("${${PROJECT_NAME}_TRIBITS_DIR}/TriBITS.cmake")

  # Do all of the processing for this Tribits project
  tribits_project()



TriBITS Repositories
--------------------

In the simplest case your project will use packages that are in the
same repository as your project but his does not have to be the case.
Packages may be defined in a TriBITS repository that can then be used
by your project by adding the extra repositories in a file
"<projectDir>/cmake/ExtraRepositoriesList.cmake" which sets up the
repositories with a call to::

  tribits_project_define_extra_repositories()

such as ::

  tribits_project_define_extra_repositories(
   <REPO_NAME> <REPO_DIR> <REPO_VCTYPE> <REPO_URL> <REPO_PACKSTAT> <REPO_CLASSIFICATION>
    ...
  )

where each line is one repo and

- **REPO_NAME** is the name of the repo

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

