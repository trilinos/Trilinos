=====================================
Tribits Dependencies 
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


TriBITS Example Project
========================

The previous tutorial had you create all the files nesseasary for a
very simple tribits project.  This time we will look at an example of
a more complicated example project that is included when you clone the
tribits repository. If you go to::

  TriBITS/tribits/examples/

you will see several examples including one very similar to what we
constructed in the last tutorial.  In this tutorial we will be using
the TribitsExampleProject.  This project has multiple packages as well
as TPL dependencies so we will be able to see how tribits deals with
both types of dependencies. In TribitsExampleProject, look at::

  PackagesList.cmake
  TPLsList.cmake

and you will see that we have a several packages defined for the
project in PackageList.cmake and a couple TPLs defined in
TPLsList.cmake.  Recall that the package definitions looks like the
following::

  tribits_repository_define_packages(
    NameOfPackageA     location/of/packageA        <options>
    NameOfPackageB     path/to/packageB            <options>
  )

you must specify the name of the package and where it is located.
Each package must contain a CmakeLists.cmake file and a
Dependencies.cmake file.  The CmakeLists file defines targets for this
package and must begin with a call to tribits_package() and end with
a call to tribits_package_postprocess().  

Tribits TPLs are defined similarly in TPLsList.cmake::

  tribits_repository_define_tpls(
    NameOfTPL-1     location/of/TPL-1        <options>
    NameOfTPL-2     path/to/TPL-2            <options>
  )

Tribits Packages
-----------------

A Package is software that will be built at the same time as the rest
of your project.  It is a modular piece of your software that you may
want to enable or disable.  A tribits package can have dependencies on
other packages or on TPLs.

Enabling Packages
+++++++++++++++++++++++

If you enable a package Then TriBITS will enable
all of that projects dependencies.  You do not have to explicitly
enable the whole chain of dependencies.  Just enable the packages that
you directly depend on and TriBITS will make sure that everything you
need for those packages is enabled as well

Package Options
++++++++++++++++

For <options> you will see PT, ST, or EX.  

- *PT (Primary Tested)* - This tpl is essential to developer
  productivity and would adversly effect customers if broken.
- *ST (Secondary Tested)* - This tpl is important to the project but
  may be difficult to install or the TPL is not available on all
  development platforms.
- *EX (Experimental)* - TPL is experimental, unstable and/or difficult to
  maintain.

Most packages belong in the PT category and you will be fine if you
only use that categroy for your packages to begin with


Important things to know when listing packages
++++++++++++++++++++++++++++++++++++++++++++++++

**Important - Read this section**

* TriBITS packages must be listed in dependency order in
  the PackageList.cmake file.  A package must be listed after all
  packages on which it depends. 
* Circular dependencies are not allowed in TriBITS.  If you have
  packages with circular dependencies then they really should be one
  package



TriBITS TPLs
-------------

A TPL is software you depend on that is built and installed prior to
building your project from source.  For example, many projects dependo
nthe boost libraries but it is rare atha anyone wants to build them as
part of building their project.  Projects assume that they areinstalled
already and we just need to point tribits to them.  TPLs are assumed to
be self contained and will not trigger additional dependencies in your
build.  TPLs are always leaf nodes in the dependency tree.
