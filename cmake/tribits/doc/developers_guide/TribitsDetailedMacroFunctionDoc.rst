
TRIBITS_DEFINE_PACKAGE_DEPENDENCIES()
-------------------------------------

Define the dependenices for a given TriBITS SE package (i.e. a top-level
package or a subpackage).

Usage::

  TRIBITS_DEFINE_PACKAGE_DEPENDENCIES(
     [LIB_REQUIRED_PACKAGES <pkg1> <pkg2> ...]
     [LIB_OPTIONAL_PACKAGES <pkg1> <pkg2> ...]
     [TEST_REQUIRED_PACKAGES <pkg1> <pkg2> ...]
     [TEST_OPTIONAL_PACKAGES <pkg1> <pkg2> ...]
     [LIB_REQUIRED_TPLS <tpl1> <tpl2> ...]
     [LIB_OPTIONAL_TPLS <tpl1> <tpl2> ...]
     [TEST_REQUIRED_TPLS <tpl1> <tpl2> ...]
     [TEST_OPTIONAL_TPLS <tpl1> <tpl2> ...]
     [REGRESSION_EMAIL_LIST  <regression-email-address>
     [SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
       <spkg1_name>  <spkg1_dir>  <spkg1_classifications>  <spkg1_optreq>
       <spkg2_name>  <spkg2_dir>  <spkg2_classifications>  <spkg2_optreq>
       ...
       ]
     )

Every argument in this macro is optional.  The arguments that apply a package
itself are:

* **LIB_REQUIRED_PACKAGES:** List of upstream packages that must be enabled
  in order to build and use the libraries (or capabilities) in this
  package.

* **LIB_OPTIONAL_PACKAGES:** List of additional optional upstream packages
  that can be used in this package if enabled.  These upstream packages need
  not be enabled in order to use this package but not enabling one or more
  of these optional upstream packages will result in diminished capabilities
  of this package.

* **TEST_REQUIRED_PACKAGES:** List of additional upstream packages that must
  be enabled in order to build and/or run the tests and/or examples in this
  packages.  If any of these upstream packages is not enabled, then there
  will be no tests or examples defined or run for this package.

* **TEST_OPTIONAL_PACKAGES:** List of additional optional upstream packages
  that can be used by the tests in this package.  These upstream packages
  need not be enabled in order to run basic tests for this package.
  Typically, extra tests that depend on optional test packages involve
  integration testing of some type.

* **LIB_REQUIRED_TPLS:** List of upstream TPLs that must be enabled in order
  to build and use the libraries (or capabilities) in this package.

* **LIB_OPTIONAL_TPLS:** List of additional optional upstream TPLs that can
  be used in this package if enabled.  These upstream TPLs need not be
  enabled in order to use this package but not enabling one or more of these
  optional upstream TPLs will result in diminished capabilities of this
  package.

* **TEST_REQUIRED_TPLS:** List of additional upstream TPLs that must
  be enabled in order to build and/or run the tests and/or examples in this
  packages.  If any of these upstream TPLs is not enabled, then there
  will be no tests or examples defined or run for this package.

* **TEST_OPTIONAL_TPLS:** List of additional optional upstream TPLs
  that can be used by the tests in this package.  These upstream TPLs
  need not be enabled in order to run basic tests for this package.
  Typically, extra tests that depend on optional test TPLs involve
  integration testing of some type.

Only direct package dependenices need to be listed.  Indirect package
dependencies are automatically handled.  For example, if this SE package
directly depends on PKG2 which depends on PKG1 (but this SE package does not
directly depend on anything in PKG1) then this package only needs to list a
dependency on PKG2, not PKG1.  The dependnecy on PKG1 will be taken care of
automatically by the TriBITS dependency tracking system.

However, currently, all TPL dependendies must be listed, even the indirect
ones.  This is a requirement that will be dropped in the future.

The packages listed in LIB_REQUIRED_PACKAGES are implicitly also
dependenices in TEST_REQUIRED_PACKAGES.  Likewise LIB_OPTIONAL_PACKAGES are
implicitly also dependenices in TEST_OPTIONAL_PACKAGES.  Same goes for TPL
dependencies.

The dependencies within a single list do not need to be listed in any order.
For example if PKG2 depends on PKG1, and this given SE package depends on
both, one can list "LIB_REQUIRED_PACKAGES PKG2 PKG1" or
"LIB_REQUIRED_PACKAGES PKG1 PKG2".  Likewise the listing of TPLs order is
not important.

If some upstream packages are allowed to be missing, this can be specified
by calling the macro `TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES()`_.

A top-level package can also have subpackages.  In this case, the following
varible must be set:

* **SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS:** 2D array with rows listing
  the subpackages and the columns:

  * **SUBPACKAGE:** The name of the subpackage <spkg_name>.  The full SE
    package name is "${PARENT_PACKAGE_NAME}<spkg_name>".  The full SE
    package name is what is used in listing dependenices in other SE
    packages.

  * **DIRS:** The subdirectory <spkg_dir> relative to the parent package's
    base directory.  All of the contents of the subpackage should be under
    this subdirectory.  This is assumed by the TriBITS testing support
    software when mapping modified files to SE packages that need to be
    tested.

  * **CLASSIFICATIONS***: The test group PT, ST, EX and the maturity level
    EP, RS, PG, PM, GRS, GPG, GPM, and UM, separated by a coma ',' with no
    spaces in between (e.g. "PT,GPM").  These have exactly the name meaning
    as for full packages (see
    `TRIBITS_DEFINE_REPOSITORY_PACKAGES_DIRS_CLASSIFICATIONS()`_).

  * **OPTREQ:** Determines if the outer parent package has an OPTIONAL or
    REQUIRED dependence on this subpackage.

Other variables that this macro handles:

* **REGRESSION_EMAIL_LIST:** The email list that is used to send CDash error
  messages.  If this is missing, then the email list that CDash errors go to
  is determined by other means (see ???).

NOTE: All this macro really does is to just define the variables:

* LIB_REQUIRED_DEP_PACKAGES
* LIB_OPTIONAL_DEP_PACKAGES
* TEST_REQUIRED_DEP_PACKAGES
* TEST_OPTIONAL_DEP_PACKAGES
* LIB_REQUIRED_DEP_TPLS
* LIB_OPTIONAL_DEP_TPLS
* TEST_REQUIRED_DEP_TPLS
* TEST_OPTIONAL_DEP_TPLS
* REGRESSION_EMAIL_LIST
* SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS

which are then read by the TriBITS cmake code to build the package
dependency graph.  The advantage of using this macro instead of just
directly setting the varibles is that you only need to list the dependencies
you have.  Otherwise, you need to set all of these varibles, even those that
are empty.  This is a error checking property of the TriBITS system to avoid
misspelling the names of these variables.

TRIBITS_DEFINE_REPOSITORY_PACKAGES_DIRS_CLASSIFICATIONS()
---------------------------------------------------------

Define the set of packages for a given TriBIT repo.  This macro is typically
called from inside of a PackagesList.cmake file for a given TriBITS repo.

Usage::

   TRIBITS_DEFINE_REPOSITORY_PACKAGES_DIRS_CLASSIFICATIONS(
      <pkg0>  <pkg0_dir>  <pkg0_classifications>
      ...
      <pkgnm1>  <pkgnm1_dir>  <pkgnm1_classifications>
      )

This macro sets up a 2D array of NumPackages by NumColumns listing out the
packages for a TriBITS repository.  Each row (with 3 entries) specifies a
package which contains the three columns:

* **PACKAGE**: The name of the TriBITS package.  This name must be unique
  across all other TriBITS packages in this or any other TriBITS repo that
  might be combined into a single TriBITS project meta-build.  The name
  should be a valid identifier (e.g. matches the regex
  ``[a-zA-Z_][a-zA-Z0-9_]*``).

* **DIR**: The relative directory for the package.  This is relative to the
  TriBITS repository base directory.  Under this directory will be a
  package-specific 'cmake/' directory with file 'cmake/Dependencies.cmake'
  and a base-level CMakeLists.txt file.  The entire contents of the package
  including all of the source code and all of the tests should be contained
  under this directory.  The TriBITS testing infrastructure relies on the
  mapping of changed files to these base directories when deciding what
  packages are modified and need to be retested (along with downstream
  packages).

* **CLASSIFICATION**: Gives the testing group PT, ST, EX and
  the maturity level EP, RS, PG, PM, GRS, GPG, GPM, UM.  These are seprated
  by a coma with no space in between such as "RS,PT" for a "Research
  Stable", "Primary Tested" package.  No spaces are allowed so that CMake
  treats this a one field in the array.  The maturity level can be left off
  in which case it is assumed to be UM for "Unspecified Maturity".

 NOTE: This macro just sets the varaible
 ${REPOSITORY_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS in the current
 scope.  The advantages of using this macro instead of directly setting this
 varible include:

 * Asserts that REPOSITORY_NAME is defined and set

 * Avoids having to hard-code the assumed repository name
   ${REPOSITORY_NAME}.  This provides more flexibility for how other TriBITS
   project name a given TriBITS repo (i.e. the name of repo subdirs).

 * Avoid mispelling the name of the varible
   ${REPOSITORY_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS.  If you misspell
   the name of the macro, it is an immediate error in CMake.

TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES()
-----------------------------------------

Macro used in Dependencies.cmake files to allow some upstream dependent packages
to be missing.

Usage::

  TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(<pack_1> <pack_2> ...)

If the missing upstream SE package <pack_i> is optional, then the effect
will be to simply ignore the missing package and remove it from the
dependency list.  However, if the missing upstream SE package <pack_i> is
required, then in addition to ignoring the missing package, the current SE
(sub)package will also ee hard disabled,
i.e. ${PROJECT_NAME}_ENABLE_{CURRENT_PACKAGE}=OFF.

This function is typically used in packages in external TriBITS repos that
are depend on other packages in other exteral TriBITS repos that might be
missing.

NOTE: Using this function effectively turns off error checking for
misspelled package names so it is important to only use it when it
absolutely is needed.

.. @FUNCTION: TRIBITS_TPL_DECLARE_LIBRARIES -
.. @MACRO: TRIBITS_PACKAGE_DECL() -
.. @MACRO: TRIBITS_PACKAGE_DEF() -
.. @MACRO: TRIBITS_ADD_TEST_DIRECTORIES() -
.. @MACRO: TRIBITS_ADD_EXAMPLE_DIRECTORIES() -
.. @MACRO: TRIBITS_PACKAGE_POSTPROCESS() -
.. @MACRO: TRIBITS_PROCESS_SUBPACKAGES() -
.. @MACRO: TRIBITS_DEFINE_REPOSITORY_TPLS_FINDMODS_CLASSIFICATIONS() -
.. @FUNCTION: TRIBITS_SET_ST_FOR_DEV_MODE() -

.. @FUNCTION: TRIBITS_ADD_LIBRARY() -

.. @FUNCTION: TRIBITS_ADD_EXECUTABLE() -
