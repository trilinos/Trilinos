# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

tribits_package_define_dependencies(
  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    CoreLibs      src                  PT  REQUIRED
    GoodStuff     src/good_stuff       ST  OPTIONAL
    CrazyStuff    src/crazy_stuff      EX  OPTIONAL
    Epetra        adapters/epetra      PT  OPTIONAL
    EpetraExt     adapters/epetraext   PT  OPTIONAL
    Tpetra        adapters/tpetra      PT  OPTIONAL
  LIB_OPTIONAL_PACKAGES MissingPackage
  REGRESSION_EMAIL_LIST thyra-boneheads@gmail.com
  )

# NOTE: The above subpackages automatically become required and optional LIB
# package dependencies (prefixed by 'Thyra') added to the variables below.
# There is no need to add them again!

tribits_allow_missing_external_packages(MissingPackage)

