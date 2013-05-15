# @HEADER
# ************************************************************************
#
#            TriBITS: Tribial Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
# @HEADER

IF (EXTRA_REPO_INCLUDE_MISSING_OPTIONAL_DEP_PACKAGE)
   SET(MISSING_OPTIONAL_DEP_PACKAGE MissingUpstreamPackage)
  TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(MissingUpstreamPackage)
ELSE()
   SET(MISSING_OPTIONAL_DEP_PACKAGE)
ENDIF()

IF (EXTRA_REPO_INCLUDE_MISSING_REQUIRED_DEP_PACKAGE)
   SET(MISSING_REQUIRED_DEP_PACKAGE MissingUpstreamPackage)
  TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(MissingUpstreamPackage)
ELSE()
   SET(MISSING_REQUIRED_DEP_PACKAGE)
ENDIF()

SET(LIB_REQUIRED_DEP_PACKAGES Teuchos ${MISSING_REQUIRED_DEP_PACKAGE})
SET(LIB_OPTIONAL_DEP_PACKAGES ${MISSING_OPTIONAL_DEP_PACKAGE})
SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES)
SET(LIB_REQUIRED_DEP_TPLS Boost)
SET(LIB_OPTIONAL_DEP_TPLS)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)

SET(REGRESSION_EMAIL_LIST ex2-package1-override@some.ornl.gov)

