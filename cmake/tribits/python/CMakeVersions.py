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

# This file designates the officially supported versions of CMake for Trilinos.
# It is updated periodically as new builds of CMake are available.
#
# The variables defined at the bottom of this file are as follows:
#
# cmake_version_min:
#   Minimum required -- updated manually, synched with cmake_minimum_required
#   Should match Trilinos/CMakeLists.txt cmake_minimum_required
#
# cmake_version_release:
#   Latest official release -- updated manually
#   May be the same as cmake_version_min.
#
# cmake_version_rc:
#   Latest available release candidate for vdir "v2.8" -- detected automatically
#   May be the same as the latest official release if there are no release
#   candidates currently available.
#
# cmake_version_dev:
#   Latest available build for vdir "vCVS" -- detected automatically
#

cmake_version_min = "2.8.0" # manual_update

cmake_version_release = "2.8.8" # manual_update

cmake_version_rc = "2.8.10.2" # auto_update v2.8

cmake_version_dev = "2.8.10.20130226-gd10c2" # auto_update vCVS
