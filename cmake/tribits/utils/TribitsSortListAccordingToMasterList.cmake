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

INCLUDE(PrintVar)
INCLUDE(AppendSet)

#
# Function that does an in-place sort of a list of items according to the
# ordering in a master list
#
# NOTE: This function has wost-case N^2 complexity as the number of packages N
# or TPLs increases.  It actually has N * n complexity where N is the total
# number of packages/TPLs and n is the number of passed-in packages/TPLs.
# However, since N is not likely to ever be more than a few hundred, this is
# likely not going to be a big performance problem.  If this does become a
# performance problem, LIST(SORT ...) could be used but would require some
# work to build up the datastructures to make this very efficient.
#

FUNCTION(TRIBITS_SORT_LIST_ACCORDING_TO_MASTER_LIST  MASTER_LIST  LIST_VAR_INOUT)

  #MESSAGE("TRIBITS_SORT_LIST_ACCORDING_TO_MASTER_LIST:")
  #PRINT_VAR(MASTER_LIST)
  #PRINT_VAR(LIST_VAR_INOUT)
  #PRINT_VAR(${LIST_VAR_INOUT})

  SET(SORTED_LIST)

  FOREACH(ITEM ${MASTER_LIST})
    LIST(FIND ${LIST_VAR_INOUT} ${ITEM} ITEM_IDX)
     IF (NOT ITEM_IDX EQUAL -1)
      APPEND_SET(SORTED_LIST ${ITEM})
    ENDIF()
  ENDFOREACH()

  #PRINT_VAR(SORTED_LIST)

  SET(${LIST_VAR_INOUT} ${SORTED_LIST} PARENT_SCOPE)

ENDFUNCTION()
