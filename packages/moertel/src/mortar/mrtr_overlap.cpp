/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2006) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
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
# Questions? Contact Glen Hansen (gahanse@sandia.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "Moertel_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_MOERTEL_EXPLICIT_INSTANTIATION
#include "mrtr_overlap_Def.hpp"
#include "mrtr_overlap_utils_Def.hpp"
#include "mrtr_convexhull_Def.hpp"

#include "mrtr_interface.H"
#ifdef HAVE_MOERTEL_TPETRA
#include "Moertel_InterfaceT.hpp"
#endif

namespace MOERTEL {

  #ifdef HAVE_MOERTEL_TPETRA
    #ifdef HAVE_MOERTEL_INST_DOUBLE_INT_INT
      MOERTEL_INSTANTIATE_NESTED_TEMPLATE_CLASS_ST_LO_GO_N(MOERTEL::Overlap, MoertelT::InterfaceT, double, int, int, KokkosNode)
    #endif
    #ifdef HAVE_MOERTEL_INST_DOUBLE_INT_LONGLONGINT
      MOERTEL_INSTANTIATE_NESTED_TEMPLATE_CLASS_ST_LO_GO_N(MOERTEL::Overlap, MoertelT::InterfaceT, double, int, long long, KokkosNode)
    #endif
  #endif

  MOERTEL_INSTANTIATE_TEMPLATE_CLASS_ON_NAME_ORD(Overlap, Interface)

} // namespace Moertel

#endif
