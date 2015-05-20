/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2015) Sandia Corporation
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



#include "mrk_default_kokkos_device_type.hpp"
#include "mrk_data_types.hpp"
#include "mrk_api_classes.hpp"
#include "mrk_interface_impl.hpp"
#include "mrk_manager_impl.hpp"

namespace morkon_exp {

typedef Morkon_Manager<default_kokkos_device_t, 3, MRK_QUAD4>    default_manager_3d_t;
typedef Teuchos::RCP< default_manager_3d_t >        default_manager_3d_ptr;

typedef Interface<default_kokkos_device_t, 3>     default_interface_3d_t;
typedef Teuchos::RCP< default_interface_3d_t >  default_interface_3d_ptr;

void just_check_if_it_compiles()
{
  default_manager_3d_ptr manager_0 = default_manager_3d_t::MakeInstance(0, 0);

  default_interface_3d_ptr interface_0 = manager_0->create_interface(0,0);

  manager_0->commit_interfaces();

  Tpetra::CrsMatrix<> *dummy_D = 0;
  Tpetra::CrsMatrix<> *dummy_M = 0;
  manager_0->mortar_integrate(dummy_D, dummy_M);
}

} // namespace morkon_exp
