// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   ETI_CellGeometry.cpp
    \brief  Explicit Template Instantiation for CellGeometry.  Each definition here should be paired with a declaration in ETI.hpp.
    \author Created by N.V. Roberts.
 */

#include "Intrepid2_CellGeometry.hpp"

#include "Intrepid2_ConfigDefs.hpp"

#include <Kokkos_Core.hpp>

template class Intrepid2::CellGeometry<double,1,Kokkos::DefaultExecutionSpace>;
template class Intrepid2::CellGeometry<double,2,Kokkos::DefaultExecutionSpace>;
template class Intrepid2::CellGeometry<double,3,Kokkos::DefaultExecutionSpace>;

#ifdef HAVE_INTREPID2_SACADO
// CellGeometry - DFad, up to 3D
using Sacado_Fad_DFadType = Sacado::Fad::DFad<double>; // Sacado type used in Intrepid2 tests
template class Intrepid2::CellGeometry<Sacado_Fad_DFadType,1,Kokkos::DefaultExecutionSpace>;
template class Intrepid2::CellGeometry<Sacado_Fad_DFadType,2,Kokkos::DefaultExecutionSpace>;
template class Intrepid2::CellGeometry<Sacado_Fad_DFadType,3,Kokkos::DefaultExecutionSpace>;
#endif
