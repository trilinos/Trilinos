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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_CubaturePolylibDef.hpp
    \brief  Definition file for the Intrepid2::CubaturePolylib class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

namespace Intrepid2 {

  template <typename SpT, typename PT, typename WT>
  CubaturePolylib<SpT,PT,WT>::
  CubaturePolylib(const ordinal_type degree,
                  const EPolyType polytype,
                  const double alpha,
                  const double beta) 
    : CubatureDirect<SpT,PT,WT>(degree, 1) {

    INTREPID2_TEST_FOR_EXCEPTION( degree < 0 ||
                                  degree > static_cast<ordinal_type>(Parameters::MaxCubatureDegreeEdge), std::out_of_range,
                                  ">>> ERROR (CubaturePolylib): No cubature rule implemented for the desired polynomial degree.");
    INTREPID2_TEST_FOR_EXCEPTION( !isValidPolyType(polytype), std::invalid_argument,
                                  ">>> ERROR (CubaturePolylib): Invalid poly type.");

    ordinal_type npts = 0;
    switch (polytype) {
    case POLYTYPE_GAUSS:             npts = (degree+2)/2; break;
    case POLYTYPE_GAUSS_RADAU_LEFT:
    case POLYTYPE_GAUSS_RADAU_RIGHT: npts = (degree == 0 ? 2 : (degree+3)/2); break;
    case POLYTYPE_GAUSS_LOBATTO:     npts = (degree+4)/2; break;
    case POLYTYPE_MAX: break;
    }
    this->cubatureData_.numPoints_ = npts;

    auto points  = Kokkos::DynRankView<PT,Kokkos::HostSpace>("CubaturePolylib::points",  npts);
    auto weights = Kokkos::DynRankView<WT,Kokkos::HostSpace>("CubaturePolylib::weights", npts);

    Polylib::Serial::getCubature( points,
                                  weights,
                                  npts,
                                  alpha,
                                  beta,
                                  polytype );
    
    this->cubatureData_.points_  = Kokkos::DynRankView<PT,SpT>("CubaturePolylib::cubatureData_::points_",  npts, 1);
    this->cubatureData_.weights_ = Kokkos::DynRankView<WT,SpT>("CubaturePolylib::cubatureData_::weights_", npts);

    Kokkos::deep_copy(Kokkos::subdynrankview(this->cubatureData_.points_, Kokkos::ALL(), 0), points );  
    Kokkos::deep_copy(this->cubatureData_.weights_, weights );  
  }


}

