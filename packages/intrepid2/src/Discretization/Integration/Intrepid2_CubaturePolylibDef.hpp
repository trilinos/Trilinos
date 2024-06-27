// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CubaturePolylibDef.hpp
    \brief  Definition file for the Intrepid2::CubaturePolylib class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

namespace Intrepid2 {

  template <typename DT, typename PT, typename WT>
  CubaturePolylib<DT,PT,WT>::
  CubaturePolylib(const ordinal_type degree,
                  const EPolyType polytype,
                  const double alpha,
                  const double beta) 
    : CubatureDirect<DT,PT,WT>(degree, 1) {

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
    
    this->cubatureData_.points_  = Kokkos::DynRankView<PT,DT>("CubaturePolylib::cubatureData_::points_",  npts, 1);
    this->cubatureData_.weights_ = Kokkos::DynRankView<WT,DT>("CubaturePolylib::cubatureData_::weights_", npts);

    Kokkos::deep_copy(Kokkos::subdynrankview(this->cubatureData_.points_, Kokkos::ALL(), 0), points );  
    Kokkos::deep_copy(this->cubatureData_.weights_, weights );  
  }


}

