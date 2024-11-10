// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CubaturePolylib.hpp
    \brief  Header file for the Intrepid2::CubaturePolylib class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_POLYLIB_HPP__
#define __INTREPID2_CUBATURE_POLYLIB_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_CubatureDirect.hpp"

#include "Intrepid2_Polylib.hpp"


namespace Intrepid2 {

  /** \class Intrepid2::CubaturePolylib
      \brief Utilizes cubature (integration) rules contained in the library Polylib
      (Spencer Sherwin, Aeronautics, Imperial College London) within Intrepid.

      They are based on zeros of Jacobi polynomials,
      e.g. Legendre (alpha=beta=0, default), Chebyshev (alpha=beta=-0.5), etc.
      They are given on the interval [-1,1] and are optimal with respect to
      the following requirements, yielding 4 subclasses:
      \li Gauss             - no restrictions (all integration points are contained in the open interval (-1,1))
      \li Gauss-Radau-Left  - left-most integration point is fixed at -1
      \li Gauss-Radau-Right - right-most integration point is fixed at +1
      \li Gauss-Lobatto     - left-most and right-most integration points are fixed at -1 and +1, respectively
  */
  template<typename DeviceType = void,
           typename pointValueType = double,
           typename weightValueType = double>
  class CubaturePolylib
    : public CubatureDirect<DeviceType,pointValueType,weightValueType> {
  public:
    typedef typename CubatureDirect<DeviceType,pointValueType,weightValueType>::PointViewType  PointViewType; 
    typedef typename CubatureDirect<DeviceType,pointValueType,weightValueType>::weightViewType weightViewType;

    CubaturePolylib( const ordinal_type degree,
                     const EPolyType polytype = POLYTYPE_GAUSS,
                     const double alpha = 0.0,
                     const double beta = 0.0 );

    /** \brief Returns cubature name.
     */
    virtual
    const char*
    getName() const override {
      return "CubaturePolylib";
    }


  };

}

#include "Intrepid2_CubaturePolylibDef.hpp"

#endif
