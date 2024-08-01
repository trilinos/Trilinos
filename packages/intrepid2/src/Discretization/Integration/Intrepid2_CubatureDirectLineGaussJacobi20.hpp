// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CubatureDirectLineGaussJacobi20.hpp
    \brief  Header file for the Intrepid2::CubatureDirectLineGaussJacobi20 class.
    \author Created by P. Bochev, D. Ridzal and M. Perego
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_DIRECT_LINE_GAUSSJACOBI20_HPP__
#define __INTREPID2_CUBATURE_DIRECT_LINE_GAUSSJACOBI20_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_CubatureDirect.hpp"

namespace Intrepid2 {

  /** \class Intrepid2::CubatureDirectLineGaussJacobi20
      \brief Defines GaussJacobi20 integration rules on a line used for Pyramid only.
  */
  template<typename DeviceType = void,
           typename pointValueType = double,
           typename weightValueType = double>
  class CubatureDirectLineGaussJacobi20
    : public CubatureDirect<DeviceType,pointValueType,weightValueType> {
  public:
    typedef typename CubatureDirect<DeviceType,pointValueType,weightValueType>::CubatureDataStatic CubatureDataStatic;
    typedef typename CubatureDirect<DeviceType,pointValueType,weightValueType>::CubatureData       CubatureData;

    typedef typename CubatureDirect<DeviceType,pointValueType,weightValueType>::PointViewType  PointViewType;
    typedef typename CubatureDirect<DeviceType,pointValueType,weightValueType>::weightViewType weightViewType;

  private:
    // static data initialize upto 1 but we only support upto Parameters::MaxCubatureDegreePyr
    static constexpr ordinal_type cubatureDataStaticSize=12;

    /** \brief Complete set of data defining line Gauss(-Legendre) rules.
     */
    static const CubatureDataStatic cubatureDataStatic_[cubatureDataStaticSize]; // initialized once

  public:

    /** \brief Constructor.
        \param degree           [in]     - The degree of polynomials that are integrated
        exactly by this cubature rule. Default: 0.
    */
    CubatureDirectLineGaussJacobi20(const ordinal_type degree = 0);

    /** \brief Returns cubature name.
     */
    virtual
    const char*
    getName() const override {
      return "CubatureDirectLineGaussJacobi20";
    }
  };

}

// include templated definitions
#include <Intrepid2_CubatureDirectLineGaussJacobi20Def.hpp>

#endif
