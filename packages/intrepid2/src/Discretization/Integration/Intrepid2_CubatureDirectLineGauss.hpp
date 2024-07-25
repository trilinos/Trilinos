// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CubatureDirectLineGauss.hpp
    \brief  Header file for the Intrepid2::CubatureDirectLineGauss class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CUBATURE_DIRECT_LINE_GAUSS_HPP__
#define __INTREPID2_CUBATURE_DIRECT_LINE_GAUSS_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_CubatureDirect.hpp"

namespace Intrepid2 {

  /** \class Intrepid2::CubatureDirectLineGauss
      \brief Defines Gauss integration rules on a line.
  */
  template<typename DeviceType = void,
           typename pointValueType = double,
           typename weightValueType = double>
  class CubatureDirectLineGauss
    : public CubatureDirect<DeviceType,pointValueType,weightValueType> {
  public:
    typedef typename CubatureDirect<DeviceType,pointValueType,weightValueType>::CubatureDataStatic CubatureDataStatic;
    typedef typename CubatureDirect<DeviceType,pointValueType,weightValueType>::CubatureData       CubatureData;

    typedef typename CubatureDirect<DeviceType,pointValueType,weightValueType>::PointViewType  PointViewType;
    typedef typename CubatureDirect<DeviceType,pointValueType,weightValueType>::weightViewType weightViewType;

  private:
    // static data initialize upto 62 but we only support upto Parameters::MaxCubatureDegreeEdge
    static constexpr ordinal_type cubatureDataStaticSize=62;

    /** \brief Complete set of data defining line Gauss(-Legendre) rules.
     */
    static const CubatureDataStatic cubatureDataStatic_[cubatureDataStaticSize]; // initialized once

  public:

    /** \brief Constructor.
        \param degree           [in]     - The degree of polynomials that are integrated
        exactly by this cubature rule. Default: 0.
    */
    CubatureDirectLineGauss(const ordinal_type degree = 0);

    /** \brief Returns cubature name.
     */
    virtual
    const char*
    getName() const override {
      return "CubatureDirectLineGauss";
    }


  };

} // end namespace Intrepid2

// include templated definitions
#include <Intrepid2_CubatureDirectLineGaussDef.hpp>

#endif
