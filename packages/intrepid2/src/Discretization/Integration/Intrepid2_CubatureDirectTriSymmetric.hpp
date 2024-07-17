// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CubatureDirectTriSymmetric.hpp
    \brief  Header file for the Intrepid2::CubatureDirectTrisymPos class.
    \author Created by P. Created by M. Perego from an initial implementation by John Burkardt 
            (https://people.sc.fsu.edu/~jburkardt/py_src/triangle_symq_rule/triangle_symq_rule.html)
*/

#ifndef __INTREPID2_CUBATURE_DIRECT_TRI_SYMMETRIC_HPP__
#define __INTREPID2_CUBATURE_DIRECT_TRI_SYMMETRIC_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_CubatureDirect.hpp"

namespace Intrepid2 {

  /** \class Intrepid2::CubatureDirectTri
      \brief Defines direct integration rules on a triangle.
            These quadrature rules are symmetric (invariant under orientation mappings)
            and the quadrature weights are all positive.
            Reference: H. Xiao and Z. Gimbutas, A numerical algorithm for the construction of efficient quadrature rules in two and higher dimensions, 
            Computers and Mathematics with Applications, 59 (2009), pp. 663-676. 
  */
  template<typename DeviceType = void,
           typename pointValueType = double,
           typename weightValueType = double>
  class CubatureDirectTriSymmetric
    : public CubatureDirect<DeviceType,pointValueType,weightValueType> {
  public:
    typedef typename CubatureDirect<DeviceType,pointValueType,weightValueType>::CubatureDataStatic CubatureDataStatic;
    typedef typename CubatureDirect<DeviceType,pointValueType,weightValueType>::CubatureData       CubatureData;

    typedef typename CubatureDirect<DeviceType,pointValueType,weightValueType>::PointViewType  PointViewType;
    typedef typename CubatureDirect<DeviceType,pointValueType,weightValueType>::weightViewType weightViewType;

  private:

    // static data initialize upto 21
    static constexpr ordinal_type cubatureDataStaticSize=51;

    /** \brief Complete set of data defining default cubature rules on a triangle.
     */
    static const CubatureDataStatic cubatureDataStatic_[cubatureDataStaticSize]; // initialized once

  public:

    /** \brief Constructor.
        \param degree           [in]     - The degree of polynomials that are integrated
        exactly by this cubature rule. Default: 0.
    */
    CubatureDirectTriSymmetric(const ordinal_type degree = 0);

    /** \brief Returns cubature name.
     */
    virtual
    const char* 
    getName() const override {
      return "CubatureDirectTriSymmetric";
    }

  }; 
  
} // end namespace Intrepid2


// include templated definitions
#include <Intrepid2_CubatureDirectTriSymmetricDef.hpp>

#endif
