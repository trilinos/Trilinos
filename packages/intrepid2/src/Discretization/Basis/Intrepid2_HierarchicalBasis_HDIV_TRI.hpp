// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HierarchicalBasis_HDIV_TRI.hpp
    \brief  H(div) basis on the triangle using a construction involving Legendre and integrated Jacobi polynomials.
    \author Created by N.V. Roberts.
 
 The H(div) basis in 2D is a rotation of the H(curl) basis.  That is how this is implemented, in terms of the H(curl) basis.
 */

#ifndef Intrepid2_HierarchicalBasis_HDIV_TRI_h
#define Intrepid2_HierarchicalBasis_HDIV_TRI_h

#include <Kokkos_DynRankView.hpp>

#include <Intrepid2_config.h>

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_HierarchicalBasis_HCURL_TRI.hpp"
#include "Intrepid2_Polynomials.hpp"
#include "Intrepid2_Utils.hpp"

namespace Intrepid2
{
  
  /** \class  Intrepid2::HierarchicalBasis_HDIV_TRI
      \brief  For mathematical details of the construction, see:
   
               Federico Fuentes, Brendan Keith, Leszek Demkowicz, Sriram Nagaraj.
               "Orientation embedded high order shape functions for the exact sequence elements of all shapes."
               Computers & Mathematics with Applications, Volume 70, Issue 4, 2015, Pages 353-458, ISSN 0898-1221.
               https://doi.org/10.1016/j.camwa.2015.04.027.
  */
  template<typename DeviceType,
           typename OutputScalar = double,
           typename PointScalar  = double,
           bool useCGBasis = true> // if useCGBasis is false, all basis functions will be associated with the interior
  class HierarchicalBasis_HDIV_TRI
  : public HierarchicalBasis_HCURL_TRI<DeviceType,OutputScalar,PointScalar,useCGBasis>
  {
  public:
    using CurlBasis = HierarchicalBasis_HCURL_TRI<DeviceType,OutputScalar,PointScalar,useCGBasis>;
    
    using BasisBase = Basis<DeviceType,OutputScalar,PointScalar>;

    using HostBasis = HierarchicalBasis_HDIV_TRI<typename Kokkos::HostSpace::device_type, OutputScalar, PointScalar, useCGBasis>;

    using typename BasisBase::OrdinalTypeArray1DHost;
    using typename BasisBase::OrdinalTypeArray2DHost;

    using typename BasisBase::OutputViewType;
    using typename BasisBase::PointViewType;
    using typename BasisBase::ScalarViewType;

    using typename BasisBase::ExecutionSpace;

  public:
    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order of the basis.
        \param [in] pointType - point type for nodal basis.  Ignored here (irrelevant for hierarchical/modal basis).
     */
    HierarchicalBasis_HDIV_TRI(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    CurlBasis(polyOrder,pointType)
    {
      this->functionSpace_ = FUNCTION_SPACE_HDIV; // CurlBasis will set the other Basis member variables correctly…
    }
    
    /** \brief  Returns basis name
     
     \return the name of the basis
     */
    const char* getName() const override {
      return "Intrepid2_HierarchicalBasis_HDIV_TRI";
    }
    
    /** \brief True if orientation is required
    */
    virtual bool requireOrientation() const override {
      return (this->getDegree() > 2);
    }

    // since the getValues() below only overrides the FEM variant, we specify that
    // we use the base class's getValues(), which implements the FVD variant by throwing an exception.
    // (It's an error to use the FVD variant on this basis.)
    using BasisBase::getValues;
    
    /** \brief  Evaluation of a FEM basis on a <strong>reference cell</strong>.

        Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
        points in the <strong>reference cell</strong> for which the basis is defined.

        \param  outputValues      [out] - variable rank array with the basis values
        \param  inputPoints       [in]  - rank-2 array (P,D) with the evaluation points
        \param  operatorType      [in]  - the operator acting on the basis functions

        \remark For rank and dimension specifications of the output array see Section
        \ref basis_md_array_sec.  Dimensions of <var>ArrayScalar</var> arguments are checked
        at runtime if HAVE_INTREPID2_DEBUG is defined.

        \remark A FEM basis spans a COMPLETE or INCOMPLETE polynomial space on the reference cell
        which is a smooth function space. Thus, all operator types that are meaningful for the
        approximated function space are admissible. When the order of the operator exceeds the
        degree of the basis, the output array is filled with the appropriate number of zeros.
    */
    virtual void getValues( OutputViewType outputValues, const PointViewType inputPoints,
                           const EOperator operatorType = OPERATOR_VALUE ) const override
    {
      // H(div) functions are a rotation of the corresponding H(curl) functions
      if (operatorType == OPERATOR_DIV)
      {
        this->CurlBasis::getValues(outputValues, inputPoints, OPERATOR_CURL);
      }
      else if (operatorType == OPERATOR_VALUE)
      {
        this->CurlBasis::getValues(outputValues, inputPoints, OPERATOR_VALUE);
        
        const ordinal_type numFields = outputValues.extent_int(0);
        const ordinal_type numPoints = outputValues.extent_int(1);

        auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>(Kokkos::Array<int,2>{0,0},Kokkos::Array<int,2>{numFields,numPoints});
        
        // multiply by [[0 1] [-1 0]]
        Kokkos::parallel_for("apply rotation to H(curl) values", policy,
        KOKKOS_LAMBDA (const int &fieldOrdinal, const int &pointOrdinal) {
          const auto tempValue = outputValues(fieldOrdinal,pointOrdinal,0);
          outputValues(fieldOrdinal,pointOrdinal,0) = outputValues(fieldOrdinal,pointOrdinal,1);
          outputValues(fieldOrdinal,pointOrdinal,1) = -tempValue;
        });
      }
    }

    /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.
     
        \return Pointer to the new Basis object.
     */
    virtual BasisPtr<typename Kokkos::HostSpace::device_type, OutputScalar, PointScalar>
    getHostBasis() const override {
      using HostDeviceType = typename Kokkos::HostSpace::device_type;
      using HostBasisType  = HierarchicalBasis_HDIV_TRI<HostDeviceType, OutputScalar, PointScalar, useCGBasis>;
      return Teuchos::rcp( new HostBasisType(this->CurlBasis::polyOrder_) );
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_HierarchicalBasis_HDIV_TRI_h */
