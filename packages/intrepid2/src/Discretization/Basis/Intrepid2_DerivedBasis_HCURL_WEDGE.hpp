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
// Questions? Contact Mauro Perego  (mperego@sandia.gov) or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid2_DerivedBasis_HCURL_WEDGE.hpp
    \brief  Implementation of H(curl) basis on the wedge that is templated on H(grad,tri), H(curl,tri), and H(vol,line).
    \author Created by N.V. Roberts.
 
 This class constructs the H(curl) space as the direct sum of two families of tensor-product bases on the triangle and line:
 - family 1: H(curl,tri)  x  H(grad,line), placed in the x and y components of vector output
 - family 2: H(grad,tri) x  H(vol,line),  placed in the z component of vector output
 
 Unfortunately, the way that Family 1 decomposes into operators on the component bases cannot be expressed in
 OperatorTensorDecomposition because Intrepid2::EOperator does not include a VALUE_X, VALUE_Y operator -- the
 the CURL evaluation requires the H(curl,tri) component to be evaluated with OP_VALUE, then rotated 90 degrees and multiplied
 by the GRAD of the H(grad,line) component.
 
 Therefore, we instead avoid any use of OperatorTensorDecomposition in our implementation of Family 1, instead overriding
 TensorBasis::getValues(BasisValues,TensorPoints,EOperator) and TensorBasis::allocateBasisValues().
 
 Our Famiy 1 corresponds to the following ESEAS entities:
 - mixed edges
 - triangle faces, Family I and II
 - quadrilateral faces, Family II
 - Interior Family I
 - Interior Family II.
 
 Our Family 2 corresponds to:
 - quadrilateral edges
 - quadrilateral faces, Family I
 - Interior Family III.
 
 See p. 449 in Fuentes et al. (https://doi.org/10.1016/j.camwa.2015.04.027).
 */

#ifndef Intrepid2_DerivedBasis_HCURL_WEDGE_h
#define Intrepid2_DerivedBasis_HCURL_WEDGE_h

#include <Kokkos_View.hpp>
#include <Kokkos_DynRankView.hpp>

#include "Intrepid2_Polynomials.hpp"
#include "Intrepid2_Sacado.hpp"

#include "Intrepid2_DirectSumBasis.hpp"
#include "Intrepid2_TensorBasis.hpp"

namespace Intrepid2
{
  template<class HCURL_TRI, class HGRAD_LINE>
  class Basis_Derived_HCURL_Family1_WEDGE
  : public Basis_TensorBasis<typename HGRAD_LINE::BasisBase>
  {
  public:
    using OutputViewType = typename HGRAD_LINE::OutputViewType;
    using PointViewType  = typename HGRAD_LINE::PointViewType ;
    using ScalarViewType = typename HGRAD_LINE::ScalarViewType;
    
    using BasisBase = typename HGRAD_LINE::BasisBase;
  
    using DeviceType      = typename BasisBase::DeviceType;
    using ExecutionSpace  = typename BasisBase::ExecutionSpace;
    using OutputValueType = typename BasisBase::OutputValueType;
    using PointValueType  = typename BasisBase::PointValueType;
    
    using TriCurlBasis  = HCURL_TRI;
    using LineGradBasis = HGRAD_LINE;
    
    using TensorBasis = Basis_TensorBasis<BasisBase>;
  public:
    /** \brief  Constructor.
        \param [in] polyOrder_xy - the polynomial order in the x and y dimensions.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HCURL_Family1_WEDGE(int polyOrder_xy, int polyOrder_z, const EPointType pointType = POINTTYPE_DEFAULT)
    :
    TensorBasis(Teuchos::rcp( new TriCurlBasis(polyOrder_xy,pointType)),
                Teuchos::rcp( new LineGradBasis(polyOrder_z,pointType)))
    {
      this->functionSpace_ = FUNCTION_SPACE_HCURL;
      this->setShardsTopologyAndTags();
    }
    
    using TensorBasis::getValues;
    
    /** \brief Returns a simple decomposition of the specified operator: what operator(s) should be applied to basis1, and what operator(s) to basis2.  A one-element vector corresponds to a single TensorData entry; a multiple-element vector corresponds to a VectorData object with axialComponents = false.
    */
    virtual OperatorTensorDecomposition getSimpleOperatorDecomposition(const EOperator &operatorType) const override
    {
      const EOperator & VALUE = OPERATOR_VALUE;
      const EOperator & GRAD  = OPERATOR_GRAD;
      const EOperator & CURL  = OPERATOR_CURL;
      if (operatorType == VALUE)
      {
        // family 1 goes in x,y components
        std::vector< std::vector<EOperator> > ops(2);
        ops[0] = std::vector<EOperator>{VALUE,VALUE}; // occupies x,y
        ops[1] = std::vector<EOperator>{}; // zero z
        std::vector<double> weights {1.0, 0.0};
        return OperatorTensorDecomposition(ops, weights);
      }
      else if (operatorType == CURL)
      {
        // curl of (f_x(x,y) g(z), f_y(x,y) g(z), 0), where f is in H(curl,tri), g in H(grad,line)
        // x,y components of curl: rot(f) dg/dz, where rot(f) is a 90-degree rotation of f.
        //   z component  of curl: curl(f) g, where curl(f) is the 2D curl, a scalar.
        std::vector< std::vector<EOperator> > ops(2);
        ops[0] = std::vector<EOperator>{VALUE,GRAD}; // occupies the x,y components
        ops[1] = std::vector<EOperator>{CURL,VALUE}; // z component
        std::vector<double> weights {1.0, 1.0};
        OperatorTensorDecomposition opDecomposition(ops, weights);
        opDecomposition.setRotateXYNinetyDegrees(true);
        return opDecomposition;
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported operator type");
      }
    }
    
//    /** \brief Allocate BasisValues container suitable for passing to the getValues() variant that takes a TensorPoints container as argument.
//
//        The basic exact-sequence operators are supported (VALUE, CURL).
//     */
//    virtual BasisValues<OutputValueType,DeviceType> allocateBasisValues( TensorPoints<PointValueType,DeviceType> points, const EOperator operatorType = OPERATOR_VALUE) const override
//    {
//      if (points.numTensorComponents() == 1)
//      {
//        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "trivial tensor points structure not (yet) supported");
//      }
//      INTREPID2_TEST_FOR_EXCEPTION(points.numTensorComponents() != 2, std::invalid_argument, "points must have either 1 or 2 tensor components");
//
//      auto points_xy = points.getTensorComponent(0);
//      INTREPID2_TEST_FOR_EXCEPTION(points_xy.extent(1) != 2, std::invalid_argument, "first component of tensor points container must have spatial dimension 2 (corresponding to the triangle)");
//
//      auto points_z = points.getTensorComponent(1);
//      INTREPID2_TEST_FOR_EXCEPTION(points_z.extent(1) != 1, std::invalid_argument, "first component of tensor points container must have spatial dimension 1 (corresponding to the line)");
//
//      BasisValues<OutputValueType,DeviceType> triCurlBasisValues, lineGradBasisValues;
//      triCurlBasisValues  = this->TensorBasis::basis1_->allocateBasisValues(points_xy, OPERATOR_VALUE);
//      lineGradBasisValues = this->TensorBasis::basis2_->allocateBasisValues(points_z,  OPERATOR_VALUE);
//
//      if (operatorType == OPERATOR_VALUE)
//      {
//        auto triCurlVectorData = triCurlBasisValues.vectorData();
//        ordinal_type triCurlFamilyCount = triCurlVectorData.numFamilies();
//
//        std::vector< std::vector<TensorData<OutputValueType,DeviceType> > > vectorComponents(triCurlFamilyCount); // outer dimension: numFamilies; inner dimension: number of vector components (here, should be 1).
//
//        INTREPID2_TEST_FOR_EXCEPTION(triCurlBasisValues.vectorData().numComponents() != 1, std::invalid_argument, "Unexpected component count for tri basis");
//
//        INTREPID2_TEST_FOR_EXCEPTION(lineGradBasisValues.tensorData().numTensorComponents() != 1, std::invalid_argument, "Unexpected tensor component count for line basis");
//
//        auto lineData = lineGradBasisValues.tensorData().getTensorComponent(0);
//
//        for (ordinal_type familyOrdinal=0; familyOrdinal<triCurlFamilyCount; familyOrdinal++)
//        {
//          std::vector<TensorData<OutputValueType,DeviceType> > &componentsForFamily = vectorComponents[familyOrdinal];
//          ordinal_type triComponentCount = triCurlVectorData.numComponents();
//          for (ordinal_type componentOrdinal=0; componentOrdinal<triComponentCount; componentOrdinal++)
//          {
//            auto &triComponentData = triCurlVectorData.getComponent(familyOrdinal,componentOrdinal);
//            INTREPID2_TEST_FOR_EXCEPTION(triComponentData.numTensorComponents().numComponents() != 1, std::invalid_argument, "Unexpected component count for tri basis");
//
//            std::vector< Data<OutputValueType,DeviceType> > componentData { triComponentData.getTensorComponent(0), lineData };
//            componentsForFamily.push_back(TensorData<OutputValueType,DeviceType>(componentData));
//          }
//          // zero in z component:
//          componentsForFamily.push_back(TensorData<OutputValueType,DeviceType>());
//        }
//
//        VectorData<OutputValueType,DeviceType> vectorData(vectorComponents);
//        return BasisValues<OutputValueType,DeviceType>(vectorData);
//      }
//      else if (operatorType == OPERATOR_CURL)
//      {
//        // curl of (f_x(x,y) g(z), f_y(x,y) g(z), 0), where f is in H(curl,tri), g in H(grad,line)
//        // x,y components of curl: rot(f) dg/dz, where rot(f) is a 90-degree rotation of f.
//        //   z component  of curl: curl(f) g, where curl(f) is the 2D curl, a scalar.
//
//        BasisValues<OutputValueType,DeviceType> triCurlBasisCurls, lineGradBasisGrads;
//        triCurlBasisCurls  = this->TensorBasis::basis1_->allocateBasisValues(points_xy, OPERATOR_CURL);
//        lineGradBasisGrads = this->TensorBasis::basis2_->allocateBasisValues(points_z,  OPERATOR_GRAD);
//
//        auto triCurlVectorData = triCurlBasisValues.vectorData();
//        ordinal_type triCurlFamilyCount = triCurlVectorData.numFamilies();
//
//        std::vector< std::vector<TensorData<OutputValueType,DeviceType> > > vectorComponents(triCurlFamilyCount); // outer dimension: numFamilies; inner dimension: number of vector components (here, should be 1).
//
//        INTREPID2_TEST_FOR_EXCEPTION(triCurlBasisValues.vectorData().numComponents() != 1, std::invalid_argument, "Unexpected component count for tri basis");
//        INTREPID2_TEST_FOR_EXCEPTION(lineGradBasisGrads.tensorData().numTensorComponents() != 1, std::invalid_argument, "Unexpected tensor component count for line basis");
//
//        auto lineGrads = lineGradBasisGrads.tensorData().getTensorComponent(0);
//
//        for (ordinal_type familyOrdinal=0; familyOrdinal<triCurlFamilyCount; familyOrdinal++)
//        {
//          INTREPID2_TEST_FOR_EXCEPTION(triCurlBasisCurls.tensorData(familyOrdinal).numComponents() != 1, std::invalid_argument, "Unexpected component count for tri basis");
//
//          std::vector<TensorData<OutputValueType,DeviceType> > &componentsForFamily = vectorComponents[familyOrdinal];
//          ordinal_type triComponentCount = triCurlVectorData.numComponents();
//          for (ordinal_type componentOrdinal=0; componentOrdinal<triComponentCount; componentOrdinal++)
//          {
//            auto &triComponentData = triCurlVectorData.getComponent(familyOrdinal,componentOrdinal);
//            INTREPID2_TEST_FOR_EXCEPTION(triComponentData.numTensorComponents().numComponents() != 1, std::invalid_argument, "Unexpected component count for tri basis");
//
//            // x,y components of curl: rot(f) dg/dz, where rot(f) is a 90-degree rotation of f (here, we use the VALUE of f, which has the same size, for the allocation).
//            std::vector< Data<OutputValueType,DeviceType> > componentData { triComponentData.getTensorComponent(0), lineGrads };
//            componentsForFamily.push_back(TensorData<OutputValueType,DeviceType>(componentData));
//          }
//          // z component  of curl: curl(f) g, where curl(f) is the 2D curl, a scalar.
//          const auto &curlTensorData = triCurlBasisCurls.tensorData(familyOrdinal);
//          std::vector< Data<OutputValueType,DeviceType> > zComponentData { curlTensorData.getTensorComponent(0), lineGrads };
//          componentsForFamily.push_back(TensorData<OutputValueType,DeviceType>(zComponentData));
//        }
//
//        VectorData<OutputValueType,DeviceType> vectorData(vectorComponents);
//        return BasisValues<OutputValueType,DeviceType>(vectorData);
//      }
//      else
//      {
//        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported operator");
//      }
//    }
    
    /** \brief  Evaluation of a FEM basis on a <strong>reference cell</strong>, using point and output value containers that allow preservation of tensor-product structure.

        Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
        points in the <strong>reference cell</strong> for which the basis is defined.

        \param  outputValues      [out] - variable rank array with the basis values.  Should be allocated using Basis::allocateBasisValues().
        \param  inputPoints       [in]  - rank-2 array (P,D) with the evaluation points.  This should be allocated using Cubature::allocateCubaturePoints() and filled using Cubature::getCubature().
        \param  operatorType      [in]  - the operator acting on the basis function
     
        This is the preferred getValues() method for TensorBasis and DirectSumBasis and their subclasses.  It allows a reduced memory footprint and optimized integration, etc.
    */
//    virtual
//    void
//    getValues(       BasisValues<OutputValueType,DeviceType> outputValues,
//               const TensorPoints<PointValueType,DeviceType>  inputPoints,
//               const EOperator operatorType = OPERATOR_VALUE ) const override
//    {
//      // TODO: implement this.  (Below, the implementation from TensorBasis.)
//
////      const ordinal_type numTensorComponents = tensorComponents_.size();
////      if (inputPoints.numTensorComponents() < numTensorComponents)
////      {
////        // then we require that both inputPoints and outputValues trivial tensor structure
////        INTREPID2_TEST_FOR_EXCEPTION( inputPoints.numTensorComponents() != 1, std::invalid_argument, "If inputPoints differs from the tensor basis in component count, then inputPoints must have trivial tensor product structure" );
////        INTREPID2_TEST_FOR_EXCEPTION( outputValues.numFamilies() != 1, std::invalid_argument, "If inputPoints differs from the tensor basis in component count, outputValues must have a single family with trivial tensor product structure" );
////        INTREPID2_TEST_FOR_EXCEPTION( outputValues.tensorData().numTensorComponents() != 1, std::invalid_argument, "If inputPoints differs from the tensor basis in component count, outputValues must have a single family with trivial tensor product structure" );
////
////        OutputViewType outputView = outputValues.tensorData().getTensorComponent(0).getUnderlyingView();
////        PointViewType   pointView = inputPoints.getTensorComponent(0);
////        this->getValues(outputView, pointView, operatorType);
////        return;
////      }
////
////      OperatorTensorDecomposition operatorDecomposition = getOperatorDecomposition(operatorType);
////
////      const ordinal_type numVectorComponents = operatorDecomposition.numVectorComponents();
////      const bool               useVectorData = numVectorComponents > 1;
////      const ordinal_type  numBasisComponents = operatorDecomposition.numBasisComponents();
////
////      for (ordinal_type vectorComponentOrdinal=0; vectorComponentOrdinal<numVectorComponents; vectorComponentOrdinal++)
////      {
////        const double weight = operatorDecomposition.weight(vectorComponentOrdinal);
////        ordinal_type pointComponentOrdinal = 0;
////        for (ordinal_type basisOrdinal=0; basisOrdinal<numBasisComponents; basisOrdinal++, pointComponentOrdinal++)
////        {
////          const EOperator op = operatorDecomposition.op(vectorComponentOrdinal, basisOrdinal);
////          // by convention, op == OPERATOR_MAX signals a zero component; skip
////          if (op != OPERATOR_MAX)
////          {
////            const int vectorFamily = 0; // TensorBasis always has just a single family; multiple families arise in DirectSumBasis
////            auto tensorData = useVectorData ? outputValues.vectorData().getComponent(vectorFamily,vectorComponentOrdinal) : outputValues.tensorData();
////            INTREPID2_TEST_FOR_EXCEPTION( ! tensorData.getTensorComponent(basisOrdinal).isValid(), std::invalid_argument, "Invalid output component encountered");
////
////            const Data<OutputValueType,DeviceType> & outputData = tensorData.getTensorComponent(basisOrdinal);
////
////            auto basisValueView = outputData.getUnderlyingView();
////            PointViewType  pointView = inputPoints.getTensorComponent(pointComponentOrdinal);
////            const ordinal_type basisDomainDimension = tensorComponents_[basisOrdinal]->getDomainDimension();
////            if (pointView.extent_int(1) == basisDomainDimension)
////            {
////              tensorComponents_[basisOrdinal]->getValues(basisValueView, pointView, op);
////            }
////            else
////            {
////              // we need to wrap the basisValueView in a BasisValues container, and to wrap the point components in a TensorPoints container.
////
////              // combine point components to build up to basisDomainDimension
////              ordinal_type dimsSoFar = 0;
////              std::vector< ScalarView<PointValueType,DeviceType> > basisPointComponents;
////              while (dimsSoFar < basisDomainDimension)
////              {
////                INTREPID2_TEST_FOR_EXCEPTION(pointComponentOrdinal >= inputPoints.numTensorComponents(), std::invalid_argument, "Error in processing points container; perhaps it is mis-sized?");
////                const auto & pointComponent = inputPoints.getTensorComponent(pointComponentOrdinal);
////                const ordinal_type numComponentDims   = pointComponent.extent_int(1);
////                dimsSoFar += numComponentDims;
////                INTREPID2_TEST_FOR_EXCEPTION(dimsSoFar > inputPoints.numTensorComponents(), std::invalid_argument, "Error in processing points container; perhaps it is mis-sized?");
////                basisPointComponents.push_back(pointComponent);
////                if (dimsSoFar < basisDomainDimension)
////                {
////                  // we will pass through this loop again, so we should increment the point component ordinal
////                  pointComponentOrdinal++;
////                }
////              }
////
////              TensorPoints<PointValueType, DeviceType> basisPoints(basisPointComponents);
////
////              bool useVectorData2 = (basisValueView.rank() == 3);
////
////              BasisValues<OutputValueType,DeviceType> basisValues;
////              if (useVectorData2)
////              {
////                VectorData<OutputValueType,DeviceType> vectorData(outputData);
////                basisValues = BasisValues<OutputValueType,DeviceType>(vectorData);
////              }
////              else
////              {
////                TensorData<OutputValueType,DeviceType> tensorData2(outputData);
////                basisValues = BasisValues<OutputValueType,DeviceType>(tensorData2);
////              }
////
////              tensorComponents_[basisOrdinal]->getValues(basisValues, basisPoints, op);
////            }
////
////            // if weight is non-trivial (not 1.0), then we need to multiply one of the component views by weight.
////            // we do that for the first basisOrdinal's values
////            if ((weight != 1.0) && (basisOrdinal == 0))
////            {
////              if (basisValueView.rank() == 2)
////              {
////                auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<2>>({0,0},{basisValueView.extent_int(0),basisValueView.extent_int(1)});
////                Kokkos::parallel_for("multiply basisValueView by weight", policy,
////                KOKKOS_LAMBDA (const int &fieldOrdinal, const int &pointOrdinal) {
////                  basisValueView(fieldOrdinal,pointOrdinal) *= weight;
////                });
////              }
////              else if (basisValueView.rank() == 3)
////              {
////                auto policy = Kokkos::MDRangePolicy<ExecutionSpace,Kokkos::Rank<3>>({0,0,0},{basisValueView.extent_int(0),basisValueView.extent_int(1),basisValueView.extent_int(2)});
////                Kokkos::parallel_for("multiply basisValueView by weight", policy,
////                KOKKOS_LAMBDA (const int &fieldOrdinal, const int &pointOrdinal, const int &d) {
////                  basisValueView(fieldOrdinal,pointOrdinal,d) *= weight;
////                });
////              }
////              else
////              {
////                INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported rank for basisValueView");
////              }
////            }
////          }
////        }
////      }
//    }
    
    /** \brief  multi-component getValues() method (required/called by TensorBasis)
        \param [out] outputValues - the view into which to place the output values
        \param [in] operatorType - the operator on the basis
        \param [in] inputPoints1 - input points in the x dimension
        \param [in] inputPoints2 - input points in the y dimension
        \param [in] tensorPoints - if true, inputPoints1 and inputPoints2 should be understood as tensorial components of the points in outputValues (i.e., the evaluation points are the tensor product of inputPoints1 and inputPoints2).  If false, inputPoints1 and inputPoints2 should correspond elementwise to the evaluation points.
     */
    virtual void getValues(OutputViewType outputValues, const EOperator operatorType,
                           const PointViewType  inputPoints1, const PointViewType inputPoints2,
                           bool tensorPoints) const override
    {
      // TODO: implement this.  (Below, the implementation from HCURL_QUAD)
      
//      EOperator op1, op2;
//      if (operatorType == OPERATOR_VALUE)
//      {
//        op1 = OPERATOR_VALUE;
//        op2 = OPERATOR_VALUE;
//
//        // family 2 goes in the y component; 0 in the x component
//        auto outputValuesComponent1 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),0);
//        auto outputValuesComponent2 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),1);
//
//        // place 0 in the x component
//        Kokkos::deep_copy(outputValuesComponent1, 0.0);
//        this->TensorBasis::getValues(outputValuesComponent2,
//                                     inputPoints1, op1,
//                                     inputPoints2, op2, tensorPoints);
//
//      }
//      else if (operatorType == OPERATOR_CURL)
//      {
//        // family 2 gets a d/dx applied to the second (nonzero) vector component
//        // since this is H(GRAD)(x) * H(VOL)(y), this amounts to taking the derivative in the first tensorial component
//        op1 = OPERATOR_GRAD;
//        op2 = OPERATOR_VALUE;
//
//        this->TensorBasis::getValues(outputValues,
//                                     inputPoints1, op1,
//                                     inputPoints2, op2, tensorPoints);
//      }
//      else
//      {
//        INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"operator not yet supported");
//      }
    }

    /** \brief  Fills in coefficients of degrees of freedom for Lagrangian basis on the reference cell
        \param [out] dofCoeffs - the container into which to place the degrees of freedom.

     dofCoeffs have shape (F,D), field dimension matches the cardinality of the basis, and D is the
     basis dimension.

     Degrees of freedom coefficients are such that
     \phi_i(dofCoords_(j)) \cdot dofCoeffs_(j)  = \delta_ij,
     where \phi_i are the basis and \delta_ij the Kronecker delta.
     Note that getDofCoeffs() is supported only for Lagrangian bases.
     */
    virtual void getDofCoeffs( ScalarViewType dofCoeffs ) const override {
      auto dofCoeffs1 = Kokkos::subview(dofCoeffs,Kokkos::ALL(),0);
      auto dofCoeffs2 = Kokkos::subview(dofCoeffs,Kokkos::ALL(),1);
      this->TensorBasis::getDofCoeffs(dofCoeffs1);
      Kokkos::deep_copy(dofCoeffs2,0.0);
    }
  };

  template<class HGRAD_TRI, class HVOL_LINE>
  class Basis_Derived_HCURL_Family2_WEDGE
  : public Basis_TensorBasis<typename HVOL_LINE::BasisBase>
  {

  public:
    using OutputViewType = typename HVOL_LINE::OutputViewType;
    using PointViewType  = typename HVOL_LINE::PointViewType ;
    using ScalarViewType = typename HVOL_LINE::ScalarViewType;
    
    using TriGradBasis  = HGRAD_TRI;
    using LineHVolBasis = HVOL_LINE;
    
    using BasisBase   = typename HVOL_LINE::BasisBase;
    using TensorBasis = Basis_TensorBasis<BasisBase>;
    
    using DeviceType      = typename BasisBase::DeviceType;
    using ExecutionSpace  = typename BasisBase::ExecutionSpace;
    using OutputValueType = typename BasisBase::OutputValueType;
    using PointValueType  = typename BasisBase::PointValueType;
    
    /** \brief  Constructor.
        \param [in] polyOrder_xy - the polynomial order in the x and y dimensions.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HCURL_Family2_WEDGE(int polyOrder_xy, int polyOrder_z, const EPointType pointType = POINTTYPE_DEFAULT)
    :
    TensorBasis(Teuchos::rcp( new  TriGradBasis(polyOrder_xy,  pointType) ),
                Teuchos::rcp( new LineHVolBasis(polyOrder_z-1, pointType) ))
    {
      this->functionSpace_ = FUNCTION_SPACE_HCURL;
      this->setShardsTopologyAndTags();
    }
    
    /** \brief Returns a simple decomposition of the specified operator: what operator(s) should be applied to basis1, and what operator(s) to basis2.  One-element ops and weights vectors correspond to a single TensorData entry; multiple-element vectors correspond to a VectorData object with axialComponents = false.
    */
    OperatorTensorDecomposition getSimpleOperatorDecomposition(const EOperator &operatorType) const override
    {
      if (operatorType == OPERATOR_CURL)
      {
        // curl of (0,0,f) is (df/dy, -df/dx, 0)
        
        // this is a rotation of gradient of the triangle H^1 basis, times the H(vol) line value
        std::vector< std::vector<EOperator> > ops(2);
        ops[0] = std::vector<EOperator>{OPERATOR_GRAD,OPERATOR_VALUE}; // occupies the x,y components
        ops[1] = std::vector<EOperator>{};
        std::vector<double> weights {-1.0, 0.0}; // -1 because the rotation goes from (df/dx,df/dy) --> (-df/dy,df/dx), and we want (df/dy,-df/dx)
        OperatorTensorDecomposition opDecomposition(ops, weights);
        opDecomposition.setRotateXYNinetyDegrees(true);
        return opDecomposition;
      }
      else if (OPERATOR_VALUE == operatorType)
      {
        // family 2 goes in z component
        std::vector< std::vector<EOperator> > ops(2);
        ops[0] = std::vector<EOperator>{}; // because family I identifies this as spanning (x,y), empty op here will also span (x,y)
        ops[1] = std::vector<EOperator>{OPERATOR_VALUE,OPERATOR_VALUE}; // z component
        std::vector<double> weights {0.0, 1.0};
        return OperatorTensorDecomposition(ops, weights);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported operator type");
      }
    }
    
    using TensorBasis::getValues;
    
    /** \brief  multi-component getValues() method (required/called by TensorBasis)
        \param [out] outputValues - the view into which to place the output values
        \param [in] operatorType - the operator on the basis
        \param [in] inputPoints1 - input points in the x dimension
        \param [in] inputPoints2 - input points in the y dimension
        \param [in] tensorPoints - if true, inputPoints1 and inputPoints2 should be understood as tensorial components of the points in outputValues (i.e., the evaluation points are the tensor product of inputPoints1 and inputPoints2).  If false, inputPoints1 and inputPoints2 should correspond elementwise to the evaluation points.
     */
    virtual void getValues(OutputViewType outputValues, const EOperator operatorType,
                           const PointViewType  inputPoints1, const PointViewType inputPoints2,
                           bool tensorPoints) const override
    {
      // TODO: implement this.  (Below, the implementation from HCURL_QUAD)
//      EOperator op1, op2;
//      if (operatorType == OPERATOR_VALUE)
//      {
//        op1 = OPERATOR_VALUE;
//        op2 = OPERATOR_VALUE;
//
//        // family 2 goes in the y component; 0 in the x component
//        auto outputValuesComponent1 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),0);
//        auto outputValuesComponent2 = Kokkos::subview(outputValues,Kokkos::ALL(),Kokkos::ALL(),1);
//
//        // place 0 in the x component
//        Kokkos::deep_copy(outputValuesComponent1, 0.0);
//        this->TensorBasis::getValues(outputValuesComponent2,
//                                     inputPoints1, op1,
//                                     inputPoints2, op2, tensorPoints);
//
//      }
//      else if (operatorType == OPERATOR_CURL)
//      {
//        // family 2 gets a d/dx applied to the second (nonzero) vector component
//        // since this is H(GRAD)(x) * H(VOL)(y), this amounts to taking the derivative in the first tensorial component
//        op1 = OPERATOR_GRAD;
//        op2 = OPERATOR_VALUE;
//
//        this->TensorBasis::getValues(outputValues,
//                                     inputPoints1, op1,
//                                     inputPoints2, op2, tensorPoints);
//      }
//      else
//      {
//        INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"operator not yet supported");
//      }
    }

    /** \brief  Fills in coefficients of degrees of freedom for Lagrangian basis on the reference cell
        \param [out] dofCoeffs - the container into which to place the degrees of freedom.

     dofCoeffs have shape (F,D), field dimension matches the cardinality of the basis, and D is the
     basis dimension.

     Degrees of freedom coefficients are such that
     \phi_i(dofCoords_(j)) \cdot dofCoeffs_(j)  = \delta_ij,
     where \phi_i are the basis and \delta_ij the Kronecker delta.
     Note that getDofCoeffs() is supported only for Lagrangian bases.
     */
    virtual void getDofCoeffs( ScalarViewType dofCoeffs ) const override {
      auto dofCoeffs1 = Kokkos::subview(dofCoeffs,Kokkos::ALL(),0);
      auto dofCoeffs2 = Kokkos::subview(dofCoeffs,Kokkos::ALL(),1);
      Kokkos::deep_copy(dofCoeffs1,0.0);
      this->TensorBasis::getDofCoeffs(dofCoeffs2);
    }
  };
  
  template<class HGRAD_TRI, class HCURL_TRI, class HGRAD_LINE, class HVOL_LINE>
  class Basis_Derived_HCURL_WEDGE
  : public Basis_DirectSumBasis <typename HGRAD_LINE::BasisBase>
  {
    using Family1 = Basis_Derived_HCURL_Family1_WEDGE<HCURL_TRI, HGRAD_LINE>;
    using Family2 = Basis_Derived_HCURL_Family2_WEDGE<HGRAD_TRI, HVOL_LINE>;
    using DirectSumBasis = Basis_DirectSumBasis <typename HGRAD_LINE::BasisBase>;
  public:
    using BasisBase = typename HGRAD_LINE::BasisBase;

  protected:
    std::string name_;
    ordinal_type order_xy_;
    ordinal_type order_z_;
    EPointType pointType_;

  public:
    using ExecutionSpace  = typename HGRAD_LINE::ExecutionSpace;
    using OutputValueType = typename HGRAD_LINE::OutputValueType;
    using PointValueType  = typename HGRAD_LINE::PointValueType;
    
    /** \brief  Constructor.
        \param [in] polyOrder_xy - the polynomial order in the x and y dimensions.
        \param [in] polyOrder_z - the polynomial order in the z dimension.
        \param [in] pointType   - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HCURL_WEDGE(int polyOrder_xy, int polyOrder_z, const EPointType pointType=POINTTYPE_DEFAULT)
    :
    DirectSumBasis(Teuchos::rcp( new Family1(polyOrder_xy, polyOrder_z, pointType) ),
                   Teuchos::rcp( new Family2(polyOrder_xy, polyOrder_z, pointType) ))
    {
      this->functionSpace_ = FUNCTION_SPACE_HCURL;

      std::ostringstream basisName;
      basisName << "HCURL_WEDGE (" << this->DirectSumBasis::getName() << ")";
      name_ = basisName.str();

      order_xy_ = polyOrder_xy;
      order_z_ = polyOrder_z;
      pointType_ = pointType;
    }
    
    /** \brief  Constructor.
        \param [in] polyOrder - the polynomial order to use in all dimensions.
        \param [in] pointType - type of lattice used for creating the DoF coordinates.
     */
    Basis_Derived_HCURL_WEDGE(int polyOrder, const EPointType pointType=POINTTYPE_DEFAULT) : Basis_Derived_HCURL_WEDGE(polyOrder, polyOrder, pointType) {}

    /** \brief True if orientation is required
    */
    virtual bool requireOrientation() const override
    {
      return (this->getDofCount(1,0) > 0); //if it has edge DOFs, than it needs orientations
    }

    /** \brief  Returns basis name

     \return the name of the basis
     */
    virtual
    const char*
    getName() const override {
      return name_.c_str();
    }

    /** \brief returns the basis associated to a subCell.

        The bases of the subCell are the restriction to the subCell of the bases of the parent cell,
        projected to the subCell line.

        TODO: test this method when different orders are used in different directions
        \param [in] subCellDim - dimension of subCell
        \param [in] subCellOrd - position of the subCell among of the subCells having the same dimension
        \return pointer to the subCell basis of dimension subCellDim and position subCellOrd
     */
    Teuchos::RCP<BasisBase>
      getSubCellRefBasis(const ordinal_type subCellDim, const ordinal_type subCellOrd) const override{
      if(subCellDim == 1) {
        switch(subCellOrd) {
        case 0:
        case 2:
          return Teuchos::rcp( new HVOL_LINE(order_xy_-1, pointType_) );
        case 1:
        case 3:
          return Teuchos::rcp( new HVOL_LINE(order_z_-1, pointType_) );
        }
      }

      INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Input parameters out of bounds");
    }

    /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.
     
        \return Pointer to the new Basis object.
     */
    virtual HostBasisPtr<OutputValueType, PointValueType>
    getHostBasis() const override {
      using HostBasis  = Basis_Derived_HCURL_WEDGE<typename HGRAD_TRI::HostBasis, typename HCURL_TRI::HostBasis, typename HGRAD_LINE::HostBasis, typename HVOL_LINE::HostBasis>;
      
      auto hostBasis = Teuchos::rcp(new HostBasis(order_xy_, order_z_, pointType_));
      
      return hostBasis;
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DerivedBasis_HCURL_WEDGE_h */
