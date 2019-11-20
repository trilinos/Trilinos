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

/** \file   Intrepid2_DirectSumBasis.hpp
    \brief  Implementation of a basis that is the direct sum of two other bases.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_DirectSumBasis_h
#define Intrepid2_DirectSumBasis_h

#include <Kokkos_View.hpp>
#include <Kokkos_DynRankView.hpp>

namespace Intrepid2
{
  /**
   \class  Intrepid2::Basis_DirectSumBasis
   \brief  A basis that is the direct sum of two other bases.
   
   The direct-sum basis is ordered such that the Basis1 members come first
   (and in the same order as they exist in Basis1), followed by the members of
   Basis2, in the same order as they exist in Basis2.
   
   The two bases must agree in their BasisType (the return value of getBasisType()).
   */
  template<typename Basis1, typename Basis2>
  class Basis_DirectSumBasis
  : public Basis<typename Basis1::ExecutionSpace,
                 typename Basis1::OutputValueType,
                 typename Basis1::PointValueType>
  {
  protected:
    Basis1 basis1_;
    Basis2 basis2_;
  public:
    using OrdinalTypeArray1DHost = typename Basis1::OrdinalTypeArray1DHost;
    using OrdinalTypeArray2DHost = typename Basis1::OrdinalTypeArray2DHost;
    
    using ExecutionSpace  = typename Basis1::ExecutionSpace;
    using OutputValueType = typename Basis1::OutputValueType;
    using PointValueType  = typename Basis1::PointValueType;
    
    using OutputViewType = typename Basis1::OutputViewType;
    using PointViewType  = typename Basis1::PointViewType;
    using ScalarViewType = typename Basis1::ScalarViewType;
  public:
    /** \brief  Constructor.
        \param [in] basis1 - the instance of Basis1
        \param [in] basis2 - the instance of Basis2
     */
    Basis_DirectSumBasis(Basis1 basis1, Basis2 basis2)
    :
    basis1_(basis1),basis2_(basis2)
    {
      INTREPID2_TEST_FOR_EXCEPTION(basis1.getBasisType() != basis2.getBasisType(), std::invalid_argument, "basis1 and basis2 must agree in basis type");
      INTREPID2_TEST_FOR_EXCEPTION(basis1.getBaseCellTopology().getKey() != basis2.getBaseCellTopology().getKey(),
                                 std::invalid_argument, "basis1 and basis2 must agree in cell topology");
      INTREPID2_TEST_FOR_EXCEPTION(basis1.getCoordinateSystem() != basis2.getCoordinateSystem(),
                                 std::invalid_argument, "basis1 and basis2 must agree in coordinate system");
      
      this->basisCardinality_  = basis1.getCardinality() + basis2.getCardinality();
      this->basisDegree_       = std::max(basis1.getDegree(), basis2.getDegree());
      
      this->basisCellTopology_ = basis1.getBaseCellTopology();
      this->basisType_         = basis1.getBasisType();
      this->basisCoordinates_  = basis1.getCoordinateSystem();

      if (this->basisType_ == BASIS_FEM_HIERARCHICAL)
      {
        int degreeLength = basis1_.getPolynomialDegreeLength();
        INTREPID2_TEST_FOR_EXCEPTION(degreeLength != basis2_.getPolynomialDegreeLength(), std::invalid_argument, "Basis1 and Basis2 must agree on polynomial degree length");
        
        this->fieldOrdinalPolynomialDegree_ = OrdinalTypeArray2DHost("DirectSumBasis degree lookup",this->basisCardinality_,degreeLength);
        // our field ordinals start with basis1_; basis2_ follows
        for (int fieldOrdinal1=0; fieldOrdinal1<basis1_.getCardinality(); fieldOrdinal1++)
        {
          int fieldOrdinal = fieldOrdinal1;
          auto polynomialDegree = basis1.getPolynomialDegreeOfField(fieldOrdinal1);
          for (int d=0; d<degreeLength; d++)
          {
            this->fieldOrdinalPolynomialDegree_(fieldOrdinal,d) = polynomialDegree(d);
          }
        }
        for (int fieldOrdinal2=0; fieldOrdinal2<basis2_.getCardinality(); fieldOrdinal2++)
        {
          int fieldOrdinal = basis1.getCardinality() + fieldOrdinal2;
          
          auto polynomialDegree = basis2.getPolynomialDegreeOfField(fieldOrdinal2);
          for (int d=0; d<degreeLength; d++)
          {
            this->fieldOrdinalPolynomialDegree_(fieldOrdinal,d) = polynomialDegree(d);
          }
        }
      }
      
      // initialize tags
      {
        const auto & cardinality = this->basisCardinality_;
        
        // Basis-dependent initializations
        const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
        const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
        const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
        const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
        
        OrdinalTypeArray1DHost tagView("tag view", cardinality*tagSize);
        
        shards::CellTopology cellTopo = this->basisCellTopology_;
        
        unsigned spaceDim  = cellTopo.getDimension();
        
        ordinal_type basis2Offset = basis1_.getCardinality();
                
        for (unsigned d=0; d<=spaceDim; d++)
        {
          unsigned subcellCount = cellTopo.getSubcellCount(d);
          for (unsigned subcellOrdinal=0; subcellOrdinal<subcellCount; subcellOrdinal++)
          {
            ordinal_type subcellDofCount1 = basis1.getDofCount(d, subcellOrdinal);
            ordinal_type subcellDofCount2 = basis2.getDofCount(d, subcellOrdinal);
            
            ordinal_type subcellDofCount = subcellDofCount1 + subcellDofCount2;
            for (ordinal_type localDofID=0; localDofID<subcellDofCount; localDofID++)
            {
              ordinal_type fieldOrdinal;
              if (localDofID < subcellDofCount1)
              {
                // first basis: field ordinal matches the basis1 ordinal
                fieldOrdinal = basis1_.getDofOrdinal(d, subcellOrdinal, localDofID);
              }
              else
              {
                // second basis: field ordinal is offset by basis1 cardinality
                fieldOrdinal = basis2Offset + basis2_.getDofOrdinal(d, subcellOrdinal, localDofID - subcellDofCount1);
              }
              tagView(fieldOrdinal*tagSize+0) = d; // subcell dimension
              tagView(fieldOrdinal*tagSize+1) = subcellOrdinal;
              tagView(fieldOrdinal*tagSize+2) = localDofID;
              tagView(fieldOrdinal*tagSize+3) = subcellDofCount;
            }
          }
        }
        //        // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
        //        // tags are constructed on host
        this->setOrdinalTagData(this->tagToOrdinal_,
                                this->ordinalToTag_,
                                tagView,
                                this->basisCardinality_,
                                tagSize,
                                posScDim,
                                posScOrd,
                                posDfOrd);
      }
    }
    
    /** \brief  Fills in spatial locations (coordinates) of degrees of freedom (nodes) on the reference cell
        \param [out] dofCoords - the container into which to place the degrees of freedom.
     
     dofCoords should have shape (F,D), where the field dimension matches the cardinality of the basis, and D is the
     spatial dimension of the topology on which the basis is defined.
     
     Note that getDofCoords() is not supported by all bases; in particular, hierarchical bases do not generally support this.
     */
    virtual void getDofCoords( ScalarViewType dofCoords ) const override {
      const int basisCardinality1 = basis1_.getCardinality();
      const int basisCardinality2 = basis2_.getCardinality();
      const int basisCardinality  = basisCardinality1 + basisCardinality2;

      auto dofCoords1 = Kokkos::subview(dofCoords, std::make_pair(0,basisCardinality1),                Kokkos::ALL());
      auto dofCoords2 = Kokkos::subview(dofCoords, std::make_pair(basisCardinality1,basisCardinality), Kokkos::ALL());
      
      basis1_.getDofCoords(dofCoords1);
      basis2_.getDofCoords(dofCoords2);
    }
    
    // since the getValues() below only overrides the FEM variant, we specify that
    // we use the base class's getValues(), which implements the FVD variant by throwing an exception.
    // (It's an error to use the FVD variant on this basis.)
    using Basis<typename Basis1::ExecutionSpace, typename Basis1::OutputValueType, typename Basis1::PointValueType>::getValues;
    
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
    virtual void getValues( OutputViewType outputValues, const PointViewType  inputPoints,
                           const EOperator operatorType = OPERATOR_VALUE ) const override
    {
      int cardinality1 = basis1_.getCardinality();
      int cardinality2 = basis2_.getCardinality();
      
      auto range1 = std::make_pair(0,cardinality1);
      auto range2 = std::make_pair(cardinality1,cardinality1+cardinality2);
      if (outputValues.rank() == 2) // F,P
      {
        auto outputValues1 = Kokkos::subview(outputValues, range1, Kokkos::ALL());
        auto outputValues2 = Kokkos::subview(outputValues, range2, Kokkos::ALL());
        
        basis1_.getValues(outputValues1, inputPoints, operatorType);
        basis2_.getValues(outputValues2, inputPoints, operatorType);
      }
      else if (outputValues.rank() == 3) // F,P,D
      {
        auto outputValues1 = Kokkos::subview(outputValues, range1, Kokkos::ALL(), Kokkos::ALL());
        auto outputValues2 = Kokkos::subview(outputValues, range2, Kokkos::ALL(), Kokkos::ALL());
        
        basis1_.getValues(outputValues1, inputPoints, operatorType);
        basis2_.getValues(outputValues2, inputPoints, operatorType);
      }
      else if (outputValues.rank() == 4) // F,P,D,D
      {
        auto outputValues1 = Kokkos::subview(outputValues, range1, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto outputValues2 = Kokkos::subview(outputValues, range2, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        
        basis1_.getValues(outputValues1, inputPoints, operatorType);
        basis2_.getValues(outputValues2, inputPoints, operatorType);
      }
      else if (outputValues.rank() == 5) // F,P,D,D,D
      {
        auto outputValues1 = Kokkos::subview(outputValues, range1, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto outputValues2 = Kokkos::subview(outputValues, range2, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        
        basis1_.getValues(outputValues1, inputPoints, operatorType);
        basis2_.getValues(outputValues2, inputPoints, operatorType);
      }
      else if (outputValues.rank() == 6) // F,P,D,D,D,D
      {
        auto outputValues1 = Kokkos::subview(outputValues, range1, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto outputValues2 = Kokkos::subview(outputValues, range2, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        
        basis1_.getValues(outputValues1, inputPoints, operatorType);
        basis2_.getValues(outputValues2, inputPoints, operatorType);
      }
      else if (outputValues.rank() == 7) // F,P,D,D,D,D,D
      {
        auto outputValues1 = Kokkos::subview(outputValues, range1, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto outputValues2 = Kokkos::subview(outputValues, range2, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        
        basis1_.getValues(outputValues1, inputPoints, operatorType);
        basis2_.getValues(outputValues2, inputPoints, operatorType);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported outputValues rank");
      }
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DirectSumBasis_h */
