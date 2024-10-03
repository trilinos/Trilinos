// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_DirectSumBasis.hpp
    \brief  Implementation of a basis that is the direct sum of two other bases.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_DirectSumBasis_h
#define Intrepid2_DirectSumBasis_h

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
  template<typename BasisBaseClass>
  class Basis_DirectSumBasis : public BasisBaseClass
  {
  public:
    using BasisBase = BasisBaseClass;
    using BasisPtr  = Teuchos::RCP<BasisBase>;
    
    using DeviceType      = typename BasisBase::DeviceType;
    using ExecutionSpace  = typename BasisBase::ExecutionSpace;
    using OutputValueType = typename BasisBase::OutputValueType;
    using PointValueType  = typename BasisBase::PointValueType;
    
    using OrdinalTypeArray1DHost = typename BasisBase::OrdinalTypeArray1DHost;
    using OrdinalTypeArray2DHost = typename BasisBase::OrdinalTypeArray2DHost;
    using OutputViewType         = typename BasisBase::OutputViewType;
    using PointViewType          = typename BasisBase::PointViewType;
    using ScalarViewType         = typename BasisBase::ScalarViewType;
  protected:
    BasisPtr basis1_;
    BasisPtr basis2_;
    
    std::string name_;
  public:
    /** \brief  Constructor.
        \param [in] basis1 - the instance of Basis1
        \param [in] basis2 - the instance of Basis2
     */
    Basis_DirectSumBasis(BasisPtr basis1, BasisPtr basis2)
    :
    basis1_(basis1),basis2_(basis2)
    {
      INTREPID2_TEST_FOR_EXCEPTION(basis1->getBasisType() != basis2->getBasisType(), std::invalid_argument, "basis1 and basis2 must agree in basis type");
      INTREPID2_TEST_FOR_EXCEPTION(basis1->getBaseCellTopology().getKey() != basis2->getBaseCellTopology().getKey(),
                                 std::invalid_argument, "basis1 and basis2 must agree in cell topology");
      INTREPID2_TEST_FOR_EXCEPTION(basis1->getNumTensorialExtrusions() != basis2->getNumTensorialExtrusions(),
                                   std::invalid_argument, "basis1 and basis2 must agree in number of tensorial extrusions");
      INTREPID2_TEST_FOR_EXCEPTION(basis1->getCoordinateSystem() != basis2->getCoordinateSystem(),
                                 std::invalid_argument, "basis1 and basis2 must agree in coordinate system");
      
      this->basisCardinality_  = basis1->getCardinality() + basis2->getCardinality();
      this->basisDegree_       = std::max(basis1->getDegree(), basis2->getDegree());
      
      {
        std::ostringstream basisName;
        basisName << basis1->getName() << " + " << basis2->getName();
        name_ = basisName.str();
      }
      
      this->basisCellTopologyKey_ = basis1->getBaseCellTopology().getKey();
      this->basisType_            = basis1->getBasisType();
      this->basisCoordinates_     = basis1->getCoordinateSystem();

      if (this->basisType_ == BASIS_FEM_HIERARCHICAL)
      {
        int degreeLength = basis1_->getPolynomialDegreeLength();
        INTREPID2_TEST_FOR_EXCEPTION(degreeLength != basis2_->getPolynomialDegreeLength(), std::invalid_argument, "Basis1 and Basis2 must agree on polynomial degree length");
        
        this->fieldOrdinalPolynomialDegree_   = OrdinalTypeArray2DHost("DirectSumBasis degree lookup",    this->basisCardinality_,degreeLength);
        this->fieldOrdinalH1PolynomialDegree_ = OrdinalTypeArray2DHost("DirectSumBasis H^1 degree lookup",this->basisCardinality_,degreeLength);
        // our field ordinals start with basis1_; basis2_ follows
        for (int fieldOrdinal1=0; fieldOrdinal1<basis1_->getCardinality(); fieldOrdinal1++)
        {
          int fieldOrdinal = fieldOrdinal1;
          auto polynomialDegree   = basis1->getPolynomialDegreeOfField(fieldOrdinal1);
          auto polynomialH1Degree = basis1->getH1PolynomialDegreeOfField(fieldOrdinal1);
          for (int d=0; d<degreeLength; d++)
          {
            this->fieldOrdinalPolynomialDegree_  (fieldOrdinal,d) = polynomialDegree(d);
            this->fieldOrdinalH1PolynomialDegree_(fieldOrdinal,d) = polynomialH1Degree(d);
          }
        }
        for (int fieldOrdinal2=0; fieldOrdinal2<basis2_->getCardinality(); fieldOrdinal2++)
        {
          int fieldOrdinal = basis1->getCardinality() + fieldOrdinal2;
          
          auto polynomialDegree   = basis2->getPolynomialDegreeOfField(fieldOrdinal2);
          auto polynomialH1Degree = basis2->getH1PolynomialDegreeOfField(fieldOrdinal2);
          for (int d=0; d<degreeLength; d++)
          {
            this->fieldOrdinalPolynomialDegree_  (fieldOrdinal,d) = polynomialDegree(d);
            this->fieldOrdinalH1PolynomialDegree_(fieldOrdinal,d) = polynomialH1Degree(d);
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
        
        shards::CellTopology cellTopo(getCellTopologyData(this->basisCellTopologyKey_));
        
        unsigned spaceDim  = cellTopo.getDimension();
        
        ordinal_type basis2Offset = basis1_->getCardinality();
                
        for (unsigned d=0; d<=spaceDim; d++)
        {
          unsigned subcellCount = cellTopo.getSubcellCount(d);
          for (unsigned subcellOrdinal=0; subcellOrdinal<subcellCount; subcellOrdinal++)
          {
            ordinal_type subcellDofCount1 = basis1->getDofCount(d, subcellOrdinal);
            ordinal_type subcellDofCount2 = basis2->getDofCount(d, subcellOrdinal);
            
            ordinal_type subcellDofCount = subcellDofCount1 + subcellDofCount2;
            for (ordinal_type localDofID=0; localDofID<subcellDofCount; localDofID++)
            {
              ordinal_type fieldOrdinal;
              if (localDofID < subcellDofCount1)
              {
                // first basis: field ordinal matches the basis1 ordinal
                fieldOrdinal = basis1_->getDofOrdinal(d, subcellOrdinal, localDofID);
              }
              else
              {
                // second basis: field ordinal is offset by basis1 cardinality
                fieldOrdinal = basis2Offset + basis2_->getDofOrdinal(d, subcellOrdinal, localDofID - subcellDofCount1);
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
    
    /** \brief Allocate BasisValues container suitable for passing to the getValues() variant that takes a TensorPoints container as argument.
     
        The default implementation employs a trivial tensor-product structure, for compatibility across all bases.  Subclasses that have tensor-product structure
        should override.  Note that only the basic exact-sequence operators are supported at the moment: VALUE, GRAD, DIV, CURL.
     */
    virtual BasisValues<OutputValueType,DeviceType> allocateBasisValues( TensorPoints<PointValueType,DeviceType> points, const EOperator operatorType = OPERATOR_VALUE) const override
    {
      BasisValues<OutputValueType,DeviceType> basisValues1 = basis1_->allocateBasisValues(points, operatorType);
      BasisValues<OutputValueType,DeviceType> basisValues2 = basis2_->allocateBasisValues(points, operatorType);
      
      const int numScalarFamilies1 = basisValues1.numTensorDataFamilies();
      if (numScalarFamilies1 > 0)
      {
        // then both basis1 and basis2 should be scalar-valued; check that for basis2:
        const int numScalarFamilies2 = basisValues2.numTensorDataFamilies();
        INTREPID2_TEST_FOR_EXCEPTION(basisValues2.numTensorDataFamilies() <=0, std::invalid_argument, "When basis1 has scalar value, basis2 must also");
        std::vector< TensorData<OutputValueType,DeviceType> > scalarFamilies(numScalarFamilies1 + numScalarFamilies2);
        for (int i=0; i<numScalarFamilies1; i++)
        {
          scalarFamilies[i] = basisValues1.tensorData(i);
        }
        for (int i=0; i<numScalarFamilies2; i++)
        {
          scalarFamilies[i+numScalarFamilies1] = basisValues2.tensorData(i);
        }
        return BasisValues<OutputValueType,DeviceType>(scalarFamilies);
      }
      else
      {
        // then both basis1 and basis2 should be vector-valued; check that:
        INTREPID2_TEST_FOR_EXCEPTION(!basisValues1.vectorData().isValid(), std::invalid_argument, "When basis1 does not have tensorData() defined, it must have a valid vectorData()");
        INTREPID2_TEST_FOR_EXCEPTION(basisValues2.numTensorDataFamilies() > 0, std::invalid_argument, "When basis1 has vector value, basis2 must also");
        
        const auto & vectorData1 = basisValues1.vectorData();
        const auto & vectorData2 = basisValues2.vectorData();
        
        const int numFamilies1  = vectorData1.numFamilies();
        const int numComponents = vectorData1.numComponents();
        INTREPID2_TEST_FOR_EXCEPTION(numComponents != vectorData2.numComponents(), std::invalid_argument, "basis1 and basis2 must agree on the number of components in each vector");
        const int numFamilies2 = vectorData2.numFamilies();
        
        const int numFamilies = numFamilies1 + numFamilies2;
        std::vector< std::vector<TensorData<OutputValueType,DeviceType> > > vectorComponents(numFamilies, std::vector<TensorData<OutputValueType,DeviceType> >(numComponents));
        
        for (int i=0; i<numFamilies1; i++)
        {
          for (int j=0; j<numComponents; j++)
          {
            vectorComponents[i][j] = vectorData1.getComponent(i,j);
          }
        }
        for (int i=0; i<numFamilies2; i++)
        {
          for (int j=0; j<numComponents; j++)
          {
            vectorComponents[i+numFamilies1][j] = vectorData2.getComponent(i,j);
          }
        }
        VectorData<OutputValueType,DeviceType> vectorData(vectorComponents);
        return BasisValues<OutputValueType,DeviceType>(vectorData);
      }
    }
    
    /** \brief  Fills in spatial locations (coordinates) of degrees of freedom (nodes) on the reference cell
        \param [out] dofCoords - the container into which to place the degrees of freedom.
     
     dofCoords should have shape (F,D), where the field dimension matches the cardinality of the basis, and D is the
     spatial dimension of the topology on which the basis is defined.
     
     Note that getDofCoords() is not supported by all bases; in particular, hierarchical bases do not generally support this.
     */
    virtual void getDofCoords( ScalarViewType dofCoords ) const override {
      const int basisCardinality1 = basis1_->getCardinality();
      const int basisCardinality2 = basis2_->getCardinality();
      const int basisCardinality  = basisCardinality1 + basisCardinality2;

      auto dofCoords1 = Kokkos::subview(dofCoords, std::make_pair(0,basisCardinality1),                Kokkos::ALL());
      auto dofCoords2 = Kokkos::subview(dofCoords, std::make_pair(basisCardinality1,basisCardinality), Kokkos::ALL());
      
      basis1_->getDofCoords(dofCoords1);
      basis2_->getDofCoords(dofCoords2);
    }
    
    /** \brief  Fills in coefficients of degrees of freedom for Lagrangian basis on the reference cell
        \param [out] dofCoeffs - the container into which to place the degrees of freedom.

     dofCoeffs have shape (F,D) or (F) if D=1, field dimension matches the cardinality of the basis, and D is the
     basis dimension.

     Degrees of freedom coefficients are such that
     \phi_i(dofCoords_(j)) \cdot dofCoeffs_(j)  = \delta_ij,
     where \phi_i are the basis and \delta_ij the Kronecker delta.
     Note that getDofCoeffs() is supported only for Lagrangian bases.
     */
    virtual void getDofCoeffs( ScalarViewType dofCoeffs ) const override {
      const int basisCardinality1 = basis1_->getCardinality();
      const int basisCardinality2 = basis2_->getCardinality();
      const int basisCardinality  = basisCardinality1 + basisCardinality2;

      auto dofCoeffs1 = Kokkos::subview(dofCoeffs, std::make_pair(0,basisCardinality1), Kokkos::ALL());
      auto dofCoeffs2 = Kokkos::subview(dofCoeffs, std::make_pair(basisCardinality1,basisCardinality), Kokkos::ALL());

      basis1_->getDofCoeffs(dofCoeffs1);
      basis2_->getDofCoeffs(dofCoeffs2);
    }


    /** \brief  Returns basis name
     
     \return the name of the basis
     */
    virtual
    const char*
    getName() const override {
      return name_.c_str();
    }
    
    // since the getValues() below only overrides the FEM variants, we specify that
    // we use the base class's getValues(), which implements the FVD variant by throwing an exception.
    // (It's an error to use the FVD variant on this basis.)
    using BasisBase::getValues;
    
    /** \brief  Evaluation of a FEM basis on a <strong>reference cell</strong>, using point and output value containers that allow preservation of tensor-product structure.

        Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
        points in the <strong>reference cell</strong> for which the basis is defined.

        \param  outputValues      [out] - variable rank array with the basis values.  Should be allocated using Basis::allocateBasisValues().
        \param  inputPoints       [in]  - rank-2 array (P,D) with the evaluation points.  This should be allocated using Cubature::allocateCubaturePoints() and filled using Cubature::getCubature().
        \param  operatorType      [in]  - the operator acting on the basis function
     
        This is the preferred getValues() method for TensorBasis and DirectSumBasis and their subclasses.  It allows a reduced memory footprint and optimized integration, etc.
    */
    virtual
    void
    getValues(       BasisValues<OutputValueType,DeviceType> outputValues,
               const TensorPoints<PointValueType,DeviceType>  inputPoints,
               const EOperator operatorType = OPERATOR_VALUE ) const override
    {
      const int fieldStartOrdinal1 = 0;
      const int numFields1         = basis1_->getCardinality();
      const int fieldStartOrdinal2 = numFields1;
      const int numFields2         = basis2_->getCardinality();
      
      auto basisValues1 = outputValues.basisValuesForFields(fieldStartOrdinal1, numFields1);
      auto basisValues2 = outputValues.basisValuesForFields(fieldStartOrdinal2, numFields2);
      
      basis1_->getValues(basisValues1, inputPoints, operatorType);
      basis2_->getValues(basisValues2, inputPoints, operatorType);
    }
    
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
      int cardinality1 = basis1_->getCardinality();
      int cardinality2 = basis2_->getCardinality();
      
      auto range1 = std::make_pair(0,cardinality1);
      auto range2 = std::make_pair(cardinality1,cardinality1+cardinality2);
      if (outputValues.rank() == 2) // F,P
      {
        auto outputValues1 = Kokkos::subview(outputValues, range1, Kokkos::ALL());
        auto outputValues2 = Kokkos::subview(outputValues, range2, Kokkos::ALL());
        
        basis1_->getValues(outputValues1, inputPoints, operatorType);
        basis2_->getValues(outputValues2, inputPoints, operatorType);
      }
      else if (outputValues.rank() == 3) // F,P,D
      {
        auto outputValues1 = Kokkos::subview(outputValues, range1, Kokkos::ALL(), Kokkos::ALL());
        auto outputValues2 = Kokkos::subview(outputValues, range2, Kokkos::ALL(), Kokkos::ALL());
        
        basis1_->getValues(outputValues1, inputPoints, operatorType);
        basis2_->getValues(outputValues2, inputPoints, operatorType);
      }
      else if (outputValues.rank() == 4) // F,P,D,D
      {
        auto outputValues1 = Kokkos::subview(outputValues, range1, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto outputValues2 = Kokkos::subview(outputValues, range2, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        
        basis1_->getValues(outputValues1, inputPoints, operatorType);
        basis2_->getValues(outputValues2, inputPoints, operatorType);
      }
      else if (outputValues.rank() == 5) // F,P,D,D,D
      {
        auto outputValues1 = Kokkos::subview(outputValues, range1, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto outputValues2 = Kokkos::subview(outputValues, range2, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        
        basis1_->getValues(outputValues1, inputPoints, operatorType);
        basis2_->getValues(outputValues2, inputPoints, operatorType);
      }
      else if (outputValues.rank() == 6) // F,P,D,D,D,D
      {
        auto outputValues1 = Kokkos::subview(outputValues, range1, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto outputValues2 = Kokkos::subview(outputValues, range2, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        
        basis1_->getValues(outputValues1, inputPoints, operatorType);
        basis2_->getValues(outputValues2, inputPoints, operatorType);
      }
      else if (outputValues.rank() == 7) // F,P,D,D,D,D,D
      {
        auto outputValues1 = Kokkos::subview(outputValues, range1, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto outputValues2 = Kokkos::subview(outputValues, range2, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        
        basis1_->getValues(outputValues1, inputPoints, operatorType);
        basis2_->getValues(outputValues2, inputPoints, operatorType);
      }
      else
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported outputValues rank");
      }
    }
    
    virtual int getNumTensorialExtrusions() const override
    {
      return basis1_->getNumTensorialExtrusions();
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_DirectSumBasis_h */
