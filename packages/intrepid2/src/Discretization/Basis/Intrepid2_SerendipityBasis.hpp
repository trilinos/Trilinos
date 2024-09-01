// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
//  Intrepid2_SerendipityBasis.hpp
//  intrepid2
//
//  Created by Roberts, Nathan V on 1/4/22.
//

#ifndef Intrepid2_SerendipityBasis_h
#define Intrepid2_SerendipityBasis_h

#include <Intrepid2_CellTopology.hpp>

namespace Intrepid2
{

/** \class  Intrepid2::SerendipityBasis
    \brief  Serendipity Basis, defined as the sub-basis of a provided basis, consisting of basis elements for which tensorial component polynomial orders satisfy the Serendipity criterion.
*/
template<typename BasisBaseClass = void>
class SerendipityBasis
:
public BasisBaseClass
{
public:
  using BasisBase = BasisBaseClass;
  using BasisPtr  = Teuchos::RCP<BasisBase>;

protected:
  BasisPtr fullBasis_;
  
  std::string name_; // name of the basis
  
  int numTensorialExtrusions_; // relative to cell topo returned by getBaseCellTopology().
public:
  using DeviceType = typename BasisBase::DeviceType;
  using ExecutionSpace  = typename BasisBase::ExecutionSpace;
  using OutputValueType = typename BasisBase::OutputValueType;
  using PointValueType  = typename BasisBase::PointValueType;
  
  using OrdinalTypeArray1D     = typename BasisBase::OrdinalTypeArray1D;
  using OrdinalTypeArray1DHost = typename BasisBase::OrdinalTypeArray1DHost;
  using OrdinalTypeArray2DHost = typename BasisBase::OrdinalTypeArray2DHost;
  using OutputViewType         = typename BasisBase::OutputViewType;
  using PointViewType          = typename BasisBase::PointViewType;
protected:
  OrdinalTypeArray1D ordinalMap_; // our field ordinal --> full basis field ordinal
public:
  /** \brief  Constructor.
      \param [in] basis - the full, hierarchical basis
   */
  SerendipityBasis(BasisPtr fullBasis)
  :
  fullBasis_(fullBasis)
  {
    INTREPID2_TEST_FOR_EXCEPTION(fullBasis_->getBasisType() != BASIS_FEM_HIERARCHICAL, std::invalid_argument, "SerendipityBasis only supports full bases whose type is BASIS_FEM_HIERARCHICAL");
    this->basisType_         = fullBasis_->getBasisType(); // BASIS_FEM_HIERARCHICAL
    
    this->functionSpace_ = fullBasis_->getFunctionSpace();
    this->basisDegree_   = fullBasis_->getDegree();
    
    {
      std::ostringstream basisName;
      basisName << "Serendipity(" << fullBasis_->getName() << ")";
      name_ = basisName.str();
    }
    
    using std::vector;
    vector<ordinal_type> fullBasisOrdinals;
    vector<vector<ordinal_type>> polynomialDegreeOfField;
    vector<vector<ordinal_type>> polynomialH1DegreeOfField;
    ordinal_type fullCardinality = fullBasis_->getCardinality();
    ordinal_type basisH1Degree = (this->functionSpace_ == FUNCTION_SPACE_HVOL) ? this->basisDegree_ + 1 : this->basisDegree_;
    for (ordinal_type i=0; i<fullCardinality; i++)
    {
      vector<ordinal_type> fieldDegree   = fullBasis_->getPolynomialDegreeOfFieldAsVector(i);
      vector<ordinal_type> fieldH1Degree = fullBasis_->getH1PolynomialDegreeOfFieldAsVector(i);
      ordinal_type superlinearDegree = 0;
      for (const ordinal_type & p : fieldH1Degree)
      {
        if (p > 1)
        {
          // superlinear; contributes to superlinearDegree
          superlinearDegree += p;
        }
      }
      if (superlinearDegree <= basisH1Degree)
      {
        // satisfies serendipity criterion
        fullBasisOrdinals.push_back(i);
        polynomialDegreeOfField.push_back(fieldDegree);
        polynomialH1DegreeOfField.push_back(fieldH1Degree);
      }
    }
    this->basisCardinality_  = fullBasisOrdinals.size();
    ordinalMap_ = OrdinalTypeArray1D("serendipity ordinal map",fullBasisOrdinals.size());
    
    auto ordinalMapHost = Kokkos::create_mirror_view(ordinalMap_);
    const ordinal_type fullBasisCardinality = fullBasisOrdinals.size();
    for (ordinal_type i=0; i<fullBasisCardinality; i++)
    {
      ordinalMapHost(i) = fullBasisOrdinals[i];
    }
    Kokkos::deep_copy(ordinalMap_, ordinalMapHost);

    // set cell topology
    this->basisCellTopologyKey_ = fullBasis_->getBaseCellTopology().getKey();
    this->numTensorialExtrusions_ = fullBasis_->getNumTensorialExtrusions();
    const ordinal_type spaceDim = fullBasis_->getBaseCellTopology().getDimension() + numTensorialExtrusions_;
    
    this->basisCoordinates_  = COORDINATES_CARTESIAN;
            
    // fill in degree lookup:
    int degreeSize = fullBasis_->getPolynomialDegreeLength();
    this->fieldOrdinalPolynomialDegree_ = OrdinalTypeArray2DHost("TensorBasis - field ordinal polynomial degree", this->basisCardinality_, degreeSize);
    this->fieldOrdinalH1PolynomialDegree_ = OrdinalTypeArray2DHost("TensorBasis - field ordinal H^1 degree", this->basisCardinality_, degreeSize);
    
    for (ordinal_type fieldOrdinal1 = 0; fieldOrdinal1 < this->basisCardinality_; fieldOrdinal1++)
    {
      for (int d=0; d<degreeSize; d++)
      {
        this->fieldOrdinalPolynomialDegree_  (fieldOrdinal1,d) = polynomialDegreeOfField  [fieldOrdinal1][d];
        this->fieldOrdinalH1PolynomialDegree_(fieldOrdinal1,d) = polynomialH1DegreeOfField[fieldOrdinal1][d];
      }
    }
    
    // build tag view
    const auto & cardinality = this->basisCardinality_;

    const ordinal_type tagSize  = 4;  // size of DoF tag, i.e., number of fields in the tag
    const ordinal_type posScDim = 0;  // position in the tag, counting from 0, of the subcell dim
    const ordinal_type posScOrd = 1;  // position in the tag, counting from 0, of the subcell ordinal
    const ordinal_type posDfOrd = 2;  // position in the tag, counting from 0, of DoF ordinal relative to the subcell
    const ordinal_type posDfCnt = 3;  // position in the tag, counting from 0, of DoF count for the subcell

    OrdinalTypeArray1DHost tagView("tag view", cardinality*tagSize);
    auto cellTopo = CellTopology::cellTopology(fullBasis_->getBaseCellTopology(), numTensorialExtrusions_);
    const OrdinalTypeArray2DHost fullBasisOrdinalToTag = fullBasis->getAllDofTags();
    
    using std::map;
    vector<vector<ordinal_type>> subcellDofCount(spaceDim+1);   // outer: dimension of subcell; inner: subcell ordinal.  Entry: number of dofs.
    vector<vector<ordinal_type>> subcellDofOrdinal(spaceDim+1); // outer: dimension of subcell; inner: subcell ordinal.  Entry: ordinal relative to subcell.
    for (ordinal_type d=0; d<=spaceDim; d++)
    {
      const ordinal_type numSubcells = cellTopo->getSubcellCount(d);
      subcellDofCount[d]   = vector<ordinal_type>(numSubcells,0);
      subcellDofOrdinal[d] = vector<ordinal_type>(numSubcells,0);
    }
    for (ordinal_type fieldOrdinal=0; fieldOrdinal<fullBasisCardinality; fieldOrdinal++)
    {
      const ordinal_type fullFieldOrdinal = fullBasisOrdinals[fieldOrdinal];
      const ordinal_type subcellDim   = fullBasisOrdinalToTag(fullFieldOrdinal,posScDim);
      const ordinal_type subcellOrd   = fullBasisOrdinalToTag(fullFieldOrdinal,posScOrd);
      subcellDofCount[subcellDim][subcellOrd]++;
    }
    for (ordinal_type fieldOrdinal=0; fieldOrdinal<fullBasisCardinality; fieldOrdinal++)
    {
      const ordinal_type fullFieldOrdinal = fullBasisOrdinals[fieldOrdinal];
      const ordinal_type subcellDim   = fullBasisOrdinalToTag(fullFieldOrdinal,posScDim);
      const ordinal_type subcellOrd   = fullBasisOrdinalToTag(fullFieldOrdinal,posScOrd);
      const ordinal_type subcellDfCnt = subcellDofCount[subcellDim][subcellOrd];
      const ordinal_type subcellDfOrd = subcellDofOrdinal[subcellDim][subcellOrd]++;
      
      tagView(tagSize*fieldOrdinal + posScDim) = subcellDim;   // subcellDim
      tagView(tagSize*fieldOrdinal + posScOrd) = subcellOrd;   // subcell ordinal
      tagView(tagSize*fieldOrdinal + posDfOrd) = subcellDfOrd; // ordinal of the specified DoF relative to the subcell
      tagView(tagSize*fieldOrdinal + posDfCnt) = subcellDfCnt; // total number of DoFs associated with the subcell
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
  
  virtual int getNumTensorialExtrusions() const override
  {
    return numTensorialExtrusions_;
  }
  
  /** \brief Allocate BasisValues container suitable for passing to the getValues() variant that takes a TensorPoints container as argument.
   
      The basic exact-sequence operators are supported (VALUE, GRAD, DIV, CURL), as are the Dn operators (OPERATOR_D1 through OPERATOR_D10).
   */
  virtual BasisValues<OutputValueType,DeviceType> allocateBasisValues( TensorPoints<PointValueType,DeviceType> points, const EOperator operatorType = OPERATOR_VALUE) const override
  {
    auto basisValues = fullBasis_->allocateBasisValues(points,operatorType);
    basisValues.setOrdinalFilter(this->ordinalMap_);
    return basisValues;
  }
  
  // since the getValues() below only overrides the FEM variant, we specify that
  // we use the base class's getValues(), which implements the FVD variant by throwing an exception.
  // (It's an error to use the FVD variant on this basis.)
  using BasisBase::getValues;

  /** \brief  Returns basis name
   
   \return the name of the basis
   */
  virtual
  const char*
  getName() const override {
    return name_.c_str();
  }
  
  /** \brief  Evaluation of a FEM basis on a <strong>reference cell</strong>, using point and output value containers that allow preservation of tensor-product structure.

      Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
      points in the <strong>reference cell</strong> for which the basis is defined.

      \param  outputValues      [out] - variable rank array with the basis values.  Should be allocated using Basis::allocateBasisValues().
      \param  inputPoints       [in]  - rank-2 array (P,D) with the evaluation points.  This should be allocated using Cubature::allocateCubaturePoints() and filled using Cubature::getCubature().
      \param  operatorType      [in]  - the operator acting on the basis function
   
      This is the preferred getValues() method for SerendipityBasis.  It allows a reduced memory footprint and optimized integration, etc.
  */
  virtual
  void
  getValues(       BasisValues<OutputValueType,DeviceType> outputValues,
             const TensorPoints<PointValueType,DeviceType>  inputPoints,
             const EOperator operatorType = OPERATOR_VALUE ) const override
  {
    fullBasis_->getValues(outputValues, inputPoints, operatorType);
    outputValues.setOrdinalFilter(this->ordinalMap_);
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
  virtual
  void getValues( OutputViewType outputValues, const PointViewType  inputPoints,
                 const EOperator operatorType = OPERATOR_VALUE ) const override
  {
    INTREPID2_TEST_FOR_EXCEPTION( true, std::logic_error,
                                  ">>> ERROR (Basis::getValues): this method is not supported by SerendipityBasis (use the getValues() method that accepts a BasisValues object instead).");
  }
  
  /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.
   
      \return Pointer to the new Basis object.
   */
  virtual HostBasisPtr<OutputValueType, PointValueType>
  getHostBasis() const override {
    using HostBasisBase = Basis<typename Kokkos::HostSpace::device_type, OutputValueType, PointValueType>;
    using HostBasis = SerendipityBasis<HostBasisBase>;
    auto hostBasis = Teuchos::rcp(new HostBasis(fullBasis_->getHostBasis()));
    return hostBasis;
  }
  
  /** \brief Returns a pointer to the underlying full basis.
   
      \return Pointer to the underlying full basis.
   */
  BasisPtr getUnderlyingBasis() const
  {
    return fullBasis_;
  }
  
  /** \brief Returns the ordinal map from the Serendipity basis ordinal to the ordinal in the underlying full basis.
   
      \return Ordinal map from the Serendipity basis ordinal to the ordinal in the underlying full basis.  (Indices to the container are Serendipity basis ordinals; values are full basis ordinals.)
   */
  OrdinalTypeArray1D ordinalMap() const
  {
    return ordinalMap_;
  }
}; // SerendipityBasis

} // namespace Intrepid2
#endif /* Intrepid2_SerendipityBasis_h */
