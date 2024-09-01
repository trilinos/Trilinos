// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_Basis.hpp
    \brief  Header file for the abstract base class Intrepid2::Basis.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
 */

#ifndef __INTREPID2_BASIS_HPP__
#define __INTREPID2_BASIS_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_BasisValues.hpp"
#include "Intrepid2_CellTopologyTags.hpp"
#include "Intrepid2_TensorPoints.hpp"
#include "Shards_CellTopology.hpp"
#include "Intrepid2_CellData.hpp"
#include <Teuchos_RCPDecl.hpp>
#include <Kokkos_Core.hpp>

#include <vector>

namespace Intrepid2 {

template<typename DeviceType = void,
         typename OutputType = double,
         typename PointType = double>
class Basis;

  /**  \brief Basis Pointer
    */
template <typename DeviceType = void, typename OutputType = double, typename PointType = double>
using BasisPtr = Teuchos::RCP<Basis<DeviceType,OutputType,PointType> >;

/** \brief Pointer to a Basis whose device type is on the host (Kokkos::HostSpace::device_type), allowing host access to input and output views, and ensuring host execution of basis evaluation.
 */
template <typename OutputType = double, typename PointType = double>
using HostBasisPtr = BasisPtr<typename Kokkos::HostSpace::device_type, OutputType, PointType>;

  /** \class  Intrepid2::Basis
      \brief  An abstract base class that defines interface for concrete basis implementations for
              Finite Element (FEM) and Finite Volume/Finite Difference (FVD) discrete spaces.

              A FEM basis spans a discrete space whose type can be either COMPLETE or INCOMPLETE.
              FEM basis functions are always defined on a reference cell and are dual to a unisolvent
              set of degrees-of-freedom (DoF). FEM basis requires cell topology with a reference cell.

              An FVD basis spans a discrete space whose type is typically BROKEN. The basis functions
              are defined directly on the physical cell and are dual to a set of DoFs on that cell.
              As a result, FVD bases require the vertex coordinates of the physical cell but the cell
              itself is not required to have a reference cell.

              Every DoF and its corresponding basis function from a given FEM or FVD basis set is
              assigned an ordinal number which which specifies its numerical position in the DoF set,
              and a 4-field DoF tag whose first 3 fields establish association between the DoF and a
              subcell of particular dimension, and the last field gives the total number of basis
              functions associated with that subcell; see Section \ref basis_dof_tag_ord_sec for details.

              Note the use of Kokkos::LayoutStride as the layout for various View definitions (including the
              argument for \ref getValues() ).  A method expecting a LayoutStride view can accept, thanks to
              some conversion mechanisms defined elsewhere, a LayoutLeft view (the default for Cuda views in
              Kokkos) or a LayoutRight view (the default on most other platforms).  This does introduce some
              additional complexity when Views need to be allocated for temporary storage; see the method
              \ref getMatchingViewWithLabel() provided in \ref Intrepid_Utils.hpp.
   
      \remark To limit memory use by factory-type objects (basis factories will be included in future
              releases of Intrepid2), tag data is not initialized by basis ctors,
              instead, whenever a function that requires tag data is invoked for a first time, it calls
              initializeTags() to fill <var>ordinalToTag_</var> and <var>tagToOrdinal_</var>. Because
              tag data is basis-specific, every concrete basis class requires its own implementation
              of initializeTags().

      \todo  restore test for inclusion of reference points in their resective reference cells in
             getValues_HGRAD_Args, getValues_CURL_Args, getValues_DIV_Args
  */
  template<typename Device,
           typename outputValueType,
           typename pointValueType>
  class Basis {
  public:
    /**  \brief (Kokkos) Device type on which Basis is templated.  Does not necessarily return true for Kokkos::is_device (may be Kokkos::Serial, for example).
     */
    using DeviceType = Device;
    
    /**  \brief (Kokkos) Execution space for basis.
     */
    using ExecutionSpace  = typename DeviceType::execution_space;
    
    
    /**  \brief Output value type for basis; default is double.
     */
    using OutputValueType = outputValueType;
    
    /**  \brief Point value type for basis; default is double.
     */
    using PointValueType  = pointValueType;
    
    /**  \brief View type for ordinal
    */
    using OrdinalViewType = Kokkos::View<ordinal_type,DeviceType>;

    /**  \brief View for basis type
    */
    using EBasisViewType = Kokkos::View<EBasis,DeviceType>;

    /**  \brief View for coordinate system type
    */
    using ECoordinatesViewType = Kokkos::View<ECoordinates,DeviceType>;

    // ** tag interface
    //  - tag interface is not decorated with Kokkos inline so it should be allocated on hostspace

    /**  \brief View type for 1d host array
    */
    using OrdinalTypeArray1DHost = Kokkos::View<ordinal_type*,typename ExecutionSpace::array_layout,Kokkos::HostSpace>;

    /**  \brief View type for 2d host array
    */
    using OrdinalTypeArray2DHost = Kokkos::View<ordinal_type**,typename ExecutionSpace::array_layout,Kokkos::HostSpace>;

    /**  \brief View type for 3d host array
    */
    using OrdinalTypeArray3DHost = Kokkos::View<ordinal_type***,typename ExecutionSpace::array_layout,Kokkos::HostSpace>;

    /**  \brief View type for 1d host array
    */
    using OrdinalTypeArrayStride1DHost = Kokkos::View<ordinal_type*, Kokkos::LayoutStride, Kokkos::HostSpace>;

    /**  \brief View type for 1d device array
    */
    using OrdinalTypeArray1D = Kokkos::View<ordinal_type*,DeviceType>;

    /**  \brief View type for 2d device array
    */
    using OrdinalTypeArray2D = Kokkos::View<ordinal_type**,DeviceType>;

    /**  \brief View type for 3d device array
    */
    using OrdinalTypeArray3D = Kokkos::View<ordinal_type***,DeviceType>;

    /**  \brief View type for 1d device array 
    */
    using OrdinalTypeArrayStride1D = Kokkos::View<ordinal_type*, Kokkos::LayoutStride, DeviceType>;

    /**  \brief Scalar type for point values
    */
    typedef typename ScalarTraits<pointValueType>::scalar_type scalarType;
  protected:

    /** \brief  Cardinality of the basis, i.e., the number of basis functions/degrees-of-freedom
     */
    ordinal_type basisCardinality_;

    /** \brief  Degree of the largest complete polynomial space that can be represented by the basis
     */
    ordinal_type basisDegree_;

    /** \brief Identifier of the base topology of the cells for which the basis is defined. See
         the <a href="https://trilinos.org/packages/shards/">Shards</a> package
         for definition of base cell topology.  For TensorBasis subclasses, by default this the cell topology that is extruded (i.e., it is a lower-dimensional CellTopology than
         the space on which the tensor basis is defined).  This allows tensor bases to be defined in higher dimensions than shards::CellTopology supports.  TensorBasis subclasses can
         opt to use an equivalent shards CellTopology for base topology, as well as using Intrepid2's tagging for tensor bases in dimensions up to 3, by calling
         TensorBasis::setShardsTopologyAndTags().
    */
    unsigned basisCellTopologyKey_;

    /** \brief  Type of the basis
     */
    EBasis basisType_;

    /** \brief  The coordinate system for which the basis is defined
     */
    ECoordinates basisCoordinates_;
    
    /** \brief  The function space in which the basis is defined
     */
    EFunctionSpace functionSpace_ = FUNCTION_SPACE_MAX;
    
    /** \brief  "true" if <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> have been initialized
     */
    //Kokkos::View<bool,DeviceType> basisTagsAreSet_;

    /** \brief  DoF ordinal to tag lookup table.

        Rank-2 array with dimensions (basisCardinality_, 4) containing the DoF tags. This array
        is left empty at instantiation and filled by initializeTags() only when tag data is
        requested.

        \li     ordinalToTag_[DofOrd][0] = dim. of the subcell associated with the specified DoF
        \li     ordinalToTag_[DofOrd][1] = ordinal of the subcell defined in the cell topology
        \li     ordinalToTag_[DodOrd][2] = ordinal of the specified DoF relative to the subcell
        \li     ordinalToTag_[DofOrd][3] = total number of DoFs associated with the subcell
    */
    OrdinalTypeArray2DHost ordinalToTag_;

    /** \brief  DoF tag to ordinal lookup table.

        Rank-3 array with dimensions (maxScDim + 1, maxScOrd + 1, maxDfOrd + 1), i.e., the
        columnwise maximums of the 1st three columns in the DoF tag table for the basis plus 1.
        For every triple (subscDim, subcOrd, subcDofOrd) that is valid DoF tag data this array
        stores the corresponding DoF ordinal. If the triple does not correspond to tag data,
        the array stores -1. This array is left empty at instantiation and filled by
        initializeTags() only when tag data is requested.

        \li     tagToOrdinal_[subcDim][subcOrd][subcDofOrd] = Degree-of-freedom ordinal
    */
    OrdinalTypeArray3DHost tagToOrdinal_;

    /** \brief  Fills <var>ordinalToTag_</var> and <var>tagToOrdinal_</var> by basis-specific tag data

      \param  tagToOrdinal     [out]  - Lookup table for the DoF's ordinal by its tag
      \param  ordinalToTag     [out]  - Lookup table for the DoF's tag by its ordinal
      \param  tags             [in]   - a set of basis-dependent tags in flat (rank-1) array format.
      \param  basisCard        [in]   - cardinality of the basis
      \param  tagSize          [in]   - number of fields in a DoF tag
      \param  posScDim         [in]   - position in the tag, counting from 0, of the subcell dim
      \param  posScOrd         [in]   - position in the tag, counting from 0, of the subcell ordinal
      \param  posDfOrd         [in]   - position in the tag, counting from 0, of DoF ordinal relative to the subcell
    */
    template<typename OrdinalTypeView3D,
             typename OrdinalTypeView2D,
             typename OrdinalTypeView1D>
    void setOrdinalTagData(       OrdinalTypeView3D &tagToOrdinal,
                                  OrdinalTypeView2D &ordinalToTag,
                            const OrdinalTypeView1D  tags,
                            const ordinal_type       basisCard,
                            const ordinal_type       tagSize,
                            const ordinal_type       posScDim,
                            const ordinal_type       posScOrd,
                            const ordinal_type       posDfOrd ) {
      // Create ordinalToTag
      ordinalToTag = OrdinalTypeView2D("ordinalToTag", basisCard, tagSize);

      // Initialize with -1
      Kokkos::deep_copy( ordinalToTag, -1 );

      // Copy tags
      for (ordinal_type i=0;i<basisCard;++i)
        for (ordinal_type j=0;j<tagSize;++j)
          ordinalToTag(i, j) = tags(i*tagSize + j);

      // Find out dimension of tagToOrdinal
      auto maxScDim = 0;  // first dimension of tagToOrdinal
      for (ordinal_type i=0;i<basisCard;++i)
        if (maxScDim < tags(i*tagSize + posScDim))
          maxScDim = tags(i*tagSize + posScDim);
      ++maxScDim;

      auto maxScOrd = 0; // second dimension of tagToOrdinal
      for (ordinal_type i=0;i<basisCard;++i)
        if (maxScOrd < tags(i*tagSize + posScOrd))
          maxScOrd = tags(i*tagSize + posScOrd);
      ++maxScOrd;

      auto maxDfOrd = 0;  // third dimension of tagToOrdinal
      for (ordinal_type i=0;i<basisCard;++i)
        if (maxDfOrd < tags(i*tagSize + posDfOrd))
          maxDfOrd = tags(i*tagSize + posDfOrd);
      ++maxDfOrd;

      // Create tagToOrdinal
      tagToOrdinal = OrdinalTypeView3D("tagToOrdinal", maxScDim, maxScOrd, maxDfOrd);

      // Initialize with -1
      Kokkos::deep_copy( tagToOrdinal, -1 );

      // Overwrite elements of the array corresponding to tags with local DoF Id's, leave all other = -1
      for (ordinal_type i=0;i<basisCard;++i)
        tagToOrdinal(tags(i*tagSize), tags(i*tagSize+1), tags(i*tagSize+2)) = i;
    }

    // dof coords
    /** \brief Coordinates of degrees-of-freedom for basis functions defined in physical space.
     */
    Kokkos::DynRankView<scalarType,DeviceType> dofCoords_;

    // dof coeffs
    /** \brief Coefficients for computing degrees of freedom for Lagrangian basis
        If P is an element of the space spanned by the basis,
        \alpha_i := P(dofCoords_(i)) \cdot dofCoeffs_(i) are the nodal coefficients associated to basis functions i.

        Rank-1 array for scalar basis with dimension (cardinality)
        Rank-2 array for vector basis with dimensions (cardinality, cell dimension)
     */
    Kokkos::DynRankView<scalarType,DeviceType> dofCoeffs_;
    
    /** \brief Polynomial degree for each degree of freedom.  Only defined for hierarchical bases right now.
     The number of entries per degree of freedom in this table depends on the basis type.  For hypercubes,
     this will be the spatial dimension.  We have not yet determined what this will be for simplices beyond 1D;
     there are not yet hierarchical simplicial bases beyond 1D in Intrepid2.
     
     Rank-2 array with dimensions (cardinality, cell dimension)
     */
    OrdinalTypeArray2DHost fieldOrdinalPolynomialDegree_;
    
    /** \brief H^1 polynomial degree for each degree of freedom.  Only defined for hierarchical bases right now.
     The number of entries per degree of freedom in this table depends on the basis type.  For hypercubes,
     this will be the spatial dimension.  We have not yet determined what this will be for simplices beyond 1D;
     there are not yet hierarchical simplicial bases beyond 1D in Intrepid2.
     
     The H^1 polynomial degree is identical to the polynomial degree for H(grad) bases.  For H(vol) bases, it is one
     higher than the polynomial degree.  Since H(div) and H(curl) bases are constructed as products of H(vol) and H(grad)
     bases, the H^1 degree in a given dimension is the H^1 degree for the multiplicand in that dimension.
     
     Rank-2 array with dimensions (cardinality, cell dimension)
     */
    OrdinalTypeArray2DHost fieldOrdinalH1PolynomialDegree_;
  public:

    Basis() = default;
    virtual~Basis() = default;

    // receives input arguments
    /** \brief Dummy array to receive input arguments
    */
    OutputValueType getDummyOutputValue() { return outputValueType(); }

    /** \brief Dummy array to receive input arguments
    */
    PointValueType getDummyPointValue() { return pointValueType(); }

    /** \brief View type for basis value output
    */
    using OutputViewType = Kokkos::DynRankView<OutputValueType,Kokkos::LayoutStride,DeviceType>;

    /** \brief View type for input points
    */
    using PointViewType = Kokkos::DynRankView<PointValueType,Kokkos::LayoutStride,DeviceType>;

    /** \brief View type for scalars 
    */
    using ScalarViewType = Kokkos::DynRankView<scalarType,Kokkos::LayoutStride,DeviceType>;

    /** \brief Allocate a View container suitable for passing to the getValues() variant that accepts Kokkos DynRankViews as arguments (as opposed to the Intrepid2 BasisValues and PointValues containers).
     
        Note that only the basic exact-sequence operators are supported at the moment: VALUE, GRAD, DIV, CURL.
     */
    Kokkos::DynRankView<OutputValueType,DeviceType> allocateOutputView( const int numPoints, const EOperator operatorType = OPERATOR_VALUE) const;
    
    /** \brief Allocate BasisValues container suitable for passing to the getValues() variant that takes a TensorPoints container as argument.
     
        The default implementation employs a trivial tensor-product structure, for compatibility across all bases.  Subclasses that have non-trivial tensor-product structure
        should override.  The basic exact-sequence operators are supported (VALUE, GRAD, DIV, CURL), as are the Dn operators (OPERATOR_D1 through OPERATOR_D10).
     */
    virtual BasisValues<OutputValueType,DeviceType> allocateBasisValues( TensorPoints<PointValueType,DeviceType> points, const EOperator operatorType = OPERATOR_VALUE) const
    {
      const bool operatorIsDk = (operatorType >= OPERATOR_D1) && (operatorType <= OPERATOR_D10);
      const bool operatorSupported = (operatorType == OPERATOR_VALUE) || (operatorType == OPERATOR_GRAD) || (operatorType == OPERATOR_CURL) || (operatorType == OPERATOR_DIV) || operatorIsDk;
      INTREPID2_TEST_FOR_EXCEPTION(!operatorSupported, std::invalid_argument, "operator is not supported by allocateBasisValues");
      
      const int numPoints = points.extent_int(0);
      
      using Scalar = OutputValueType;
      
      auto dataView = allocateOutputView(numPoints, operatorType);
      Data<Scalar,DeviceType> data(dataView);
      
      bool useVectorData = (rank(dataView) == 3);
      
      if (useVectorData)
      {
        VectorData<Scalar,DeviceType> vectorData(data);
        return BasisValues<Scalar,DeviceType>(vectorData);
      }
      else
      {
        TensorData<Scalar,DeviceType> tensorData(data);
        return BasisValues<Scalar,DeviceType>(tensorData);
      }
    }

    /** \brief  Evaluation of a FEM basis on a <strong>reference cell</strong>.

        Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
        points in the <strong>reference cell</strong> for which the basis is defined.

        \param  space             [in]  - execution space instance
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
    void
    getValues( const ExecutionSpace& /* space */,
                     OutputViewType /* outputValues */,
               const PointViewType  /* inputPoints */,
               const EOperator /* operatorType */ = OPERATOR_VALUE ) const {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::logic_error,
                                    ">>> ERROR (Basis::getValues): this method (FEM) is not supported or should be overridden accordingly by derived classes.");
    }

    //! @overload For backward compatibility.
    virtual void getValues(
            OutputViewType outputValues,
      const PointViewType  inputPoints,
      const EOperator      operatorType  = OPERATOR_VALUE
    ) const {
      this->getValues(ExecutionSpace{}, outputValues, inputPoints, operatorType);
    }

    /** \brief  Evaluation of a FEM basis on a <strong>reference cell</strong>, using point and output value containers that allow preservation of tensor-product structure.

        Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
        points in the <strong>reference cell</strong> for which the basis is defined.

        \param  outputValues      [out] - variable rank array with the basis values.  Should be allocated using Basis::allocateBasisValues().
        \param  inputPoints       [in]  - rank-2 array (P,D) with the evaluation points.  This should be allocated using Cubature::allocateCubaturePoints() and filled using Cubature::getCubature().
        \param  operatorType      [in]  - the operator acting on the basis function
     
        The default implementation does not take advantage of any tensor-product structure; subclasses with tensor-product support must override allocateBasisValues() and this getValues() method.
    */
    virtual
    void
    getValues(       BasisValues<OutputValueType,DeviceType> outputValues,
               const TensorPoints<PointValueType,DeviceType>  inputPoints,
               const EOperator operatorType = OPERATOR_VALUE ) const {
      // note the extra allocation/copy here (this is one reason, among several, to override this method):
      auto rawExpandedPoints = inputPoints.allocateAndFillExpandedRawPointView();
      
      OutputViewType rawOutputView;
      Data<OutputValueType,DeviceType> outputData;
      if (outputValues.numTensorDataFamilies() > 0)
      {
        INTREPID2_TEST_FOR_EXCEPTION(outputValues.tensorData(0).numTensorComponents() != 1, std::invalid_argument, "default implementation of getValues() only supports outputValues with trivial tensor-product structure");
        outputData = outputValues.tensorData().getTensorComponent(0);
      }
      else if (outputValues.vectorData().isValid())
      {
        INTREPID2_TEST_FOR_EXCEPTION(outputValues.vectorData().numComponents() != 1, std::invalid_argument, "default implementation of getValues() only supports outputValues with trivial tensor-product structure");
        INTREPID2_TEST_FOR_EXCEPTION(outputValues.vectorData().getComponent(0).numTensorComponents() != 1, std::invalid_argument, "default implementation of getValues() only supports outputValues with trivial tensor-product structure");
        outputData = outputValues.vectorData().getComponent(0).getTensorComponent(0);
      }
      
      this->getValues(outputData.getUnderlyingView(), rawExpandedPoints, operatorType);
    }

    /** \brief  Evaluation of an FVD basis evaluation on a <strong>physical cell</strong>.

        Returns values of <var>operatorType</var> acting on FVD basis functions for a set of
        points in the <strong>physical cell</strong> for which the FVD basis is defined.

        \param  outputValues      [out] - variable rank array with the basis values
        \param  inputPoints       [in]  - rank-2 array (P,D) with the evaluation points
        \param  cellVertices      [in]  - rank-2 array (V,D) with the vertices of the physical cell
        \param  operatorType      [in]  - the operator acting on the basis functions

        \remark For rank and dimension specifications of the output array see Section
        \ref basis_md_array_sec. Dimensions of <var>ArrayScalar</var> arguments are checked
        at runtime if HAVE_INTREPID2_DEBUG is defined.

        \remarks A typical FVD basis spans a BROKEN discrete space which is only piecewise smooth. For
        example, it could be a piecewise constant space defined with respect to a partition
        of the cell into simplices. Because differential operators are not meaningful for such
        spaces, the default operator type in this method is set to OPERATOR_VALUE.
    */
    virtual
    void
    getValues(        OutputViewType /* outputValues */,
                const PointViewType  /* inputPoints */,
                const PointViewType  /* cellVertices */,
                const EOperator /* operatorType */ = OPERATOR_VALUE ) const {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::logic_error,
                                    ">>> ERROR (Basis::getValues): this method (FVM) is not supported or should be overridden accordingly by derived classes.");
    }


    /** \brief  Returns spatial locations (coordinates) of degrees of freedom on the reference cell
     */

    virtual
    void
    getDofCoords( ScalarViewType /* dofCoords */ ) const {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::logic_error,
                                    ">>> ERROR (Basis::getDofCoords): this method is not supported or should be overridden accordingly by derived classes.");
    }

    /** \brief Coefficients for computing degrees of freedom for Lagrangian basis
        If P is an element of the space spanned by the basis,
        \alpha_i := P(dofCoords(i)) \cdot dofCoeffs(i) are the nodal coefficients associated to basis function i.

        Rank-1 array for scalar basis with dimension (cardinality)
        Rank-2 array for vector basis with dimensions (cardinality, cell dimension)
     */

    virtual
    void
    getDofCoeffs( ScalarViewType /* dofCoeffs */ ) const {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::logic_error,
                                    ">>> ERROR (Basis::getDofCoeffs): this method is not supported or should be overridden accordingly by derived classes.");
    }

    /** \brief For hierarchical bases, returns the field ordinals that have at most the specified degree in each dimension.
     Assuming that these are less than or equal to the polynomial orders provided at Basis construction, the corresponding polynomials will form a superset of the Basis of the same type constructed with polynomial orders corresponding to the specified degrees.
     
     \param  degrees      [in] - 1D host ordinal array of length specified by getPolynomialDegreeLength(), indicating what the maximum degree in each dimension should be
     
     \return a 1D host ordinal array containing the ordinals of matching basis functions
     */
    OrdinalTypeArray1DHost getFieldOrdinalsForDegree(OrdinalTypeArray1DHost &degrees) const
    {
      INTREPID2_TEST_FOR_EXCEPTION( basisType_ != BASIS_FEM_HIERARCHICAL, std::logic_error,
                                   ">>> ERROR (Basis::getFieldOrdinalsForDegree): this method is not supported for non-hierarchical bases.");
      int degreeEntryLength     = fieldOrdinalPolynomialDegree_.extent_int(1);
      int requestedDegreeLength = degrees.extent_int(0);
      INTREPID2_TEST_FOR_EXCEPTION(degreeEntryLength != requestedDegreeLength, std::invalid_argument, "length of degrees does not match the entries in fieldOrdinalPolynomialDegree_");
      std::vector<int> fieldOrdinalsVector;
      for (int basisOrdinal=0; basisOrdinal<fieldOrdinalPolynomialDegree_.extent_int(0); basisOrdinal++)
      {
        bool matches = true;
        for (int d=0; d<degreeEntryLength; d++)
        {
          if (fieldOrdinalPolynomialDegree_(basisOrdinal,d) > degrees(d)) matches = false;
        }
        if (matches) fieldOrdinalsVector.push_back(basisOrdinal);
      }
      OrdinalTypeArray1DHost fieldOrdinals("fieldOrdinalsForDegree",fieldOrdinalsVector.size());
      for (unsigned i=0; i<fieldOrdinalsVector.size(); i++)
      {
        fieldOrdinals(i) = fieldOrdinalsVector[i];
      }
      return fieldOrdinals;
    }
    
    /** \brief For hierarchical bases, returns the field ordinals that have at most the specified H^1 degree in each dimension.
     Assuming that these are less than or equal to the polynomial orders provided at Basis construction, the corresponding polynomials will form a superset of the Basis of the same type constructed with polynomial orders corresponding to the specified degrees.
     
     \param  degrees      [in] - 1D host ordinal array of length specified by getPolynomialDegreeLength(), indicating what the maximum degree in each dimension should be
     
     \return a 1D host ordinal array containing the ordinals of matching basis functions
     */
    OrdinalTypeArray1DHost getFieldOrdinalsForH1Degree(OrdinalTypeArray1DHost &degrees) const
    {
      INTREPID2_TEST_FOR_EXCEPTION( basisType_ != BASIS_FEM_HIERARCHICAL, std::logic_error,
                                   ">>> ERROR (Basis::getFieldOrdinalsForDegree): this method is not supported for non-hierarchical bases.");
      int degreeEntryLength     = fieldOrdinalH1PolynomialDegree_.extent_int(1);
      int requestedDegreeLength = degrees.extent_int(0);
      INTREPID2_TEST_FOR_EXCEPTION(degreeEntryLength != requestedDegreeLength, std::invalid_argument, "length of degrees does not match the entries in fieldOrdinalPolynomialDegree_");
      std::vector<int> fieldOrdinalsVector;
      for (int basisOrdinal=0; basisOrdinal<fieldOrdinalH1PolynomialDegree_.extent_int(0); basisOrdinal++)
      {
        bool matches = true;
        for (int d=0; d<degreeEntryLength; d++)
        {
          if (fieldOrdinalH1PolynomialDegree_(basisOrdinal,d) > degrees(d)) matches = false;
        }
        if (matches) fieldOrdinalsVector.push_back(basisOrdinal);
      }
      OrdinalTypeArray1DHost fieldOrdinals("fieldOrdinalsForH1Degree",fieldOrdinalsVector.size());
      for (unsigned i=0; i<fieldOrdinalsVector.size(); i++)
      {
        fieldOrdinals(i) = fieldOrdinalsVector[i];
      }
      return fieldOrdinals;
    }
    
    /** \brief For hierarchical bases, returns the field ordinals that have at most the specified degree in each dimension.
     Assuming that these are less than or equal to the polynomial orders provided at Basis construction, the corresponding polynomials will form a superset of the Basis of the same type constructed with polynomial orders corresponding to the specified degrees.
     

     This variant takes a std::vector of polynomial degrees and returns a std::vector of field ordinals.  It calls the other variant, which uses Kokkos Views on the host.
     
     \param  degrees      [in] - std::vector<int> of length specified by getPolynomialDegreeLength(), indicating what the maximum degree in each dimension should be
     
     \return a std::vector<int> containing the ordinals of matching basis functions
     
     */
    std::vector<int> getFieldOrdinalsForDegree(std::vector<int> &degrees) const
    {
      INTREPID2_TEST_FOR_EXCEPTION( basisType_ != BASIS_FEM_HIERARCHICAL, std::logic_error,
                                   ">>> ERROR (Basis::getFieldOrdinalsForDegree): this method is not supported for non-hierarchical bases.");
      OrdinalTypeArray1DHost degreesView("degrees",degrees.size());
      for (unsigned d=0; d<degrees.size(); d++)
      {
        degreesView(d) = degrees[d];
      }
      auto fieldOrdinalsView = getFieldOrdinalsForDegree(degreesView);
      std::vector<int> fieldOrdinalsVector(fieldOrdinalsView.extent_int(0));
      for (int i=0; i<fieldOrdinalsView.extent_int(0); i++)
      {
        fieldOrdinalsVector[i] = fieldOrdinalsView(i);
      }
      return fieldOrdinalsVector;
    }
    
    /** \brief For hierarchical bases, returns the field ordinals that have at most the specified H^1 degree in each dimension.
     Assuming that these are less than or equal to the polynomial orders provided at Basis construction, the corresponding polynomials will form a superset of the Basis of the same type constructed with polynomial orders corresponding to the specified degrees.

     This variant takes a std::vector of polynomial degrees and returns a std::vector of field ordinals.  It calls the other variant, which uses Kokkos Views on the host.
     
     \param  degrees      [in] - std::vector<int> of length specified by getPolynomialDegreeLength(), indicating what the maximum degree in each dimension should be
     
     \return a std::vector<int> containing the ordinals of matching basis functions
     
     */
    std::vector<int> getFieldOrdinalsForH1Degree(std::vector<int> &degrees) const
    {
      INTREPID2_TEST_FOR_EXCEPTION( basisType_ != BASIS_FEM_HIERARCHICAL, std::logic_error,
                                   ">>> ERROR (Basis::getFieldOrdinalsForDegree): this method is not supported for non-hierarchical bases.");
      OrdinalTypeArray1DHost degreesView("degrees",degrees.size());
      for (unsigned d=0; d<degrees.size(); d++)
      {
        degreesView(d) = degrees[d];
      }
      auto fieldOrdinalsView = getFieldOrdinalsForH1Degree(degreesView);
      std::vector<int> fieldOrdinalsVector(fieldOrdinalsView.extent_int(0));
      for (int i=0; i<fieldOrdinalsView.extent_int(0); i++)
      {
        fieldOrdinalsVector[i] = fieldOrdinalsView(i);
      }
      return fieldOrdinalsVector;
    }
    
    /** \brief For hierarchical bases, returns the polynomial degree (which may have multiple values in higher spatial dimensions) for the specified basis ordinal as a host array.
     
        \param fieldOrdinal     [in] - ordinal of the basis function whose polynomial degree is requested.
     
        \return a 1D host array of length matching getPolynomialDegreeLength(), with the polynomial degree of the basis function in each dimension.
     */
    OrdinalTypeArray1DHost getPolynomialDegreeOfField(int fieldOrdinal) const
    {
      INTREPID2_TEST_FOR_EXCEPTION( basisType_ != BASIS_FEM_HIERARCHICAL, std::logic_error,
                                   ">>> ERROR (Basis::getPolynomialDegreeOfField): this method is not supported for non-hierarchical bases.");
      INTREPID2_TEST_FOR_EXCEPTION(fieldOrdinal < 0, std::invalid_argument, "field ordinal must be non-negative");
      INTREPID2_TEST_FOR_EXCEPTION(fieldOrdinal >= fieldOrdinalPolynomialDegree_.extent_int(0), std::invalid_argument, "field ordinal out of bounds");
      
      int polyDegreeLength = getPolynomialDegreeLength(); 
      OrdinalTypeArray1DHost polyDegree("polynomial degree", polyDegreeLength);
      for (int d=0; d<polyDegreeLength; d++)
      {
        polyDegree(d) = fieldOrdinalPolynomialDegree_(fieldOrdinal,d);
      }
      return polyDegree;
    }
    
    /** \brief For hierarchical bases, returns the polynomial degree (which may have multiple values in higher spatial dimensions) for the specified basis ordinal as a host array.
     
        \param fieldOrdinal     [in] - ordinal of the basis function whose polynomial degree is requested.
     
        \return a 1D host array of length matching getPolynomialDegreeLength(), with the H^1 polynomial degree of the basis function in each dimension.
     */
    OrdinalTypeArray1DHost getH1PolynomialDegreeOfField(int fieldOrdinal) const
    {
      INTREPID2_TEST_FOR_EXCEPTION( basisType_ != BASIS_FEM_HIERARCHICAL, std::logic_error,
                                   ">>> ERROR (Basis::getPolynomialDegreeOfField): this method is not supported for non-hierarchical bases.");
      INTREPID2_TEST_FOR_EXCEPTION(fieldOrdinal < 0, std::invalid_argument, "field ordinal must be non-negative");
      INTREPID2_TEST_FOR_EXCEPTION(fieldOrdinal >= fieldOrdinalH1PolynomialDegree_.extent_int(0), std::invalid_argument, "field ordinal out of bounds");
      
      int polyDegreeLength = getPolynomialDegreeLength();
      OrdinalTypeArray1DHost polyDegree("polynomial degree", polyDegreeLength);
      for (int d=0; d<polyDegreeLength; d++)
      {
        polyDegree(d) = fieldOrdinalH1PolynomialDegree_(fieldOrdinal,d);
      }
      return polyDegree;
    }
    
    /**
     \brief For hierarchical bases, returns the polynomial degree (which may have multiple values in higher spatial dimensions) for the specified basis ordinal as a host array.
     
     \param fieldOrdinal     [in] - ordinal of the basis function whose polynomial degree is requested.
     
     \return a std::vector<int> of length matching getPolynomialDegreeLength(), with the polynomial degree of the basis function in each dimension.
     */
    std::vector<int> getPolynomialDegreeOfFieldAsVector(int fieldOrdinal) const
    {
      INTREPID2_TEST_FOR_EXCEPTION( basisType_ != BASIS_FEM_HIERARCHICAL, std::logic_error,
                                   ">>> ERROR (Basis::getPolynomialDegreeOfFieldAsVector): this method is not supported for non-hierarchical bases.");
      auto polynomialDegreeView = getPolynomialDegreeOfField(fieldOrdinal);
      std::vector<int> polynomialDegree(polynomialDegreeView.extent_int(0));
      
      for (unsigned d=0; d<polynomialDegree.size(); d++)
      {
        polynomialDegree[d] = polynomialDegreeView(d);
      }
      return polynomialDegree;
    }
    
    /**
     \brief For hierarchical bases, returns the polynomial degree (which may have multiple values in higher spatial dimensions) for the specified basis ordinal as a host array.
     
     \param fieldOrdinal     [in] - ordinal of the basis function whose polynomial degree is requested.
     
     \return a std::vector<int> of length matching getPolynomialDegreeLength(), with the polynomial degree of the basis function in each dimension.
     */
    std::vector<int> getH1PolynomialDegreeOfFieldAsVector(int fieldOrdinal) const
    {
      INTREPID2_TEST_FOR_EXCEPTION( basisType_ != BASIS_FEM_HIERARCHICAL, std::logic_error,
                                   ">>> ERROR (Basis::getPolynomialDegreeOfFieldAsVector): this method is not supported for non-hierarchical bases.");
      auto polynomialDegreeView = getH1PolynomialDegreeOfField(fieldOrdinal);
      std::vector<int> polynomialDegree(polynomialDegreeView.extent_int(0));
      
      for (unsigned d=0; d<polynomialDegree.size(); d++)
      {
        polynomialDegree[d] = polynomialDegreeView(d);
      }
      return polynomialDegree;
    }
    
    /** \brief For hierarchical bases, returns the number of entries required to specify the polynomial degree of a basis function.
     */
    int getPolynomialDegreeLength() const
    {
      INTREPID2_TEST_FOR_EXCEPTION( basisType_ != BASIS_FEM_HIERARCHICAL, std::logic_error,
                                   ">>> ERROR (Basis::getPolynomialDegreeLength): this method is not supported for non-hierarchical bases.");
      return fieldOrdinalPolynomialDegree_.extent_int(1);
    }
    
    /** \brief  Returns basis name

        \return the name of the basis
    */
    virtual
    const char*
    getName() const {
      return "Intrepid2_Basis";
    }

    /** \brief True if orientation is required
    */
    virtual
    bool
    requireOrientation() const {
      return false;
    }

    /** \brief  Returns cardinality of the basis

        \return the number of basis functions in the basis
    */
    ordinal_type
    getCardinality() const {
      return basisCardinality_;
    }


    /** \brief  Returns the degree of the basis.

        \return max. degree of the complete polynomials that can be represented by the basis.
    */
    ordinal_type
    getDegree() const {
      return basisDegree_;
    }
    
    /** \brief  Returns the function space for the basis.
     
        \return the function space.
     */
    EFunctionSpace
    getFunctionSpace() const {
      return functionSpace_;
    }
    
    /** \brief  Returns the base cell topology for which the basis is defined. See Shards documentation
        https://trilinos.org/packages/shards for definition of base cell topology.

        \return Base cell topology
    */
    shards::CellTopology
    getBaseCellTopology() const {
      return getCellTopologyData(basisCellTopologyKey_);
    }


    /** \brief  Returns the basis type.

        \return Basis type
    */
    EBasis
    getBasisType() const {
      return basisType_;
    }


    /** \brief  Returns the type of coordinate system for which the basis is defined

        \return Type of the coordinate system (Cartesian, polar, R-Z, etc.).
    */
    ECoordinates
    getCoordinateSystem() const {
      return basisCoordinates_;
    }

    /** \brief  DoF count for specified subcell
     
     \param  subcDim           [in]  - tag field 0: dimension of the subcell associated with the DoFs
     \param  subcOrd           [in]  - tag field 1: ordinal of the subcell defined by cell topology
     
     \return the number of DoFs associated with the specified subcell.
     */
    ordinal_type
    getDofCount( const ordinal_type subcDim,
                 const ordinal_type subcOrd ) const {
      if ( subcDim >= 0   &&   subcDim < static_cast<ordinal_type>(tagToOrdinal_.extent(0)) &&
           subcOrd >= 0   &&   subcOrd < static_cast<ordinal_type>(tagToOrdinal_.extent(1)) )
      {
        int firstDofOrdinal = tagToOrdinal_(subcDim, subcOrd, 0); // will be -1 if no dofs for subcell
        if (firstDofOrdinal == -1) return static_cast<ordinal_type>(0);
        // otherwise, lookup its tag and return the dof count stored there
        return static_cast<ordinal_type>(this->getDofTag(firstDofOrdinal)[3]);
      }
      else
      {
        // subcDim and/or subcOrd out of bounds -> no dofs associated with subcell
        return static_cast<ordinal_type>(0);
      }
    }

    /** \brief  DoF tag to ordinal lookup.

        \param  subcDim           [in]  - tag field 0: dimension of the subcell associated with the DoF
        \param  subcOrd           [in]  - tag field 1: ordinal of the subcell defined by cell topology
        \param  subcDofOrd        [in]  - tag field 2: ordinal of the DoF relative to the subcell.

        \return the DoF ordinal corresponding to the specified DoF tag data.
    */
    ordinal_type
    getDofOrdinal( const ordinal_type subcDim,
                   const ordinal_type subcOrd,
                   const ordinal_type subcDofOrd ) const {
      // this should be abort and able to be called as a device function
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_EXCEPTION( subcDim < 0 || subcDim >= static_cast<ordinal_type>(tagToOrdinal_.extent(0)), std::out_of_range,
                                    ">>> ERROR (Basis::getDofOrdinal): subcDim is out of range");
      INTREPID2_TEST_FOR_EXCEPTION( subcOrd < 0 || subcOrd >= static_cast<ordinal_type>(tagToOrdinal_.extent(1)), std::out_of_range,
                                    ">>> ERROR (Basis::getDofOrdinal): subcOrd is out of range");
      INTREPID2_TEST_FOR_EXCEPTION( subcDofOrd < 0 || subcDofOrd >= static_cast<ordinal_type>(tagToOrdinal_.extent(2)), std::out_of_range,
                                    ">>> ERROR (Basis::getDofOrdinal): subcDofOrd is out of range");
#endif
      ordinal_type r_val = -1;
      if ( subcDim    < static_cast<ordinal_type>(tagToOrdinal_.extent(0)) &&
           subcOrd    < static_cast<ordinal_type>(tagToOrdinal_.extent(1)) &&
           subcDofOrd < static_cast<ordinal_type>(tagToOrdinal_.extent(2)) )
        r_val = tagToOrdinal_(subcDim, subcOrd, subcDofOrd);
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_EXCEPTION( r_val == -1, std::runtime_error,
                                    ">>> ERROR (Basis::getDofOrdinal): Invalid DoF tag is found.");
#endif
      return r_val;
    }
    
    /** \brief returns the number of tensorial extrusions relative to the cell topology returned by getBaseCellTopology().  Base class returns 0; overridden by TensorBasis.
     */
    virtual int getNumTensorialExtrusions() const
    {
      return 0;
    }

    /** \brief DoF tag to ordinal data structure */
    const OrdinalTypeArray3DHost
    getAllDofOrdinal() const {
      return tagToOrdinal_;
    }


    /** \brief  DoF ordinal to DoF tag lookup.

        \param  dofOrd            [in]  - ordinal of the DoF whose tag is being retrieved

        \return reference to a vector with dimension (4) such that \n
        \li     element [0] = tag field 0  ->  dim. of the subcell associated with the specified DoF
        \li     element [1] = tag field 1  ->  ordinal of the subcell defined by cell topology
        \li     element [2] = tag field 2  ->  ordinal of the specified DoF relative to the subcell
        \li     element [3] = tag field 3  ->  total number of DoFs associated with the subcell
    */
    const OrdinalTypeArrayStride1DHost
    getDofTag( const ordinal_type dofOrd ) const {
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_EXCEPTION( dofOrd < 0 || dofOrd >= static_cast<ordinal_type>(ordinalToTag_.extent(0)), std::out_of_range,
                                    ">>> ERROR (Basis::getDofTag): dofOrd is out of range");
#endif
      return Kokkos::subview(ordinalToTag_, dofOrd, Kokkos::ALL());
    }


    /** \brief  Retrieves all DoF tags.

        \return reference to a vector of vectors with dimensions (basisCardinality_, 4) such that \n
        \li     element [DofOrd][0] = tag field 0 for the DoF with the specified ordinal
        \li     element [DofOrd][1] = tag field 1 for the DoF with the specified ordinal
        \li     element [DofOrd][2] = tag field 2 for the DoF with the specified ordinal
        \li     element [DofOrd][3] = tag field 3 for the DoF with the specified ordinal
    */
    const OrdinalTypeArray2DHost
    getAllDofTags() const {
      return ordinalToTag_;
    }

    /** \brief returns the basis associated to a subCell.

        HGRAD case: The bases of the subCell are the restriction to the subCell
        of the bases of the parent cell.
        HCURL case: The bases of the subCell are the restriction to the subCell
        of the bases of the parent cell, projected onto the manifold tangent to the subCell
        HDIV case: The bases of the subCell are the restriction to the subCell
        of the bases of the parent cell, projected along the normal of the subCell

        This method is not supported by all bases (e.g. bases defined on a line and HVOL bases).
        \param [in] subCellDim - dimension of subCell
        \param [in] subCellOrd - position of the subCell among of the subCells having the same dimension
        \return pointer to the subCell basis of dimension subCellDim and position subCellOrd
     */
    virtual BasisPtr<DeviceType, OutputValueType, PointValueType>
      getSubCellRefBasis(const ordinal_type subCellDim, const ordinal_type subCellOrd) const {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::logic_error,
                                    ">>> ERROR (Basis::getSubCellRefBasis): this method is not supported or should be overridden accordingly by derived classes.");
    }

    /** \brief Returns the spatial dimension of the domain of the basis; this is equal to getBaseCellTopology().getDimension() + getNumTensorialExtrusions().
    
       \return The spatial dimension of the domain.
    */
    ordinal_type getDomainDimension() const
    {
      return this->getBaseCellTopology().getDimension() + this->getNumTensorialExtrusions();
    }
    
    /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.
    
       \return Pointer to the new Basis object.
    */
    virtual HostBasisPtr<OutputValueType, PointValueType>
    getHostBasis() const {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::logic_error,
                                    ">>> ERROR (Basis::getHostBasis): this method is not supported or should be overridden accordingly by derived classes.");
    }

  }; // class Basis

  //--------------------------------------------------------------------------------------------//
  //                                                                                            //
  //                            Helper functions of the Basis class                             //
  //                                                                                            //
  //--------------------------------------------------------------------------------------------//

  //
  // functions for orders, cardinality, etc.
  //


  /** \brief  Returns the rank of fields in a function space of the specified type.

      Field rank is defined as the number of indices needed to specify function value and
      equals 0, 1,or 2 for scalars, vectors and tensors, respectively. The scalar field
      spaces in Intrepid are FUNCTION_SPACE_HGRAD and FUNCTION_SPACE_HVOL. The vector field
      spaces are FUNCTION_SPACE_HCURL, FUNCTION_SPACE_HDIV and FUNCTION_SPACE_VECTOR_HGRAD.
      FUNCTION_SPACE_TENSOR_HGRAD contains rank-2 tensors.

      \param  spaceType         [in]     -  function space type
      \return rank of the fields in the specified function space
  */
  KOKKOS_INLINE_FUNCTION
  ordinal_type getFieldRank(const EFunctionSpace spaceType);

  /** \brief  Returns rank of an operator.

      When an operator acts on a field of a certain rank, the result can be a field with the
      same or a different rank. Operator rank is defined the difference between the ranks of
      the output field and the input field:
      \verbatim
      Rank(OPERATOR) = Rank(OPERATOR(FIELD)) - Rank(FIELD)
      \endverbatim
      Therefore, operator rank allows us to figure out the rank of the result:
      \verbatim
      Rank(OPERATOR(FIELD)) = Rank(FIELD) + Rank(OPERATOR)
      \endverbatim
      and provides means to size properly arrays for output results. The following table
      summarizes operator ranks (~ denotes undefined, below slash means 3D).
      By default, in 1D any operator other than VALUE has rank 1, i.e., GRAD, CURL and DIV
      reduce to d/dx and Dk are the higher-order derivatives d^k/dx^k. Only scalar functions
      are allowed in 1D.

      \verbatim
      |========|======|============================|=========|==========|==========|==========|
      | field  | rank |  FUNCTION_SPACE_[type]     |  VALUE  | GRAD, Dk |   CURL   |    DIV   |
      |--------|------|----------------------------|---------|----------|----------|----------|
      | scalar |   0  |  HGRAD, HVOL               |    0    |     1    | 3-dim/~  |     ~    |
      | vector |   1  |  HCURL, HDIV, VECTOR_HGRAD |    0    |     1    | dim - 3  |    -1    |
      | tensor |   2  |  TENSOR_HGRAD              |    0    |     1    | dim - 3  |    -1    |
      |--------|------|----------------------------|---------|----------|----------|----------|
      |   1D   |   0  |  HGRAD, HVOL only          |    0    |     1    |     1    |     1    |
      |=======================================================================================|
      \endverbatim

      \param  spaceType        [in]    - function space type
      \param  operatorType     [in]    - the operator acting on the specified function space
      \param  spaceDim         [in]    - spatial dimension
      \return rank of the operator as defined in the table
  */
  KOKKOS_INLINE_FUNCTION
  ordinal_type getOperatorRank(const EFunctionSpace spaceType,
                               const EOperator      operatorType,
                               const ordinal_type   spaceDim);

  /** \brief  Returns order of an operator.

      \param  operatorType       [in]    - type of the operator whose order we want to know
      \return result ranges from 0 (for OPERATOR_VALUE) to 10 (OPERATOR_D10)
  */
  KOKKOS_INLINE_FUNCTION
  ordinal_type getOperatorOrder(const EOperator operatorType);

  template<EOperator operatorType>
  KOKKOS_INLINE_FUNCTION
  constexpr ordinal_type getOperatorOrder();

  /** \brief  Returns the ordinal of a partial derivative of order k based on the multiplicities of
      the partials dx, dy, and dz.

      By default, any implementation of Intrepid2::Basis method returns partials of order k
      (specified by OPERATOR_Dk) as a multiset ordered by the lexicographical order of the
      partial derivatives multiplicities. For example, the 10 derivatives of order 3 in 3D
      are enumerated as:
      \verbatim
      D3={(3,0,0),(2,1,0),(2,0,1),(1,2,0),(1,1,1),(1,0,2),(0,3,0),(0,2,1),(0,1,2),(0,0,3)}
      \endverbatim
      The enumeration formula for this lexicographical order is
      <table>
      <tr> <td>\f$i()            = 0\f$                            </td> <td>in 1D (only 1 derivative)</td> </tr>
      <tr> <td>\f$i(yMult)      =yMult\f$                         </td> <td>in 2D</td> </tr>
      <tr> <td>\f$i(yMult,zMult)=zMult+\sum_{r = 0}^{k-xMult} r\f$</td> <td>in 3D</td> </tr>
      </table>
      where the order k of Dk is defined by xMult + yMult + zMult. However, xMult is not really needed needed.
      Enumeration goes from 0 to ( (k+spaceDim-1) choose (spaceDim-1) )-1.

      \param  yMult            [in]    - multiplicity of dy (default = -1)
      \param  zMult            [in]    - multiplicity of dz (default = -1)
      \return the ordinal of partial derivative of order k as function of dx, dy, dz multiplicities
  */
  template<ordinal_type spaceDim>
  KOKKOS_INLINE_FUNCTION
  ordinal_type getDkEnumeration(const ordinal_type xMult,
                                const ordinal_type yMult = -1,
                                const ordinal_type zMult = -1);  


  /** \brief  Returns the index of the term x^p y^q z^r of a polynomial of degree n  (p+q+r <= n).
      In 2D, the terms of a polynomial of degree 2 are ordered as  1,  x,  y, x^2,  xy, y^2.
      So if p=q=1, the term x^p y^q has index 4 (counting from 0), while p=2, q=0 has index 3.
      Enumeration goes from 0 to ( (n+spaceDim) choose (spaceDim) )-1.

      \param  p            [in]    - exponent of x in the polynomial term x^p y^q z^r
      \param  q            [in]    - exponent of y in the polynomial term x^p y^q z^r (default = 0)
      \param  r            [in]    - exponent of z in the polynomial term x^p y^q z^r (default = 0)
      \return the index of the term x^p y^q z^r of a polynomial of degree k  (p+q+r <= n)
  */
  template<ordinal_type spaceDim>
  KOKKOS_INLINE_FUNCTION
  ordinal_type getPnEnumeration(const ordinal_type p,
                                const ordinal_type q = 0,
                                const ordinal_type r = 0);



  /** \brief function for computing the Jacobi recurrence coefficients so that
    \param an    [out] - the a weight for recurrence
    \param bn    [out] - the b weight for recurrence
    \param cn    [out] - the c weight for recurrence
    \param alpha [in] - the first Jacobi weight
    \param beta  [in] - the second Jacobi weight
    \param n     [n]  - the polynomial degree

    The recurrence is
    \f[
    P^{\alpha,\beta}_{n+1} = \left( a_n + b_n x\right) P^{\alpha,\beta}_n - c_n P^{\alpha,\beta}_{n-1}
    \f],
    where
    \f[
    P^{\alpha,\beta}_0 = 1
    \f]
*/

template<typename value_type>
KOKKOS_INLINE_FUNCTION
void getJacobyRecurrenceCoeffs (
          value_type  &an,
          value_type  &bn,
          value_type  &cn,
    const ordinal_type alpha,
    const ordinal_type beta ,
    const ordinal_type n);





  // /** \brief  Returns multiplicities of dx, dy, and dz based on the enumeration of the partial
  //     derivative, its order and the space dimension. Inverse of the getDkEnumeration() method.

  //     \param  partialMult      [out]    - array with the multiplicities f dx, dy and dz
  //     \param  derivativeEnum   [in]     - enumeration of the partial derivative
  //     \param  operatorType     [in]     - k-th partial derivative Dk
  //     \param  spaceDim         [in]     - space dimension
  // */
  // template<typename OrdinalArrayType>
  // KOKKOS_INLINE_FUNCTION
  // void getDkMultiplicities(OrdinalArrayType   partialMult,
  //                          const ordinal_type derivativeEnum,
  //                          const EOperator    operatorType,
  //                          const ordinal_type spaceDim);

  /** \brief  Returns cardinality of Dk, i.e., the number of all derivatives of order k.

      The set of all partial derivatives of order k is isomorphic to the set of all multisets
      of cardinality k with elements taken from the sets {x}, {x,y}, and {x,y,z} in 1D, 2D,
      and 3D respectively. For example, the partial derivative
      \f$\displaystyle D\{1,2,1\}f = \frac{d^4 f}{dx dy^2 dz}\f$  maps to the multiset
      \f$\{x, y, z\}\f$ with multiplicities \f$\{1,2,1\}\f$. The number of all such multisets
      is given by the binomial coefficient
      \f[       \begin{pmatrix} spaceDim + k - 1 \\ spaceDim - 1 \end{pmatrix}              \f]
      Therefore:
      \li     in 1D: cardinality = 1\n
      \li     in 2D: cardinality = k + 1\n
      \li     in 3D: cardinality = (k + 1)*(k + 2)/2

      \param  operatorType     [in]     - k-th derivative operator Dk
      \param  spaceDim         [in]     - space dimension
      \return the number of all partial derivatives of order k
  */
  KOKKOS_INLINE_FUNCTION
  ordinal_type getDkCardinality(const EOperator    operatorType,
                                const ordinal_type spaceDim);

  template<EOperator operatorType, ordinal_type spaceDim>
  KOKKOS_INLINE_FUNCTION
  constexpr ordinal_type getDkCardinality();



  /** \brief  Returns cardinality of Polynomials of order n (P^n).

       \li     in 1D: cardinality = n+1
       \li     in 2D: cardinality = (n+1)*(n+2)/2
       \li     in 3D: cardinality = (n+1)*(n+2)*(n+3)/6

       \param  n     [in]     - polynomial order
       \return dimension of polynomial space
   */
   template<ordinal_type spaceDim>
   KOKKOS_INLINE_FUNCTION
   ordinal_type getPnCardinality (ordinal_type n);

   template<ordinal_type spaceDim, ordinal_type n>
   KOKKOS_INLINE_FUNCTION
   constexpr ordinal_type getPnCardinality ();



  //
  // Argument check
  //


  /** \brief  Runtime check of the arguments for the getValues method in an HGRAD-conforming
      FEM basis. Verifies that ranks and dimensions of <var>ViewType</var> input and output
      arrays are consistent with the specified <var>operatorType</var>.

      \param  outputValues     [in]  - array of variable rank for the output basis values
      \param  inputPoints      [in]  - rank-2 array with dimensions (P,D) containing the points
      \param  operatorType     [in]  - operator applied to basis functions
      \param  cellTopo         [in]  - base cell topology on which the basis is defined
      \param  basisCard        [in]  - cardinality of the basis
  */
  template<typename outputValueViewType,
           typename inputPointViewType>
  void getValues_HGRAD_Args( const outputValueViewType   outputValues,
                             const inputPointViewType    inputPoints,
                             const EOperator             operatorType,
                             const shards::CellTopology  cellTopo,
                             const ordinal_type          basisCard );

  /** \brief  Runtime check of the arguments for the getValues method in an HCURL-conforming
      FEM basis. Verifies that ranks and dimensions of <var>ViewType</var> input and output
      arrays are consistent with the specified <var>operatorType</var>.

      \param  outputValues     [in]  - array of variable rank for the output basis values
      \param  inputPoints      [in]  - rank-2 array with dimensions (P,D) containing the points
      \param  operatorType     [in]  - operator applied to basis functions
      \param  cellTopo         [in]  - base cell topology on which the basis is defined
      \param  basisCard        [in]  - cardinality of the basis
  */
  template<typename outputValueViewType,
           typename inputPointViewType>
  void getValues_HCURL_Args( const outputValueViewType   outputValues,
                             const inputPointViewType    inputPoints,
                             const EOperator             operatorType,
                             const shards::CellTopology  cellTopo,
                             const ordinal_type          basisCard );

  /** \brief  Runtime check of the arguments for the getValues method in an HDIV-conforming
      FEM basis. Verifies that ranks and dimensions of <var>ViewType</var> input and output
      arrays are consistent with the specified <var>operatorType</var>.

      \param  outputValues     [in]  - array of variable rank for the output basis values
      \param  inputPoints      [in]  - rank-2 array with dimensions (P,D) containing the points
      \param  operatorType     [in]  - operator applied to basis functions
      \param  cellTopo         [in]  - base cell topology on which the basis is defined
      \param  basisCard        [in]  - cardinality of the basis
  */
  template<typename outputValueViewType,
           typename inputPointViewType>
  void getValues_HDIV_Args( const  outputValueViewType  outputValues,
                            const inputPointViewType    inputPoints,
                            const EOperator             operatorType,
                            const shards::CellTopology  cellTopo,
                            const ordinal_type          basisCard );

  /** \brief  Runtime check of the arguments for the getValues method in an HVOL-conforming
      FEM basis. Verifies that ranks and dimensions of <var>ViewType</var> input and output
      arrays are consistent with the specified <var>operatorType</var>.

      \param  outputValues     [in]  - array of variable rank for the output basis values
      \param  inputPoints      [in]  - rank-2 array with dimensions (P,D) containing the points
      \param  operatorType     [in]  - operator applied to basis functions
      \param  cellTopo         [in]  - base cell topology on which the basis is defined
      \param  basisCard        [in]  - cardinality of the basis
  */
  template<typename outputValueViewType,
           typename inputPointViewType>
  void getValues_HVOL_Args( const outputValueViewType   outputValues,
                             const inputPointViewType    inputPoints,
                             const EOperator             operatorType,
                             const shards::CellTopology  cellTopo,
                             const ordinal_type          basisCard );

}// namespace Intrepid2

#include <Intrepid2_BasisDef.hpp>

//--------------------------------------------------------------------------------------------//
//                                                                                            //
//                            D O C U M E N T A T I O N   P A G E S                           //
//                                                                                            //
//--------------------------------------------------------------------------------------------//
/**
 \page basis_page                       Intrepid2 basis class



 \section basis_dof_tag_ord_sec         Degree of freedom ordinals and tags

 Regardless of the basis type, i.e., FEM or FVD, each DoF is assigned an ordinal number which specifies
 its numerical position in the DoF set, and a 4-field DoF tag whose first 3 fields establish association
 between the DoF and a subcell of particular dimension. The last field in the DoF tag is for convenience
 and stores the total number of DoFs associated with the specified subcell. In summary, the tag contains
 the following information about a DoF with a given ordinal:

 \li field 0: dimension of the subcell associated with the specified DoF ordinal;
 \li field 1: ordinal of the subcell relative to its parent cell;
 \li field 2: ordinal of the DoF relative to the subcell;
 \li field 3: cardinality of the DoF set associated with this subcell.

 DoF definition, DoF ordinals and DoF tags are basis-dependent and are documemented in the concrete
 basis implementation. A typical entry in a DoF tag table has the following format:
 \verbatim
 |-------------------------------------------------------------------------------------------------|
 |         |           degree-of-freedom-tag table                    |                            |
 |   DoF   |----------------------------------------------------------|      DoF definition        |
 | ordinal |  subc dim    | subc ordinal | subc DoF ord |subc num DoF |                            |
 |-------------------------------------------------------------------------------------------------|
 |---------|--------------|--------------|--------------|-------------|----------------------------|
 |    k    |       1      |       2      |       1      |      3      |   L_k(u) = (definition)    |
 |---------|--------------|--------------|--------------|-------------|----------------------------|
 |-------------------------------------------------------------------------------------------------|
 \endverbatim

 The tag in this example establishes an association between the DoF with ordinal <var>k</var> and the 3rd
 edge of the parent cell on which the basis is defined. Furthermore, the tag specifies that relative
 to that edge, this DoF has ordinal 1, i.e., it is the second DoF on that edge. The last field in the
 tag indicates that there are a total of 3 DoFs associated with that subcell.



 \section basis_md_array_sec            MD array template arguments for basis methods

 FEM and FVD basis evaluation methods in Intrepid2 use <a href="https://github.com/kokkos/">Kokkos</a> 
 dynamic rank views as the default MD arrays to pass the evaluation points (and cell vertices for FVD 
 evaluation) and to return the basis values.
 The ranks and the dimensions of the MD array arguments for these methods are as follows.



 \subsection basis_md_array_out_sec     Rank and dimensions of the output MD array

 Rank and dimensions of the output array depend on the field rank of the basis functions, which can be 0
 (scalar fields), 1 (vector fields), or 2 (tensor fields), the space  dimension, and the <var>operatorType</var>.
 The following table summarizes all admissible combinations:
 \verbatim
 |-------------------------------------------------------------------------------------------------|
 |            Rank and multi-dimensions of the output MD array in getValues methods                |
 |--------------------|-------------------------|-------------------------|------------------------|
 |operator/field rank |  rank 0                 |  rank 1 2D/3D           |  rank 2 2D/3D          |
 |--------------------|-------------------------|-------------------------|------------------------|
 |       VALUE        |  (F,P)                  | (F,P,D)                 | (F,P,D,D)              |
 |--------------------|-------------------------|-------------------------|------------------------|
 |     GRAD, D1       |  (F,P,D)                | (F,P,D,D)               | (F,P,D,D,D)            |
 |--------------------|-------------------------|-------------------------|------------------------|
 |        CURL        |  (F,P,D) (undef. in 3D) | (F,P)/(F,P,D)           | (F,P,D)/(F,P,D,D)      |
 |--------------------|-------------------------|-------------------------|------------------------|
 |        DIV         |  (F,P,D) (only in 1D)   | (F,P)                   | (F,P,D)                |
 |--------------------|-------------------------|-------------------------|------------------------|
 |    D1,D2,..,D10    |  (F,P,K)                | (F,P,D,K)               | (F,P,D,D,K)            |
 |-------------------------------------------------------------------------------------------------|
 \endverbatim

 \remarks
 \li The totality of all derivatives whose order equals k (OPERATOR_Dk in Intrepid) forms a multiset;
 see http://mathworld.wolfram.com/Multiset.html In Intrepid this multiset is enumerated using the
 lexicographical order of the partial derivatives; see getDkEnumeration() for details.

 \li The last dimension of the output array for D1,...,D10 is the cardinality of the Dk multiset
 (computed by DkCardinality). The array is filled with zeroes whenever the order of the derivative
 Dk exceeed the polynomial degree.



 \subsection basis_md_array_in_sec     Rank and dimensions of the input MD arrays

 The FEM evaluation method has one MD array input argument which is used to pass the coordinates of
 P evaluation points in the reference cell for which the concrete basis is defined. The FVD method
 has two MD array input arguments. The first one passes the coordinates of P evaluation points in the
 physical cell for which the concrete basis is defined. The second MD array passes the vertices of the
 physical cell. Ranks and dimensions of these arrays are summarized in the following table:
 \verbatim
 |-------------------------------------------------------------------------------------------------|
 |             Rank and multi-dimensions of the input MD arrays in getValues methods               |
 |--------------------|------------------------|---------------------------------------------------|
 | MD array           | rank | multi-dimension |  Description                                      |
 |--------------------|------------------------|---------------------------------------------------|
 | evaluation points  |   2  |  (P,D)          |  Coordinates of P points in D-dimensions          |
 |--------------------|------------------------|---------------------------------------------------|
 | cell vertices      |   2  |  (V,D)          |  Coordinates of V vertices of D-dimensional cell  |
 |-------------------------------------------------------------------------------------------------|
 \endverbatim
 */
#endif
