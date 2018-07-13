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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
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

#include "Intrepid2_CellTopologyTags.hpp"
#include "Shards_CellTopology.hpp"

namespace Intrepid2 {

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

      \remark To limit memory use by factory-type objects (basis factories will be included in future
              releases of Intrepid2), tag data is not initialized by basis ctors,
              instead, whenever a function that requires tag data is invoked for a first time, it calls
              initializeTags() to fill <var>ordinalToTag_</var> and <var>tagToOrdinal_</var>. Because
              tag data is basis specific, every concrete basis class requires its own implementation
              of initializeTags().

      \todo  restore test for inclusion of reference points in their resective reference cells in
             getValues_HGRAD_Args, getValues_CURL_Args, getValues_DIV_Args
  */
  template<typename ExecSpaceType = void,
           typename outputValueType = double,
           typename pointValueType = double>
  class Basis {
  public:

    /**  \brief View type for ordinal
    */
    typedef Kokkos::View<ordinal_type,ExecSpaceType> ordinal_view_type;

    /**  \brief View for basis type
    */
    typedef Kokkos::View<EBasis,ExecSpaceType> ebasis_view_type;

    /**  \brief View for coordinate system type
    */
    typedef Kokkos::View<ECoordinates,ExecSpaceType> ecoordinates_view_type;

    // ** tag interface
    //  - tag interface is not decorated with Kokkos inline so it should be allocated on hostspace

    /**  \brief View type for 1d host array
    */
    typedef Kokkos::View<ordinal_type*  ,typename ExecSpaceType::array_layout,Kokkos::HostSpace> ordinal_type_array_1d_host;

    /**  \brief View type for 2d host array
    */
    typedef Kokkos::View<ordinal_type** ,typename ExecSpaceType::array_layout,Kokkos::HostSpace> ordinal_type_array_2d_host;

    /**  \brief View type for 3d host array
    */
    typedef Kokkos::View<ordinal_type***,typename ExecSpaceType::array_layout,Kokkos::HostSpace> ordinal_type_array_3d_host;

    /**  \brief View type for 1d host array
    */
    typedef Kokkos::View<ordinal_type*  , Kokkos::LayoutStride, Kokkos::HostSpace> ordinal_type_array_stride_1d_host;

    /**  \brief View type for 1d device array
    */
    typedef Kokkos::View<ordinal_type*  ,ExecSpaceType> ordinal_type_array_1d;

    /**  \brief View type for 2d device array
    */
    typedef Kokkos::View<ordinal_type** ,ExecSpaceType> ordinal_type_array_2d;

    /**  \brief View type for 3d device array
    */
    typedef Kokkos::View<ordinal_type***,ExecSpaceType> ordinal_type_array_3d;

    /**  \brief View type for 1d device array 
    */
    typedef Kokkos::View<ordinal_type*  , Kokkos::LayoutStride, ExecSpaceType> ordinal_type_array_stride_1d;

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

    /** \brief  Base topology of the cells for which the basis is defined. See
         the <a href="https://trilinos.org/packages/shards/">Shards</a> package
         for definition of base cell topology.
    */
    shards::CellTopology basisCellTopology_;

    /** \brief  Type of the basis
     */
    EBasis basisType_;

    /** \brief  The coordinate system for which the basis is defined
     */
    ECoordinates basisCoordinates_;

    /** \brief  "true" if <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> have been initialized
     */
    //Kokkos::View<bool,ExecSpaceType> basisTagsAreSet_;

    /** \brief  DoF ordinal to tag lookup table.

        Rank-2 array with dimensions (basisCardinality_, 4) containing the DoF tags. This array
        is left empty at instantiation and filled by initializeTags() only when tag data is
        requested.

        \li     ordinalToTag_[DofOrd][0] = dim. of the subcell associated with the specified DoF
        \li     ordinalToTag_[DofOrd][1] = ordinal of the subcell defined in the cell topology
        \li     ordinalToTag_[DodOrd][2] = ordinal of the specified DoF relative to the subcell
        \li     ordinalToTag_[DofOrd][3] = total number of DoFs associated with the subcell
    */
    ordinal_type_array_2d_host ordinalToTag_;

    /** \brief  DoF tag to ordinal lookup table.

        Rank-3 array with dimensions (maxScDim + 1, maxScOrd + 1, maxDfOrd + 1), i.e., the
        columnwise maximums of the 1st three columns in the DoF tag table for the basis plus 1.
        For every triple (subscDim, subcOrd, subcDofOrd) that is valid DoF tag data this array
        stores the corresponding DoF ordinal. If the triple does not correspond to tag data,
        the array stores -1. This array is left empty at instantiation and filled by
        initializeTags() only when tag data is requested.

        \li     tagToOrdinal_[subcDim][subcOrd][subcDofOrd] = Degree-of-freedom ordinal
    */
    ordinal_type_array_3d_host tagToOrdinal_;

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
    Kokkos::DynRankView<scalarType,ExecSpaceType> dofCoords_;

    // dof coeffs
    /** \brief Coefficients for computing degrees of freedom for Lagrangian basis
        If P is an element of the space spanned by the basis,
        \alpha_i := P(dofCoords_(i)) \cdot dofCoeffs_(i) are the nodal coefficients associated to basis functions i.

        Rank-1 array for scalar basis with dimension (cardinality)
        Rank-2 array for vector basis with dimensions (cardinality, cell dimension)
     */
    Kokkos::DynRankView<scalarType,ExecSpaceType> dofCoeffs_;

  public:

    Basis() = default;
    virtual~Basis() = default;

    // receives input arguments
    /** \brief Dummy array to receive input arguments
    */
    outputValueType getDummyOutputValue() { return outputValueType(); }

    /** \brief Dummy array to receive input arguments
    */
    pointValueType getDummyPointValue() { return pointValueType(); }

    /** \brief View type for basis value output
    */
    typedef Kokkos::DynRankView<outputValueType,Kokkos::LayoutStride,ExecSpaceType> outputViewType;

    /** \brief View type for input points
    */
    typedef Kokkos::DynRankView<pointValueType,Kokkos::LayoutStride,ExecSpaceType>  pointViewType;

    /** \brief View type for scalars 
    */
    typedef Kokkos::DynRankView<scalarType,Kokkos::LayoutStride,ExecSpaceType>      scalarViewType;

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
    void
    getValues(       outputViewType outputValues,
               const pointViewType  inputPoints,
               const EOperator operatorType = OPERATOR_VALUE ) const {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::logic_error,
                                    ">>> ERROR (Basis::getValues): this method (FEM) is not supported or should be over-riden accordingly by derived classes.");
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
    getValues(        outputViewType outputValues,
                const pointViewType  inputPoints,
                const pointViewType  cellVertices,
                const EOperator operatorType = OPERATOR_VALUE ) const {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::logic_error,
                                    ">>> ERROR (Basis::getValues): this method (FVM) is not supported or should be over-riden accordingly by derived classes.");
    }


    /** \brief  Returns spatial locations (coordinates) of degrees of freedom on the reference cell
     */

    virtual
    void
    getDofCoords( scalarViewType dofCoords ) const {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::logic_error,
                                    ">>> ERROR (Basis::getDofCoords): this method is not supported or should be over-riden accordingly by derived classes.");
    }

    /** \brief Coefficients for computing degrees of freedom for Lagrangian basis
        If P is an element of the space spanned by the basis,
        \alpha_i := P(dofCoords(i)) \cdot dofCoeffs(i) are the nodal coefficients associated to basis function i.

        Rank-1 array for scalar basis with dimension (cardinality)
        Rank-2 array for vector basis with dimensions (cardinality, cell dimension)
     */

    virtual
    void
    getDofCoeffs( scalarViewType dofCoeffs ) const {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::logic_error,
                                    ">>> ERROR (Basis::getDofCoeffs): this method is not supported or should be over-riden accordingly by derived classes.");
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


    /** \brief  Returns the base cell topology for which the basis is defined. See Shards documentation
        https://trilinos.org/packages/shards for definition of base cell topology.

        \return Base cell topology
    */
    shards::CellTopology
    getBaseCellTopology() const {
      return basisCellTopology_;
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

    /** \brief DoF tag to ordinal data structure */
    const ordinal_type_array_3d_host
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
    const ordinal_type_array_stride_1d_host
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
    const ordinal_type_array_2d_host
    getAllDofTags() const {
      return ordinalToTag_;
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
