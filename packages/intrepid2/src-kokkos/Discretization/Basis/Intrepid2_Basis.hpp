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

/** \file   Intrepid_Basis.hpp
    \brief  Header file for the abstract base class Intrepid2::Basis.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
 */

#ifndef __INTREPID2_BASIS_HPP__
#define __INTREPID2_BASIS_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

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
              releases of Intrepid), tag data is not initialized by basis ctors,
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
    typedef Kokkos::View<ordinal_type,ExecSpaceType> ordinal_view_type;
    typedef Kokkos::View<EBasis,ExecSpaceType> ebasis_view_type;
    typedef Kokkos::View<ECoordinates,ExecSpaceType> ecoordiantes_view_type;

    typedef Kokkos::View<ordinal_type*  ,ExecSpaceType> ordinal_type_array_1d;
    typedef Kokkos::View<ordinal_type** ,ExecSpaceType> ordinal_type_array_2d;
    typedef Kokkos::View<ordinal_type***,ExecSpaceType> ordinal_type_array_3d;

  protected:

    /** \brief  Cardinality of the basis, i.e., the number of basis functions/degrees-of-freedom
     */
    ordinal_type basisCardinality_;

    /** \brief  Degree of the largest complete polynomial space that can be represented by the basis
     */
    ordinal_type basisDegree_;

    /** \brief  Base topology of the cells for which the basis is defined. See  the Shards package
        http://trilinos.sandia.gov/packages/shards for definition of base cell topology.
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
    ordinal_type_array_2d ordinalToTag_;

    /** \brief  DoF tag to ordinal lookup table.

        Rank-3 array with dimensions (maxScDim + 1, maxScOrd + 1, maxDfOrd + 1), i.e., the
        columnwise maximums of the 1st three columns in the DoF tag table for the basis plus 1.
        For every triple (subscDim, subcOrd, subcDofOrd) that is valid DoF tag data this array
        stores the corresponding DoF ordinal. If the triple does not correspond to tag data,
        the array stores -1. This array is left empty at instantiation and filled by
        initializeTags() only when tag data is requested.

        \li     tagToOrdinal_[subcDim][subcOrd][subcDofOrd] = Degree-of-freedom ordinal
    */
    ordinal_type_array_3d tagToOrdinal_;

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
    void setOrdinalTagData( /**/  ordinal_type_array_3d &tagToOrdinal,
                            /**/  ordinal_type_array_2d &ordinalToTag,
                            const ordinal_type_array_1d  tags,
                            const ordinal_type           basisCard,
                            const ordinal_type           tagSize,
                            const ordinal_type           posScDim,
                            const ordinal_type           posScOrd,
                            const ordinal_type           posDfOrd ) {
      // Create ordinalToTag
      ordinalToTag = ordinal_type_array_2d("ordinalToTag", basisCard, 4);

      // Initialize with -1
      Kokkos::deep_copy( ordinalToTag, -1 );

      // Copy tags
      for (auto i=0;i<basisCard;++i)
        for (auto j=0;j<tagSize;++j)
          ordinalToTag(i, j) = tags(i*tagSize + j);

      // Find out dimension of tagToOrdinal
      auto maxScDim = 0;  // first dimension of tagToOrdinal
      for (auto i=0;i<basisCard;++i)
        if (maxScDim < tags(i*tagSize + posScDim))
          maxScDim = tags(i*tagSize + posScDim);
      ++maxScDim;

      auto maxScOrd = 0; // second dimension of tagToOrdinal
      for (auto i=0;i<basisCard;++i)
        if (maxScOrd < tags(i*tagSize + posScOrd))
          maxScOrd = tags(i*tagSize + posScOrd);
      ++maxScOrd;

      auto maxDfOrd = 0;  // third dimension of tagToOrdinal
      for (auto i=0;i<basisCard;++i)
        if (maxDfOrd < tags(i*tagSize + posDfOrd))
          maxDfOrd = tags[i*tagSize + posDfOrd];
      ++maxDfOrd;

      // Create tagToOrdinal
      tagToOrdinal = ordinal_type_array_3d("tagToOrdinal", maxScDim, maxScOrd, maxDfOrd);

      // Initialize with -1
      Kokkos::deep_copy( tagToOrdinal, -1 );

      // Overwrite elements of the array corresponding to tags with local DoF Id's, leave all other = -1
      for (auto i=0;i<basisCard;++i)
        tagToOrdinal(tags(i*tagSize), tags(i*tagSize+1), tags(i*tagSize+2)) = i;
    }



  public:

    Basis() = default;
    virtual~Basis() = default;

    typedef Kokkos::DynRankView<outputValueType,ExecSpaceType> outputViewType;
    typedef Kokkos::DynRankView<pointValueType,ExecSpaceType>  pointViewType;

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
    getValues( /**/  outputViewType outputValues,
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
    getValues(  /**/  outputViewType outputValues,
                const pointViewType  inputPoints,
                const pointViewType  cellVertices,
                const EOperator operatorType = OPERATOR_VALUE ) const {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::logic_error,
                                    ">>> ERROR (Basis::getValues): this method (FVM) is not supported or should be over-riden accordingly by derived classes.");
    }


    /** \brief  Returns spatial locations (coordinates) of degrees of freedom on a
        <strong>reference Quadrilateral</strong>.

        \param  DofCoords      [out] - array with the coordinates of degrees of freedom,
        dimensioned (F,D)
    */
    virtual
    void
    getDofCoords( pointViewType dofCoords ) const {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::logic_error,
                                    ">>> ERROR (Basis::getDofCoords): this method is not supported or should be over-riden accordingly by derived classes.");
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
        http://trilinos.sandia.gov/packages/shards for definition of base cell topology.

        \return Base cell topology
    */
    const shards::CellTopology
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
      INTREPID2_TEST_FOR_EXCEPTION( subcDim < 0 || subcDim >= static_cast<ordinal_type>(tagToOrdinal_.dimension(0)), std::out_of_range,
                                    ">>> ERROR (Basis::getDofOrdinal): subcDim is out of range");
      INTREPID2_TEST_FOR_EXCEPTION( subcOrd < 0 || subcOrd >= static_cast<ordinal_type>(tagToOrdinal_.dimension(1)), std::out_of_range,
                                    ">>> ERROR (Basis::getDofOrdinal): subcOrd is out of range");
      INTREPID2_TEST_FOR_EXCEPTION( subcDofOrd < 0 || subcDofOrd >= static_cast<ordinal_type>(tagToOrdinal_.dimension(2)), std::out_of_range,
                                    ">>> ERROR (Basis::getDofOrdinal): subcDofOrd is out of range");
#endif
      const auto r_val = tagToOrdinal_(subcDim, subcOrd, subcDofOrd);
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_EXCEPTION( r_val == -1, std::runtime_error,
                                    ">>> ERROR (Basis::getDofOrdinal): Invalid DoF tag is found.");
#endif
      return r_val;
    }

    /** \brief DoF tag to ordinal data structure */
    const ordinal_type_array_3d
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
    const ordinal_type_array_1d
    getDofTag( const ordinal_type dofOrd ) const {
#ifdef HAVE_INTREPID2_DEBUG
      INTREPID2_TEST_FOR_EXCEPTION( dofOrd < 0 || dofOrd >= static_cast<ordinal_type>(ordinalToTag_.dimension(0)), std::out_of_range,
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
    const ordinal_type_array_2d
    getAllDofTags() const {
      return ordinalToTag_;
    }

  }; // class Basis


  //--------------------------------------------------------------------------------------------//
  //                                                                                            //
  //                            Helper functions of the Basis class                             //
  //                                                                                            //
  //--------------------------------------------------------------------------------------------//

  ///
  /// functions for orders, cardinality, etc.
  ///


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
      <tr> <td>\f$i(xMult)            = 0\f$                            </td> <td>in 1D (only 1 derivative)</td> </tr>
      <tr> <td>\f$i(xMult,yMult)      =yMult\f$                         </td> <td>in 2D</td> </tr>
      <tr> <td>\f$i(xMult,yMult,zMult)=zMult+\sum_{r = 0}^{k-xMult} r\f$</td> <td>in 3D</td> </tr>
      </table>
      where the order k of Dk is implicitly defined by xMult + yMult + zMult. Space dimension is
      implicitly defined by the default values of the multiplicities of y and z derivatives.

      \param  xMult            [in]    - multiplicity of dx
      \param  yMult            [in]    - multiplicity of dy (default = -1)
      \param  zMult            [in]    - multiplicity of dz (default = -1)
      \return the ordinal of partial derivative of order k as function of dx, dy, dz multiplicities
  */
  KOKKOS_INLINE_FUNCTION
  ordinal_type getDkEnumeration(const ordinal_type xMult,
                                const ordinal_type yMult = -1,
                                const ordinal_type zMult = -1);

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



  ///
  /// Argument check
  ///


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

}// namespace Intrepid2

#include <Intrepid2_BasisDef.hpp>

//--------------------------------------------------------------------------------------------//
//                                                                                            //
//                            D O C U M E N T A T I O N   P A G E S                           //
//                                                                                            //
//--------------------------------------------------------------------------------------------//
/**
 \page basis_page                       Intrepid basis class



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

 FEM and FVD basis evaluation methods use generic MD arrays (see \ref md_array_page for details) to
 pass the  evaluation points (and cell vertices for FVD evaluation) and to return the basis values.
 The ranks  and the dimensions of the MD array arguments for these methods are as follows.



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
