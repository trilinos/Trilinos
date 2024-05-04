// @HEADER
// ************************************************************************
//
//                           Intrepid Package
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_Basis.hpp
    \brief  Header file for the abstract base class Intrepid::Basis.
    \author Created by P. Bochev and D. Ridzal.
 */

#ifndef INTREPID_BASIS_HPP
#define INTREPID_BASIS_HPP
#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Types.hpp"
#include "Intrepid_Utils.hpp"
#include "Shards_CellTopology.hpp"

namespace Intrepid {
    
/** \class  Intrepid::Basis
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
  
  \remark   To limit memory use by factory-type objects (basis factories will be included in future
            releases of Intrepid), tag data is not initialized by basis ctors,
            instead, whenever a function that requires tag data is invoked for a first time, it calls
            initializeTags() to fill <var>ordinalToTag_</var> and <var>tagToOrdinal_</var>. Because  
            tag data is basis specific, every concrete basis class requires its own implementation
            of initializeTags(). 
  
  \todo     restore test for inclusion of reference points in their resective reference cells in 
            getValues_HGRAD_Args, getValues_CURL_Args, getValues_DIV_Args
  
 */
template<class Scalar, class ArrayScalar>
class Basis {
private:

  /** \brief  Initializes <var>tagToOrdinal_</var> and <var>ordinalToTag_</var> lookup arrays.
   */
  virtual void initializeTags() = 0;
  
protected: 
  
  /** \brief  Cardinality of the basis, i.e., the number of basis functions/degrees-of-freedom
   */
  int basisCardinality_;

  /** \brief  Degree of the largest complete polynomial space that can be represented by the basis
   */
  int basisDegree_;
  
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
  bool basisTagsAreSet_;
  
  /** \brief  DoF ordinal to tag lookup table.
              
              Rank-2 array with dimensions (basisCardinality_, 4) containing the DoF tags. This array 
              is left empty at instantiation and filled by initializeTags() only when tag data is
              requested. 
  
      \li     ordinalToTag_[DofOrd][0] = dim. of the subcell associated with the specified DoF 
      \li     ordinalToTag_[DofOrd][1] = ordinal of the subcell defined in the cell topology
      \li     ordinalToTag_[DodOrd][2] = ordinal of the specified DoF relative to the subcell
      \li     ordinalToTag_[DofOrd][3] = total number of DoFs associated with the subcell
   */
  std::vector<std::vector<int> > ordinalToTag_;
  
  /** \brief  DoF tag to ordinal lookup table.
  
              Rank-3 array with dimensions (maxScDim + 1, maxScOrd + 1, maxDfOrd + 1), i.e., the 
              columnwise maximums of the 1st three columns in the DoF tag table for the basis plus 1. 
              For every triple (subscDim, subcOrd, subcDofOrd) that is valid DoF tag data this array 
              stores the corresponding DoF ordinal. If the triple does not correspond to tag data, 
              the array stores -1. This array is left empty at instantiation and filled by 
              initializeTags() only when tag data is requested. 
  
      \li     tagToOrdinal_[subcDim][subcOrd][subcDofOrd] = Degree-of-freedom ordinal
   */
  std::vector<std::vector<std::vector<int> > > tagToOrdinal_;  
  
public:
    
  /** \brief  Destructor
   */
  virtual ~Basis() {}
  
        
  /** \brief  Evaluation of a FEM basis on a <strong>reference cell</strong>. 
       
              Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
              points in the <strong>reference cell</strong> for which the basis is defined. 

      \param  outputValues      [out] - variable rank array with the basis values
      \param  inputPoints       [in]  - rank-2 array (P,D) with the evaluation points
      \param  operatorType      [in]  - the operator acting on the basis functions    
  
      \remark For rank and dimension specifications of the output array see Section 
              \ref basis_md_array_sec.  Dimensions of <var>ArrayScalar</var> arguments are checked  
              at runtime if HAVE_INTREPID_DEBUG is defined.

      \remark A FEM basis spans a COMPLETE or INCOMPLETE polynomial space on the reference cell 
              which is a smooth function space. Thus, all operator types that are meaningful for the
              approximated function space are admissible. When the order of the operator exceeds the 
              degree of the basis, the output array is filled with the appropriate number of zeros.  
   */
  virtual void getValues(ArrayScalar &          outputValues,
                         const ArrayScalar &    inputPoints,
                         const EOperator        operatorType) const = 0;
  
    
  /** \brief  Evaluation of an FVD basis evaluation on a <strong>physical cell</strong>.  
  
              Returns values of <var>operatorType</var> acting on FVD basis functions for a set of
              points in the <strong>physical cell</strong> for which the FVD basis is defined. 
  
      \param  outputValues      [out] - variable rank array with the basis values
      \param  inputPoints       [in]  - rank-2 array (P,D) with the evaluation points
      \param  cellVertices      [in]  - rank-2 array (V,D) with the vertices of the physical cell
      \param  operatorType      [in]  - the operator acting on the basis functions
  
      \remark For rank and dimension specifications of the output array see Section  
              \ref basis_md_array_sec. Dimensions of <var>ArrayScalar</var> arguments are checked 
               at runtime if HAVE_INTREPID_DEBUG is defined.

      \remarks A typical FVD basis spans a BROKEN discrete space which is only piecewise smooth. For
               example, it could be a piecewise constant space defined with respect to a partition
               of the cell into simplices. Because differential operators are not meaningful for such
               spaces, the default operator type in this method is set to OPERATOR_VALUE.
   */    
  virtual void getValues(ArrayScalar &          outputValues,
                         const ArrayScalar &    inputPoints,
                         const ArrayScalar &    cellVertices,
                         const EOperator        operatorType = OPERATOR_VALUE) const = 0;
  
  /** \brief  Returns cardinality of the basis
    
      \return the number of basis functions in the basis
   */  
  virtual int getCardinality() const;
  
  
  /** \brief  Returns the degree of the basis.
    
      \return max. degree of the complete polynomials that can be represented by the basis.
   */
  virtual int getDegree() const;
  
  
  /** \brief  Returns the base cell topology for which the basis is defined. See Shards documentation
              http://trilinos.sandia.gov/packages/shards for definition of base cell topology.
    
      \return Base cell topology
   */
  virtual const shards::CellTopology getBaseCellTopology() const;

  
  /** \brief  Returns the basis type.
    
      \return Basis type
   */  
  virtual EBasis getBasisType() const;
  
  
  /** \brief  Returns the type of coordinate system for which the basis is defined
    
      \return Type of the coordinate system (Cartesian, polar, R-Z, etc.).
   */
  virtual ECoordinates getCoordinateSystem() const;
  
  
  /** \brief  DoF tag to ordinal lookup. 
    
      \param  subcDim           [in]  - tag field 0: dimension of the subcell associated with the DoF 
      \param  subcOrd           [in]  - tag field 1: ordinal of the subcell defined by cell topology
      \param  subcDofOrd        [in]  - tag field 2: ordinal of the DoF relative to the subcell.
    
      \return the DoF ordinal corresponding to the specified DoF tag data.
   */
  virtual int getDofOrdinal(const int subcDim,
                            const int subcOrd,
                            const int subcDofOrd);  

  /** \brief DoF tag to ordinal data structure */
  virtual const std::vector<std::vector<std::vector<int> > > &getDofOrdinalData( );

  
  /** \brief  DoF ordinal to DoF tag lookup.
    
      \param  dofOrd            [in]  - ordinal of the DoF whose tag is being retrieved
    
      \return reference to a vector with dimension (4) such that \n
      \li     element [0] = tag field 0  ->  dim. of the subcell associated with the specified DoF 
      \li     element [1] = tag field 1  ->  ordinal of the subcell defined by cell topology
      \li     element [2] = tag field 2  ->  ordinal of the specified DoF relative to the subcell
      \li     element [3] = tag field 3  ->  total number of DoFs associated with the subcell
   */
  virtual const std::vector<int>&  getDofTag(const int dofOrd);
  
  
  /** \brief  Retrieves all DoF tags. 
    
      \return reference to a vector of vectors with dimensions (basisCardinality_, 4) such that \n
      \li     element [DofOrd][0] = tag field 0 for the DoF with the specified ordinal
      \li     element [DofOrd][1] = tag field 1 for the DoF with the specified ordinal
      \li     element [DofOrd][2] = tag field 2 for the DoF with the specified ordinal
      \li     element [DofOrd][3] = tag field 3 for the DoF with the specified ordinal
   */
  virtual const std::vector<std::vector<int> >& getAllDofTags();
  
  


}; // class Basis


//--------------------------------------------------------------------------------------------//
//                                                                                            //
//                            Helper functions of the Basis class                             //
//                                                                                            //
//--------------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------------//
//                                                                                            //            
//                                        Argument checks                                     //
//                                                                                            //            
//--------------------------------------------------------------------------------------------//

/** \brief  Runtime check of the arguments for the getValues method in an HGRAD-conforming 
FEM basis. Verifies that ranks and dimensions of <var>ArrayScalar</var> input and output
arrays are consistent with the specified <var>operatorType</var>.

    \param  outputValues     [in]  - array of variable rank for the output basis values
    \param  inputPoints      [in]  - rank-2 array with dimensions (P,D) containing the points
    \param  operatorType     [in]  - operator applied to basis functions  
    \param  cellTopo         [in]  - base cell topology on which the basis is defined
    \param  basisCard        [in]  - cardinality of the basis
 */
template<class Scalar, class ArrayScalar>
void getValues_HGRAD_Args(ArrayScalar &                outputValues,
                          const ArrayScalar &          inputPoints,
                          const EOperator              operatorType,
                          const shards::CellTopology&  cellTopo,
                          const int                    basisCard);

/** \brief  Runtime check of the arguments for the getValues method in an HCURL-conforming 
FEM basis. Verifies that ranks and dimensions of <var>ArrayScalar</var> input and output
arrays are consistent with the specified <var>operatorType</var>.

    \param  outputValues     [in]  - array of variable rank for the output basis values
    \param  inputPoints      [in]  - rank-2 array with dimensions (P,D) containing the points
    \param  operatorType     [in]  - operator applied to basis functions  
    \param  cellTopo         [in]  - base cell topology on which the basis is defined
    \param  basisCard        [in]  - cardinality of the basis
 */
template<class Scalar, class ArrayScalar>
void getValues_HCURL_Args(ArrayScalar &                outputValues,
                          const ArrayScalar &          inputPoints,
                          const EOperator              operatorType,
                          const shards::CellTopology&  cellTopo,
                          const int                    basisCard);

/** \brief  Runtime check of the arguments for the getValues method in an HDIV-conforming 
FEM basis. Verifies that ranks and dimensions of <var>ArrayScalar</var> input and output
arrays are consistent with the specified <var>operatorType</var>.

    \param  outputValues     [in]  - array of variable rank for the output basis values
    \param  inputPoints      [in]  - rank-2 array with dimensions (P,D) containing the points
    \param  operatorType     [in]  - operator applied to basis functions  
    \param  cellTopo         [in]  - base cell topology on which the basis is defined
    \param  basisCard        [in]  - cardinality of the basis
 */
template<class Scalar, class ArrayScalar>
void getValues_HDIV_Args(ArrayScalar &                outputValues,
                          const ArrayScalar &          inputPoints,
                          const EOperator              operatorType,
                          const shards::CellTopology&  cellTopo,
                          const int                    basisCard);


/** \brief  This is an interface class for bases whose degrees of freedom
            can be associated with spatial locations in a reference element
            (typically interpolation points for interpolatory bases).
*/
template<class ArrayScalar>
class DofCoordsInterface {
public:
   /** \brief Pure virtual destructor (gives warnings if not included).
     *        Following "Effective C++: 3rd Ed." item 7 the implementation
     *        is included in the definition file. 
     */
   virtual ~DofCoordsInterface() = 0;

  /** \brief  Returns spatial locations (coordinates) of degrees of freedom on a
              <strong>reference cell</strong>; defined for interpolatory bases.

      \param  DofCoords      [out] - array with the coordinates of degrees of freedom,
                                     dimensioned (F,D)
   */
   virtual void getDofCoords(ArrayScalar & DofCoords) const = 0;
};


// include templated definitions
#include <Intrepid_BasisDef.hpp>
  

}// namespace Intrepid


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

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

