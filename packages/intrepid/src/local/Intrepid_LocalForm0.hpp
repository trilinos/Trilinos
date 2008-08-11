// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copytest (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_LocalForm0.hpp
    \brief  Header file for the Intrepid::LocalForm0 class.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_LOCAL_FORM_0_HPP
#define INTREPID_LOCAL_FORM_0_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_LocalField.hpp"
#include "Intrepid_Basis.hpp"
#include "Intrepid_Cubature.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_BLAS.hpp"


namespace Intrepid {

/** \class Intrepid::LocalForm0
    \brief Defines the base class for the representation of differential 0 forms in Intrepid.

    The local field interface relies on two member functions, <var>getOperator</var>
    and <var>getFunctional</var>, with varying signatures. Detailed info here ...
*/
template<class Scalar, class ArrayType = FieldContainer<Scalar> >
class LocalForm0 : public LocalField<Scalar> {
  private:

  /** \brief Type of the local field (the <strong>native</strong> field)
  */
  const static EField                                                      fieldType_ = FIELD_FORM_0;

  /** \brief Refcount pointer to a basis class for the native field (the <strong>native</strong> basis)
  */
  Teuchos::RCP<Basis<Scalar> >                                             basis_;
  
  /** \brief Type of the cell for which the native basis is defined, i.e, the native cell type
  */
  ECell                                                                    basisCell_;
  
  /** \brief The coordinate system for the native basis functions
  */
  ECoordinates                                                             basisCoordSystem_;
  
  /** \brief Type of the native basis (FEM_DEFAULT, FVD_DEFAULT, etc)
  */
  EBasis                                                                   basisType_;
  
  /** \brief Degree of the native basis functions
  */
  int                                                                      basisDegree_;
  
  /** \brief Array, indexed by basis function Id, containing the DoF tags assigned to each basis function 
  */
  Teuchos::Array<LocalDofTag>                                              basisDofTags_;
  
  /** \brief Dimension of the native basis, i.e., the number of basis functions
  */
  int                                                                      basisNumDofs_;
  
  /** \brief A 3-dimensional ragged array of FieldContainer objects whose leading dimension equals the 
    number of operator types in Intrepid and the 2nd and 3rd dimensions are the same as the 1st and 
    2nd dimensions of the cubature_ array. The first dimension of this array is resized by the
    constructor, the 2nd and 3rd are resized when the private getOperator(EOperator,int,int) is called
    for a first time, by copying the dimensions from cubature_ array. In this way, basisVals_ has the
    minimal dimensions needed to store basis values at the specified cubature poits. The FieldContainer
    at <var>basisVals_[operatorType][dimIndex][subcellId]</var> is always empty after the above
    resizings take place and is sized and filled with data only when this data is requested for a first
    time. After such an event, this FieldContainer will contain values of the <strong>native</strong> 
    basis functions evaluated at a set of cubature points where
    
    \li <var>operatorType</var> is the operator applied to the basis functions;
    \li <var>dimIndex</var> is cell dim. minus dimension of the subcell with the cubature points;
    \li <var>subcellId</var> is the local order of the subcell (remember, subcells of a given
        dimension are implicitely ordered in the cell template)                                                        
    
    Depending on <var>operatorType</var> the rank of the FieldContainer at 
    <var>basisVals_[operatorType][dimIndex][subcellId]</var> is 2 or 3
  */
  Teuchos::Array<Teuchos::Array<Teuchos::Array<FieldContainer<Scalar> > > >  basisVals_;

  
  /** \brief Two-dimensional ragged array of refcount pointers to subcell cubature classes. This array is
    an input argument to the LocalForm0 constructor and should be filled according to the following rules:
    
     - every subcell in the native cell can have its own cubature;
     - if cell dimension is 3 and cubatures are provided for all subcells, including the native cell itself,  
      -#  <var>cubature_[0][0]</var> is the cell cubature;
      -# <var>cubature_[1][faceId]</var> is the cubature on the face with local <var>faceId</var>;
      -# <var>cubature_[2][edgeId]</var> is the cubature on the edge with local <var>edgeId</var>;
    
    - if cell dimension is 2 and cubatures are provided for all subcells, including the native cell itself,  
      -# <var>cubature_[0][0]</var> is the cell cubature;
      -# <var>cubature_[1][subcellId]</var> is the cubature on 1-subcell with local <var>subcellId</var>;
  
    - cubatures can be provided only for selected subcell dimensions, e.g., faces. In this case, 
      cubature_ array is always shaped using the <strong>minimal</strong> dimensions needed to store them.
    
    Examples: dimensions of cubature_ assuming  CELL_TET and different cubature allocations
    
    - if only native cell cubature is specified: <var>cubature_</var> has dimensions (1,1)
      -#  <var>cubature_[0][0]</var> contains the cell cubature
    
    - if only edge cubatures are specified: <var>cubature_</var> has dimension 3;
      -# <var>cubature_[0]</var> has dimension 0 (is empty);
      -# <var>cubature_[1]</var> has dimension 0 (is empty)
      -# <var>cubature_[2]</var> has dimension 6 and contains RCP's to 6 possibly different edge cubatures: \n
         <var>cubature_[2][edgeId]</var> is a valid CELL_EDGE cubature for the edge with the specified edgeId 
    
    - if cell and edge cubatures are specified: <var>cubature_</var> has dimension 3;
      -# <var>cubature_[0]</var> has dimension 1 and contains the cell cubature: \n
         <var>cubature_[0][0]</var> is any valid CEL_TRI cubature;
      -# <var>cubature_[1]</var> has dimension 0 (is empty)
      -# <var>cubature_[2]</var> has dimension 6 and contains RCP's to 6 possibly different edge cubatures: \n
         <var>cubature_[2][edgeId]</var> is a valid CELL_EDGE cubature for the edge with the specified edgeId  
  */
  Teuchos::Array<Teuchos::Array<Teuchos::RCP<Cubature<Scalar> > > >        cubature_;
  
  
  /** \brief Three-dimensional array of cubature point sets, one set per subcell, i.e.,
    
    \li <var>cubPoints_[0][0][i]</var> is the ith cubature point on the cell;
    \li <var>cubPoints_[1][faceId][i]</var> is the ith cubature point on the face with local <var>faceId</var>;
    \li <var>cubPoints_[2][edgeId][i]</var> is the ith cubature point on the edge with local <var>edgeId</var>;
  
    For two-dimensional cells only the first two sets of cubature points are initialized because edges
    and faces are the same 1-subcells.
  */
  Teuchos::Array<Teuchos::Array<Teuchos::Array<Point<Scalar> > > >         cubPoints_; 
  
  
  /** \brief Three-dimensional array of cubature point weight sets, one set per subcell, i.e.,
    
    \li <var>cubWeights_[0][0][i]</var> is the ith cubature weight on the cell;
    \li <var>cubWeights_[1][faceId][i]</var> is the ith cubature weight on the face with local <var>faceId</var>;
    \li <var>cubWeights_[2][edgeId][i]</var> is the ith cubature weight on the edge with local <var>edgeId</var>;
  
  For two-dimensional cells only the first two sets of cubature weights are initialized because edges
  and faces are the same 1-subcells.
  */
  Teuchos::Array<Teuchos::Array<Teuchos::Array<Scalar> > >                 cubWeights_;
  
  
  /** \brief Two-dimensional array of the number of cubature points per subcell, i.e., 
    
    \li <var>numCubPoints_[0][0]</var> is the number of cubature points in the cell;
    \li <var>numCubPoints_[1][faceId]</var> is the number of cubature points in the face with local <var>faceId</var>
    \li <var>numCubPoints_[2][edgeId]</var> is the number of cubature points in the edge with local <var>edgeId</var>
    
    For two-dimensional cells only the first two sets of values are initialized because edges and
    faces are the same 1-subcells.
  */
  Teuchos::Array<Teuchos::Array<int> >                                     numCubPoints_;

  
  /** \brief Specifies how to carry out low-level linear algebra operations (BLAS or native)
  */
  ECompEngine                                                              compEngine_;
  
  //===========================================================================//
  //                                                                           //
  //                         Private methods of LocalForm0                     //
  //                                                                           //
  //===========================================================================//
  
  /** \brief Returns const reference to the FieldContainer at <var>basisVals_[primOp][dimIndex][subcellId]</var>. 
    This is a helper function which implements the logic for the on-demand computation of basis function
    values at (preset) subcell cubature points. If the FieldContainer at <var>basisVals_[primOp][dimIndex][subcellId]</var>
    is empty it is filled with the appropriate values, otherwise, a const reference to that
    FieldContainer is returned directly without computing the values again. Note that <var>basisVals_</var>
    stores values using <var>dimIndex</var> defined as cell dimension - subcell dimension.
    
    \param primOp           [in]     - Input operator (primitive).
    \param subDim           [in]     - Dimension of subcells at whose cubature points the values are computed.
    \param subcellId        [in]     - Subcell id.
    
    \return
    - Output container of rank 2 or 3. For formatting details see documentation of the public method
      getOperator(FieldContainer<Scalar> &, const Teuchos::Array<Point<Scalar> > &, const EOperator)
    */
  const FieldContainer<Scalar> & getOperator(const EOperator  primOp,
                                             const int        subDim,
                                             const int        subcellId);
  
  
  /** \brief Returns a FieldContainer of rank 3,4 or 5 with properly transformed (to the physical cells of the 
    <var>mCell</var> argument) reference values from <var>basisVals_[primOp][subcellDim][subcellId]</var>.
    Recall that <var>basisVals_[primOp][subcellDim][subcellId]</var> stores values of the specified operator 
    applied to the <strong>native</strong> basis functions and evaluated at the cubature point set
    associated with the subcell with <var>subcellDim</var> and <var>subcellId</var>.
    The rank of the output container, its dimensions and its multi-index are as follows:
    \verbatim
    |------------------|----------|-------------------|------------------------------|
    |   output field   | its rank |   container rank  |     container multi-index    |
    |------------------|----------|-------------------|------------------------------|
    |      scalar      |     0    |         3         |     (C, P, F)                |
    |------------------|----------|-------------------|------------------------------|
    |      vector      |     1    |         4         |     (C, P, D, F)             |
    |------------------|----------|-------------------|------------------------------|
    |      tensor      |     2    |         5         |     (C, P, D, D, F)          |
    |------------------|----------|-------------------|------------------------------|

    |------|----------------------|--------------------------------------------------|
    |      |         Index        |                        Dimension                 |
    |------|----------------------|--------------------------------------------------|
    |   C  |         cell         |  0 <= C < num. integration domains               |
    |   P  |         point        |  0 <= P < num. cubature points on the subcell    |
    |   D  |   field coordinate   |  0 <= D < dim. of the cell of the LocalField0    |
    |   F  |         field        |  0 <= F < dim. of the native basis               |
    |------|----------------------|--------------------------------------------------|
    \endverbatim
    \remarks
    \li The type of <var>intDomain</var> implicitly defines the value of <var>subcellDim</var>. 
    \li Number of cubature points for the specified subcell is  <var>numCubPoints_[dimIndex][subcellId]</var>
    \li The field coordinate index (D) in this container preceedes the field index (F) in order to enable
        utilization of BLAS. This index order corresponds to a column-major storage, per cell, of
        the basis function values needed to form the integrals.
    
    \param transVals        [out]    - Output array.
    \param primOp           [in]     - Input operator (primitive).
    \param mCell            [in]     - Multicell.
    \param reuseJacobians   [in]     - if true forces method to store Jacobian and integration measure values
                                       used to compute integrals in the MultiCell argument. Otherwise
                                       all necessary values for this method are computed on the fly.
    \param intDomain        [in]     - Integration domain (line, surface, cell).
  */
  void transformBasisVals(FieldContainer<Scalar> &    transVals,
                          const EOperator             primOp,
                          MultiCell<Scalar> &         mCell,
                          const bool                  reuseJacobians,
                          const EIntegrationDomain    intDomain);
  
  
  /** \brief Returns a FieldContainer with the same type of values as the previous method but using 
    reference values of the <strong>auxiliary</strong> basis from the specified <var>primOpField</var>.
    
    \param transVals        [out]    - Output array.
    \param primOp           [in]     - Input operator (primitive).
    \param primOpField      [in]     - Domain of primOp - specifies the <strong>auxiliary</strong> basis 
    \param mCell            [in]     - Multicell.
    \param reuseJacobians   [in]     - if true forces method to store Jacobian and measure values
                                        used to compute integrals in the MultiCell argument. Otherwise
                                        all necessary for this method values are computed on the fly.
    \param intDomain        [in]     - Integration domain (line, surface, cell)
  */
  void transformBasisVals(FieldContainer<Scalar> &    transVals,
                          const EOperator             primOp,
                          const LocalField<Scalar> &  primOpField,
                          MultiCell<Scalar> &         mCell, 
                          const bool                  reuseJacobians,
                          const EIntegrationDomain    intDomain);
  
  
  /** \brief Returns FieldContainer of the same rank as the input <var>transValues</var> container
    with the weighted integration measure applied to the already transformed reference basis values.
    
    \param finalVals        [out]    - Output array with weighted measure applied to transformed values
    \param transVals        [in]     - Array with transformed basis function values.
    \param primOp           [in]     - Input operator (primitive).
    \param mCell            [in]     - Multicell.
    \param reuseJacobians   [in]     - if true forces method to store Jacobian and measure values
                                       used to compute integrals in the MultiCell argument. Otherwise
                                       all necessary for this method values are computed on the fly.
    \param intDomain        [in]     - Integration domain (line, surface, cell).
  */  
  void applyWeightedMeasure(FieldContainer<Scalar> &         finalVals,
                            const FieldContainer<Scalar> &   transVals,
                            const EOperator                  primOp,
                            MultiCell<Scalar> &              mCell,
                            const bool                       reuseJacobians,
                            const EIntegrationDomain         intDomain);
  
  
  /** \brief Obsolete
    \param outputValues     [out]     - Output array: rank-3 container with multi-index (C,LBF,RBF)
    \param trialValues       [in]      - Trial array.
    \param testValues      [in]      - Test array.
 */
  void integrate(FieldContainer<Scalar> &        outputValues,
                 const FieldContainer<Scalar> &  trialValues,
                 const FieldContainer<Scalar> &  testValues) const;
  
  public:

  /** brief  Creates a LocalForm0 from a basis, cubature and computational engine. 

      \param basis         [in]     - Refcount pointer to an Intrepid::Basis object.
      \param cubature      [in]     - Two-dimensional array of refcount pointers to an Intrepid::Cubature 
                                       object, one for each subcell. See below for formatting rules.
      \param compEngine    [in]     - Computational engine (manual, BLAS, etc.).
    
    The <var>cubature</var> array has to be filled with refcount pointers to cubature classses as follows:
    
    \li <var>cubature[0][0]</var> - RCP for the cell cubature;
    \li <var>cubature[1][faceId]</var> - RCP for the cubature on the face with local <var>faceId</var>;
    \li <var>cubature[2][edgeId]</var> - RCP for the cubature on the edge with local <var>edgeId</var>.
    
    For two-dimensional cells only the first two sets of cubatures have to be initialized because
    edges and faces are the same 1-subcells.
  */
  LocalForm0(const Teuchos::RCP<Basis<Scalar> > basis,
             const Teuchos::Array<Teuchos::Array<Teuchos::RCP<Cubature<Scalar> > > > cubature, 
             const ECompEngine compEngine = COMP_CPP);


  /** \brief Destructor
  */
  ~LocalForm0() {}

  /** \brief  Fills <var>outputValues</var> with values of <var>primOp</var>, acting on a set of FEM 
    basis functions, at a set of <strong>reference cell</strong> points <var>inputPoints</var>.     
    The context in which this interface function can be called is limited to FEM reconstructions 
    (i.e. via a reference cell). The function returns operator values <strong>decoupled from the 
    geometry of a particular cell</strong> (untransformed values).  
    
    Admissible <var>primOp</var> arguments, rank, format and dimensions of <var>outputValues</var>
    are as follows (see FieldContainer::resize(int,int,EField,EOperator,int) for a complete list)
    \verbatim
    |--------------------|----------------------|----------------|
    |       primOp       |  outputValues format | container rank |
    |--------------------|----------------------|----------------|
    |       VALUE        |      (P,F)           |       2        | 
    |--------------------|----------------------|----------------|
    |     GRAD, D1       |      (P,F,D)         |       3        |
    |--------------------|----------------------|----------------|
    |        CURL        |      (P,F,D)         |       3        |
    |--------------------|----------------------|----------------|
    |        DIV         |      undefined       |       -        | 
    |--------------------|----------------------|----------------|
    |    D1,D2,..,D10    |      (P,F,K)         |       3        | 
    |--------------------|----------------------|----------------|
    
    |------|----------------------|---------------------------|
    |      |         Index        |         Dimension         |
    |------|----------------------|---------------------------|
    |   P  |         point        |  0 <= P < numPoints       |
    |   F  |         field        |  0 <= F < numFields       |
    |   D  |   field coordinate   |  0 <= D < spaceDim (= 2)  |
    |   K  |   enumeration of Dk  |  0 <= K < DkCardinality   |
    |------|----------------------|---------------------------|
    \endverbatim
    
    \remarks 
    \li The field index (F) in this container preceeded the field coordinate index (D). This choice
    places field coordinates in a contiguous order and facilitates application of transformations.
    \li Enumeration of Dk (derivatives of total order k) follows the lexicographical order of 
    the partial derivatives; see getDkEnumeration() for details.
    
    \li For a particular FEM basis all Dk values may be zero for all k greater than some fixed
    integer value. Nevertheless, the output container for such Dk is still shaped using DkCardinality
    as an upper bound for the last index, i.e., the output container is filled with as many zeroes as 
    there are partial derivatives of order k; see getDkCardinality.     
    
    \param outputValues   [out]         - FieldContainer of rank 2 or 3 with the computed values
    \param inputPoints     [in]         - evaluation points on the reference cell  
    \param primOp          [in]         - the (primitive) operator being applied to the basis function    
  */
  void getOperator(ArrayType &                             outputValues,
                   const Teuchos::Array<Point<Scalar> > &  inputPoints,
                   const EOperator                         primOp);


  /** \brief  Fills <var>outputValues</var> with values of a (primitive) operator <var>primOp</var> 
    acting on FEM or FVD basis functions, at a set of <strong>physical cell</strong> points 
    <var>inputPoints</var>. Admissible <var>primOp</var> arguments and the format of the output container 
    are as in  getOperator(ArrayType &, const Teuchos::Array<Point<Scalar> > &, const EOperator).
    
    The context in which this interface function can be called is two-fold. For FEM reconstructions 
    (i.e. those based on a reference cell), the function returns the values of the operator applied 
    to basis functions in <strong>physical</strong> space, i.e. relevant geometric transformations 
    will be incorporated into the evaluation. For FV/D reconstructions, the function computes operator 
    values <strong>directly</strong> on physical cells.

      \param outputValues    [out]     - Output array.
      \param inputPoints      [in]     - Array of input (physical) points.
      \param primOp           [in]     - Input operator (primitive).
      \param cell             [in]     - Physical cell.
  */
  void getOperator(ArrayType &                            outputValues,
                   const Teuchos::Array<Point<Scalar> > & inputPoints,
                   const EOperator                        primOp,
                   const Cell<Scalar> &                   cell);


  void getOperator(ArrayType &                  outputValues,
                   const EOperator              trialOp,
                   const EOperator              testOp,
                   MultiCell<Scalar> &          mCell,
                   const ArrayType &            inputData,
                   const bool                   reuseJacobians = false,
                   const EIntegrationDomain     intDomain = INTEGRATION_DOMAIN_CELL);
  
  
  void getOperator(ArrayType &                  outputValues,
                   const EOperator              trialOp,
                   const EOperator              testOp,
                   MultiCell <Scalar> &         mCell,
                   const bool                   reuseJacobians = false,
                   const EIntegrationDomain     intDomain = INTEGRATION_DOMAIN_CELL);

  
  void getOperator(ArrayType &                  outputValues,
                   const EOperator              trialOp,
                   const EOperator              testOp,
                   const LocalField<Scalar> &   testOpField,
                   MultiCell<Scalar> &          mCell,
                   const ArrayType &            inputData,
                   const bool                   reuseJacobians = false,
                   const EIntegrationDomain     intDomain = INTEGRATION_DOMAIN_CELL);

  
  void getOperator(ArrayType &                  outputValues,
                   const EOperator              trialOp,
                   const EOperator              testOp,
                   const LocalField<Scalar> &   testOpField,
                   MultiCell<Scalar> &          mCell,
                    const bool                  reuseJacobians = false,
                   const EIntegrationDomain     intDomain = INTEGRATION_DOMAIN_CELL);
  
  
  void getFunctional(ArrayType &              outputValues,
                     const ArrayType &        trialData,
                     const EOperator          testOp,
                     MultiCell<Scalar> &      mCell,
                     const bool               reuseJacobians = false,
                     const EIntegrationDomain intDomain = INTEGRATION_DOMAIN_CELL);
  
  
  EField getFieldType() const {return fieldType_;}
  
  
  ECell  getCellType() const  {return basisCell_;}
  
  
  int    getNumCubPoints(const int subcellDim,
                         const int subcellId) const;

}; // class LocalForm0

}// end namespace Intrepid

#include "Intrepid_LocalForm0Def.hpp"

#endif
