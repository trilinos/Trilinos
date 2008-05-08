// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
template<class Scalar, class ArrayType = LexContainer<Scalar> >
class LocalForm0 : public LocalField<Scalar> {
  private:

  /** \brief Type of the local field (the <strong>native</strong> field)
  */
  const static EField                                                      fieldType_ = FIELD_FORM_0;

  
  /** \brief Refcount pointer to a basis class for the native field (the <strong>native</strong> basis)
  */
  Teuchos::RCP<Basis<Scalar> >                                             basis_;
  
  
  /** \brief Type of the cell for which the native basis is defined
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
  
  /** \brief A three-dimensional array of VarContainers whose leading dimension equals the number of
    operator types in Intrepid. Initially the array is empty and is filled by values on demand, i.e.,
    only when a particular set of values on particular subcells is requested by invocation of a method.
    Specifically: <var>basisVals_[operatorType][subcellDim][subcellId]</var>  is a VarContainer that
    stores basis functions evaluated at a set of cubature points such that
    
    \li <var>operatorType</var> is the operator applied to the basis functions;
    \li <var>subcellDim</var> is subcell dimension where the cubature points are located;
    \li <var>subcellId</var> is the local order of the subcell (remember, subcells of a given
        dimension are implicitely ordered in the cell template)                                                        
    
    Depending on <var>operatorType</var> the rank of the VarContainer at 
    <var>basisVals_[operatorType][subcellDim][subcellId]</var> is 2 or 3
  */
  Teuchos::Array<Teuchos::Array<Teuchos::Array<VarContainer<Scalar> > > >  basisVals_;

  
  /** \brief Two-dimensional array of refcount pointers to cubature classes, one per subcell, i.e.,
    
    \li <var>cubature_[0][0]</var> is the cell cubature class;
    \li <var>cubature_[1][faceId]</var> is the cubature on the face with local <var>faceId</var>;
    \li <var>cubature_[2][edgeId]</var> is the cubature on the edge with local <var>edgeId</var>.

    For two-dimensional cells only the first two sets of cubatures are initialized because edges and
    faces are the same 1-subcells.
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
  
    
  /** \brief Returns const reference to the VarContainer at <var>basisVals_[primOp][subDim][subcellId]</var>. 
    This is a helper function which implements the logic for the on-demand computation of basis function
    values at (preset) subcell cubature points. If the VarContainer at <var>basisVals_[primOp][subDim][subcellId]</var>
    is empty it is filled with the appropriate values, otherwise, a const reference to that
    VarContainer is returned directly without computing the values again.
    
    \param primOp           [in]     - Input operator (primitive).
    \param subDim           [in]     - Dimension of subcells at whose cubature points the values are computed.
    \param subcellId        [in]     - Subcell id.
    
    \return
    - Output container of rank 2 or 3. For formatting details see documentation of the public method
      getOperator(VarContainer<Scalar> &, const Teuchos::Array<Point<Scalar> > &, const EOperator)
    */
  const VarContainer<Scalar> & getOperator(const EOperator  primOp,
                                           const int        subDim,
                                           const int        subcellId);
  
  
  /** \brief Fills the Intrepid::LexContainer <var>values</var> with appropriately transformed 
    (to physical cells) values of the specified operator applied to the <strong>native</strong> basis functions. 
    
    \param transVals        [out]    - Output array.
    \param primOp           [in]     - Input operator (primitive).
    \param mCell            [in]     - Multicell.
    \param intDomain        [in]     - Integration domain (line, surface, cell).
    \param reuseJacobians   [in]     - if true forces method to store Jacobian and measure values
                                       used to compute integrals in the MultiCell argument. Otherwise
                                       all necessary for this method values are computed on the fly.
  */
  void transformBasisVals(LexContainer<Scalar> &      transVals,
                          const EOperator             primOp,
                          MultiCell<Scalar> &         mCell,
                          const bool                  reuseJacobians,
                          const EIntegrationDomain    intDomain);
  
  
  /** \briefFills the Intrepid::LexContainer <var>values</var> with appropriately transformed 
    (to physical cells) values of the specified operator applied to <strong>auxiliary</strong> basis 
    functions.
    Acts on auxiliary basis.
    
    \param transVals        [out]    - Output array.
    \param primOp           [in]     - Input operator (primitive).
    \param primOpField      [in]     - Domain of primOp - specifies the <strong>auxiliary</strong> basis 
    \param mCell            [in]     - Multicell.
    \param intDomain        [in]     - Integration domain (line, surface, cell).
    \param reuseJacobians   [in]     - if true forces method to store Jacobian and measure values
                                        used to compute integrals in the MultiCell argument. Otherwise
                                        all necessary for this method values are computed on the fly.
  */
  void transformBasisVals(LexContainer<Scalar> &  transVals,
                 const EOperator                  primOp,
                 const LocalField<Scalar> &       primOpField,
                 MultiCell<Scalar> &              mCell, 
                 const bool                       reuseJacobians,
                 const EIntegrationDomain         intDomain);
  
  
  /** \brief Fills the Intrepid::LexContainer <var>values</var> with appropriately transformed 
    (to physical cells) values of the specified operator applied to the <strong>native</strong> basis functions. 
    
    \param finalVals        [out]    - Output array with weighted measure applied to transformed values
    \param transVals        [in]     - Array with transformed basis function values.
    \param primOp           [in]     - Input operator (primitive).
    \param mCell            [in]     - Multicell.
    \param intDomain        [in]     - Integration domain (line, surface, cell).
    \param reuseJacobians   [in]     - if true forces method to store Jacobian and measure values
                                       used to compute integrals in the MultiCell argument. Otherwise
                                       all necessary for this method values are computed on the fly.
  */  
  void applyWeightedMeasure(LexContainer<Scalar> &         finalVals,
                            const LexContainer<Scalar> &   transVals,
                            const EOperator                primOp,
                            MultiCell<Scalar> &            mCell,
                            const bool                     reuseJacobians,
                            const EIntegrationDomain       intDomain);
  
  
  /** \brief Performs numerical integration (note that geometric transformations
    and integration weights have already been included in <var>leftValues</var>
    and <var>rightValues</var>).

  \param outputValues     [out]     - Output array.
  \param leftValues       [out]     - Left array.
  \param rightValues      [out]     - Right array.
*/
  void integrate(LexContainer<Scalar> &        outputValues,
                 const LexContainer<Scalar> &  leftValues,
                 const LexContainer<Scalar> &  rightValues) const;
  
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


  /** \brief Returns a VarContainer with the values of <var>primOp</var> applied to the FEM basis
    functions, evaluated at the array <var>inputPoints</var> of <strong>reference</strong> points.
    The context in which this interface function can be called is limited to FEM reconstructions 
    (i.e. via a reference cell). The function returns operator values <strong>decoupled from the 
    geometry of a particular cell</strong> (untransformed values). 
    
    Admissible <var>primOp</var> arguments and the format of the output container are as follows 
    (see also VarContainer::getEnumeration for a detailed list of VarContainer shapes)
    \verbatim
    |--------------------|----------------------|----------------|
    |       primOp       |  outputValues format | container rank |
    |--------------------|----------------------|----------------|
    |       VALUE        |    [P][F]            |       2        | 
    |--------------------|----------------------|----------------|
    |     GRAD, D1       |    [P][F][D]         |       3        |
    |--------------------|----------------------|----------------|
    |        CURL        |    [P][F][D]         |       3        |
    |--------------------|----------------------|----------------|
    |        DIV         |    undefined         |       -        | 
    |--------------------|----------------------|----------------|
    |    D1,D2,..,D10    |    [P][F][K]         |       3        | 
    |--------------------|----------------------|----------------|
    
Legend:
    P -> point index            range: 0 <= P < numPoints = inputPoints.size()
    F -> field index            range: 0 <= F < numFields = 4
    D -> field component index  range: 0 <= D < spaceDim  = 2
    K -> enumeration of Dk      range: 0 <= K < DkCardinality 
    \endverbatim
    
    \param outputValues   [out]         - VarContainer of rank 2 or 3 with the computed values
    \param inputPoints     [in]         - evaluation points on the reference cell  
    \param operatorType    [in]         - the operator being applied to the basis function    
    
    \remarks 
    \li Enumeration of Dk (derivatives of total order k) follows the lexicographical order of 
    the partial derivatives; see getDkEnumeration() for details.
    
    \li For a particular FEM basis all Dk values may be zero for all k greater than some fixed
    integer value. Nevertheless, the output container for such Dk is still shaped using DkCardinality
    as an upper bound for the last index, i.e., the output container is filled with as many zeroes as 
    there are partial derivatives of order k; see getDkCardinality. 
  */
  void getOperator(VarContainer<Scalar> &                  outputValues,
                   const Teuchos::Array<Point<Scalar> > &  inputPoints,
                   const EOperator                         primOp);


  /** \brief  Returns a VarContainer with the values of <var>primOp</var>  applied to FEM or FVD 
    basis functions, evaluated at the array <var>inputPoints</var> of <strong>physical</strong> points.
    Admissible <var>primOp</var> arguments and the format of the output container are as in 
    getOperator(VarContainer<Scalar> &, const Teuchos::Array<Point<Scalar> > &, const EOperator).
    
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
  void getOperator(VarContainer<Scalar> &                  outputValues,
                   const Teuchos::Array<Point<Scalar> > &  inputPoints,
                   const EOperator                         primOp,
                   const Cell<Scalar> &                    cell);


  void getOperator(LexContainer<Scalar> &          outputValues,
                   const EOperator                 leftOp,
                   const EOperator                 rightOp,
                   MultiCell<Scalar> &             mCell,
                   const Teuchos::Array<Scalar> &  inputData,
                   const EDataFormat               inputFormat,
                   const bool                      reuseJacobians = false,
                   const EIntegrationDomain        intDomain = INTEGRATION_DOMAIN_CELL);
  
  
  void getOperator(ArrayType &                 outputValues,
                   const EOperator             leftOp,
                   const EOperator             rightOp,
                   MultiCell <Scalar> &        mCell,
                   const bool                  reuseJacobians = false,
                   const EIntegrationDomain    intDomain = INTEGRATION_DOMAIN_CELL);

  
  void getOperator(LexContainer<Scalar> &           outputValues,
                   const EOperator                  leftOp,
                   const EOperator                  rightOp,
                   const LocalField<Scalar> &       rightOpField,
                   MultiCell<Scalar> &              mCell,
                   const Teuchos::Array<Scalar> &   inputData,
                   const EDataFormat                inputFormat,
                   const bool                       reuseJacobians = false,
                   const EIntegrationDomain         intDomain = INTEGRATION_DOMAIN_CELL);

  
  void getOperator(LexContainer<Scalar> &           outputValues,
                   const EOperator                  leftOp,
                   const EOperator                  rightOp,
                   const LocalField<Scalar> &       rightOpField,
                   MultiCell<Scalar> &              mCell,
                    const bool                      reuseJacobians = false,
                   const EIntegrationDomain         intDomain = INTEGRATION_DOMAIN_CELL);
  
  
  EField getFieldType() const {return fieldType_;}
  
  
  ECell  getCellType() const  {return basisCell_;}
  
  
  int    getNumCubPoints(const int subcellDim,
                         const int subcellId) const;

}; // class LocalForm0

}// end namespace Intrepid

#include "Intrepid_LocalForm0Def.hpp"

#endif
