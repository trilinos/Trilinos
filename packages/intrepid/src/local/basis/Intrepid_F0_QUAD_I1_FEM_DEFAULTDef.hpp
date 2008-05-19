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

/** \file   Intrepid_F0_QUAD_I1_FEM_DEFAULTDef.hpp
    \brief  Definition file for bi-linear FEM basis functions for 0-forms on QUAD cells.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {


template<class Scalar>
void Basis_F0_QUAD_I1_FEM_DEFAULT<Scalar>::initialize() {
  
  // Basis-dependent intializations
  int tagSize  = 4;         // size of DoF tag
  int posScDim = 0;         // poisition in the tag, counting from 0, of the subcell dim 
  int posScId  = 1;         // poisition in the tag, counting from 0, of the subcell id
  int posBfId  = 2;         // position in the tag, counting from 0, of DoF Id relative to the subcell

  // An array with local DoF tags assigned to basis functions, in the order of their local enumeration 
  int tags[]  = {
                  0, 0, 0, 1,
                  0, 1, 0, 1,
                  0, 2, 0, 1,
                  0, 3, 0, 1,
                };
  
  // Basis-independent function sets tag and enum data in the static arrays:
  setEnumTagData(tagToEnum_,
                 enumToTag_,
                 tags,
                 numDof_,
                 tagSize,
                 posScDim,
                 posScId,
                 posBfId);
}



template<class Scalar> 
void Basis_F0_QUAD_I1_FEM_DEFAULT<Scalar>::getValues(FieldContainer<Scalar>&               outputValues,
                                                    const Teuchos::Array< Point<Scalar> >& inputPoints,
                                                    const EOperator                        operatorType) const {

  // Number of evaluation points =  size of outputValues
  int numPoints = inputPoints.size();             
  
  // Shape the FieldContainer for the output:
  outputValues.resize(numPoints,                   // number of evaluation points
                      numDof_,                     // number of fields = number of DoFs in the basis
                      FIELD_FORM_0,                // field type of the basis functions
                      operatorType,                // operator type that is applied to basis functions
                      2);                          // space dimension for QUAD cell is 2
  
  // Temporaries: point counter and (x,y) coordinates of the evaluation point
  int countPt  = 0;                               
  Scalar x(0);                                    
  Scalar y(0);                                    
  
  switch (operatorType) {
    case OPERATOR_VALUE:
      for (countPt=0; countPt<numPoints; countPt++) {
#ifdef HAVE_INTREPID_DEBUG 
        // Verify argument: check if input point is inside the reference QUAD
        TEST_FOR_EXCEPTION( !MultiCell<Scalar>::inReferenceCell(CELL_QUAD, inputPoints[countPt]),
                            std::invalid_argument,
                            ">>> ERROR (Basis_F0_QUAD_I1_FEM_DEFAULT): Evaluation point is outside the QUAD reference cell");
#endif
        x = (inputPoints[countPt])[0];
        y = (inputPoints[countPt])[1];
        
        // Output container has rank 2. The indices are (P,F)
        outputValues(countPt, 0) = (1.0 - x)*(1.0 - y)/4.0;
        outputValues(countPt, 1) = (1.0 + x)*(1.0 - y)/4.0;
        outputValues(countPt, 2) = (1.0 + x)*(1.0 + y)/4.0;
        outputValues(countPt, 3) = (1.0 - x)*(1.0 + y)/4.0;
      }
      break;
      
    case OPERATOR_GRAD:
    case OPERATOR_D1:
      for (countPt=0; countPt<numPoints; countPt++) {
#ifdef HAVE_INTREPID_DEBUG 
        // Verify argument: check if input point is inside the reference QUAD
        TEST_FOR_EXCEPTION( !MultiCell<Scalar>::inReferenceCell(CELL_QUAD, inputPoints[countPt]),
                            std::invalid_argument,
                            ">>> ERROR (Basis_F0_QUAD_I1_FEM_DEFAULT): Evaluation point is outside the QUAD reference cell");
#endif
        x = (inputPoints[countPt])[0];
        y = (inputPoints[countPt])[1];
        
        // Output container has rank 3. The indices are (P,F,D)
        outputValues(countPt, 0, 0) = -(1.0 - y)/4.0;
        outputValues(countPt, 0, 1) = -(1.0 - x)/4.0;
        
        outputValues(countPt, 1, 0) =  (1.0 - y)/4.0;
        outputValues(countPt, 1, 1) = -(1.0 + x)/4.0;
        
        outputValues(countPt, 2, 0) =  (1.0 + y)/4.0;
        outputValues(countPt, 2, 1) =  (1.0 + x)/4.0;
  
        outputValues(countPt, 3, 0) = -(1.0 + y)/4.0;
        outputValues(countPt, 3, 1) =  (1.0 - x)/4.0;
      }
      break;
      
    case OPERATOR_CURL:
      for (countPt=0; countPt<numPoints; countPt++) {
#ifdef HAVE_INTREPID_DEBUG 
        // Verify argument: check if input point is inside the reference QUAD
        TEST_FOR_EXCEPTION( !MultiCell<Scalar>::inReferenceCell(CELL_QUAD, inputPoints[countPt]),
                            std::invalid_argument,
                            ">>> ERROR (Basis_F0_QUAD_I1_FEM_DEFAULT): Evaluation point is outside the QUAD reference cell");
#endif
        x = (inputPoints[countPt])[0];
        y = (inputPoints[countPt])[1];
        
        // Output container has rank 3. The indices are (P,F,D)
        outputValues(countPt, 0, 0) = -(1.0 - x)/4.0;
        outputValues(countPt, 0, 1) =  (1.0 - y)/4.0;
        
        outputValues(countPt, 1, 0) = -(1.0 + x)/4.0;
        outputValues(countPt, 1, 1) = -(1.0 - y)/4.0;
        
        outputValues(countPt, 2, 0) =  (1.0 + x)/4.0;
        outputValues(countPt, 2, 1) = -(1.0 + y)/4.0;
        
        outputValues(countPt, 3, 0) =  (1.0 - x)/4.0;
        outputValues(countPt, 3, 1) =  (1.0 + y)/4.0;
      }
      break;
      
    case OPERATOR_DIV:
      // should have been caught already as an exception by outputValues.reset(...)
      break;
      
    case OPERATOR_D2:
      for (countPt=0; countPt<numPoints; countPt++) {
#ifdef HAVE_INTREPID_DEBUG 
        // Verify argument: check if input point is inside the reference QUAD
        TEST_FOR_EXCEPTION( !MultiCell<Scalar>::inReferenceCell(CELL_QUAD, inputPoints[countPt]),
                            std::invalid_argument,
                            ">>> ERROR (Basis_F0_QUAD_I1_FEM_DEFAULT): Evaluation point is outside the QUAD reference cell");
#endif
        // Output container has rank 3. The indices are (P,F,K) 
        outputValues(countPt, 0, 0) =  0.0;
        outputValues(countPt, 0, 1) =  0.25;
        outputValues(countPt, 0, 2) =  0.0;

        outputValues(countPt, 1, 0) =  0.0;
        outputValues(countPt, 1, 1) = -0.25;
        outputValues(countPt, 1, 2) =  0.0;
        
        outputValues(countPt, 2, 0) =  0.0;
        outputValues(countPt, 2, 1) =  0.25;
        outputValues(countPt, 2, 2) =  0.0;
        
        outputValues(countPt, 3, 0) =  0.0;
        outputValues(countPt, 3, 1) = -0.25;
        outputValues(countPt, 3, 2) =  0.0;
      }
      break;
      
    case OPERATOR_D3:
    case OPERATOR_D4:
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10:
      
      // Output container has rank 3. The indices are (P,F,K). All values are zero
      outputValues.storeZero();
      break;
      
    default:
      TEST_FOR_EXCEPTION( ( (operatorType != OPERATOR_VALUE) &&
                            (operatorType != OPERATOR_GRAD)  &&
                            (operatorType != OPERATOR_CURL)  &&
                            (operatorType != OPERATOR_DIV)   &&
                            (operatorType != OPERATOR_D1)    &&
                            (operatorType != OPERATOR_D2)    &&
                            (operatorType != OPERATOR_D3)    &&
                            (operatorType != OPERATOR_D4)    &&
                            (operatorType != OPERATOR_D5)    &&
                            (operatorType != OPERATOR_D6)    &&
                            (operatorType != OPERATOR_D7)    &&
                            (operatorType != OPERATOR_D8)    &&
                            (operatorType != OPERATOR_D9)    &&
                            (operatorType != OPERATOR_D10) ),
                          std::invalid_argument,
                          ">>> ERROR (Basis_F0_QUAD_I1_FEM_DEFAULT): Invalid operator type");
  }
}


  
template<class Scalar>
void Basis_F0_QUAD_I1_FEM_DEFAULT<Scalar>::getValues(FieldContainer<Scalar>&                  outputValues,
                                                    const Teuchos::Array< Point<Scalar> >& inputPoints,
                                                    const Cell<Scalar>&                    cell) const {
  TEST_FOR_EXCEPTION( (true),
                      std::logic_error,
                      ">>> ERROR (Basis_F0_QUAD_I1_FEM_DEFAULT): FEM Basis calling an FVD member function");
}


template<class Scalar>
inline int Basis_F0_QUAD_I1_FEM_DEFAULT<Scalar>::getNumLocalDof() const {
    return numDof_;   
}


template<class Scalar>
int Basis_F0_QUAD_I1_FEM_DEFAULT<Scalar>::getLocalDofEnumeration(const LocalDofTag dofTag) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return tagToEnum_[dofTag.tag_[0]][dofTag.tag_[1]][dofTag.tag_[2]];
}



template<class Scalar>
LocalDofTag Basis_F0_QUAD_I1_FEM_DEFAULT<Scalar>::getLocalDofTag(int dofId) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_[dofId];
}



template<class Scalar>
const Teuchos::Array<LocalDofTag> & Basis_F0_QUAD_I1_FEM_DEFAULT<Scalar>::getAllLocalDofTags() {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_;
}



template<class Scalar>
inline ECell Basis_F0_QUAD_I1_FEM_DEFAULT<Scalar>::getCellType() const {
  return CELL_QUAD;
}



template<class Scalar>
inline EBasis Basis_F0_QUAD_I1_FEM_DEFAULT<Scalar>::getBasisType() const {
  return BASIS_FEM_DEFAULT;
}



template<class Scalar>
inline ECoordinates Basis_F0_QUAD_I1_FEM_DEFAULT<Scalar>::getCoordinateSystem() const {
  return COORDINATES_CARTESIAN;
}



template<class Scalar>
inline int Basis_F0_QUAD_I1_FEM_DEFAULT<Scalar>::getDegree() const {
  return 1;
}


}// namespace Intrepid
