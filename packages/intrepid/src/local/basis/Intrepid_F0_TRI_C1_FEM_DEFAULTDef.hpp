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

/** \file   Intrepid_F0_TRI_C1_FEM_DEFAULTDef.hpp
    \brief  Definition file for FEM basis functions of degree 1 for 0-forms on TRI cells.
    \author Created by P. Bochev and D. Ridzal.
*/

namespace Intrepid {
    
template<class Scalar>
int Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::numDof_ = 3;
    
    
template<class Scalar>
Teuchos::Array<Teuchos::Array<Teuchos::Array<int> > > Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::tagToEnum_;



template<class Scalar>
Teuchos::Array<LocalDofTag> Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::enumToTag_;



template<class Scalar>
bool Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::isSet_ = false;



template<class Scalar>
Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::Basis_F0_TRI_C1_FEM_DEFAULT() {}



template<class Scalar>
void Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::initialize() {
  
  // Basis-dependent initializations
  int numBf    = 3;         // number of basis functions
  int tagSize  = 4;         // size of DoF tag
  int posScDim = 0;         // position in the tag, counting from 0, of the subcell dim 
  int posScId  = 1;         // position in the tag, counting from 0, of the subcell id
  int posBfId  = 2;         // position in the tag, counting from 0, of DoF Id relative to the subcell

  // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 
  int tags[]  = {
                  0, 0, 0, 1,
                  0, 1, 0, 1,
                  0, 2, 0, 1
                };
  
  // Basis-independent function sets tag and enum data in the static arrays:
  Intrepid::setEnumTagData(tagToEnum_,
                           enumToTag_,
                           tags,
                           numBf,
                           tagSize,
                           posScDim,
                           posScId,
                           posBfId);
}



template<class Scalar> 
void Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::getValues(VarContainer<Scalar>&                  outputValues,
                                                    const Teuchos::Array< Point<Scalar> >& inputPoints,
                                                    const EOperator                        operatorType) const {

  // Determine parameters to shape outputValues: number of points = size of input array
  int numPoints = inputPoints.size();
	
  // Complete polynomial basis of degree 1 (C1) has 3 basis functions on a triangle that are 0-forms in 2D
  int    numFields = 3;
  EField fieldType = FIELD_FORM_0;
  int    spaceDim  = 2;

  // temporaries
  int countPt  = 0;               // point counter
  Scalar x(0);                    // x coord
  Scalar y(0);                    // y coord
  Teuchos::Array<int> indexV(2);  // multi-index for values
  Teuchos::Array<int> indexD(3);  // multi-index for gradients and D1's

  // Shape the VarContainer for the output values using these values:
  outputValues.reset(numPoints,
                     numFields,
                     fieldType,
                     operatorType,
                     spaceDim);

  switch (operatorType) {
      
    case OPERATOR_VALUE:
      for (countPt = 0; countPt < numPoints; countPt++) {
          
          // Verify that all points are inside the TRI reference cell
#ifdef HAVE_INTREPID_DEBUG
          TEST_FOR_EXCEPTION( !MultiCell<Scalar>::inReferenceCell(CELL_TRI, inputPoints[countPt]),
                              std::invalid_argument,
                              ">>> ERROR (Basis_F0_TRI_C1_FEM_DEFAULT): Evaluation point is outside the TRI reference cell");
#endif
        x = (inputPoints[countPt])[0];
        y = (inputPoints[countPt])[1];
        indexV[0] = countPt;
          
        indexV[1] = 0;
        outputValues.setValue(1.0-x-y, indexV);
        
        indexV[1] = 1;
        outputValues.setValue(x, indexV);
        
        indexV[1] = 2;
        outputValues.setValue(y, indexV);
      }
      break;

    case OPERATOR_GRAD:
    case OPERATOR_D1:
      for (countPt=0; countPt<numPoints; countPt++) {
        indexD[0] = countPt;

        indexD[1] = 0;
        indexD[2] = 0;
        outputValues.setValue(-1.0, indexD);
        indexD[2] = 1;
        outputValues.setValue(-1.0, indexD);

        indexD[1] = 1;
        indexD[2] = 0;
        outputValues.setValue(1.0, indexD);
        indexD[2] = 1;
        outputValues.setValue(0.0, indexD);

        indexD[1] = 2;
        indexD[2] = 0;
        outputValues.setValue(0.0, indexD);
        indexD[2] = 1;
        outputValues.setValue(1.0, indexD);
      }
      break;

    case OPERATOR_CURL:
      for (countPt=0; countPt<numPoints; countPt++) {
        indexD[0] = countPt;

        indexD[1] = 0;
        indexD[2] = 0;
        outputValues.setValue(-1.0, indexD);
        indexD[2] = 1;
        outputValues.setValue(1.0, indexD);

        indexD[1] = 1;
        indexD[2] = 0;
        outputValues.setValue(0.0, indexD);
        indexD[2] = 1;
        outputValues.setValue(-1.0, indexD);

        indexD[1] = 2;
        indexD[2] = 0;
        outputValues.setValue(1.0, indexD);
        indexD[2] = 1;
        outputValues.setValue(0.0, indexD);
      }
      break;

    case OPERATOR_DIV:
      // should have been caught already as an exception
      break;

    case OPERATOR_D2:
    case OPERATOR_D3:
    case OPERATOR_D4:
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10:
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
                          ">>> ERROR (Basis_F0_TRI_C1_FEM_DEFAULT): Invalid operator type");
  }
}
  

  
template<class Scalar>
void Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::getValues(VarContainer<Scalar>&                  outputValues,
                                                    const Teuchos::Array< Point<Scalar> >& inputPoints,
                                                    const Cell<Scalar>&                    cell) const {
  TEST_FOR_EXCEPTION( (true),
                      std::logic_error,
                      ">>> ERROR (Basis_F0_TRI_C1_FEM_DEFAULT): FEM Basis calling an FV/D member function");
}


template<class Scalar>
int Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::getNumLocalDof() const {
    return numDof_;   
}


template<class Scalar>
int Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::getLocalDofEnumeration(const LocalDofTag dofTag) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return tagToEnum_[dofTag.tag_[0]][dofTag.tag_[1]][dofTag.tag_[2]];
}



template<class Scalar>
LocalDofTag Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::getLocalDofTag(int id) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_[id];
}



template<class Scalar>
const Teuchos::Array<LocalDofTag> & Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::getAllLocalDofTags() {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_;
}



template<class Scalar>
ECell Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::getCellType() const {
  return CELL_TRI;
}



template<class Scalar>
EBasis Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::getBasisType() const {
  return BASIS_FEM_DEFAULT;
}



template<class Scalar>
ECoordinates Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::getCoordinateSystem() const {
  return COORDINATES_CARTESIAN;
}



template<class Scalar>
int Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::getDegree() const {
  return 1;
}


}// namespace Intrepid
