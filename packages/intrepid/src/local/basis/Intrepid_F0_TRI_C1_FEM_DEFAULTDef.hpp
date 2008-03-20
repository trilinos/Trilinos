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
Teuchos::Array<Teuchos::Array<Teuchos::Array<int> > > Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::tagToEnum_;


template<class Scalar>
Teuchos::Array<LocalDofTag> Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::enumToTag_;


template<class Scalar>
bool Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::isSet_ = false;


template<class Scalar>
Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::Basis_F0_TRI_C1_FEM_DEFAULT() {
  if (!isSet_) {

    /************ This section will differ depending on the basis *************/
    ECell mycell = CELL_TRI;  // cell type
    int numBf    = 3;         // number of basis functions
    int tagSize  = 4;         // size of DoF tag
    int posScId  = 1;         // position, counting from 0, of subcell id in the tag
    int posBfId  = 2;         // position, counting from 0, of local Bf id in the tag

    int tags[]  = {
                    0, 0, 0, 1,
                    0, 1, 0, 1,
                    0, 2, 0, 1
                  };


    /************ This section is generic *************/
    // build enumToTag array
    enumToTag_.resize(numBf);
    for (int i=0; i<numBf; i++)
      for (int j=0; j<tagSize; j++)
        enumToTag_[i].tag_[j] = tags[i*tagSize+j];
    // build tagToEnum array
    int maxBfId      = 0;
    int maxScId = 0;
    int maxDim       = MultiCell<Scalar>::getTopologicalDim(mycell);  // max dimension for this cell type
    for (int i=0; i<numBf; i++)  // find max local Bf id
      if (maxBfId < tags[i*tagSize+posBfId])
        maxBfId = tags[i*tagSize+posBfId];
    for (int i=0; i<numBf; i++)  // find max subcell id
      if (maxBfId < tags[i*tagSize+posScId])
        maxScId = tags[i*tagSize+posScId];
    Teuchos::Array<int> level1(maxBfId+1, -1);
    Teuchos::Array<Teuchos::Array<int> > level2(maxScId+1, level1);
    tagToEnum_.assign(maxDim+1, level2);
    for (int i=0; i<numBf; i++)
      tagToEnum_[tags[i*tagSize]][tags[i*tagSize+1]][tags[i*tagSize+2]] = i;
  }
  isSet_ = true;
}

  
template<class Scalar> 
void Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::getValues(VarContainer<Scalar>&                  outputValues,
                                                    const Teuchos::Array< Point<Scalar> >& inputPoints,
                                                    const EOperator                        operatorType) {

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
      for (countPt=0; countPt<numPoints; countPt++) {
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
                          ">>> ERROR (VarContainer): Invalid operator type");
  }

}
  
  
template<class Scalar>
void Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::getValues(VarContainer<Scalar>&                  outputValues,
                                                    const Teuchos::Array< Point<Scalar> >& inputPoints,
                                                    const MultiCell<Scalar>&               mCell) {
  TEST_FOR_EXCEPTION( (true),
                      std::logic_error,
                      ">>> ERROR (VarContainer): FEM Basis calling an FV/D member function");
}


template<class Scalar>
int Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::getLocalDofEnumeration(const LocalDofTag dofTag) {
  return tagToEnum_[dofTag.tag_[0]][dofTag.tag_[1]][dofTag.tag_[2]];
}


template<class Scalar>
LocalDofTag Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::getLocalDofTag(int id) {
  return enumToTag_[id];
}


template<class Scalar>
void Basis_F0_TRI_C1_FEM_DEFAULT<Scalar>::getAllLocalDofTags(Teuchos::Array<LocalDofTag>& dofTags) {
  dofTags = enumToTag_;
}

  
}// namespace Intrepid
