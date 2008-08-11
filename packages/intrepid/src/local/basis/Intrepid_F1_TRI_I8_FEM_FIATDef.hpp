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
//                    Denis Ridzal (dridzal@sandia.gov) or
//                    Robert Kirby (robert.c.kirby@ttu.edu)
//                    FIAT (fiat-dev@fenics.org) (must join mailing list first, see www.fenics.org)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_F1_TRI_I8_FEM_FIATDef.hpp
    \brief  Definition file for FEM basis functions of degree 7 for 1-forms on TRI cells.
    \author Created by R. Kirby via the FIAT project
*/

namespace Intrepid {

template<class Scalar>
void Basis_F1_TRI_I8_FEM_FIAT<Scalar>::initialize() {

  // Basis-dependent initializations
  int tagSize  = 4;         // size of DoF tag
  int posScDim = 0;         // position in the tag, counting from 0, of the subcell dim 
  int posScId  = 1;         // position in the tag, counting from 0, of the subcell id
  int posBfId  = 2;         // position in the tag, counting from 0, of DoF Id relative to the subcell

  // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 
  int tags[] = { 1 , 0 , 0 , 8 , 1 , 0 , 1 , 8 , 1 , 0 , 2 , 8 , 1 , 0 , 3 , 8 , 1 , 0 , 4 , 8 , 1 , 0 , 5 , 8 , 1 , 0 , 6 , 8 , 1 , 0 , 7 , 8 , 1 , 1 , 0 , 8 , 1 , 1 , 1 , 8 , 1 , 1 , 2 , 8 , 1 , 1 , 3 , 8 , 1 , 1 , 4 , 8 , 1 , 1 , 5 , 8 , 1 , 1 , 6 , 8 , 1 , 1 , 7 , 8 , 1 , 2 , 0 , 8 , 1 , 2 , 1 , 8 , 1 , 2 , 2 , 8 , 1 , 2 , 3 , 8 , 1 , 2 , 4 , 8 , 1 , 2 , 5 , 8 , 1 , 2 , 6 , 8 , 1 , 2 , 7 , 8 , 2 , 0 , 0 , 56 , 2 , 0 , 1 , 56 , 2 , 0 , 2 , 56 , 2 , 0 , 3 , 56 , 2 , 0 , 4 , 56 , 2 , 0 , 5 , 56 , 2 , 0 , 6 , 56 , 2 , 0 , 7 , 56 , 2 , 0 , 8 , 56 , 2 , 0 , 9 , 56 , 2 , 0 , 10 , 56 , 2 , 0 , 11 , 56 , 2 , 0 , 12 , 56 , 2 , 0 , 13 , 56 , 2 , 0 , 14 , 56 , 2 , 0 , 15 , 56 , 2 , 0 , 16 , 56 , 2 , 0 , 17 , 56 , 2 , 0 , 18 , 56 , 2 , 0 , 19 , 56 , 2 , 0 , 20 , 56 , 2 , 0 , 21 , 56 , 2 , 0 , 22 , 56 , 2 , 0 , 23 , 56 , 2 , 0 , 24 , 56 , 2 , 0 , 25 , 56 , 2 , 0 , 26 , 56 , 2 , 0 , 27 , 56 , 2 , 0 , 28 , 56 , 2 , 0 , 29 , 56 , 2 , 0 , 30 , 56 , 2 , 0 , 31 , 56 , 2 , 0 , 32 , 56 , 2 , 0 , 33 , 56 , 2 , 0 , 34 , 56 , 2 , 0 , 35 , 56 , 2 , 0 , 36 , 56 , 2 , 0 , 37 , 56 , 2 , 0 , 38 , 56 , 2 , 0 , 39 , 56 , 2 , 0 , 40 , 56 , 2 , 0 , 41 , 56 , 2 , 0 , 42 , 56 , 2 , 0 , 43 , 56 , 2 , 0 , 44 , 56 , 2 , 0 , 45 , 56 , 2 , 0 , 46 , 56 , 2 , 0 , 47 , 56 , 2 , 0 , 48 , 56 , 2 , 0 , 49 , 56 , 2 , 0 , 50 , 56 , 2 , 0 , 51 , 56 , 2 , 0 , 52 , 56 , 2 , 0 , 53 , 56 , 2 , 0 , 54 , 56 , 2 , 0 , 55 , 56 };
  
  // Basis-independent function sets tag and enum data in the static arrays:
  Intrepid::setEnumTagData(tagToEnum_,
                           enumToTag_,
                           tags,
                           numDof_,
                           tagSize,
                           posScDim,
                           posScId,
                           posBfId);
}

template<class Scalar> 
void Basis_F1_TRI_I8_FEM_FIAT<Scalar>::getValues(FieldContainer<Scalar>&                  outputValues,
          const Teuchos::Array< Point<Scalar> >& inputPoints,
          const EOperator                        operatorType) const {

  // Determine parameters to shape outputValues: number of points = size of input array
  int numPoints = inputPoints.size();
	
  // Incomplete polynomial basis of degree 7 (I7) has 80 basis functions on a triangle that are 1-forms in 2D
  int    numFields = 80;
  EField fieldType = FIELD_FORM_1;
  int    spaceDim  = 2;

  // temporaries
  int countPt  = 0;               // point counter
  Teuchos::Array<int> indexV(3);  // multi-index for values
  Teuchos::Array<int> indexD(2);  // multi-index for curl

  // Shape the FieldContainer for the output values using these values:
  outputValues.resize(numPoints,
                     numFields,
                     fieldType,
                     operatorType,
                     spaceDim);

#ifdef HAVE_INTREPID_DEBUG
  for (countPt=0; countPt<numPoints; countPt++) {
    // Verify that all points are inside the TRI reference cell
    TEST_FOR_EXCEPTION( !MultiCell<Scalar>::inReferenceCell(CELL_TRI, inputPoints[countPt]),
                        std::invalid_argument,
                        ">>> ERROR (Basis_F1_TRI_I8_FEM_FIAT): Evaluation point is outside the TRI reference cell");
  }
#endif
  switch(operatorType) {
    case OPERATOR_VALUE:   
    {
      Teuchos::SerialDenseMatrix<int,Scalar> expansions(45,numPoints);
      Teuchos::SerialDenseMatrix<int,Scalar> result(80,numPoints);
      OrthogonalExpansions<Scalar>::tabulate(CELL_TRI,8,inputPoints,expansions);
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm0_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<80;i++) {
          outputValues(countPt,i,0) = result(i,countPt);
        }
      }
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm1_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<80;i++) {
          outputValues(countPt,i,1) = result(i,countPt);
        }
      }
    }
    break;

    case OPERATOR_D1:
    case OPERATOR_D2:
    case OPERATOR_D3:
    case OPERATOR_D4:
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10:
    {
      TEST_FOR_EXCEPTION( true , std::invalid_argument, ">>>ERROR F1_TRI_I8_FEM_FIAT: operator not implemented" );
    }
    break; 
    case OPERATOR_CURL:
    {
      Teuchos::SerialDenseMatrix<int,Scalar> expansions(45,numPoints);
      Teuchos::SerialDenseMatrix<int,Scalar> tmp(45,45);
      Teuchos::SerialDenseMatrix<int,Scalar> curl(80,numPoints);
      OrthogonalExpansions<Scalar>::tabulate(CELL_TRI,8,inputPoints,expansions);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats0_,expansions,0.0);
      curl.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm1_,tmp,0.0);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats1_,expansions,0.0);
      curl.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm0_,tmp,-1.0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<80;i++) {
          outputValues(countPt,i) = curl(i,countPt);
        }
      }
    }
    break;

    case OPERATOR_DIV:
    {
      TEST_FOR_EXCEPTION( true , std::invalid_argument, ">>>ERROR F1_TRI_I8_FEM_FIAT: operator not implemented" );
    }
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
                          ">>> ERROR (Basis_F1_TRI_I8_FEM_DEFAULT): Invalid operator type");
  }


}


template<class Scalar>
void Basis_F1_TRI_I8_FEM_FIAT<Scalar>::getValues(FieldContainer<Scalar>&                  outputValues,
                                                    const Teuchos::Array< Point<Scalar> >& inputPoints,
                                                    const Cell<Scalar>&                    cell) const {
  TEST_FOR_EXCEPTION( (true),
                      std::logic_error,
                      ">>> ERROR (Basis_F1_TRI_I8_FEM_FIAT): FEM Basis calling an FV/D member function");
}



template<class Scalar>
int Basis_F1_TRI_I8_FEM_FIAT<Scalar>::getNumLocalDof() const {
    return numDof_;   
}



template<class Scalar>
int Basis_F1_TRI_I8_FEM_FIAT<Scalar>::getLocalDofEnumeration(const LocalDofTag dofTag) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return tagToEnum_[dofTag.tag_[0]][dofTag.tag_[1]][dofTag.tag_[2]];
}



template<class Scalar>
LocalDofTag Basis_F1_TRI_I8_FEM_FIAT<Scalar>::getLocalDofTag(int id) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_[id];
}



template<class Scalar>
const Teuchos::Array<LocalDofTag> & Basis_F1_TRI_I8_FEM_FIAT<Scalar>::getAllLocalDofTags() {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_;
}



template<class Scalar>
ECell Basis_F1_TRI_I8_FEM_FIAT<Scalar>::getCellType() const {
  return CELL_TRI;
}



template<class Scalar>
EBasis Basis_F1_TRI_I8_FEM_FIAT<Scalar>::getBasisType() const {
  return BASIS_FEM_FIAT;
}



template<class Scalar>
ECoordinates Basis_F1_TRI_I8_FEM_FIAT<Scalar>::getCoordinateSystem() const {
  return COORDINATES_CARTESIAN;
}



template<class Scalar>
int Basis_F1_TRI_I8_FEM_FIAT<Scalar>::getDegree() const {
  return 1;
}


}// namespace Intrepid

