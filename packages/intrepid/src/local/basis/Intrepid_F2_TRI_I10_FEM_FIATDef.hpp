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

/** \file   Intrepid_F2_TRI_I10_FEM_FIATDef.hpp
    \brief  Definition file for FEM basis functions of degree 9 for 2-forms on TRI cells.
    \author Created by R. Kirby via the FIAT project
*/

namespace Intrepid {

template<class Scalar>
void Basis_F2_TRI_I10_FEM_FIAT<Scalar>::initialize() {

  // Basis-dependent initializations
  int tagSize  = 4;         // size of DoF tag
  int posScDim = 0;         // position in the tag, counting from 0, of the subcell dim 
  int posScId  = 1;         // position in the tag, counting from 0, of the subcell id
  int posBfId  = 2;         // position in the tag, counting from 0, of DoF Id relative to the subcell

  // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 
  int tags[] = { 1 , 0 , 0 , 10 , 1 , 0 , 1 , 10 , 1 , 0 , 2 , 10 , 1 , 0 , 3 , 10 , 1 , 0 , 4 , 10 , 1 , 0 , 5 , 10 , 1 , 0 , 6 , 10 , 1 , 0 , 7 , 10 , 1 , 0 , 8 , 10 , 1 , 0 , 9 , 10 , 1 , 1 , 0 , 10 , 1 , 1 , 1 , 10 , 1 , 1 , 2 , 10 , 1 , 1 , 3 , 10 , 1 , 1 , 4 , 10 , 1 , 1 , 5 , 10 , 1 , 1 , 6 , 10 , 1 , 1 , 7 , 10 , 1 , 1 , 8 , 10 , 1 , 1 , 9 , 10 , 1 , 2 , 0 , 10 , 1 , 2 , 1 , 10 , 1 , 2 , 2 , 10 , 1 , 2 , 3 , 10 , 1 , 2 , 4 , 10 , 1 , 2 , 5 , 10 , 1 , 2 , 6 , 10 , 1 , 2 , 7 , 10 , 1 , 2 , 8 , 10 , 1 , 2 , 9 , 10 , 2 , 0 , 0 , 90 , 2 , 0 , 1 , 90 , 2 , 0 , 2 , 90 , 2 , 0 , 3 , 90 , 2 , 0 , 4 , 90 , 2 , 0 , 5 , 90 , 2 , 0 , 6 , 90 , 2 , 0 , 7 , 90 , 2 , 0 , 8 , 90 , 2 , 0 , 9 , 90 , 2 , 0 , 10 , 90 , 2 , 0 , 11 , 90 , 2 , 0 , 12 , 90 , 2 , 0 , 13 , 90 , 2 , 0 , 14 , 90 , 2 , 0 , 15 , 90 , 2 , 0 , 16 , 90 , 2 , 0 , 17 , 90 , 2 , 0 , 18 , 90 , 2 , 0 , 19 , 90 , 2 , 0 , 20 , 90 , 2 , 0 , 21 , 90 , 2 , 0 , 22 , 90 , 2 , 0 , 23 , 90 , 2 , 0 , 24 , 90 , 2 , 0 , 25 , 90 , 2 , 0 , 26 , 90 , 2 , 0 , 27 , 90 , 2 , 0 , 28 , 90 , 2 , 0 , 29 , 90 , 2 , 0 , 30 , 90 , 2 , 0 , 31 , 90 , 2 , 0 , 32 , 90 , 2 , 0 , 33 , 90 , 2 , 0 , 34 , 90 , 2 , 0 , 35 , 90 , 2 , 0 , 36 , 90 , 2 , 0 , 37 , 90 , 2 , 0 , 38 , 90 , 2 , 0 , 39 , 90 , 2 , 0 , 40 , 90 , 2 , 0 , 41 , 90 , 2 , 0 , 42 , 90 , 2 , 0 , 43 , 90 , 2 , 0 , 44 , 90 , 2 , 0 , 45 , 90 , 2 , 0 , 46 , 90 , 2 , 0 , 47 , 90 , 2 , 0 , 48 , 90 , 2 , 0 , 49 , 90 , 2 , 0 , 50 , 90 , 2 , 0 , 51 , 90 , 2 , 0 , 52 , 90 , 2 , 0 , 53 , 90 , 2 , 0 , 54 , 90 , 2 , 0 , 55 , 90 , 2 , 0 , 56 , 90 , 2 , 0 , 57 , 90 , 2 , 0 , 58 , 90 , 2 , 0 , 59 , 90 , 2 , 0 , 60 , 90 , 2 , 0 , 61 , 90 , 2 , 0 , 62 , 90 , 2 , 0 , 63 , 90 , 2 , 0 , 64 , 90 , 2 , 0 , 65 , 90 , 2 , 0 , 66 , 90 , 2 , 0 , 67 , 90 , 2 , 0 , 68 , 90 , 2 , 0 , 69 , 90 , 2 , 0 , 70 , 90 , 2 , 0 , 71 , 90 , 2 , 0 , 72 , 90 , 2 , 0 , 73 , 90 , 2 , 0 , 74 , 90 , 2 , 0 , 75 , 90 , 2 , 0 , 76 , 90 , 2 , 0 , 77 , 90 , 2 , 0 , 78 , 90 , 2 , 0 , 79 , 90 , 2 , 0 , 80 , 90 , 2 , 0 , 81 , 90 , 2 , 0 , 82 , 90 , 2 , 0 , 83 , 90 , 2 , 0 , 84 , 90 , 2 , 0 , 85 , 90 , 2 , 0 , 86 , 90 , 2 , 0 , 87 , 90 , 2 , 0 , 88 , 90 , 2 , 0 , 89 , 90 };
  
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
void Basis_F2_TRI_I10_FEM_FIAT<Scalar>::getValues(FieldContainer<Scalar>&                  outputValues,
          const Teuchos::Array< Point<Scalar> >& inputPoints,
          const EOperator                        operatorType) const {

  // Determine parameters to shape outputValues: number of points = size of input array
  int numPoints = inputPoints.size();
	
  // Incomplete polynomial basis of degree 9 (I9) has 120 basis functions on a triangle that are 2-forms in 2D
  int    numFields = 120;
  EField fieldType = FIELD_FORM_2;
  int    spaceDim  = 2;

  // temporaries
  int countPt  = 0;               // point counter
  Teuchos::Array<int> indexV(3);  // multi-index for values
  Teuchos::Array<int> indexD(2);  // multi-index for divergences

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
                        ">>> ERROR (Basis_F2_TRI_I10_FEM_FIAT): Evaluation point is outside the TRI reference cell");
  }
#endif
  switch(operatorType) {
    case OPERATOR_VALUE:   
    {
      Teuchos::SerialDenseMatrix<int,Scalar> expansions(66,numPoints);
      Teuchos::SerialDenseMatrix<int,Scalar> result(120,numPoints);
      OrthogonalExpansions<Scalar>::tabulate(CELL_TRI,10,inputPoints,expansions);
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm0_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<120;i++) {
          outputValues(countPt,i,0) = result(i,countPt);
        }
      }
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm1_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<120;i++) {
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
      TEST_FOR_EXCEPTION( true , std::invalid_argument, ">>>ERROR F2_TRI_I10_FEM_FIAT: operator not implemented" );
    }
    break; 
    case OPERATOR_CURL:
    {
      TEST_FOR_EXCEPTION( true , std::invalid_argument, ">>>ERROR F2_TRI_I10_FEM_FIAT: operator not implemented" );
    }
    break; 
    case OPERATOR_DIV:
    {
      Teuchos::SerialDenseMatrix<int,Scalar> expansions(66,numPoints);
      Teuchos::SerialDenseMatrix<int,Scalar> tmp(66,66);
      Teuchos::SerialDenseMatrix<int,Scalar> div(120,numPoints);
      div.putScalar(0.0);
      OrthogonalExpansions<Scalar>::tabulate(CELL_TRI,10,inputPoints,expansions);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats0_,expansions,0.0);
      div.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm0_,tmp,1.0);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats1_,expansions,0.0);
      div.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm1_,tmp,1.0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<120;i++) {
          outputValues(countPt,i) = div(i,countPt);
        }
      }
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
                          ">>> ERROR (Basis_F2_TRI_I10_FEM_DEFAULT): Invalid operator type");
  }


}


template<class Scalar>
void Basis_F2_TRI_I10_FEM_FIAT<Scalar>::getValues(FieldContainer<Scalar>&                  outputValues,
                                                    const Teuchos::Array< Point<Scalar> >& inputPoints,
                                                    const Cell<Scalar>&                    cell) const {
  TEST_FOR_EXCEPTION( (true),
                      std::logic_error,
                      ">>> ERROR (Basis_F2_TRI_I10_FEM_FIAT): FEM Basis calling an FV/D member function");
}



template<class Scalar>
int Basis_F2_TRI_I10_FEM_FIAT<Scalar>::getNumLocalDof() const {
    return numDof_;   
}



template<class Scalar>
int Basis_F2_TRI_I10_FEM_FIAT<Scalar>::getLocalDofEnumeration(const LocalDofTag dofTag) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return tagToEnum_[dofTag.tag_[0]][dofTag.tag_[1]][dofTag.tag_[2]];
}



template<class Scalar>
LocalDofTag Basis_F2_TRI_I10_FEM_FIAT<Scalar>::getLocalDofTag(int id) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_[id];
}



template<class Scalar>
const Teuchos::Array<LocalDofTag> & Basis_F2_TRI_I10_FEM_FIAT<Scalar>::getAllLocalDofTags() {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_;
}



template<class Scalar>
ECell Basis_F2_TRI_I10_FEM_FIAT<Scalar>::getCellType() const {
  return CELL_TRI;
}



template<class Scalar>
EBasis Basis_F2_TRI_I10_FEM_FIAT<Scalar>::getBasisType() const {
  return BASIS_FEM_FIAT;
}



template<class Scalar>
ECoordinates Basis_F2_TRI_I10_FEM_FIAT<Scalar>::getCoordinateSystem() const {
  return COORDINATES_CARTESIAN;
}



template<class Scalar>
int Basis_F2_TRI_I10_FEM_FIAT<Scalar>::getDegree() const {
  return 1;
}


}// namespace Intrepid

