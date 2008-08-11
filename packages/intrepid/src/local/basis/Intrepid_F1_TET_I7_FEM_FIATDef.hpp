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

/** \file   Intrepid_F1_TET_I7_FEM_FIATDef.hpp
    \brief  Definition file for FEM basis functions of degree 6 for 1-forms on TET cells.
    \author Created by R. Kirby via the FIAT project
*/

namespace Intrepid {

template<class Scalar>
void Basis_F1_TET_I7_FEM_FIAT<Scalar>::initialize() {

  // Basis-dependent initializations
  int tagSize  = 4;         // size of DoF tag
  int posScDim = 0;         // position in the tag, counting from 0, of the subcell dim 
  int posScId  = 1;         // position in the tag, counting from 0, of the subcell id
  int posBfId  = 2;         // position in the tag, counting from 0, of DoF Id relative to the subcell

  // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 
  int tags[] = { 1 , 0 , 0 , 7 , 1 , 0 , 1 , 7 , 1 , 0 , 2 , 7 , 1 , 0 , 3 , 7 , 1 , 0 , 4 , 7 , 1 , 0 , 5 , 7 , 1 , 0 , 6 , 7 , 1 , 1 , 0 , 7 , 1 , 1 , 1 , 7 , 1 , 1 , 2 , 7 , 1 , 1 , 3 , 7 , 1 , 1 , 4 , 7 , 1 , 1 , 5 , 7 , 1 , 1 , 6 , 7 , 1 , 2 , 0 , 7 , 1 , 2 , 1 , 7 , 1 , 2 , 2 , 7 , 1 , 2 , 3 , 7 , 1 , 2 , 4 , 7 , 1 , 2 , 5 , 7 , 1 , 2 , 6 , 7 , 1 , 3 , 0 , 7 , 1 , 3 , 1 , 7 , 1 , 3 , 2 , 7 , 1 , 3 , 3 , 7 , 1 , 3 , 4 , 7 , 1 , 3 , 5 , 7 , 1 , 3 , 6 , 7 , 1 , 4 , 0 , 7 , 1 , 4 , 1 , 7 , 1 , 4 , 2 , 7 , 1 , 4 , 3 , 7 , 1 , 4 , 4 , 7 , 1 , 4 , 5 , 7 , 1 , 4 , 6 , 7 , 1 , 5 , 0 , 7 , 1 , 5 , 1 , 7 , 1 , 5 , 2 , 7 , 1 , 5 , 3 , 7 , 1 , 5 , 4 , 7 , 1 , 5 , 5 , 7 , 1 , 5 , 6 , 7 , 2 , 0 , 0 , 42 , 2 , 0 , 1 , 42 , 2 , 0 , 2 , 42 , 2 , 0 , 3 , 42 , 2 , 0 , 4 , 42 , 2 , 0 , 5 , 42 , 2 , 0 , 6 , 42 , 2 , 0 , 7 , 42 , 2 , 0 , 8 , 42 , 2 , 0 , 9 , 42 , 2 , 0 , 10 , 42 , 2 , 0 , 11 , 42 , 2 , 0 , 12 , 42 , 2 , 0 , 13 , 42 , 2 , 0 , 14 , 42 , 2 , 0 , 15 , 42 , 2 , 0 , 16 , 42 , 2 , 0 , 17 , 42 , 2 , 0 , 18 , 42 , 2 , 0 , 19 , 42 , 2 , 0 , 20 , 42 , 2 , 0 , 21 , 42 , 2 , 0 , 22 , 42 , 2 , 0 , 23 , 42 , 2 , 0 , 24 , 42 , 2 , 0 , 25 , 42 , 2 , 0 , 26 , 42 , 2 , 0 , 27 , 42 , 2 , 0 , 28 , 42 , 2 , 0 , 29 , 42 , 2 , 0 , 30 , 42 , 2 , 0 , 31 , 42 , 2 , 0 , 32 , 42 , 2 , 0 , 33 , 42 , 2 , 0 , 34 , 42 , 2 , 0 , 35 , 42 , 2 , 0 , 36 , 42 , 2 , 0 , 37 , 42 , 2 , 0 , 38 , 42 , 2 , 0 , 39 , 42 , 2 , 0 , 40 , 42 , 2 , 0 , 41 , 42 , 2 , 1 , 0 , 42 , 2 , 1 , 1 , 42 , 2 , 1 , 2 , 42 , 2 , 1 , 3 , 42 , 2 , 1 , 4 , 42 , 2 , 1 , 5 , 42 , 2 , 1 , 6 , 42 , 2 , 1 , 7 , 42 , 2 , 1 , 8 , 42 , 2 , 1 , 9 , 42 , 2 , 1 , 10 , 42 , 2 , 1 , 11 , 42 , 2 , 1 , 12 , 42 , 2 , 1 , 13 , 42 , 2 , 1 , 14 , 42 , 2 , 1 , 15 , 42 , 2 , 1 , 16 , 42 , 2 , 1 , 17 , 42 , 2 , 1 , 18 , 42 , 2 , 1 , 19 , 42 , 2 , 1 , 20 , 42 , 2 , 1 , 21 , 42 , 2 , 1 , 22 , 42 , 2 , 1 , 23 , 42 , 2 , 1 , 24 , 42 , 2 , 1 , 25 , 42 , 2 , 1 , 26 , 42 , 2 , 1 , 27 , 42 , 2 , 1 , 28 , 42 , 2 , 1 , 29 , 42 , 2 , 1 , 30 , 42 , 2 , 1 , 31 , 42 , 2 , 1 , 32 , 42 , 2 , 1 , 33 , 42 , 2 , 1 , 34 , 42 , 2 , 1 , 35 , 42 , 2 , 1 , 36 , 42 , 2 , 1 , 37 , 42 , 2 , 1 , 38 , 42 , 2 , 1 , 39 , 42 , 2 , 1 , 40 , 42 , 2 , 1 , 41 , 42 , 2 , 2 , 0 , 42 , 2 , 2 , 1 , 42 , 2 , 2 , 2 , 42 , 2 , 2 , 3 , 42 , 2 , 2 , 4 , 42 , 2 , 2 , 5 , 42 , 2 , 2 , 6 , 42 , 2 , 2 , 7 , 42 , 2 , 2 , 8 , 42 , 2 , 2 , 9 , 42 , 2 , 2 , 10 , 42 , 2 , 2 , 11 , 42 , 2 , 2 , 12 , 42 , 2 , 2 , 13 , 42 , 2 , 2 , 14 , 42 , 2 , 2 , 15 , 42 , 2 , 2 , 16 , 42 , 2 , 2 , 17 , 42 , 2 , 2 , 18 , 42 , 2 , 2 , 19 , 42 , 2 , 2 , 20 , 42 , 2 , 2 , 21 , 42 , 2 , 2 , 22 , 42 , 2 , 2 , 23 , 42 , 2 , 2 , 24 , 42 , 2 , 2 , 25 , 42 , 2 , 2 , 26 , 42 , 2 , 2 , 27 , 42 , 2 , 2 , 28 , 42 , 2 , 2 , 29 , 42 , 2 , 2 , 30 , 42 , 2 , 2 , 31 , 42 , 2 , 2 , 32 , 42 , 2 , 2 , 33 , 42 , 2 , 2 , 34 , 42 , 2 , 2 , 35 , 42 , 2 , 2 , 36 , 42 , 2 , 2 , 37 , 42 , 2 , 2 , 38 , 42 , 2 , 2 , 39 , 42 , 2 , 2 , 40 , 42 , 2 , 2 , 41 , 42 , 2 , 3 , 0 , 42 , 2 , 3 , 1 , 42 , 2 , 3 , 2 , 42 , 2 , 3 , 3 , 42 , 2 , 3 , 4 , 42 , 2 , 3 , 5 , 42 , 2 , 3 , 6 , 42 , 2 , 3 , 7 , 42 , 2 , 3 , 8 , 42 , 2 , 3 , 9 , 42 , 2 , 3 , 10 , 42 , 2 , 3 , 11 , 42 , 2 , 3 , 12 , 42 , 2 , 3 , 13 , 42 , 2 , 3 , 14 , 42 , 2 , 3 , 15 , 42 , 2 , 3 , 16 , 42 , 2 , 3 , 17 , 42 , 2 , 3 , 18 , 42 , 2 , 3 , 19 , 42 , 2 , 3 , 20 , 42 , 2 , 3 , 21 , 42 , 2 , 3 , 22 , 42 , 2 , 3 , 23 , 42 , 2 , 3 , 24 , 42 , 2 , 3 , 25 , 42 , 2 , 3 , 26 , 42 , 2 , 3 , 27 , 42 , 2 , 3 , 28 , 42 , 2 , 3 , 29 , 42 , 2 , 3 , 30 , 42 , 2 , 3 , 31 , 42 , 2 , 3 , 32 , 42 , 2 , 3 , 33 , 42 , 2 , 3 , 34 , 42 , 2 , 3 , 35 , 42 , 2 , 3 , 36 , 42 , 2 , 3 , 37 , 42 , 2 , 3 , 38 , 42 , 2 , 3 , 39 , 42 , 2 , 3 , 40 , 42 , 2 , 3 , 41 , 42 , 3 , 0 , 0 , 105 , 3 , 0 , 1 , 105 , 3 , 0 , 2 , 105 , 3 , 0 , 3 , 105 , 3 , 0 , 4 , 105 , 3 , 0 , 5 , 105 , 3 , 0 , 6 , 105 , 3 , 0 , 7 , 105 , 3 , 0 , 8 , 105 , 3 , 0 , 9 , 105 , 3 , 0 , 10 , 105 , 3 , 0 , 11 , 105 , 3 , 0 , 12 , 105 , 3 , 0 , 13 , 105 , 3 , 0 , 14 , 105 , 3 , 0 , 15 , 105 , 3 , 0 , 16 , 105 , 3 , 0 , 17 , 105 , 3 , 0 , 18 , 105 , 3 , 0 , 19 , 105 , 3 , 0 , 20 , 105 , 3 , 0 , 21 , 105 , 3 , 0 , 22 , 105 , 3 , 0 , 23 , 105 , 3 , 0 , 24 , 105 , 3 , 0 , 25 , 105 , 3 , 0 , 26 , 105 , 3 , 0 , 27 , 105 , 3 , 0 , 28 , 105 , 3 , 0 , 29 , 105 , 3 , 0 , 30 , 105 , 3 , 0 , 31 , 105 , 3 , 0 , 32 , 105 , 3 , 0 , 33 , 105 , 3 , 0 , 34 , 105 , 3 , 0 , 35 , 105 , 3 , 0 , 36 , 105 , 3 , 0 , 37 , 105 , 3 , 0 , 38 , 105 , 3 , 0 , 39 , 105 , 3 , 0 , 40 , 105 , 3 , 0 , 41 , 105 , 3 , 0 , 42 , 105 , 3 , 0 , 43 , 105 , 3 , 0 , 44 , 105 , 3 , 0 , 45 , 105 , 3 , 0 , 46 , 105 , 3 , 0 , 47 , 105 , 3 , 0 , 48 , 105 , 3 , 0 , 49 , 105 , 3 , 0 , 50 , 105 , 3 , 0 , 51 , 105 , 3 , 0 , 52 , 105 , 3 , 0 , 53 , 105 , 3 , 0 , 54 , 105 , 3 , 0 , 55 , 105 , 3 , 0 , 56 , 105 , 3 , 0 , 57 , 105 , 3 , 0 , 58 , 105 , 3 , 0 , 59 , 105 , 3 , 0 , 60 , 105 , 3 , 0 , 61 , 105 , 3 , 0 , 62 , 105 , 3 , 0 , 63 , 105 , 3 , 0 , 64 , 105 , 3 , 0 , 65 , 105 , 3 , 0 , 66 , 105 , 3 , 0 , 67 , 105 , 3 , 0 , 68 , 105 , 3 , 0 , 69 , 105 , 3 , 0 , 70 , 105 , 3 , 0 , 71 , 105 , 3 , 0 , 72 , 105 , 3 , 0 , 73 , 105 , 3 , 0 , 74 , 105 , 3 , 0 , 75 , 105 , 3 , 0 , 76 , 105 , 3 , 0 , 77 , 105 , 3 , 0 , 78 , 105 , 3 , 0 , 79 , 105 , 3 , 0 , 80 , 105 , 3 , 0 , 81 , 105 , 3 , 0 , 82 , 105 , 3 , 0 , 83 , 105 , 3 , 0 , 84 , 105 , 3 , 0 , 85 , 105 , 3 , 0 , 86 , 105 , 3 , 0 , 87 , 105 , 3 , 0 , 88 , 105 , 3 , 0 , 89 , 105 , 3 , 0 , 90 , 105 , 3 , 0 , 91 , 105 , 3 , 0 , 92 , 105 , 3 , 0 , 93 , 105 , 3 , 0 , 94 , 105 , 3 , 0 , 95 , 105 , 3 , 0 , 96 , 105 , 3 , 0 , 97 , 105 , 3 , 0 , 98 , 105 , 3 , 0 , 99 , 105 , 3 , 0 , 100 , 105 , 3 , 0 , 101 , 105 , 3 , 0 , 102 , 105 , 3 , 0 , 103 , 105 , 3 , 0 , 104 , 105 };
  
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
void Basis_F1_TET_I7_FEM_FIAT<Scalar>::getValues(FieldContainer<Scalar>&                  outputValues,
          const Teuchos::Array< Point<Scalar> >& inputPoints,
          const EOperator                        operatorType) const {

  // Determine parameters to shape outputValues: number of points = size of input array
  int numPoints = inputPoints.size();
	
  // Incomplete polynomial basis of degree 6 (I6) has 315 basis functions on a tetrahedron that are 1-forms in 3D
  int    numFields = 315;
  EField fieldType = FIELD_FORM_1;
  int    spaceDim  = 3;

  // temporaries
  int countPt  = 0;               // point counter
  Teuchos::Array<int> indexV(3);  // multi-index for values
  Teuchos::Array<int> indexD(3);  // multi-index for curl

  // Shape the FieldContainer for the output values using these values:
  outputValues.resize(numPoints,
                     numFields,
                     fieldType,
                     operatorType,
                     spaceDim);

#ifdef HAVE_INTREPID_DEBUG
  for (countPt=0; countPt<numPoints; countPt++) {
    // Verify that all points are inside the TET reference cell
    TEST_FOR_EXCEPTION( !MultiCell<Scalar>::inReferenceCell(CELL_TET, inputPoints[countPt]),
                        std::invalid_argument,
                        ">>> ERROR (Basis_F1_TET_I7_FEM_FIAT): Evaluation point is outside the TET reference cell");
  }
#endif
  switch(operatorType) {
    case OPERATOR_VALUE:   
    {
      Teuchos::SerialDenseMatrix<int,Scalar> expansions(120,numPoints);
      Teuchos::SerialDenseMatrix<int,Scalar> result(315,numPoints);
      OrthogonalExpansions<Scalar>::tabulate(CELL_TET,7,inputPoints,expansions);
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm0_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<315;i++) {
          outputValues(countPt,i,0) = result(i,countPt);
        }
      }
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm1_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<315;i++) {
          outputValues(countPt,i,1) = result(i,countPt);
        }
      }
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm2_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<315;i++) {
          outputValues(countPt,i,2) = result(i,countPt);
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
      TEST_FOR_EXCEPTION( true , std::invalid_argument, ">>>ERROR F1_TET_I7_FEM_FIAT: operator not implemented" );
    }
    break; 
    case OPERATOR_CURL:
    {
      Teuchos::SerialDenseMatrix<int,Scalar> expansions(120,numPoints);
      Teuchos::SerialDenseMatrix<int,Scalar> tmp(120,120);
      Teuchos::SerialDenseMatrix<int,Scalar> curlcomp(315,numPoints);
      OrthogonalExpansions<Scalar>::tabulate(CELL_TET,7,inputPoints,expansions);
      curlcomp.putScalar(0.0);
      // x component of curl
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats1_,expansions,0.0);
      curlcomp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm2_,tmp,0.0);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats2_,expansions,0.0);
      curlcomp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,-1.0,*vdm1_,tmp,1.0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<315;i++) {
          outputValues(countPt,i,0) = curlcomp(i,countPt);
        }
      }
      // y component of curl
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats2_,expansions,0.0);
      curlcomp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm0_,tmp,0.0);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats0_,expansions,0.0);
      curlcomp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,-1.0,*vdm2_,tmp,1.0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<315;i++) {
          outputValues(countPt,i,1) = curlcomp(i,countPt);
        }
      }
      // z component of curl
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats0_,expansions,0.0);
      curlcomp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm1_,tmp,0.0);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats1_,expansions,0.0);
      curlcomp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,-1.0,*vdm0_,tmp,1.0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<315;i++) {
          outputValues(countPt,i,2) = curlcomp(i,countPt);
        }
      }
    }
    break;

    case OPERATOR_DIV:
    {
      TEST_FOR_EXCEPTION( true , std::invalid_argument, ">>>ERROR F1_TET_I7_FEM_FIAT: operator not implemented" );
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
                          ">>> ERROR (Basis_F1_TET_I7_FEM_DEFAULT): Invalid operator type");
  }


}


template<class Scalar>
void Basis_F1_TET_I7_FEM_FIAT<Scalar>::getValues(FieldContainer<Scalar>&                  outputValues,
                                                    const Teuchos::Array< Point<Scalar> >& inputPoints,
                                                    const Cell<Scalar>&                    cell) const {
  TEST_FOR_EXCEPTION( (true),
                      std::logic_error,
                      ">>> ERROR (Basis_F1_TET_I7_FEM_FIAT): FEM Basis calling an FV/D member function");
}



template<class Scalar>
int Basis_F1_TET_I7_FEM_FIAT<Scalar>::getNumLocalDof() const {
    return numDof_;   
}



template<class Scalar>
int Basis_F1_TET_I7_FEM_FIAT<Scalar>::getLocalDofEnumeration(const LocalDofTag dofTag) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return tagToEnum_[dofTag.tag_[0]][dofTag.tag_[1]][dofTag.tag_[2]];
}



template<class Scalar>
LocalDofTag Basis_F1_TET_I7_FEM_FIAT<Scalar>::getLocalDofTag(int id) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_[id];
}



template<class Scalar>
const Teuchos::Array<LocalDofTag> & Basis_F1_TET_I7_FEM_FIAT<Scalar>::getAllLocalDofTags() {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_;
}



template<class Scalar>
ECell Basis_F1_TET_I7_FEM_FIAT<Scalar>::getCellType() const {
  return CELL_TET;
}



template<class Scalar>
EBasis Basis_F1_TET_I7_FEM_FIAT<Scalar>::getBasisType() const {
  return BASIS_FEM_FIAT;
}



template<class Scalar>
ECoordinates Basis_F1_TET_I7_FEM_FIAT<Scalar>::getCoordinateSystem() const {
  return COORDINATES_CARTESIAN;
}



template<class Scalar>
int Basis_F1_TET_I7_FEM_FIAT<Scalar>::getDegree() const {
  return 1;
}


}// namespace Intrepid

