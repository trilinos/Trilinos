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

/** \file   Intrepid_F1_TET_I2_FEM_FIATDef.hpp
    \brief  Definition file for FEM basis functions of degree 1 for 1-forms on TET cells.
    \author Created by R. Kirby via the FIAT project
*/

namespace Intrepid {

template<class Scalar>
void Basis_F1_TET_I2_FEM_FIAT<Scalar>::initialize() {

  // Basis-dependent initializations
  int tagSize  = 4;         // size of DoF tag
  int posScDim = 0;         // position in the tag, counting from 0, of the subcell dim 
  int posScId  = 1;         // position in the tag, counting from 0, of the subcell id
  int posBfId  = 2;         // position in the tag, counting from 0, of DoF Id relative to the subcell

  // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 
  int tags[] = { 1 , 0 , 0 , 2 , 1 , 0 , 1 , 2 , 1 , 1 , 0 , 2 , 1 , 1 , 1 , 2 , 1 , 2 , 0 , 2 , 1 , 2 , 1 , 2 , 1 , 3 , 0 , 2 , 1 , 3 , 1 , 2 , 1 , 4 , 0 , 2 , 1 , 4 , 1 , 2 , 1 , 5 , 0 , 2 , 1 , 5 , 1 , 2 , 2 , 0 , 0 , 2 , 2 , 0 , 1 , 2 , 2 , 1 , 0 , 2 , 2 , 1 , 1 , 2 , 2 , 2 , 0 , 2 , 2 , 2 , 1 , 2 , 2 , 3 , 0 , 2 , 2 , 3 , 1 , 2 };
  
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
void Basis_F1_TET_I2_FEM_FIAT<Scalar>::getValues(FieldContainer<Scalar>&                  outputValues,
          const Teuchos::Array< Point<Scalar> >& inputPoints,
          const EOperator                        operatorType) const {

  // Determine parameters to shape outputValues: number of points = size of input array
  int numPoints = inputPoints.size();
	
  // Incomplete polynomial basis of degree 1 (I1) has 20 basis functions on a tetrahedron that are 1-forms in 3D
  int    numFields = 20;
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
                        ">>> ERROR (Basis_F1_TET_I2_FEM_FIAT): Evaluation point is outside the TET reference cell");
  }
#endif
  switch(operatorType) {
    case OPERATOR_VALUE:   
    {
      Teuchos::SerialDenseMatrix<int,Scalar> expansions(10,numPoints);
      Teuchos::SerialDenseMatrix<int,Scalar> result(20,numPoints);
      OrthogonalExpansions<Scalar>::tabulate(CELL_TET,2,inputPoints,expansions);
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm0_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<20;i++) {
          outputValues(countPt,i,0) = result(i,countPt);
        }
      }
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm1_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<20;i++) {
          outputValues(countPt,i,1) = result(i,countPt);
        }
      }
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm2_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<20;i++) {
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
      TEST_FOR_EXCEPTION( true , std::invalid_argument, ">>>ERROR F1_TET_I2_FEM_FIAT: operator not implemented" );
    }
    break; 
    case OPERATOR_CURL:
    {
      Teuchos::SerialDenseMatrix<int,Scalar> expansions(10,numPoints);
      Teuchos::SerialDenseMatrix<int,Scalar> tmp(10,10);
      Teuchos::SerialDenseMatrix<int,Scalar> curlcomp(20,numPoints);
      OrthogonalExpansions<Scalar>::tabulate(CELL_TET,2,inputPoints,expansions);
      curlcomp.putScalar(0.0);
      // x component of curl
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats1_,expansions,0.0);
      curlcomp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm2_,tmp,0.0);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats2_,expansions,0.0);
      curlcomp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,-1.0,*vdm1_,tmp,1.0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<20;i++) {
          outputValues(countPt,i,0) = curlcomp(i,countPt);
        }
      }
      // y component of curl
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats2_,expansions,0.0);
      curlcomp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm0_,tmp,0.0);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats0_,expansions,0.0);
      curlcomp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,-1.0,*vdm2_,tmp,1.0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<20;i++) {
          outputValues(countPt,i,1) = curlcomp(i,countPt);
        }
      }
      // z component of curl
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats0_,expansions,0.0);
      curlcomp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm1_,tmp,0.0);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats1_,expansions,0.0);
      curlcomp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,-1.0,*vdm0_,tmp,1.0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<20;i++) {
          outputValues(countPt,i,2) = curlcomp(i,countPt);
        }
      }
    }
    break;

    case OPERATOR_DIV:
    {
      TEST_FOR_EXCEPTION( true , std::invalid_argument, ">>>ERROR F1_TET_I2_FEM_FIAT: operator not implemented" );
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
                          ">>> ERROR (Basis_F1_TET_I2_FEM_DEFAULT): Invalid operator type");
  }


}


template<class Scalar>
void Basis_F1_TET_I2_FEM_FIAT<Scalar>::getValues(FieldContainer<Scalar>&                  outputValues,
                                                    const Teuchos::Array< Point<Scalar> >& inputPoints,
                                                    const Cell<Scalar>&                    cell) const {
  TEST_FOR_EXCEPTION( (true),
                      std::logic_error,
                      ">>> ERROR (Basis_F1_TET_I2_FEM_FIAT): FEM Basis calling an FV/D member function");
}



template<class Scalar>
int Basis_F1_TET_I2_FEM_FIAT<Scalar>::getNumLocalDof() const {
    return numDof_;   
}



template<class Scalar>
int Basis_F1_TET_I2_FEM_FIAT<Scalar>::getLocalDofEnumeration(const LocalDofTag dofTag) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return tagToEnum_[dofTag.tag_[0]][dofTag.tag_[1]][dofTag.tag_[2]];
}



template<class Scalar>
LocalDofTag Basis_F1_TET_I2_FEM_FIAT<Scalar>::getLocalDofTag(int id) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_[id];
}



template<class Scalar>
const Teuchos::Array<LocalDofTag> & Basis_F1_TET_I2_FEM_FIAT<Scalar>::getAllLocalDofTags() {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_;
}



template<class Scalar>
ECell Basis_F1_TET_I2_FEM_FIAT<Scalar>::getCellType() const {
  return CELL_TET;
}



template<class Scalar>
EBasis Basis_F1_TET_I2_FEM_FIAT<Scalar>::getBasisType() const {
  return BASIS_FEM_FIAT;
}



template<class Scalar>
ECoordinates Basis_F1_TET_I2_FEM_FIAT<Scalar>::getCoordinateSystem() const {
  return COORDINATES_CARTESIAN;
}



template<class Scalar>
int Basis_F1_TET_I2_FEM_FIAT<Scalar>::getDegree() const {
  return 1;
}


}// namespace Intrepid

