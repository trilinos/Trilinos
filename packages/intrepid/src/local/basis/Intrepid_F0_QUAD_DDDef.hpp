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
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_F0_QUAD_DDDef.hpp
    \brief  Definition file for FEM basis functions on quads using divided difference
            representations of Lagrange bases.
    \author Created by R. Kirby.
*/

namespace Intrepid {


template<class Scalar>
void Basis_F0_QUAD_DD<Scalar>::initialize() {
  // initialize the Lagrange object
  vector<Scalar> equispacedPoints( degree_ + 1 );
  Lagrange::equispacedPoints( degree_ , -1.0 , 1.0 , equispacedPoints );
  poly_ = Teuchos::rcp( new Lagrange::Lagrange<Scalar>( equispacedPoints ) );

  // Basis-dependent intializations
  int tagSize  = 4;         // size of DoF tag
  int posScDim = 0;         // poisition in the tag, counting from 0, of the subcell dim 
  int posScId  = 1;         // poisition in the tag, counting from 0, of the subcell id
  int posBfId  = 2;         // position in the tag, counting from 0, of DoF Id relative to the subcell

  // An array with local DoF tags assigned to basis functions, in the order of their local enumeration 
  // I need to fill this!
  int *tags = new int[tagSize*numDof_];

  // order dof lexicographically from left to right (inner) and bottom to top (outer).
  // this means I have to get the boundary right for what subcell goes to what

  // get the bottom
  // bottom left is on vertex 0 (subcdim = 0, subcid = 0, dofid = 0)

  int *tagsCur = tags;

  tagsCur[0] = 0;
  tagsCur[1] = 0;
  tagsCur[2] = 0;
  tagsCur[3] = 1;

  tagsCur += tagSize;


  // get interior nodes, leftmost edge (edge 3)
  for (int j=1;j<degree_;j++) {
    tagsCur[0] = 1;
    tagsCur[1] = 3;
    tagsCur[2] = j-1;
    tagsCur[3] = degree_ - 1;
    tagsCur += tagSize;
  }

  // get top left vertex
  tagsCur[0] = 0;
  tagsCur[1] = 3;
  tagsCur[2] = 0;
  tagsCur[3] = 1;
  tagsCur += tagSize;

  // now loop over interior rows
  for (int i=1;i<degree_;i++) {
    // first dof is on edge 0
    tagsCur[0] = 1;
    tagsCur[1] = 0;
    tagsCur[2] = i - 1;
    tagsCur[3] = degree_ - 1;
    tagsCur += tagSize;
        
    // interior dof 
    for (int j=1;j<degree_;j++) {
      tagsCur[0] = 2;
      tagsCur[1] = 0;
      tagsCur[2] = (i-1)*(degree_-1) + (j-1);
      tagsCur[3] = (degree_ - 1) * (degree_ - 1);
      tagsCur += tagSize;
    }

    // last dof is on edge 1
    tagsCur[0] = 1;
    tagsCur[1] = 2;
    tagsCur[2] = i-1;
    tagsCur[3] = degree_ - i;
    tagsCur += tagSize;
  
  }

  // now get right row.
  // bottom right vertex number 1
  tagsCur[0] = 0;
  tagsCur[1] = 1;
  tagsCur[2] = 0;
  tagsCur[3] = 1;
  tagsCur += tagSize;

  // interior of right edge (#1)
  for (int j=1;j<degree_;j++) {
    tagsCur[0] = 1;
    tagsCur[1] = 1;
    tagsCur[2] = j - 1;
    tagsCur[3] = degree_ - 1;
    tagsCur += tagSize;
  }

  // top right vertex is vertex number 2
  tagsCur[0] = 0;
  tagsCur[1] = 2;
  tagsCur[2] = 0;
  tagsCur[3] = 1;


  // Basis-independent function sets tag and enum data in the static arrays:
  setEnumTagData(tagToEnum_,
                 enumToTag_,
                 tags,
                 numDof_,
                 tagSize,
                 posScDim,
                 posScId,
                 posBfId);

  delete []tags;
}



template<class Scalar> 
void Basis_F0_QUAD_DD<Scalar>::getValues(FieldContainer<Scalar>&               outputValues,
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
  
  switch (operatorType) {
    case OPERATOR_VALUE:
      for (countPt=0; countPt<numPoints; countPt++) {
#ifdef HAVE_INTREPID_DEBUG 
        // Verify argument: check if input point is inside the reference QUAD
        TEST_FOR_EXCEPTION( !MultiCell<Scalar>::inReferenceCell(CELL_QUAD, inputPoints[countPt]),
                            std::invalid_argument,
                            ">>> ERROR (Basis_F0_QUAD_DD): Evaluation point is outside the QUAD reference cell");
#endif
        Scalar x((inputPoints[countPt])[0]);
        Scalar y((inputPoints[countPt])[1]);
        
        // Output container has rank 2. The indices are (P,F)
        int bfCur = 0;
        for (int i=0;i<degree_+1;i++) {
          for (int j=0;j<degree_+1;j++) {
	          outputValues(countPt, bfCur) = poly_->eval(i,x) * poly_->eval(j,y);
	          bfCur++;
	      }
	    }
      }
      break;
      
    case OPERATOR_GRAD:
    case OPERATOR_D1:
      for (countPt=0; countPt<numPoints; countPt++) {
#ifdef HAVE_INTREPID_DEBUG 
        // Verify argument: check if input point is inside the reference QUAD
        TEST_FOR_EXCEPTION( !MultiCell<Scalar>::inReferenceCell(CELL_QUAD, inputPoints[countPt]),
                            std::invalid_argument,
                            ">>> ERROR (Basis_F0_QUAD_DD): Evaluation point is outside the QUAD reference cell");
#endif
	    // need to turn these into AD types as needed
        
        // Output container has rank 3. The indices are (P,F,D)
        int bfCur = 0;
        for (int i=0;i<degree_+1;i++) {
          for (int j=0;j<degree_+1;j++) {
            Scalar x((inputPoints[countPt])[0]);
            Scalar y((inputPoints[countPt])[1]);

            Sacado::Fad::SFad<Scalar,1> xfad(1,0,x);
            Sacado::Fad::SFad<Scalar,1> yfad(1,0,y);
            Sacado::Fad::SFad<Scalar,1> f1_x;
            Sacado::Fad::SFad<Scalar,1> f2_y;
            f1_x = poly_->eval( i , xfad );
            f2_y = poly_->eval( j , yfad );
            outputValues(countPt,bfCur,0) = f1_x.dx(0) * f2_y.val(); 
            outputValues(countPt,bfCur,1) = f1_x.val() * f2_y.dx(0);
            bfCur++;
          }
        }
      }
      break;
      
    case OPERATOR_CURL:
    case OPERATOR_DIV:
    case OPERATOR_D2:
    case OPERATOR_D3:
    case OPERATOR_D4:
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10:
      TEST_FOR_EXCEPTION( (true) ,
			  std::invalid_argument ,
			  ">>> ERROR (Basis_F0_QUAD_DD: Operator not implemented" );
      
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
                          ">>> ERROR (Basis_F0_QUAD_DD): Invalid operator type");
  }
}


  
template<class Scalar>
void Basis_F0_QUAD_DD<Scalar>::getValues(FieldContainer<Scalar>&                  outputValues,
                                                    const Teuchos::Array< Point<Scalar> >& inputPoints,
                                                    const Cell<Scalar>&                    cell) const {
  TEST_FOR_EXCEPTION( (true),
                      std::logic_error,
                      ">>> ERROR (Basis_F0_QUAD_DD): FEM Basis calling an FVD member function");
}



template<class Scalar>
int Basis_F0_QUAD_DD<Scalar>::getLocalDofEnumeration(const LocalDofTag dofTag) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return tagToEnum_[dofTag.tag_[0]][dofTag.tag_[1]][dofTag.tag_[2]];
}



template<class Scalar>
LocalDofTag Basis_F0_QUAD_DD<Scalar>::getLocalDofTag(int dofId) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_[dofId];
}



template<class Scalar>
const Teuchos::Array<LocalDofTag> & Basis_F0_QUAD_DD<Scalar>::getAllLocalDofTags() {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_;
}



template<class Scalar>
inline ECell Basis_F0_QUAD_DD<Scalar>::getCellType() const {
  return CELL_QUAD;
}



template<class Scalar>
inline EBasis Basis_F0_QUAD_DD<Scalar>::getBasisType() const {
  return BASIS_FEM_DEFAULT;
}

template<class Scalar>
inline ECoordinates Basis_F0_QUAD_DD<Scalar>::getCoordinateSystem() const {
  return COORDINATES_CARTESIAN;
}

}// namespace Intrepid
