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

/** \file   Intrepid_F0_HEX_DDDef.hpp
    \brief  Definition file for FEM basis functions on quads using divided difference
            representations of Lagrange bases.
    \author Created by R. Kirby.
*/

namespace Intrepid {


template<class Scalar>
void Basis_F0_HEX_DD<Scalar>::initialize() {
  // initialize the Lagrange object
  vector<Scalar> equispacedPoints( degree_ + 1 );
  Lagrange::equispacedPoints( degree_ , -1.0 , 1.0 , equispacedPoints );
  poly_ = Teuchos::rcp( new Lagrange::Lagrange<Scalar>( equispacedPoints ) );

  // Basis-dependent intializations
  int tagSize  = 4;         // size of DoF tag
  int posScDim = 0;         // poisition in the tag, counting from 0, of the subcell dim 
  int posScId  = 1;         // poisition in the tag, counting from 0, of the subcell id
  int posBfId  = 2;         // position in the tag, counting from 0, of DoF Id relative to the subcell


  int degm1sq = (degree_-1)*(degree_-1);

  // An array with local DoF tags assigned to basis functions, in the order of their local enumeration 
  // I need to fill this!
  int *tags = new int[tagSize*numDof_];

  int *tagsCur = tags;

  // indices: i runs along points on an edge
  //          j runs over edges in a face
  //          k runs over faces in the cube

  // start at (-1,-1,-1), which is dimension 0, facet 0, dof 0, total dof on vertex, 1
  // it is also the case of i == j == k == 0;

  tagsCur[0] = 0;
  tagsCur[1] = 0;
  tagsCur[2] = 0;
  tagsCur[3] = 1;
  tagsCur += tagSize;

  // now vary k (going up in the z direction).
  // I stay on the edge 4, running from vertex 0 to vertex 4
  for (int k=1;k<degree_;k++) {
    tagsCur[0] = 1;
    tagsCur[1] = 4;
    tagsCur[2] = k-1;
    tagsCur[3] = degree_-1;
    tagsCur += tagSize;
  }

  // now get vertex (-1,-1,1), which is dimension 0, facet 4, dof 0, dof on cell: 1
  tagsCur[0] = 0;
  tagsCur[1] = 4;
  tagsCur[2] = 0;
  tagsCur[3] = 1;
  tagsCur += tagSize;


  // now keep x=-1, increment y 
  for (int j=1;j<degree_;j++) {
    // first node is down on bottom edge running from (-1,-1,-1) (v0) to (-1,1,-1) (v3)
    // this is edge 3
    tagsCur[0] = 1;
    tagsCur[1] = 3;
    tagsCur[2] = j;
    tagsCur[3] = degree_-1;
    tagsCur += tagSize;
    
    // internal nodes to this line will lie on face 3.
    for (int k=1;k<degree_;k++) {
      tagsCur[0] = 2;
      tagsCur[1] = 3;
      tagsCur[2] = (j-1)*(degree_-1)+k-1;
      tagsCur[3] = degm1sq;
      tagsCur += tagSize;
    }

    // last node is on top edge running from (-1,-1,1) (v4) to (-1,1,1) (v7)
    // this is edge 11
    tagsCur[0] = 1;
    tagsCur[1] = 11;
    tagsCur[2] = j-1;
    tagsCur[3] = degree_-1;
    tagsCur += tagSize;
  }

  // now go up line (-1,1,z)
  // first vertex is (-1,1,-1), vertex 3
  tagsCur[0] = 0;
  tagsCur[1] = 3;
  tagsCur[2] = 0;
  tagsCur[3] = 1;
  tagsCur += tagSize;

  // middle of line, which is edge 7
  for (int k=1;k<degree_;k++) {
    tagsCur[0] = 1;
    tagsCur[1] = 7;
    tagsCur[2] = k-1;
    tagsCur[3] = degree_-1;
    tagsCur += tagSize;
  }

  // top of line is vertex 7
  tagsCur[0] = 0;
  tagsCur[1] = 7;
  tagsCur[2] = 0;
  tagsCur[3] = 1;
  tagsCur += tagSize;

  // now do middle planes for each fixed i (x value)
  for (int i=1;i<degree_;i++) {
    // j == 0 \\ y == -1 and 

    // k == 0 // z == -1, live on e0
    tagsCur[0] = 1;
    tagsCur[1] = 0;
    tagsCur[2] = i-1;
    tagsCur[3] = degree_-1;
    tagsCur += tagSize;

    // now go up the z axis (increase k along face 0)
    for (int k=1;k<degree_;k++) {
      tagsCur[0] = 2;
      tagsCur[1] = 0;
      tagsCur[2] = (i-1)*(degree_-1)+(k-1);
      tagsCur[3] = degm1sq;
      tagsCur += tagSize;
    }

    // now z == 1, edge 8
    tagsCur[0] = 1;
    tagsCur[1] = 8;
    tagsCur[2] = i-1;
    tagsCur[3] = degree_-1;
    tagsCur += tagSize;

    for (int j=1;j<degree_;j++) {
      // k == 0, z == -1, live on face 4
      tagsCur[0] = 2;
      tagsCur[1] = 4;
      tagsCur[2] = (i-1)*(degree_-1)+(j-1);
      tagsCur[3] = degm1sq;
      tagsCur += tagSize;

      // loop over interior k
      for (int k=1;k<degree_;k++) {
	tagsCur[0] = 3;
	tagsCur[1] = 0;
	tagsCur[2] = (i-1)*(degm1sq) + (j-1)*(degree_-1) + (k-1);
	tagsCur[3] = degm1sq * (degree_-1);
	tagsCur += tagSize;
      }

      // k == degree_-1, z == 1, live on face 5
      tagsCur[0] = 2;
      tagsCur[1] = 5;
      tagsCur[2] = (i-1)*(degree_-1)+(j-1);
      tagsCur[3] = degm1sq;
      tagsCur += tagSize;
    }

    // now have to handle j == degree_-1, y == 1.
    // k == 0 // z == -1, live on e2
    tagsCur[0] = 1;
    tagsCur[1] = 2;
    tagsCur[2] = i-1;
    tagsCur[3] = degree_-1;
    tagsCur += tagSize;

    // now go up the z axis (increase k along face 2)
    for (int k=1;k<degree_;k++) {
      tagsCur[0] = 2;
      tagsCur[1] = 2;
      tagsCur[2] = (i-1)*(degree_-1)+(k-1);
      tagsCur[3] = degm1sq;
      tagsCur += tagSize;
    }

    // now z == 1, edge 10
    tagsCur[0] = 1;
    tagsCur[1] = 10;
    tagsCur[2] = i-1;
    tagsCur[3] = degree_-1;
    tagsCur += tagSize;    

  }

  // now do final i (x=1 plane)
  // first fix y==-1 (j=0)
  // this is vertex 1
  tagsCur[0] = 0;
  tagsCur[1] = 1;
  tagsCur[2] = 0;
  tagsCur[3] = 1;
  tagsCur += tagSize;

  // edge dof along edge 5
  for (int k=1;k<degree_;k++) {
    tagsCur[0] = 1;
    tagsCur[1] = 5;
    tagsCur[2] = k-1;
    tagsCur[3] = degree_-1;
    tagsCur += tagSize;
  }

  // vertex 5
  tagsCur[0] = 0;
  tagsCur[1] = 5;
  tagsCur[2] = 0;
  tagsCur[3] = 1;
  tagsCur += tagSize;

  // interior lines of x=1 plane
  for (int j=1;j<degree_;j++) {
    // first vertex is on bottom (z=-1), edge 1
    tagsCur[0] = 1;
    tagsCur[1] = 1;
    tagsCur[2] = j-1;
    tagsCur[3] = degree_-1;
    tagsCur += tagSize;

    // go up the line in interior of face 1
    for (int k=1;k<degree_;k++) {
      tagsCur[0] = 2;
      tagsCur[1] = 1;
      tagsCur[2] = (j-1)*(degree_-1)+(k-1);
      tagsCur[3] = degm1sq;
      tagsCur += tagSize;
    }

    // last vertex is on top (z=1), edge 9
    tagsCur[0] = 1;
    tagsCur[1] = 9;
    tagsCur[2] = j-1;
    tagsCur[3] = degree_-1;
    tagsCur += tagSize;
  }

  // edge 6 runs from v2 to v6 (1,1,-1) to (1,1,1)
  tagsCur[0] = 0;
  tagsCur[1] = 2;
  tagsCur[2] = 0;
  tagsCur[3] = 1;
  tagsCur += tagSize;

  // go up edge 6
  for (int k=1;k<degree_;k++) {
    tagsCur[0] = 1;
    tagsCur[1] = 6;
    tagsCur[2] = k-1;
    tagsCur[3] = degree_-1;
    tagsCur += 1;
  }

  // vertex 6
  tagsCur[0] = 0;
  tagsCur[1] = 6;
  tagsCur[2] = 0;
  tagsCur[3] = 1;
  tagsCur += tagSize;


//   // indices: k runs along points on an edge
//   //          j runs over edges in a face
//   //          i runs over faces in the cube

//   // start bottom face k == 0

//   // start at (-1,-1,-1), which is dimension 0, facet 0, dof 0, total dof on vertex, 1
//   // it is also the case of i == j == k == 0;

//   tagsCur[0] = 0;
//   tagsCur[1] = 0;
//   tagsCur[2] = 0;
//   tagsCur[3] = 1;
//   tagsCur += tagSize;

//   // run along bottom edge from (-1,-1,-1) to (1,-1,-1), which is dimension 1, facet 0, degree_-1 dof on edge
  
//   for (int k=1;k<degree_;k++) {
//     tagsCur[0] = 1;
//     tagsCur[1] = 0;
//     tagsCur[2] = k-1;
//     tagsCur[3] = degree_-1;
//     tagsCur += tagSize;
//   }

//   // now get vertex (1,-1,-1), which is dimension 0, facet 1, dof 0, dof on cell: 1
//   tagsCur[0] = 0;
//   tagsCur[1] = 1;
//   tagsCur[2] = 0;
//   tagsCur[3] = 1;
//   tagsCur += tagSize;

//   // now get middle rows on bottom face
//   for (int j=1;j<degree_;j++) {
//     // first node is on the edge from vertex 3 to vertex 0 (edge 3)
//     tagsCur[0] = 1;
//     tagsCur[1] = 3;
//     tagsCur[2] = j - 1;
//     tagsCur[3] = (degree_-1);
//     tagsCur += tagSize;

//     // now do the interior vertices on the i:th row of the bottom face, which is dim 2, facet 4
//     //
//     for (int k=1;k<degree_;k++) {
//       tagsCur[0] = 2;
//       tagsCur[1] = 4;
//       tagsCur[2] = (j-1)*(degree_-1) + (k-1);
//       tagsCur[3] = (degree_-1)*(degree_-1);
//       tagsCur += tagSize;
//     }

//     // the last node is on edge (1,2) (facet 1 of dimension 0), traversed forward
//     tagsCur[0] = 1;
//     tagsCur[1] = 1;
//     tagsCur[2] = j - 1;
//     tagsCur[3] = (degree_-1);
//     tagsCur += tagSize;
//   }

//   // now get vertex 3 (beginning of last row on bottom face)
//   tagsCur[0] = 0;  // dim 0
//   tagsCur[1] = 3;  // id: 3
//   tagsCur[2] = 0;  // first dof
//   tagsCur[3] = 1;  // one dof on the vertex
//   tagsCur += tagSize;

//   // get last edge on bottom face  This is edge 2, running from vertex 2 to vertex 3.

//   for (int k=1;k<degree_;k++) {
//     tagsCur[0] = 1; // dim 1
//     tagsCur[1] = 2; // id: 2
//     tagsCur[2] = k - 1; // dof on edge
//     tagsCur[3] = degree_ - 1; // total number of nodes on edge
//     tagsCur += tagSize;
//   }

//   // now get vertex 2 (end of last row on bottom face)
//   tagsCur[0] = 0; // dim 0
//   tagsCur[1] = 2; // id: 2
//   tagsCur[2] = 0; // first dof on vertex
//   tagsCur[3] = 1; // one dof on the vertex
//   tagsCur += tagSize;

//   // end bottom face

//   // begin interior faces
//   for (int i=1;i<degree_;i++) {
//     // j == 0 lives on face y==-1, which is face 0

//     // k == 0, node lives on edge 4
//     tagsCur[0] = 1;
//     tagsCur[1] = 4;
//     tagsCur[2] = i-1;
//     tagsCur[3] = (degree_-1);
//     tagsCur += tagSize;

//     // interior nodes
//     for (int k=1;k<degree_;k++) {
//       tagsCur[0] = 2;
//       tagsCur[1] = 0;
//       tagsCur[2] = (i-1)*(degree_-1) + (k-1);
//       tagsCur[3] = (degree_-1)*(degree_-1);
//       tagsCur += tagSize;
//     }    

//     // k == degree, node lives on edge 5
//     tagsCur[0] = 1;
//     tagsCur[1] = 5;
//     tagsCur[2] = i-1;
//     tagsCur[3] = (degree_-1);
//     tagsCur += tagSize;

//     // interior edges in i:th face
//     for (int j=1;j<degree_;j++) { 
//       // k == 0, node lives on face 3
//       tagsCur[0] = 2;
//       tagsCur[1] = 3;
//       tagsCur[2] = (i-1)*(degree_-1) + (j-1);
//       tagsCur[3] = (degree_-1)*(degree_-1);
//       tagsCur += tagSize;

//       // interior nodes on j:th edge
//       for (int k=1;k<degree_;k++) {
// 	tagsCur[0] = 3;
// 	tagsCur[1] = 0;
// 	tagsCur[2] = (i-1)*(degree_-1)*(degree_-1) + (j-1)*(degree_-1) + (k-1);
// 	tagsCur[3] = (degree_-1)*(degree_-1)*(degree_-1);
// 	tagsCur += tagSize;
//       }

//       // k == degree, node lives on face 1
//       tagsCur[0] = 2;
//       tagsCur[1] = 1;
//       tagsCur[2] = (i-1)*(degree_-1)+(j-1);
//       tagsCur[3] = (degree_-1)*(degree_-1);
//       tagsCur += tagSize;
      
//     }
//     // end interior j

//     // j == degree lives on face y==1, which is number 2
    
//     // first node (k==0), lives on edge 7
//     tagsCur[0] = 1;
//     tagsCur[1] = 7;
//     tagsCur[2] = i-1;
//     tagsCur[3] = (degree_-1);
//     tagsCur += tagSize;

//     // interior nodes live on face 2
//     for (int k=1;k<degree_;k++) {
//       tagsCur[0] = 2;
//       tagsCur[1] = 2;
//       tagsCur[2] = (i-1)*(degree_-1) + (k-1);
//       tagsCur[3] = (degree_-1)*(degree_-1);
//       tagsCur += tagSize;
//     }

//     // last node (k==degree) lives on edge 6
//     tagsCur[0] = 1;
//     tagsCur[1] = 6;
//     tagsCur[2] = i-1;
//     tagsCur[3] = (degree_-1);
//     tagsCur += tagSize;

//     // end j
//   }  
//   // end interior faces

//   // begin top face
//   // start with vertex 4
//   tagsCur[0] = 0;
//   tagsCur[1] = 4;
//   tagsCur[2] = 0;
//   tagsCur[3] = 1;
//   tagsCur += tagSize;

//   // now do edge 8
//   for (int k=1;k<degree_;k++) {
//     tagsCur[0] = 1;
//     tagsCur[1] = 8;
//     tagsCur[2] = (k-1);
//     tagsCur[3] = (degree_-1);
//     tagsCur += tagSize;
//   }
  
//   // now do vertex 5
//   tagsCur[0] = 0;
//   tagsCur[1] = 5;
//   tagsCur[2] = 0;
//   tagsCur[3] = 1;
//   tagsCur += tagSize;

//   // now do interior edges of top face
//   for (int j=1;j<degree_;j++) {
//     // k==0, lives on edge 11
//     tagsCur[0] = 1;
//     tagsCur[1] = 11;
//     tagsCur[2] = j-1;
//     tagsCur[3] = (degree_-1);
//     tagsCur += tagSize;

//     // interior nodes of edges
//     for (int k=1;k<degree_;k++) {
//       tagsCur[0] = 2;
//       tagsCur[1] = 5;
//       tagsCur[2] = (j-1)*(degree_-1) + (k-1);
//       tagsCur[3] = (degree_-1)*(degree_-1);
//       tagsCur += tagSize;
//     }
    
//     // k == degree, node lives on edge 9
//     tagsCur[0] = 1;
//     tagsCur[1] = 9;
//     tagsCur[2] = j-1;
//     tagsCur[3] = (degree_-1);
//     tagsCur += tagSize;
//   }
//   // done with interior edges of top face.
//   // do last edge of top face
//   // vertex 7
//   tagsCur[0] = 0;
//   tagsCur[1] = 7;
//   tagsCur[2] = 0;
//   tagsCur[3] = 1;
//   tagsCur += tagSize;

//   // interior nodes on edge 10
//   for (int k=1;k<degree_;k++) {
//     tagsCur[0] = 1;
//     tagsCur[1] = 10;
//     tagsCur[2] = k-1;
//     tagsCur[3] = (degree_-1);
//     tagsCur += tagSize;
//   }

//   // finally do verex 6
//   tagsCur[0] = 0;
//   tagsCur[1] = 6;
//   tagsCur[2] = 0;
//   tagsCur[3] = 1;
  
//   // end top face (finally!)


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
void Basis_F0_HEX_DD<Scalar>::getValues(FieldContainer<Scalar>&               outputValues,
					 const Teuchos::Array< Point<Scalar> >& inputPoints,
					 const EOperator                        operatorType) const {

  // Number of evaluation points =  size of outputValues
  int numPoints = inputPoints.size();       
  
  // Shape the FieldContainer for the output:
  outputValues.resize(numPoints,                   // number of evaluation points
                      numDof_,                     // number of fields = number of DoFs in the basis
                      FIELD_FORM_0,                // field type of the basis functions
                      operatorType,                // operator type that is applied to basis functions
                      3);                          // space dimension for HEX cell is 3
  
  // Temporaries: point counter and (x,y) coordinates of the evaluation point
  int countPt  = 0;                               
  
  switch (operatorType) {
  case OPERATOR_VALUE:
    for (countPt=0; countPt<numPoints; countPt++) {
#ifdef HAVE_INTREPID_DEBUG 
      // Verify argument: check if input point is inside the reference HEX
      TEST_FOR_EXCEPTION( !MultiCell<Scalar>::inReferenceCell(CELL_HEX, inputPoints[countPt]),
			  std::invalid_argument,
			  ">>> ERROR (Basis_F0_HEX_DD): Evaluation point is outside the HEX reference cell");
#endif
      Scalar x((inputPoints[countPt])[0]);
      Scalar y((inputPoints[countPt])[1]);
      Scalar z((inputPoints[countPt])[2]);
      
      // Output container has rank 2. The indices are (P,F)
      int bfCur = 0;
      for (int i=0;i<degree_+1;i++) {
	for (int j=0;j<degree_+1;j++) {
	  for (int k=0;k<degree_+1;k++) {
	    outputValues(countPt, bfCur) = poly_->eval(i,x) * poly_->eval(j,y) * poly_->eval(k,z);
	    bfCur++;
	  }
	}
      }
    }
    break;
      
  case OPERATOR_GRAD:
  case OPERATOR_D1:
    for (countPt=0; countPt<numPoints; countPt++) {
#ifdef HAVE_INTREPID_DEBUG 
      // Verify argument: check if input point is inside the reference HEX
      TEST_FOR_EXCEPTION( !MultiCell<Scalar>::inReferenceCell(CELL_HEX, inputPoints[countPt]),
			  std::invalid_argument,
			  ">>> ERROR (Basis_F0_HEX_DD): Evaluation point is outside the HEX reference cell");
#endif
      // need to turn these into AD types as needed
      
      // Output container has rank 3. The indices are (P,F,D)
      int bfCur = 0;
      for (int i=0;i<degree_+1;i++) {
	for (int j=0;j<degree_+1;j++) {
	  for (int k=0;k<degree_+1;k++) {
	    Scalar x((inputPoints[countPt])[0]);
	    Scalar y((inputPoints[countPt])[1]);
	    Scalar z((inputPoints[countPt])[2]);
	    Sacado::Fad::SFad<Scalar,1> xfad(1,0,x);
	
	    Sacado::Fad::SFad<Scalar,1> yfad(1,0,y);
	    Sacado::Fad::SFad<Scalar,1> zfad(1,0,z);
	    Sacado::Fad::SFad<Scalar,1> f1_x;
	    Sacado::Fad::SFad<Scalar,1> f2_y;
	    Sacado::Fad::SFad<Scalar,1> f3_z;
	    f1_x = poly_->eval( i , xfad );
	    f2_y = poly_->eval( j , yfad );
	    f3_z = poly_->eval( k , zfad );
	    outputValues(countPt,bfCur,0) = f1_x.dx(0) * f2_y.val() * f3_z.val(); 
	    outputValues(countPt,bfCur,1) = f1_x.val() * f2_y.dx(0) * f3_z.val();
	    outputValues(countPt,bfCur,2) = f1_x.val() * f2_y.val() * f3_z.dx(0);
	    bfCur++;
	  }
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
			">>> ERROR (Basis_F0_HEX_DD: Operator not implemented" );
    
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
			">>> ERROR (Basis_F0_HEX_DD): Invalid operator type");
  }
}
  

  
template<class Scalar>
void Basis_F0_HEX_DD<Scalar>::getValues(FieldContainer<Scalar>&                  outputValues,
                                                    const Teuchos::Array< Point<Scalar> >& inputPoints,
                                                    const Cell<Scalar>&                    cell) const {
  TEST_FOR_EXCEPTION( (true),
                      std::logic_error,
                      ">>> ERROR (Basis_F0_HEX_DD): FEM Basis calling an FVD member function");
}



template<class Scalar>
int Basis_F0_HEX_DD<Scalar>::getLocalDofEnumeration(const LocalDofTag dofTag) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return tagToEnum_[dofTag.tag_[0]][dofTag.tag_[1]][dofTag.tag_[2]];
}



template<class Scalar>
LocalDofTag Basis_F0_HEX_DD<Scalar>::getLocalDofTag(int dofId) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_[dofId];
}



template<class Scalar>
const Teuchos::Array<LocalDofTag> & Basis_F0_HEX_DD<Scalar>::getAllLocalDofTags() {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_;
}



template<class Scalar>
inline ECell Basis_F0_HEX_DD<Scalar>::getCellType() const {
  return CELL_HEX;
}



template<class Scalar>
inline EBasis Basis_F0_HEX_DD<Scalar>::getBasisType() const {
  return BASIS_FEM_DEFAULT;
}

template<class Scalar>
inline ECoordinates Basis_F0_HEX_DD<Scalar>::getCoordinateSystem() const {
  return COORDINATES_CARTESIAN;
}

}// namespace Intrepid
