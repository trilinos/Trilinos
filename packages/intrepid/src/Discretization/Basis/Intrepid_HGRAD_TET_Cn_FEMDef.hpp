#ifndef INTREPID_HGRAD_TET_CN_FEMDEF_HPP
#define INTREPID_HGRAD_TET_CN_FEMDEF_HPP
// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_HGRAD_TET_Cn_FEM_Def.hpp
    \brief  Definition file for FEM basis functions of degree n for H(grad) functions on TET.
    \author Created by R. Kirby and P. Bochev and D. Ridzal.
*/

namespace Intrepid {

  template<class Scalar, class ArrayScalar>
  Basis_HGRAD_TET_Cn_FEM<Scalar,ArrayScalar>::Basis_HGRAD_TET_Cn_FEM( const int n ,
                                                                      const EPointType pointType ):
    Phis( n ),
    V((n+1)*(n+2)*(n+3)/6,(n+1)*(n+2)*(n+3)/6),
    Vinv((n+1)*(n+2)*(n+3)/6,(n+1)*(n+2)*(n+3)/6),
    latticePts( (n+1)*(n+2)*(n+3)/6 , 3 )
  {
    const int N = (n+1)*(n+2)*(n+3)/6;
    this -> basisCardinality_  = N;
    this -> basisDegree_       = n;
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<4> >() );
    this -> basisType_         = BASIS_FEM_FIAT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;

    // construct lattice

    PointTools::getLattice<Scalar,FieldContainer<Scalar> >( latticePts ,
                                                            this->getBaseCellTopology() ,
                                                            n ,
                                                            0 ,
                                                            pointType );

    
    // form Vandermonde matrix.  Actually, this is the transpose of the VDM,
    // so we transpose on copy below.
  
    Phis.getValues( V , latticePts , OPERATOR_VALUE );

    // now I need to copy V into a Teuchos array to do the inversion
    Teuchos::SerialDenseMatrix<int,Scalar> Vsdm(N,N);
    for (int i=0;i<N;i++) {
      for (int j=0;j<N;j++) {
        Vsdm(i,j) = V(i,j);
      }
    }

    // invert the matrix
    Teuchos::SerialDenseSolver<int,Scalar> solver;
    solver.setMatrix( rcp( &Vsdm , false ) );
    solver.invert( );

    // now I need to copy the inverse into Vinv
    for (int i=0;i<N;i++) {
      for (int j=0;j<N;j++) {
        Vinv(i,j) = Vsdm(j,i);
      }
    }

  }  
  
  
  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_TET_Cn_FEM<Scalar, ArrayScalar>::initializeTags() {
  
    // Basis-dependent initializations
    int tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
    int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
    int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
    int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
  
    // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 

    int *tags = new int[ tagSize * this->getCardinality() ];
    int *tag_cur = tags;
    const int degree = this->getDegree();
    const int numEdgeDof = degree - 1;
    const int numFaceDof = PointTools::getLatticeSize( shards::CellTopology( shards::getCellTopologyData<shards::Triangle<3> >() ) ,
                                                      degree , 
                                                      1);
    const int numCellDof = PointTools::getLatticeSize( this->getBaseCellTopology() ,
                                                      degree ,
                                                      1 );
    int edge_dof_cur[] = {0,0,0,0,0,0};
    int face_dof_cur[] = {0,0,0,0};
    int cell_dof_cur = 0;

    // this is the really big mess :(
    // BEGIN DOF ON BOTTOM FACE
    // first vertex: 0
    tag_cur[0] = 0;  tag_cur[1] = 0;  tag_cur[2] = 0;  tag_cur[3] = 1;
    tag_cur += tagSize;
    // end first vertex

    // internal points on line from vertex 0 to vertex 1.  This is edge 0
    for (int i=1;i<degree;i++) {
      tag_cur[0] = 1;  tag_cur[1] = 0;  tag_cur[2] = edge_dof_cur[0];  tag_cur[3] = numEdgeDof;
      edge_dof_cur[0]++;
      tag_cur += tagSize;
    }
    // end line from vertex 0 to vertex 1

    // begin vertex 1
    tag_cur[0] = 0;  tag_cur[1] = 1;  tag_cur[2] = 0;  tag_cur[3] = 1;
    tag_cur += tagSize;
    // end vertex 1

    // internal lines on bottom face
    for (int i=1;i<degree;i++) {
      // first dof is on edge 2
      tag_cur[0] = 1;  tag_cur[1] = 2;  tag_cur[2] = edge_dof_cur[2];  tag_cur[3] = numEdgeDof;
      edge_dof_cur[2]++;
      tag_cur += tagSize;
      // end dof on edge 2

      // internal points are on bottom face, which is face 3
      for (int j=1;j<degree-i;j++) {
        tag_cur[0] = 2;  tag_cur[1] = 3;  tag_cur[2] = face_dof_cur[3];  tag_cur[3] = numFaceDof;
        face_dof_cur[3]++;
        tag_cur += tagSize;
      }
      // end internal points on face 

      // last dof is on edge 1
      tag_cur[0] = 1;  tag_cur[1] = 1;  tag_cur[2] = edge_dof_cur[1];  tag_cur[3] = numEdgeDof;
      edge_dof_cur[1]++;
      tag_cur += tagSize;
      // end dof on edge 1
    }
    // end internal lines on bottom face

    // vertex 2 on bottom face
    tag_cur[0] = 0;  tag_cur[1] = 2;  tag_cur[2] = 0;  tag_cur[3] = 1;
    tag_cur += tagSize;
    // end vertex 2 on bottom face

    // END DOF ON BOTTOM FACE
    
    // BEGIN DOF ON INTERNAL FACE SLICES (ascending z)
    for (int i=1;i<degree;i++) {
      //   bottom line of internal face
      //     first point is on edge 3 (from vertex 0 to 3)
      tag_cur[0] = 1;  tag_cur[1] = 3;  tag_cur[2] = edge_dof_cur[3];  tag_cur[3] = numEdgeDof;
      edge_dof_cur[3]++;
      tag_cur += tagSize;
      //     end first point
      //     points internal to face of vertices (0,1,3), which is face 0
      for (int j=1;j<degree-i;j++) {
        tag_cur[0] = 2;  tag_cur[1] = 0;  tag_cur[2] = face_dof_cur[0];  tag_cur[3] = numFaceDof;
        face_dof_cur[0]++;
        tag_cur += tagSize;
      }
      //     end points internal to face 0
      //     last point on bottom line is on edge 4
      tag_cur[0] = 1;  tag_cur[1] = 4;  tag_cur[2] = edge_dof_cur[4];  tag_cur[3] = numEdgeDof;
      edge_dof_cur[4]++;
      tag_cur += tagSize;
      //     end last point on bottom edge
      //  end bottom line of internal face

      //  begin internal lines of internal face
      for (int j=1;j<degree-i;j++) {
        //    first point on line is on face of vertices (0,3,2), which is face 2
        tag_cur[0] = 2;  tag_cur[1] = 2;  tag_cur[2] = face_dof_cur[2];  tag_cur[3] = numFaceDof;
        face_dof_cur[2]++;
        tag_cur += tagSize;
        //    end first point of line
        //    begin internal points on the cell
        for (int k=1;k<degree-i-j;k++) {
          tag_cur[0] = 3;  tag_cur[1] = 0;  tag_cur[2] = cell_dof_cur;  tag_cur[3] = numCellDof;
          cell_dof_cur++;
          tag_cur += tagSize;
        }
        //    end internal points on the cell
        //    last point on the line is on face with vertices (1,2,3) , which is face 1
        tag_cur[0] = 2;  tag_cur[1] = 1;  tag_cur[2] = face_dof_cur[1];  tag_cur[3] = numFaceDof;
        face_dof_cur[1]++;
        tag_cur += tagSize;
        //    end last point of line
      }
      //  end internal lines of internal face
      // begin top point on current face slice:  on edge 5
      tag_cur[0] = 1;  tag_cur[1] = 5;  tag_cur[2] = edge_dof_cur[5];  tag_cur[3] = numEdgeDof;
      edge_dof_cur[5]++;
      tag_cur += 4;
      // end top point on current face slice
    }
    // END DOF ON INTERNAL FACE SLICES

    // TOP VERTEX: 3
    tag_cur[0] = 0;  tag_cur[1] = 3;  tag_cur[2] = 0;  tag_cur[3] = 1;
    // END TOP VERTEX:3

    // end of really big mess :)
    

    
    Intrepid::setOrdinalTagData(this -> tagToOrdinal_,
                                this -> ordinalToTag_,
                                tags,
                                this -> basisCardinality_,
                                tagSize,
                                posScDim,
                                posScOrd,
                                posDfOrd);

    delete []tags;
  
  }  



  template<class Scalar, class ArrayScalar> 
  void Basis_HGRAD_TET_Cn_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
                                                              const ArrayScalar &  inputPoints,
                                                              const EOperator      operatorType) const {
  
    // Verify arguments
#ifdef HAVE_INTREPID_DEBUG
    Intrepid::getValues_HGRAD_Args<Scalar, ArrayScalar>(outputValues,
                                                        inputPoints,
                                                        operatorType,
                                                        this -> getBaseCellTopology(),
                                                        this -> getCardinality() );
#endif
    const int numPts = inputPoints.dimension(0);
    const int numBf = this->getCardinality();

    try {
      switch (operatorType) {
      case OPERATOR_VALUE:
        {
          FieldContainer<Scalar> phisCur( numBf , numPts );
          Phis.getValues( phisCur , inputPoints , operatorType );
          for (int i=0;i<outputValues.dimension(0);i++) {
            for (int j=0;j<outputValues.dimension(1);j++) {
              outputValues(i,j) = 0.0;
              for (int k=0;k<this->getCardinality();k++) {
                outputValues(i,j) += this->Vinv(k,i) * phisCur(k,j);
              }
            }
          }
        }
        break;
      case OPERATOR_GRAD:
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
          const int dkcard = 
            (operatorType == OPERATOR_GRAD)? getDkCardinality(OPERATOR_D1,3): getDkCardinality(operatorType,3);
          
          FieldContainer<Scalar> phisCur( numBf , numPts , dkcard );
          Phis.getValues( phisCur , inputPoints , operatorType );

          for (int i=0;i<outputValues.dimension(0);i++) {
            for (int j=0;j<outputValues.dimension(1);j++) {
              for (int k=0;k<outputValues.dimension(2);k++) {
                outputValues(i,j,k) = 0.0;
                for (int l=0;l<this->getCardinality();l++) {
                  outputValues(i,j,k) += this->Vinv(l,i) * phisCur(l,j,k);
                }
              }
            }
          }

        }
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument,
                            ">>> ERROR (Basis_HGRAD_TET_Cn_FEM): Operator type not implemented");
        break;
      }
    }
    catch (std::invalid_argument &exception){
      TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_TET_Cn_FEM): Operator type not implemented");    
    }

  }
  

  
  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_TET_Cn_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                              const ArrayScalar &    inputPoints,
                                                              const ArrayScalar &    cellVertices,
                                                              const EOperator        operatorType) const {
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                        ">>> ERROR (Basis_HGRAD_TET_Cn_FEM): FEM Basis calling an FVD member function");
  }


}// namespace Intrepid
#endif

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

