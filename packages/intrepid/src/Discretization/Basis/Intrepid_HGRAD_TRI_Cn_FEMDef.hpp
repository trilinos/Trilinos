#ifndef INTREPID_HGRAD_TRI_CN_FEMDEF_HPP
#define INTREPID_HGRAD_TRI_CN_FEMDEF_HPP
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

/** \file   Intrepid_HGRAD_TRI_Cn_FEM_Def.hpp
    \brief  Definition file for FEM basis functions of degree n for H(grad) functions on TRI.
    \author Created by R. Kirby and P. Bochev and D. Ridzal.
*/

namespace Intrepid {

  template<class Scalar, class ArrayScalar>
  Basis_HGRAD_TRI_Cn_FEM<Scalar,ArrayScalar>::Basis_HGRAD_TRI_Cn_FEM( const int n ,
                                                                      const EPointType pointType ):
    Phis( n ),
    V((n+1)*(n+2)/2,(n+1)*(n+2)/2),
    Vinv((n+1)*(n+2)/2,(n+1)*(n+2)/2),
    latticePts( (n+1)*(n+2)/2 , 2 )
  {
    const int N = (n+1)*(n+2)/2;
    this -> basisCardinality_  = N;
    this -> basisDegree_       = n;
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3> >() );
    this -> basisType_         = BASIS_FEM_FIAT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;

    // construct lattice

    shards::CellTopology myTri_3( shards::getCellTopologyData< shards::Triangle<3> >() );  
    PointTools::getLattice<Scalar,FieldContainer<Scalar> >( latticePts ,
                                                            myTri_3 ,
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
  void Basis_HGRAD_TRI_Cn_FEM<Scalar, ArrayScalar>::initializeTags() {
  
    // Basis-dependent initializations
    int tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
    int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
    int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
    int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
  
    // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 

    int *tags = new int[ tagSize * this->getCardinality() ];
    int *tag_cur = tags;
    const int degree = this->getDegree();

    // BEGIN DOF ALONG BOTTOM EDGE

    // the first dof is on vertex 0
    tag_cur[0] = 0;  tag_cur[1] = 0;  tag_cur[2] = 0;  tag_cur[3] = 1;
    tag_cur += tagSize;

    // next degree-1 dof are on edge 0
    for (int i=1;i<degree;i++) {
      tag_cur[0] = 1;  tag_cur[1] = 0; tag_cur[2] = i-1;  tag_cur[3] = degree-1;
      tag_cur += tagSize;
    }

    // last dof is on vertex 1
    tag_cur[0] = 0;  tag_cur[1] = 1;  tag_cur[2] = 0;  tag_cur[3] = 1;
    tag_cur += tagSize;

    // END DOF ALONG BOTTOM EDGE

    int num_internal_dof = PointTools::getLatticeSize( this->getBaseCellTopology() ,
                                                      this->getDegree() ,
                                                      1 );

    int internal_dof_cur = 0;
  
    // BEGIN DOF ALONG INTERNAL HORIZONTAL LINES
    for (int i=1;i<degree;i++) {      
      // first dof along line is on edge #2
      tag_cur[0] = 1;  tag_cur[1] = 2;  tag_cur[2] = i-1;  tag_cur[3] = degree-1;
      tag_cur += tagSize;

      // next dof are internal
      for (int j=1;j<degree-i;j++) {
        tag_cur[0] = 2;  tag_cur[1] = 0;  tag_cur[2] = internal_dof_cur;  tag_cur[3] = num_internal_dof;
        internal_dof_cur++;
        tag_cur += tagSize;
      }

      // last dof along line is on edge 1
      tag_cur[0] = 1;  tag_cur[1] = 1;  tag_cur[2] = i-1;  tag_cur[3] = degree-1;
      tag_cur += tagSize;

    }
    // END DOF ALONG INTERNAL HORIZONTAL LINES
  
    // LAST DOF IS ON VERTEX 2
    tag_cur[0] = 0;  tag_cur[1] = 2;  tag_cur[2] = 0;  tag_cur[3] = 1;
    // END LAST DOF

    
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
  void Basis_HGRAD_TRI_Cn_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
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
            (operatorType == OPERATOR_GRAD)? getDkCardinality(OPERATOR_D1,2): getDkCardinality(operatorType,2);
          
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
      case OPERATOR_CURL:  // only works in 2d. first component is -d/dy, second is d/dx
        {
          FieldContainer<Scalar> phisCur( numBf , numPts , getDkCardinality( OPERATOR_D1 , 2 ) );
          Phis.getValues( phisCur , inputPoints , OPERATOR_D1 );

          for (int i=0;i<outputValues.dimension(0);i++) {
            for (int j=0;j<outputValues.dimension(1);j++) {
              outputValues(i,j,0) = 0.0;
              outputValues(i,j,1) = 0.0;
              for (int k=0;k<this->getCardinality();k++) {
                outputValues(i,j,0) = this->Vinv(k,i) * phisCur(k,j,1);
              }
              for (int k=0;k<this->getCardinality();k++) {
                outputValues(i,j,1) -= this->Vinv(k,i) * phisCur(k,j,0);
              }
            }
          }
        }
        break;
      default:
        TEST_FOR_EXCEPTION( true , std::invalid_argument,
                            ">>> ERROR (Basis_HGRAD_TRI_Cn_FEM): Operator type not implemented");
        break;
      }
    }
    catch (std::invalid_argument &exception){
      TEST_FOR_EXCEPTION( true , std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_TRI_Cn_FEM): Operator type not implemented");    
    }

  }
  

  
  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_TRI_Cn_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                              const ArrayScalar &    inputPoints,
                                                              const ArrayScalar &    cellVertices,
                                                              const EOperator        operatorType) const {
    TEST_FOR_EXCEPTION( (true), std::logic_error,
                        ">>> ERROR (Basis_HGRAD_TRI_Cn_FEM): FEM Basis calling an FVD member function");
  }


}// namespace Intrepid
#endif
