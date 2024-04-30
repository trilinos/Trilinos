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

/** \file   Intrepid_HDIV_TET_In_FEM_Def.hpp
    \brief  Definition file for FEM basis functions (Raviart-Thomas)
    of degree n for H(div) functions on TET.  
    \author Created by R. Kirby.
*/

namespace Intrepid {

  template<class Scalar, class ArrayScalar>
  Basis_HDIV_TET_In_FEM<Scalar,ArrayScalar>::Basis_HDIV_TET_In_FEM( const int n ,
                                                                    const EPointType pointType ):
    Phis_( n ),
    coeffs_( (n+1)*(n+2)*(n+3)/2 , n*(n+1)*(n+3)/2 )
  {
    const int N = n*(n+1)*(n+3)/2;
    this -> basisCardinality_  = N;
    this -> basisDegree_       = n;
    this -> basisCellTopology_ 
      = shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<4> >() );
    this -> basisType_         = BASIS_FEM_FIAT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;


    const int littleN = n*(n+1)*(n+2)/2;   // dim of (P_{n-1})^3 -- smaller space
    const int bigN = (n+1)*(n+2)*(n+3)/2;  // dim of (P_{n})^2 -- larger space
    const int start_PkH = (n-1)*n*(n+1)/6; // dim of P({n-2}), offset into 
    const int dim_PkH = n*(n+1)*(n+2)/6 - start_PkH;
    const int scalarLittleN = littleN/3;
    const int scalarBigN = bigN/3;

    // first, need to project the basis for RT space onto the
    // orthogonal basis of degree n
    // get coefficients of PkHx

    Teuchos::SerialDenseMatrix<int,Scalar> V1(bigN, N);

    // basis for the space is 
    // { (phi_i,0,0) }_{i=0}^{scalarLittleN-1} ,
    // { (0,phi_i,0) }_{i=0}^{scalarLittleN-1} ,
    // { (0,0,phi_i) }_{i=0}^{scalarLittleN-1} ,
    // { (x,y) . phi_i}_{i=startPKH}^{scalarLittleN-1}
    // columns of V1 are expansion of this basis in terms of the orthogonal basis
    // for P_{n}^3


    // these two loops get the first three sets of basis functions
    for (int i=0;i<scalarLittleN;i++) {
      for (int k=0;k<3;k++) {
        V1(i+k*scalarBigN,i+k*scalarLittleN) = 1.0;
      }
    }

    // now I need to integrate { (x,y,z) phi } against the big basis
    // first, get a cubature rule.
    CubatureDirectTetDefault<Scalar,FieldContainer<Scalar> > myCub( 2 * n );
    FieldContainer<Scalar> cubPoints( myCub.getNumPoints() , 3 );
    FieldContainer<Scalar> cubWeights( myCub.getNumPoints() );
    myCub.getCubature( cubPoints , cubWeights );

    // tabulate the scalar orthonormal basis at cubature points
    FieldContainer<Scalar> phisAtCubPoints( scalarBigN , myCub.getNumPoints() );
    Phis_.getValues( phisAtCubPoints , cubPoints , OPERATOR_VALUE );


    // now do the integration
    for (int i=0;i<dim_PkH;i++) {
      for (int j=0;j<scalarBigN;j++) {  // int (x,y,z) phi_i \cdot (phi_j,0,0)
        V1(j,littleN+i) = 0.0;
        for (int d=0;d<3;d++) {
          for (int k=0;k<myCub.getNumPoints();k++) {
            V1(j+d*scalarBigN,littleN+i) += 
              cubWeights(k) * cubPoints(k,d) 
              * phisAtCubPoints(start_PkH+i,k) 
              * phisAtCubPoints(j,k);
          }
        }
      }
    }


    // next, apply the RT nodes (rows) to the basis for (P_n)^3 (columns)
    Teuchos::SerialDenseMatrix<int,Scalar> V2(N , bigN);

    shards::CellTopology faceTop(shards::getCellTopologyData<shards::Triangle<3> >() );
    const int numPtsPerFace = PointTools::getLatticeSize( faceTop ,
                                                          n+2 ,
                                                          1 );

    FieldContainer<Scalar> twoDPts( numPtsPerFace , 2 );
    PointTools::getLattice<Scalar,FieldContainer<Scalar> >( twoDPts ,
                                                            faceTop ,
                                                            n+2 ,
                                                            1 ,
                                                            pointType );

    // holds the image of the triangle points on each face.
    FieldContainer<Scalar> facePts( numPtsPerFace , 3 );
    FieldContainer<Scalar> phisAtFacePoints( scalarBigN , 
                                            numPtsPerFace );
    
    

    // these are scaled by the appropriate face areas.
    // area of faces 0,2,3 are 0.5
    // area of face 1 is sqrt(3)/2

    Scalar normal[][4] = { {0.0,0.5,-0.5,0.0},
                          {-0.5,0.5,0.0,0.0},
                          {0.0,0.5,0.0,-0.5} };

    for (int i=0;i<4;i++) {  // loop over faces
      CellTools<Scalar>::mapToReferenceSubcell( facePts ,
                                                twoDPts ,
                                                2 ,
                                                i ,
                                                this->basisCellTopology_ );

      Phis_.getValues( phisAtFacePoints , facePts , OPERATOR_VALUE );

      // loop over points (rows of V2)
      for (int j=0;j<numPtsPerFace;j++) {
        // loop over orthonormal basis functions (columns of V2)
        for (int k=0;k<scalarBigN;k++) {
          for (int l=0;l<3;l++) {
            V2(numPtsPerFace*i+j,k+l*scalarBigN) = normal[l][i] * phisAtFacePoints(k,j);
          }
        }
      }
    }

    // remaining nodes point values of each vector component on interior
    // points of a lattice of degree+2
    // This way, RT0 --> degree = 1 and internal lattice has no points
    // RT1 --> degree = 2, and internal lattice has one point (inside of quartic)
    if (n > 1) {
      const int numInternalPoints = PointTools::getLatticeSize( this->getBaseCellTopology() ,
                                                                n + 2 ,
                                                                1 );

      FieldContainer<Scalar> internalPoints( numInternalPoints , 3 );
      PointTools::getLattice<Scalar,FieldContainer<Scalar> >( internalPoints ,
                                                              this->getBaseCellTopology() , 
                                                              n + 2 ,
                                                              1 ,
                                                              pointType );
    
      FieldContainer<Scalar> phisAtInternalPoints( scalarBigN , numInternalPoints );
      Phis_.getValues( phisAtInternalPoints , internalPoints , OPERATOR_VALUE );

      // copy values into right positions of V2
      for (int i=0;i<numInternalPoints;i++) {
        for (int j=0;j<scalarBigN;j++) {
          for (int k=0;k<3;k++) {
            V2(4*numPtsPerFace+k*numInternalPoints+i,k*scalarBigN+j) = phisAtInternalPoints(j,i);
          }
        }
      }
    }
    
    Teuchos::SerialDenseMatrix<int,Scalar> Vsdm( N , N );

    // multiply V2 * V1 --> V
    Vsdm.multiply( Teuchos::NO_TRANS , Teuchos::NO_TRANS , 1.0 , V2 , V1 , 0.0 );

    //     std::cout << "Vandermonde:\n";
    //     std::cout << Vsdm << "\n";
    //     std::cout << "End Vandermonde\n";
    
    Teuchos::SerialDenseSolver<int,Scalar> solver;
    solver.setMatrix( rcp( &Vsdm , false ) );
    solver.invert( );


    Teuchos::SerialDenseMatrix<int,Scalar> Csdm( bigN , N );
    Csdm.multiply( Teuchos::NO_TRANS , Teuchos::NO_TRANS , 1.0 , V1 , Vsdm , 0.0 );

    //std::cout << Csdm << "\n";

    for (int i=0;i<bigN;i++) {
      for (int j=0;j<N;j++) {
        coeffs_(i,j) = Csdm(i,j);
      }
    }
  }  
    
  template<class Scalar, class ArrayScalar>
  void Basis_HDIV_TET_In_FEM<Scalar, ArrayScalar>::initializeTags() {
  
    // Basis-dependent initializations
    int tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
    int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
    int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
    int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
  
    // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 

    int *tags = new int[ tagSize * this->getCardinality() ];
    int *tag_cur = tags;
    const int deg = this->getDegree();

    const int numPtsPerFace = deg*(deg+1)/2;

    // there are degree internal dofs on each edge -- normals.  Let's do them
    for (int f=0;f<4;f++) {
      for (int i=0;i<numPtsPerFace;i++) {
        tag_cur[0] = 2;  tag_cur[1] = f;  tag_cur[2] = i;  tag_cur[3] = numPtsPerFace;
        tag_cur += tagSize;
      }
    }
    // end face dofs

    // the rest of the dofs are internal
    const int numInternalDof = this->getCardinality() - 4 * numPtsPerFace;
    int internalDofCur=0;
    for (int i=4*numPtsPerFace;i<this->getCardinality();i++) {
      tag_cur[0] = 3;  tag_cur[1] = 0;  tag_cur[2] = internalDofCur; tag_cur[3] = numInternalDof;
      tag_cur += tagSize;
      internalDofCur++;
    }
    
    
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
  void Basis_HDIV_TET_In_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
                                                            const ArrayScalar &  inputPoints,
                                                            const EOperator      operatorType) const {
  
    // Verify arguments
#ifdef HAVE_INTREPID_DEBUG
    Intrepid::getValues_HDIV_Args<Scalar, ArrayScalar>(outputValues,
                                                      inputPoints,
                                                      operatorType,
                                                      this -> getBaseCellTopology(),
                                                      this -> getCardinality() );
#endif
    const int numPts = inputPoints.dimension(0);
    const int deg = this -> getDegree();
    const int scalarBigN = (deg+1)*(deg+2)*(deg+3)/6;

    try {
      switch (operatorType) {
      case OPERATOR_VALUE:
        {
          FieldContainer<Scalar> phisCur( scalarBigN , numPts );
          Phis_.getValues( phisCur , inputPoints , OPERATOR_VALUE );

          for (int i=0;i<outputValues.dimension(0);i++) { // RT bf
            for (int j=0;j<outputValues.dimension(1);j++) {  // point
              for (int l=0;l<3;l++) {
                outputValues(i,j,l) = 0.0;
              }
              for (int k=0;k<scalarBigN;k++) { // Dubiner bf
                for (int l=0;l<3;l++) { // vector components
                  outputValues(i,j,l) += coeffs_(k+l*scalarBigN,i) * phisCur(k,j);
                }
              }
            }
          }
        }
        break;
      case OPERATOR_DIV:
        {
          FieldContainer<Scalar> phisCur( scalarBigN , numPts , 3 );
          Phis_.getValues( phisCur , inputPoints , OPERATOR_GRAD );
          for (int i=0;i<outputValues.dimension(0);i++) { // bf loop
            for (int j=0;j<outputValues.dimension(1);j++) { // point loop
              outputValues(i,j) = 0.0;
              for (int k=0;k<scalarBigN;k++) {
                outputValues(i,j) += coeffs_(k,i) * phisCur(k,j,0);
                outputValues(i,j) += coeffs_(k+scalarBigN,i) * phisCur(k,j,1);
                outputValues(i,j) += coeffs_(k+2*scalarBigN,i) * phisCur(k,j,2);
              }
            }
          }
        }
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument,
                            ">>> ERROR (Basis_HDIV_TET_In_FEM): Operator type not implemented");
        break;
      }
    }
    catch (std::invalid_argument &exception){
      TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument,
                          ">>> ERROR (Basis_HDIV_TET_In_FEM): Operator type not implemented");    
    }

  }
  

  
  template<class Scalar, class ArrayScalar>
  void Basis_HDIV_TET_In_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                            const ArrayScalar &    inputPoints,
                                                            const ArrayScalar &    cellVertices,
                                                            const EOperator        operatorType) const {
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                        ">>> ERROR (Basis_HDIV_TET_In_FEM): FEM Basis calling an FVD member function");
  }


}// namespace Intrepid

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

