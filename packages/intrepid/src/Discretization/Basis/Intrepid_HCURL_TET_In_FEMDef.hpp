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

/** \file   Intrepid_HCURL_TET_In_FEM_Def.hpp
    \brief  Definition file for FEM basis functions (Raviart-Thomas)
            of degree n for H(div) functions on TET.  
    \author Created by R. Kirby.
*/

namespace Intrepid {

  template<class Scalar, class ArrayScalar>
  Basis_HCURL_TET_In_FEM<Scalar,ArrayScalar>::Basis_HCURL_TET_In_FEM( const int n ,
                                                                      const EPointType pointType ):
    Phis_( n ),
    coeffs_( (n+1)*(n+2)*(n+3)/2 , n*(n+2)*(n+3)/2  )
  {
     const int N = n*(n+2)*(n+3)/2;
     this -> basisCardinality_  = N;
     this -> basisDegree_       = n;
     this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<4> >() );
     this -> basisType_         = BASIS_FEM_FIAT;
     this -> basisCoordinates_  = COORDINATES_CARTESIAN;
     this -> basisTagsAreSet_   = false;

     const int littleN = n*(n+1)*(n+2)/2;    // dim of (P_{n-1})^3 -- smaller space
     const int bigN = (n+1)*(n+2)*(n+3)/2;   // dim of (P_{n})^3 -- larger space
     const int start_PkH = (n-1)*n*(n+1)/6;  // dim of P({n-2}), offset into 
     const int dim_PkH = n*(n+1)*(n+2)/6 - start_PkH;
     const int scalarLittleN = littleN/3;
     const int scalarBigN = bigN/3;

     // first, need to project the basis for Nedelec space onto the
     // orthogonal basis of degree n
     // get coefficients of PkHx

     Teuchos::SerialDenseMatrix<int,Scalar> V1(bigN, littleN + 3 * dim_PkH);

     // these two loops get the first three sets of basis functions
     for (int i=0;i<scalarLittleN;i++) {
       for (int k=0;k<3;k++) {
         V1(i+k*scalarBigN,i+k*scalarLittleN) = 1.0;
       }
     }

     // first 3*scalarLittleN columns are (P_{n-1})^3 space


     // now I need to integrate { (x,y,z) \times } against the big basis
     // first, get a cubature rule.
     CubatureDirectTetDefault<Scalar,FieldContainer<Scalar> > myCub( 2 * n );
     FieldContainer<Scalar> cubPoints( myCub.getNumPoints() , 3 );
     FieldContainer<Scalar> cubWeights( myCub.getNumPoints() );
     myCub.getCubature( cubPoints , cubWeights );

     // tabulate the scalar orthonormal basis at cubature points
     FieldContainer<Scalar> phisAtCubPoints( scalarBigN , myCub.getNumPoints() );
     Phis_.getValues( phisAtCubPoints , cubPoints , OPERATOR_VALUE );



     // first set of these functions will write into the first dimPkH columns of remainder

     for (int j=0;j<dim_PkH;j++) { // loop over homogeneous polynomials
       // write into second spatial component, where
       // I integrate z phi_j phi_i
       for (int i=0;i<scalarBigN;i++) {
         V1(scalarBigN+i,littleN+j) = 0.0;
         for (int k=0;k<myCub.getNumPoints();k++) {
           V1(scalarBigN+i,littleN+j) -= cubWeights(k) * cubPoints(k,2) 
             * phisAtCubPoints(start_PkH+j,k)
             * phisAtCubPoints(i,k);
         }
       }  
       // write into third spatial component (-y phi_j, phi_i)
       for (int i=0;i<scalarBigN;i++) {
         V1(2*scalarBigN+i,littleN+j) = 0.0;
         for (int k=0;k<myCub.getNumPoints();k++) {
           V1(2*scalarBigN+i,littleN+j) += cubWeights(k) * cubPoints(k,1) 
             * phisAtCubPoints(start_PkH+j,k)
             * phisAtCubPoints(i,k);
         }
       }  
     }

     // second set of basis functions, write into second set of dimPkH columns
     for (int j=0;j<dim_PkH;j++) { // loop over homogeneous polynomials
       // write into first spatial component, where
       // I integrate -z phi_j phi_i
       for (int i=0;i<scalarBigN;i++) {
         V1(i,littleN+dim_PkH+j) = 0.0;
         for (int k=0;k<myCub.getNumPoints();k++) {
           V1(i,littleN+dim_PkH+j) += cubWeights(k) * cubPoints(k,2) 
             * phisAtCubPoints(start_PkH+j,k)
             * phisAtCubPoints(i,k);
         }
       } 
     
       // third spatial component, x phi_j phi_i
       for (int i=0;i<scalarBigN;i++) {
         V1(2*scalarBigN+i,littleN+dim_PkH+j) = 0.0;
         for (int k=0;k<myCub.getNumPoints();k++) {
           V1(2*scalarBigN+i,littleN+dim_PkH+j) -= cubWeights(k) * cubPoints(k,0) 
             * phisAtCubPoints(start_PkH+j,k)
             * phisAtCubPoints(i,k);
         }
       }  
     }
    
     // third clump of dimPkH columns
     for (int j=0;j<dim_PkH;j++) { // loop over homogeneous polynomials
       // write into first spatial component, where
       // I integrate y phi_j phi_i
       for (int i=0;i<scalarBigN;i++) {
         V1(i,littleN+2*dim_PkH+j) = 0.0;
         for (int k=0;k<myCub.getNumPoints();k++) {
           V1(i,littleN+2*dim_PkH+j) -= cubWeights(k) * cubPoints(k,1) 
             * phisAtCubPoints(start_PkH+j,k)
             * phisAtCubPoints(i,k);
         }
       }  
       // second spatial component, -x phi_j phi_i
       for (int i=0;i<scalarBigN;i++) {
         V1(scalarBigN+i,littleN+2*dim_PkH+j) = 0.0;
         for (int k=0;k<myCub.getNumPoints();k++) {
           V1(scalarBigN+i,littleN+2*dim_PkH+j) += cubWeights(k) * cubPoints(k,0) 
             * phisAtCubPoints(start_PkH+j,k)
             * phisAtCubPoints(i,k);
         }
       }  
     }

     // now I need to set up an SVD to get a basis for the space
     Teuchos::SerialDenseMatrix<int,Scalar> S(bigN,1);
     Teuchos::SerialDenseMatrix<int,Scalar> U(bigN, bigN);
     Teuchos::SerialDenseMatrix<int,Scalar> Vt(bigN,bigN);
     Teuchos::SerialDenseMatrix<int,Scalar> work(5*bigN,1);
     Teuchos::SerialDenseMatrix<int,Scalar> rWork(1,1);
     int info;

     Teuchos::LAPACK<int,Scalar> lala;

     lala.GESVD( 'A',  
                 'N',
                 V1.numRows() ,
                 V1.numCols() ,
                 V1.values() ,
                 V1.stride() ,
                 S.values() ,
                 U.values() ,
                 U.stride() ,
                 Vt.values() ,
                 Vt.stride() ,
                 work.values() ,
                 5*bigN ,
                 rWork.values() ,
                 &info );
                        
     int num_nonzero_sv = 0;
     for (int i=0;i<bigN;i++) {
       if (S(i,0) > INTREPID_TOL) {
	 num_nonzero_sv++;
       }
     }
     
     Teuchos::SerialDenseMatrix<int,Scalar> Uslender(bigN, num_nonzero_sv);
     for (int j=0;j<num_nonzero_sv;j++) {
       for (int i=0;i<bigN;i++) {
	 Uslender(i,j) = U(i,j);
       }
     }

     // apply nodes to big space
     Teuchos::SerialDenseMatrix<int,Scalar> V2(N, bigN);

     shards::CellTopology edgeTop(shards::getCellTopologyData<shards::Line<2> >() );
     shards::CellTopology faceTop(shards::getCellTopologyData<shards::Triangle<3> >() );


     const int numPtsPerEdge = PointTools::getLatticeSize( edgeTop ,
                                                           n+1 ,
                                                           1 );

     const int numPtsPerFace = PointTools::getLatticeSize( faceTop ,
                                                           n+1 ,
                                                           1 );

     const int numPtsPerCell = PointTools::getLatticeSize( this->basisCellTopology_ ,
                                                           n+1 ,
                                                           1 );
    
     // these hold the reference domain points that will be mapped to each edge or face
     FieldContainer<Scalar> oneDPts( numPtsPerEdge , 1 );
     FieldContainer<Scalar> twoDPts( numPtsPerFace , 2 );

     if (pointType == POINTTYPE_WARPBLEND) {
       CubatureDirectLineGauss<Scalar> edgeRule( numPtsPerEdge );
       FieldContainer<Scalar> edgeCubWts( numPtsPerEdge );
       edgeRule.getCubature( oneDPts , edgeCubWts );
     }
     else if (pointType == POINTTYPE_EQUISPACED ) {
       PointTools::getLattice<Scalar,FieldContainer<Scalar> >( oneDPts , 
                                                               edgeTop ,
                                                               n+1 , 
                                                               1 ,
                                                               pointType );
     }

     PointTools::getLattice<Scalar,FieldContainer<Scalar> >( twoDPts ,
                                                             faceTop ,
                                                             n+1 ,
                                                             1 ,
                                                             pointType );

     FieldContainer<Scalar> edgePts( numPtsPerEdge , 3 );
     FieldContainer<Scalar> phisAtEdgePoints( scalarBigN , numPtsPerEdge );

     FieldContainer<Scalar> facePts( numPtsPerFace , 3 );
     FieldContainer<Scalar> phisAtFacePoints( scalarBigN , 
                                             numPtsPerFace );   

     FieldContainer<Scalar> edgeTan( 3 );
    
     // loop over the edges
     for (int edge=0;edge<6;edge++) {
       CellTools<Scalar>::getReferenceEdgeTangent( edgeTan ,
                                                   edge ,
                                                   this->basisCellTopology_ );
       /* multiply by 2.0 to account for a scaling in Pavel's definition */
       for (int j=0;j<3;j++) {
         edgeTan(j) *= 2.0;
       }

       CellTools<Scalar>::mapToReferenceSubcell( edgePts ,
                                                 oneDPts ,
                                                 1 ,
                                                 edge ,
                                                 this->basisCellTopology_ );
      
       Phis_.getValues( phisAtEdgePoints , edgePts , OPERATOR_VALUE );
   
       // loop over points (rows of V2)
       for (int j=0;j<numPtsPerEdge;j++) {
	 // loop over orthonormal basis functions (columns of V2)
         for (int k=0;k<scalarBigN;k++) {
           for (int d=0;d<3;d++) {
             V2(edge*numPtsPerEdge+j,k+scalarBigN*d) = edgeTan(d) * phisAtEdgePoints(k,j);
           }
         }
       }
     }

   // handle the faces, if needed
    if (n > 1) {
      FieldContainer<Scalar> refFaceTanU(3);
      FieldContainer<Scalar> refFaceTanV(3);
      for (int face=0;face<4;face++) {
        CellTools<Scalar>::getReferenceFaceTangents( refFaceTanU ,
                                                    refFaceTanV ,
                                                    face ,
                                                    this->basisCellTopology_ );
        CellTools<Scalar>::mapToReferenceSubcell( facePts ,
                                                  twoDPts ,
                                                  2 ,
                                                  face ,
                                                  this->basisCellTopology_ );
        Phis_.getValues( phisAtFacePoints , facePts , OPERATOR_VALUE );
        for (int j=0;j<numPtsPerFace;j++) {
          for (int k=0;k<scalarBigN;k++) {
            for (int d=0;d<3;d++) {
              V2(6*numPtsPerEdge+2*face*numPtsPerFace+2*j,k+scalarBigN*d) =
                refFaceTanU(d) * phisAtFacePoints(k,j);
              V2(6*numPtsPerEdge+2*face*numPtsPerFace+2*j+1,k+scalarBigN*d) =
                refFaceTanV(d) * phisAtFacePoints(k,j);
            }
          }
        }
      }
   }

     // internal dof, if needed
     if (n > 2) {
       FieldContainer<Scalar> cellPoints( numPtsPerCell , 3 );
       PointTools::getLattice<Scalar,FieldContainer<Scalar> >( cellPoints ,
							       this->getBaseCellTopology() , 
							       n + 1 ,
							       1 ,
							       pointType );
       FieldContainer<Scalar> phisAtCellPoints( scalarBigN , numPtsPerCell );
       Phis_.getValues( phisAtCellPoints , cellPoints , OPERATOR_VALUE );
       for (int i=0;i<numPtsPerCell;i++) {
	 for (int j=0;j<scalarBigN;j++) {
	   for (int k=0;k<3;k++) {
	     V2(6*numPtsPerEdge+8*numPtsPerFace+k*numPtsPerCell+i,k*scalarBigN+j) = phisAtCellPoints(j,i);
	   }
	 }
       }
     }

    Teuchos::SerialDenseMatrix<int,Scalar> Vsdm( N , N );
    
    // multiply V2 * U --> V
    Vsdm.multiply( Teuchos::NO_TRANS , Teuchos::NO_TRANS , 1.0 , V2 , Uslender , 0.0 );

    Teuchos::SerialDenseSolver<int,Scalar> solver;
    solver.setMatrix( rcp( &Vsdm , false ) );

    solver.invert( );


    Teuchos::SerialDenseMatrix<int,Scalar> Csdm( bigN , N );
    Csdm.multiply( Teuchos::NO_TRANS , Teuchos::NO_TRANS , 1.0 , Uslender , Vsdm , 0.0 );

    //std::cout << Csdm << "\n";

    for (int i=0;i<bigN;i++) {
      for (int j=0;j<N;j++) {
        coeffs_(i,j) = Csdm(i,j);
      }
    }

    //std::cout << coeffs_ << std::endl;

  }
    
  template<class Scalar, class ArrayScalar>
  void Basis_HCURL_TET_In_FEM<Scalar, ArrayScalar>::initializeTags() {
    // Basis-dependent initializations
    int tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
    int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
    int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
    int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
    
    // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 
    
    int *tags = new int[ tagSize * this->getCardinality() ];
    int *tag_cur = tags;
    const int deg = this->getDegree();

    shards::CellTopology edgeTop(shards::getCellTopologyData<shards::Line<2> >() );
    shards::CellTopology faceTop(shards::getCellTopologyData<shards::Triangle<3> >() );


    const int numPtsPerEdge = PointTools::getLatticeSize( edgeTop ,
                                                          deg+1 ,
                                                          1 );

    const int numPtsPerFace = PointTools::getLatticeSize( faceTop ,
                                                          deg+1 ,
                                                          1 );

    const int numPtsPerCell = PointTools::getLatticeSize( this->basisCellTopology_ ,
                                                          deg+1 ,
                                                          1 );

    // edge dof first
    for (int e=0;e<6;e++) {
      for (int i=0;i<numPtsPerEdge;i++) {
        tag_cur[0] = 1;  tag_cur[1] = e;  tag_cur[2] = i;  tag_cur[3] = numPtsPerEdge;
        tag_cur += tagSize;
      }
    }

    // face dof, 2 * numPtsPerFace dof per face
    for (int f=0;f<4;f++) {
      for (int i=0;i<2*numPtsPerFace;i++) {
        tag_cur[0] = 2;  tag_cur[1] = f;  tag_cur[2] = i;  tag_cur[3] = 2*numPtsPerFace;
        tag_cur+= tagSize;
      }
    }
    
    // internal dof, 3 * numPtsPerCell
    for (int i=0;i<3*numPtsPerCell;i++) {
      tag_cur[0] = 3;  tag_cur[1] = 0;  tag_cur[2] = i;  tag_cur[3] = 3*numPtsPerCell;
      tag_cur += tagSize;
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
  void Basis_HCURL_TET_In_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
                                                              const ArrayScalar &  inputPoints,
                                                              const EOperator      operatorType) const {
  
    // Verify arguments
#ifdef HAVE_INTREPID_DEBUG
    Intrepid::getValues_HCURL_Args<Scalar, ArrayScalar>(outputValues,
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
              for (int d=0;d<3;d++) {
                outputValues(i,j,d) = 0.0;
              }
              for (int k=0;k<scalarBigN;k++) { // Dubiner bf
                for (int d=0;d<3;d++) {
                  outputValues(i,j,d) += coeffs_(k+d*scalarBigN,i) * phisCur(k,j);
                }
              }
            }
          }
        }
        break;
      case OPERATOR_CURL:
        {
          FieldContainer<Scalar> phisCur( scalarBigN , numPts , 3 );
          Phis_.getValues( phisCur , inputPoints , OPERATOR_GRAD );
          for (int i=0;i<outputValues.dimension(0);i++) { // bf loop
            for (int j=0;j<outputValues.dimension(1);j++) { // point loop
              outputValues(i,j,0) = 0.0;
              for (int k=0;k<scalarBigN;k++) {
                outputValues(i,j,0) += coeffs_(k+2*scalarBigN,i) * phisCur(k,j,1);
                outputValues(i,j,0) -= coeffs_(k+scalarBigN,i) * phisCur(k,j,2);
              }
              
              outputValues(i,j,1) = 0.0;
              for (int k=0;k<scalarBigN;k++) {
                outputValues(i,j,1) += coeffs_(k,i) * phisCur(k,j,2);
                outputValues(i,j,1) -= coeffs_(k+2*scalarBigN,i) * phisCur(k,j,0);
              }

              outputValues(i,j,2) = 0.0;
              for (int k=0;k<scalarBigN;k++) {
                outputValues(i,j,2) += coeffs_(k+scalarBigN,i) * phisCur(k,j,0);
                outputValues(i,j,2) -= coeffs_(k,i) * phisCur(k,j,1);
              }
            }
          }
        }
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument,
                            ">>> ERROR (Basis_HCURL_TET_In_FEM): Operator type not implemented");
        break;
      }
    }
    catch (std::invalid_argument &exception){
      TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument,
                          ">>> ERROR (Basis_HCURL_TET_In_FEM): Operator type not implemented");    
    }

  }
  

  
  template<class Scalar, class ArrayScalar>
  void Basis_HCURL_TET_In_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                              const ArrayScalar &    inputPoints,
                                                              const ArrayScalar &    cellVertices,
                                                              const EOperator        operatorType) const {
    TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                        ">>> ERROR (Basis_HCURL_TET_In_FEM): FEM Basis calling an FVD member function");
  }


}// namespace Intrepid

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

