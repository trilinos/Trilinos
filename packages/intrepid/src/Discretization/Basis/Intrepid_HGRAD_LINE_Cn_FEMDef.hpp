#ifndef INTREPID_HGRAD_LINE_CN_FEMDEF_HPP
#define INTREPID_HGRAD_LINE_CN_FEMDEF_HPP
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
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_HGRAD_LINE_Cn_FEM_Def.hpp
    \brief  Definition file for FEM basis functions of degree n for H(grad) functions on LINE.
    \author Created by R. Kirby and P. Bochev and D. Ridzal.
*/

namespace Intrepid {

  template<class Scalar, class ArrayScalar>
  Basis_HGRAD_LINE_Cn_FEM<Scalar,ArrayScalar>::Basis_HGRAD_LINE_Cn_FEM( const int n ,
									const ArrayScalar &pts ):
    latticePts_( n+1 , 1 ),
    Phis_( n ),
    V_(n+1,n+1),
    Vinv_(n+1,n+1)
  {
    const int N = n+1;
    this -> basisCardinality_  = N;
    this -> basisDegree_       = n;
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >() );
    this -> basisType_         = BASIS_FEM_FIAT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;


    // check validity of points
    for (int i=0;i<n;i++) {
      TEST_FOR_EXCEPTION( pts(i,0) >= pts(i+1,0) ,
			  std::runtime_error ,
			  "Intrepid::Basis_HGRAD_LINE_Cn_FEM Illegal points given to constructor" );
    }

    // copy points int latticePts, correcting endpoints if needed
    if (std::abs(pts(0,0)+1.0) < INTREPID_TOL) {
      latticePts_(0,0) = -1.0;
    }
    else {
      latticePts_(0,0) = pts(0,0);
    }
    for (int i=1;i<n;i++) {
      latticePts_(i,0) = pts(i,0);
    }
    if (std::abs(pts(n,0)-1.0) < INTREPID_TOL) {
      latticePts_(n,0) = 1.0;
    }
    else {
      latticePts_(n,0) = pts(n,0);
    }
    
    // form Vandermonde matrix.  Actually, this is the transpose of the VDM,
    // so we transpose on copy below.
  
    Phis_.getValues( V_ , latticePts_ , OPERATOR_VALUE );

    // now I need to copy V into a Teuchos array to do the inversion
    Teuchos::SerialDenseMatrix<int,Scalar> Vsdm(N,N);
    for (int i=0;i<N;i++) {
      for (int j=0;j<N;j++) {
        Vsdm(i,j) = V_(i,j);
      }
    }

    // invert the matrix
    Teuchos::SerialDenseSolver<int,Scalar> solver;
    solver.setMatrix( rcp( &Vsdm , false ) );
    solver.invert( );

    // now I need to copy the inverse into Vinv
    for (int i=0;i<N;i++) {
      for (int j=0;j<N;j++) {
        Vinv_(i,j) = Vsdm(j,i);
      }
    }

  }  

  template<class Scalar, class ArrayScalar>
  Basis_HGRAD_LINE_Cn_FEM<Scalar,ArrayScalar>::Basis_HGRAD_LINE_Cn_FEM( const int n ,
									const EPointType &pointType ):
    latticePts_( n+1 , 1 ),
    Phis_( n ),
    V_(n+1,n+1),
    Vinv_(n+1,n+1)
  {
    const int N = n+1;
    this -> basisCardinality_  = N;
    this -> basisDegree_       = n;
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >() );
    this -> basisType_         = BASIS_FEM_FIAT;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;

    switch(pointType) {
    case POINTTYPE_EQUISPACED:
      PointTools::getLattice<Scalar,FieldContainer<Scalar> >( latticePts_ ,  this->basisCellTopology_ , n , 0 , POINTTYPE_EQUISPACED );
      break;
    case POINTTYPE_SPECTRAL: 
      PointTools::getLattice<Scalar,FieldContainer<Scalar> >( latticePts_ ,  this->basisCellTopology_ , n , 0 , POINTTYPE_WARPBLEND );
      break;
    case POINTTYPE_SPECTRAL_OPEN: 
      PointTools::getGaussPoints<Scalar,FieldContainer<Scalar> >( latticePts_ , n );
      break;
    default:
      TEST_FOR_EXCEPTION( true , std::invalid_argument , "Basis_HGRAD_LINE_Cn_FEM:: invalid point type" );
      break;
    }

    // form Vandermonde matrix.  Actually, this is the transpose of the VDM,
    // so we transpose on copy below.
  
    Phis_.getValues( V_ , latticePts_ , OPERATOR_VALUE );

    // now I need to copy V into a Teuchos array to do the inversion
    Teuchos::SerialDenseMatrix<int,Scalar> Vsdm(N,N);
    for (int i=0;i<N;i++) {
      for (int j=0;j<N;j++) {
        Vsdm(i,j) = V_(i,j);
      }
    }

    // invert the matrix
    Teuchos::SerialDenseSolver<int,Scalar> solver;
    solver.setMatrix( rcp( &Vsdm , false ) );
    solver.invert( );

    // now I need to copy the inverse into Vinv
    for (int i=0;i<N;i++) {
      for (int j=0;j<N;j++) {
        Vinv_(i,j) = Vsdm(j,i);
      }
    }
  }  
  
  
  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_LINE_Cn_FEM<Scalar, ArrayScalar>::initializeTags() {
  
    // Basis-dependent initializations
    int tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
    int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
    int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
    int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
  
    // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 

    int *tags = new int[ tagSize * this->getCardinality() ];

    int internal_dof;
    int edge_dof;

    const int n = this->getDegree();

    // now we check the points for association 
    if (latticePts_(0,0) == -1.0) {
      tags[0] = 0;
      tags[1] = 0;
      tags[2] = 0;
      tags[3] = 1;
      edge_dof = 1;
      internal_dof = n-1;
    }
    else {
      tags[0] = 1;
      tags[1] = 0;
      tags[2] = 0;
      tags[3] = n+1;
      edge_dof = 0;
      internal_dof = n+1;
    }
    for (int i=1;i<n;i++) {
      tags[4*i] = 1;
      tags[4*i+1] = 0;
      tags[4*i+2] = -edge_dof + i;
      tags[4*i+3] = internal_dof;
    }
    if (latticePts_(n,0) == 1.0) {
      tags[4*n] = 0;
      tags[4*n+1] = 1;
      tags[4*n+2] = 0;
      tags[4*n+3] = 1;
    }
    else {
      tags[4*n] = 1;
      tags[4*n+1] = 0;
      tags[4*n+2] = n;
      tags[4*n+3] = n;
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
  void Basis_HGRAD_LINE_Cn_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
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
          Phis_.getValues( phisCur , inputPoints , operatorType );
          for (int i=0;i<outputValues.dimension(0);i++) {
            for (int j=0;j<outputValues.dimension(1);j++) {
              outputValues(i,j) = 0.0;
              for (int k=0;k<this->getCardinality();k++) {
                outputValues(i,j) += this->Vinv_(k,i) * phisCur(k,j);
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
            (operatorType == OPERATOR_GRAD)? getDkCardinality(OPERATOR_D1,1): getDkCardinality(operatorType,1);
          
          FieldContainer<Scalar> phisCur( numBf , numPts , dkcard );
          Phis_.getValues( phisCur , inputPoints , operatorType );

          for (int i=0;i<outputValues.dimension(0);i++) {
            for (int j=0;j<outputValues.dimension(1);j++) {
              for (int k=0;k<outputValues.dimension(2);k++) {
                outputValues(i,j,k) = 0.0;
                for (int l=0;l<this->getCardinality();l++) {
                  outputValues(i,j,k) += this->Vinv_(l,i) * phisCur(l,j,k);
                }
              }
            }
          }
        }
        break;
      default:
        TEST_FOR_EXCEPTION( true , std::invalid_argument,
                            ">>> ERROR (Basis_HGRAD_LINE_Cn_FEM): Operator type not implemented" );
        break;
      }
    }
    catch (std::invalid_argument &exception){
      TEST_FOR_EXCEPTION( true , std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_LINE_Cn_FEM): Operator failed");    
    }

  }
  

  
  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_LINE_Cn_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                              const ArrayScalar &    inputPoints,
                                                              const ArrayScalar &    cellVertices,
                                                              const EOperator        operatorType) const {
    TEST_FOR_EXCEPTION( (true), std::logic_error,
                        ">>> ERROR (Basis_HGRAD_LINE_Cn_FEM): FEM Basis calling an FVD member function");
  }


  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_LINE_Cn_FEM<Scalar,ArrayScalar>::getDofCoords( ArrayScalar & dofCoords ) const
  {
    for (int i=0;i<latticePts_.dimension(0);i++)
      {
	for (int j=0;j<latticePts_.dimension(1);j++)
	  {
	    dofCoords(i,j) = latticePts_(i,j);
	  }
      }
    return;
  }

}// namespace Intrepid
#endif
