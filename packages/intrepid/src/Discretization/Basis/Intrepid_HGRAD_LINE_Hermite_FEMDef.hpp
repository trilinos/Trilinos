#ifndef INTREPID_HGRAD_LINE_HERMITE_FEMDEF_HPP
#define INTREPID_HGRAD_LINE_HERMITE_FEMDEF_HPP
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

/** \file   Intrepid_HGRAD_LINE_Hermite_FEMDef.hpp
    \brief  Definition file for Hermite FEM basis functions of degree 2n for H(grad) functions on a Line.
    \author Created by G. von Winckel
 */

#include<array>
#include<iostream>
#include<iomanip>

namespace Intrepid {


// No-arg constructor uses cubic Hermite interpolants based at the cell vertices
template<class Scalar, class ArrayScalar>
Basis_HGRAD_LINE_Hermite_FEM<Scalar,ArrayScalar>::Basis_HGRAD_LINE_Hermite_FEM() : 
    latticePts_(2,1) {
    this->basisCardinality_  = 4;    // Four basis functions
    this->basisDegree_       = 3;    // Cubic polynomials
    this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >() );
    this->basisType_         = BASIS_FEM_DEFAULT;
    this->basisCoordinates_  = COORDINATES_CARTESIAN;
    this->basisTagsAreSet_   = false;

    latticePts_(0,0) = -1.0;
    latticePts_(1,0) =  1.0;
 
    setupVandermonde();

} // no-arg constructor



// Constructor with points as argument
template<class Scalar, class ArrayScalar>
Basis_HGRAD_LINE_Hermite_FEM<Scalar,ArrayScalar>::Basis_HGRAD_LINE_Hermite_FEM( const ArrayScalar &pts) : 
    latticePts_( pts.dimension(0), 1 ) {

    int n = pts.dimension(0);

    this->basisCardinality_  = 2*n;                   
    this->basisDegree_       = 2*n-1;    
    this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >() );
    this->basisType_         = BASIS_FEM_DEFAULT;
    this->basisCoordinates_  = COORDINATES_CARTESIAN;
    this->basisTagsAreSet_   = false;

    for( int i=0; i<n-1; ++i ) {
      TEUCHOS_TEST_FOR_EXCEPTION( pts(i,0) >= pts(i+1,0), std::runtime_error ,
        "Intrepid::Basis_HGRAD_LINE_Hermite_FEM Illegal points given to constructor" );
    }  

    // copy points int latticePts, correcting endpoints if needed
    if (std::abs(pts(0,0)+1.0) < INTREPID_TOL) {
      latticePts_(0,0) = -1.0;
    }
    else {
      latticePts_(0,0) = pts(0,0);
    }
    for (int i=1;i<n-1;i++) {
      latticePts_(i,0) = pts(i,0);
    }
    if (std::abs(pts(n-1,0)-1.0) < INTREPID_TOL) {
      latticePts_(n-1,0) = 1.0;
    }
    else {
      latticePts_(n-1,0) = pts(n-1,0);
    }  

    setupVandermonde();

} // Constructor with points given
 

// Constructor with point type as argument
template<class Scalar, class ArrayScalar>
Basis_HGRAD_LINE_Hermite_FEM<Scalar,ArrayScalar>::Basis_HGRAD_LINE_Hermite_FEM( const int n, 
                                                                                const EPointType &pointType ) : 
  latticePts_(n,1) {

  TEUCHOS_TEST_FOR_EXCEPTION(n<2,std::invalid_argument,"Intrepid::Basis_HGRAD_LINE_Hermite_FEM requires the "
    "number of interpolation points to be at least 2");

  this->basisCardinality_  = 2*n;                  
  this->basisDegree_       = 2*n-1;    
  this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >() );
  this->basisType_         = BASIS_FEM_DEFAULT;
  this->basisCoordinates_  = COORDINATES_CARTESIAN;
  this->basisTagsAreSet_   = false;

  switch(pointType) {
    case POINTTYPE_EQUISPACED:
      PointTools::getLattice<Scalar,FieldContainer<Scalar> >( latticePts_ ,  
        this->basisCellTopology_ , n-1 , 0 , POINTTYPE_EQUISPACED );
    break;
  case POINTTYPE_SPECTRAL: 
    PointTools::getLattice<Scalar,FieldContainer<Scalar> >( latticePts_ ,  
        this->basisCellTopology_ , n-1 , 0 , POINTTYPE_WARPBLEND );
    break;
  case POINTTYPE_SPECTRAL_OPEN: 
    PointTools::getGaussPoints<Scalar,FieldContainer<Scalar> >( latticePts_ , n-1 );
    break;
  default:
    TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument , 
        "Basis_HGRAD_LINE_Hermite_FEM:: invalid point type" );
    break;
  }

    setupVandermonde();

} // Constructor with point type given



template<class Scalar, class ArrayScalar>
void Basis_HGRAD_LINE_Hermite_FEM<Scalar,ArrayScalar>::setupVandermonde( bool factor  ) {

    initializeTags();

    int nBf = this->getCardinality();
    int n   = nBf/2;

    V_.shape(nBf,nBf);

    // Make containers to store the Legendre polynomials and their derivatives
    // at a given point 
    ArrayScalar P ( nBf );
    ArrayScalar Px( nBf );

    // Loop over grid points
    for( int i=0; i<n; ++i ) {

      recurrence(P,Px,latticePts_(i,0));
      
      // Loop over basis functions     
      for(int j=0; j<nBf; ++j ) {
        V_(j, i  ) = P (j);
        V_(j, i+n) = Px(j); 
      }
    }

    solver_.setMatrix(Teuchos::rcpFromRef(V_));

    if(factor) {
      solver_.factorWithEquilibration(true);
      solver_.factor();
      isFactored_ = true;
    }
    else {
      isFactored_ = false;
    }

}


template<class Scalar, class ArrayScalar>
void Basis_HGRAD_LINE_Hermite_FEM<Scalar,ArrayScalar>::recurrence( ArrayScalar &P,
                                                                   ArrayScalar &Px,
                                                                   const Scalar x ) const {

  int    n = P.dimension(0);
  Scalar q = x*x-1.0;

  P (0) = 1.0;
  Px(0) = 0.0;
  
  // Loop over basis indices
  for( int j=0; j<n-1; ++j ) {
    P (j+1) =     x*P(j) + q*Px(j)/(j+1); // Compute \f$P_{j+1}(x_i)\f$
    Px(j+1) = (j+1)*P(j) + x*Px(j);       // Compute \f$P'_{j+1}(x_i)\f$
  }

} // recurrence()


// Computes the derivatives of Legendre polynomials
template<class Scalar, class ArrayScalar>
void Basis_HGRAD_LINE_Hermite_FEM<Scalar,ArrayScalar>::legendre_d( ArrayScalar &P,
                                                                   ArrayScalar &Px, 
                                                                   const int m, 
                                                                   const Scalar x ) const {
  // Compute P,P'
  recurrence(P,Px,x);

  int C = this->getCardinality();

  // Loop over derivative orders
  for( int k=1;k<m;++k) {
    P = Px;

    // Loop over polynomial indices
    for( int j=0; j<C; ++j ) {
 
     if( j<k ) {
        Px(j) = 0;         
      }

      else {
        Px(j) = (j+k)*P(j-1) + x*Px(j-1);
      }
    }
  }

} // legendre_d()


template<class Scalar, class ArrayScalar> 
void Basis_HGRAD_LINE_Hermite_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
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
  // Number of evaluation points = dim 0 of inputPoints
  int nPts = inputPoints.dimension(0);  
  int nBf  = this->getCardinality();  

  int n = nBf/2;

  // Legendre polynomials and their derivatives evaluated on inputPoints
  SerialDenseMatrix legendre(nBf,nPts); 

  // Hermite interpolants evaluated on inputPoints
  SerialDenseMatrix hermite(nBf,nPts);
    
  ArrayScalar P (nBf);
  ArrayScalar Px(nBf);
         
  int derivative_order;
  int derivative_case = static_cast<int>(operatorType);

  if( derivative_case == 0 ) {
    derivative_order = 0;
  }
  else if( derivative_case > 0 && derivative_case < 5 ) {
    derivative_order = 1;
  }
  else {
    derivative_order = derivative_case - 3;
  }
  
  try {
    // GRAD,CURL,DIV, and D1 are all the first derivative
    switch (operatorType) {
      case OPERATOR_VALUE: 
      { 
        for( int i=0; i<nPts; ++i ) { 
          recurrence( P, Px, inputPoints(i,0) );
          for( int j=0; j<nBf; ++j ) {
             legendre(j,i) = P(j);
          }
        }
        break; 
      }
      case OPERATOR_GRAD:
      case OPERATOR_DIV:
      case OPERATOR_CURL:
      case OPERATOR_D1:
      {
        for( int i=0; i<nPts; ++i ) { 
          recurrence( P, Px, inputPoints(i,0) );
          for( int j=0; j<nBf; ++j ) {
             legendre(j,i) = Px(j);
          }
        }
        break;
      }  
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
        for( int i=0; i<nPts; ++i ) {
          legendre_d( P, Px, derivative_order, inputPoints(i,0));
          for( int j=0; j<nBf; ++j ) {
            legendre(j,i) = Px(j); 
          }
        }
        break;
      }
      default:
        TEUCHOS_TEST_FOR_EXCEPTION( !( Intrepid::isValidOperator(operatorType) ), std::invalid_argument,
                            ">>> ERROR (Basis_HGRAD_LINE_Hermite_FEM): Invalid operator type");

    } // switch(operatorType)
  }
  catch (std::invalid_argument &exception){
    TEUCHOS_TEST_FOR_EXCEPTION( true , std::invalid_argument,
                        ">>> ERROR (Basis_HGRAD_LINE_Hermite_FEM): Operator failed");    
  }

  if( !isFactored_ ) {
    solver_.factorWithEquilibration(true);
    solver_.factor();
  }

  solver_.setVectors(Teuchos::rcpFromRef(hermite),Teuchos::rcpFromRef(legendre));
  solver_.solve();

  if(derivative_order > 0)
  {
    for( int i=0; i<n; ++i ) {
      for( int j=0; j<nPts; ++j ) {
        outputValues(2*i,  j,0)  = hermite(i,  j);
        outputValues(2*i+1,j,0)  = hermite(i+n,j);
      } 
    }
  }
  else {
    for( int i=0; i<n; ++i ) {
      for( int j=0; j<nPts; ++j ) {
        outputValues(2*i  ,j)   = hermite(i,  j);
        outputValues(2*i+1,j)   = hermite(i+n,j);
      } 
    }
  }
     
} // getValues()


template<class Scalar, class ArrayScalar>
void Basis_HGRAD_LINE_Hermite_FEM<Scalar, ArrayScalar>::initializeTags() {

  // Basis-dependent intializations
  int tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
  int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
  int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
  int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
    
  // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 

  int C = this->getCardinality();
  tags_.reserve( tagSize * C );

  int n = C/2;

  int hasLeftVertex  = static_cast<int>( latticePts_(0  , 0)  == -1.0 );
  int hasRightVertex = static_cast<int>( latticePts_(n-1, 0)  ==  1.0 );

  int internal_dof = C - 2*(hasLeftVertex+hasRightVertex);

  if( hasLeftVertex ) {

    // Value interpolant 
    tags_[0] = 0;               //  this is a vertex (interval end point)
    tags_[1] = 0;               //  this is the first subcell
    tags_[2] = 0;               //  this is the first DoF for this vertex
    tags_[3] = 2;               //  this vertex has 2 DoF

    // Derivative interpolant 
    tags_[4] = 0;                //  this is a vertex (interval end point)
    tags_[5] = 0;                //  this is the first subcell  
    tags_[6] = 1;                //  this is the second DoF for this vertex
    tags_[7] = 2;                //  this vertex has 2 DoF

  }
  else { // no left vertex

    // Value interpolant 
    tags_[0] = 1;                     // this point is on a line  
    tags_[1] = 0;                     // no subcells
    tags_[2] = 0;                     // this is the first DoF for this line
    tags_[3] = C-2*hasRightVertex;    // this cell has 2n DoF

    // Derivative interpolant 
    tags_[4] = 1;                     // this point is on a line  
    tags_[5] = 0;                     // no subcells
    tags_[6] = 1;                     // this is the second DoF for this line
    tags_[7] = C-2*hasRightVertex;    // this cell has 2n DoF
  } 

  if( hasRightVertex ) {
    int i0 = C-2; 
    int i1 = C-1;
    
    // Value interpolant 
    tags_[4*i0  ] = 0;                      
    tags_[4*i0+1] = hasLeftVertex;        
    tags_[4*i0+2] = 0;                          
    tags_[4*i0+3] = 2;                         

    // Derivative interpolant 
    tags_[4*i1  ] = 0;   
    tags_[4*i1+1] = hasLeftVertex; 
    tags_[4*i1+2] = 1;   
    tags_[4*i1+3] = 2;   
  }
  else { // no right vertex 
    int i0 = C-2; 
    int i1 = C-1;
    
    // Value interpolant 
    tags_[4*i0  ] = 1;                      
    tags_[4*i0+1] = 0;        
    tags_[4*i0+2] = internal_dof-2;
    tags_[4*i0+3] = internal_dof;                         

    // Derivative interpolant 
    tags_[4*i1  ] = 1;      
    tags_[4*i1+1] = 0;      
    tags_[4*i1+2] = internal_dof-1;      
    tags_[4*i1+3] = internal_dof;   
  }

  for( int i=1; i<n-1; ++i ) {
    int i0 = 2*i; int i1 = 2*i+1;

    // Value interpolant 
    tags_[4*i0  ] = 1;  // Points on a line (1 dimensional)
    tags_[4*i0+1] = 0;
    tags_[4*i0+2] = i0 - 2*hasLeftVertex;
    tags_[4*i0+3] = internal_dof;

    // Derivative interpolant 
    tags_[4*i1  ] = 1;
    tags_[4*i1+1] = 0;
    tags_[4*i1+2] = i1 - 2*hasLeftVertex;
    tags_[4*i1+3] = internal_dof;
  }

  // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
  Intrepid::setOrdinalTagData(this -> tagToOrdinal_,
                              this -> ordinalToTag_,
                              tags_.data(),
                              this -> basisCardinality_,
                              tagSize,
                              posScDim,
                              posScOrd,
                              posDfOrd);

}


template<class Scalar, class ArrayScalar> 
void Basis_HGRAD_LINE_Hermite_FEM<Scalar, ArrayScalar>::printTags( std::ostream &os ) {

  int nBf = this->getCardinality();

  os << "Tags:" << std::endl;
  os << "-----" << std::endl;


  os << "Index:        ";
  for( int i=0; i<nBf; ++i ) {
    os << std::setw(4) << i; 
  }
  os << std::endl;


  os << "Subcell dim:  ";
  for( int i=0; i<nBf; ++i ) {
    os << std::setw(4) << tags_[4*i]; 
  }
  os << std::endl;

  
  os << "Subcell ord:  ";
  for( int i=0; i<nBf; ++i ) {
    os << std::setw(4) << tags_[4*i+1]; 
  }
  os << std::endl;


  os << "Subcell DoF:  ";
  for( int i=0; i<nBf; ++i ) {
    os << std::setw(4) << tags_[4*i+2]; 
  }
  os << std::endl;


  os << "Total Sc DoF: ";
  for( int i=0; i<nBf; ++i ) {
    os << std::setw(4) << tags_[4*i+3]; 
  }
  os << std::endl;
  os << std::endl;

}




template<class Scalar, class ArrayScalar>
void Basis_HGRAD_LINE_Hermite_FEM<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                                  const ArrayScalar &    inputPoints,
                                                                  const ArrayScalar &    cellVertices,
                                                                  const EOperator        operatorType) const {
  TEUCHOS_TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (Basis_HGRAD_LINE_Hermite_FEM): FEM Basis calling an FVD member function");
}


template<class Scalar, class ArrayScalar>
void Basis_HGRAD_LINE_Hermite_FEM<Scalar,ArrayScalar>::getDofCoords( ArrayScalar & dofCoords ) const
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

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

