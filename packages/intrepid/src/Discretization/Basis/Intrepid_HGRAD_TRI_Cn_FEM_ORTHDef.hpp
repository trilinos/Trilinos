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

/** \file   Intrepid_HGRAD_TRI_Cn_FEM_ORTHDef.hpp
    \brief  Definition file for FEM orthogonal basis functions of arbitrary degree 
            for H(grad) functions on TRI.
    \author Created by R. Kirby
 */

namespace Intrepid {
  
template<class Scalar, class ArrayScalar>
Basis_HGRAD_TRI_Cn_FEM_ORTH<Scalar,ArrayScalar>::Basis_HGRAD_TRI_Cn_FEM_ORTH( int degree )
{
    this -> basisCardinality_  = (degree+1)*(degree+2)/2;
    this -> basisDegree_       = degree;
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3> >() );
    this -> basisType_         = BASIS_FEM_HIERARCHICAL;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;
}
  
  
  
template<class Scalar, class ArrayScalar>
void Basis_HGRAD_TRI_Cn_FEM_ORTH<Scalar, ArrayScalar>::initializeTags() {
  
  // Basis-dependent initializations
  int tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
  int posScDim = 0;        // poisition in the tag, counting from 0, of the subcell dim 
  int posScOrd = 1;        // poisition in the tag, counting from 0, of the subcell ordinal
  int posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
  
  // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 
  int *tags = new int[tagSize * this->getCardinality()];
  for (int i=0;i<this->getCardinality();i++) {
    tags[4*i] = 2;
    tags[4*i+1] = 0;
    tags[4*i+2] = i;
    tags[4*i+3] = this->getCardinality();
  }
  
  // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
  Intrepid::setOrdinalTagData(this -> tagToOrdinal_,
                              this -> ordinalToTag_,
                              tags,
                              this -> basisCardinality_,
                              tagSize,
                              posScDim,
                              posScOrd,
                              posDfOrd);
}  



template<class Scalar, class ArrayScalar> 
void Basis_HGRAD_TRI_Cn_FEM_ORTH<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
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
  // need to figure out AD/Sacado and put in tabulation.
}

template<class Scalar, class ArrayScalar>
void Basis_HGRAD_TRI_Cn_FEM_ORTH<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
								 const ArrayScalar &    inputPoints,
								 const ArrayScalar &    cellVertices,
								 const EOperator        operatorType) const {
  TEST_FOR_EXCEPTION( (true), std::logic_error,
                      ">>> ERROR (Basis_HGRAD_TRI_C1_FEM): FEM Basis calling an FVD member function");
}

template<class Scalar, class ArrayScalar>
void Basis_HGRAD_TRI_Cn_FEM_ORTH<Scalar,ArrayScalar>::tabulate( const ArrayScalar& z ,
								const int n ,
								ArrayScalar & poly_val )
{
  const int np = z.dimension( 0 );
  
  // each point needs to be transformed from Pavel's element
  // z(i,0) --> (2.0 * z(i,0) - 1.0)
  // z(i,1) --> (2.0 * z(i,1) - 1.0)
  
  // set up constant term
  int idx_cur = Basis_HGRAD_TRI_Cn_FEM_ORTH::idx(0,0);
  int idx_curp1,idx_curm1;
  
  // set D^{0,0} = 1.0
  for (int i=0;i<np;i++) {
    poly_val(idx_cur,i) = 1.0;
  }
  
  Teuchos::Array<Scalar> f1(np),f2(np),f3(np);
  
  for (int i=0;i<np;i++) {
    f1[i] = 0.5 * (1.0+2.0*(2.0*z(i,0)-1.0)+(2.0*z(i,1)-1.0));
    f2[i] = 0.5 * (1.0-(2.0*z(i,1)-1.0));
    f3[i] = f2[i] * f2[i];
  }
  
  // set D^{1,0} = f1
  idx_cur = Basis_HGRAD_TRI_Cn_FEM_ORTH::idx(1,0);
  for (int i=0;i<np;i++) {
    poly_val(idx_cur,i) = f1[i];
  }
  
  // recurrence in p
  for (int p=1;p<n;p++) {
    idx_cur = Basis_HGRAD_TRI_Cn_FEM_ORTH::idx(p,0);
    idx_curp1 = Basis_HGRAD_TRI_Cn_FEM_ORTH::idx(p+1,0);
    idx_curm1 = Basis_HGRAD_TRI_Cn_FEM_ORTH::idx(p-1,0);
    Scalar a = (2.0*p+1.0)/(1.0+p);
    Scalar b = p / (p+1.0);
    
    for (int i=0;i<np;i++) {
      poly_val(idx_curp1,i) = a * f1[i] * poly_val(idx_cur,i)
	- b * f3[i] * poly_val(idx_curm1,i);
    }
  }
  
  // D^{p,1}
  for (int p=0;p<n;p++) {
    int idxp0 = Basis_HGRAD_TRI_Cn_FEM_ORTH::idx(p,0);
    int idxp1 = Basis_HGRAD_TRI_Cn_FEM_ORTH::idx(p,1);
    for (int i=0;i<np;i++) {
      poly_val(idxp1,i) = poly_val(idxp0,i)
	*0.5*(1.0+2.0*p+(3.0+2.0*p)*(2.0*z(i,1)-1.0));
    }
  }
  
  
  // recurrence in q
  for (int p=0;p<n-1;p++) {
    for (int q=1;q<n-p;q++) {
      int idxpqp1=Basis_HGRAD_TRI_Cn_FEM_ORTH::idx(p,q+1);
      int idxpq=Basis_HGRAD_TRI_Cn_FEM_ORTH::idx(p,q);
      int idxpqm1=Basis_HGRAD_TRI_Cn_FEM_ORTH::idx(p,q-1);
      Scalar a,b,c;
      jrc((Scalar)(2*p+1),(Scalar)0,q,a,b,c);
      for (int i=0;i<np;i++) {
	poly_val(idxpqp1,i)
	  = (a*(2.0*z(i,1)-1.0)+b)*poly_val(idxpq,i)
	  - c*poly_val(idxpqm1,i);
      }
    }
  }
  
  return;
}
  
template<class Scalar, class ArrayScalar>
int Basis_HGRAD_TRI_Cn_FEM_ORTH<Scalar,ArrayScalar>::idx(int p , int q)
{
  return (p+q)*(p+q+1)/2+q;
}

template<class Scalar, class ArrayScalar>
void Basis_HGRAD_TRI_Cn_FEM_ORTH<Scalar,ArrayScalar>::jrc( const Scalar &alpha , const Scalar &beta , 
							   const int &n ,
							   Scalar &an , Scalar &bn, Scalar &cn )
{
  an = (2.0 * n + 1.0 + alpha + beta) * ( 2.0 * n + 2.0 + alpha + beta ) 
    / ( 2.0 * ( n + 1 ) * ( n + 1 + alpha + beta ) );
  bn = (alpha*alpha-beta*beta)*(2.0*n+1.0+alpha+beta) 
    / ( 2.0*(n+1.0)*(2.0*n+alpha+beta)*(n+1.0+alpha+beta) );
  cn = (n+alpha)*(n+beta)*(2.0*n+2.0+alpha+beta) 
    / ( (n+1.0)*(n+1.0+alpha+beta)*(2.0*n+alpha+beta) );
  
  return;
}



}// namespace Intrepid
