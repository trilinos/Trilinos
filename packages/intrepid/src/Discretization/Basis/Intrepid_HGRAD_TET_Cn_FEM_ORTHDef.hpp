#ifndef INTREPID_HGRAD_TET_CN_FEM_ORTHDEF_HPP
#define INTREPID_HGRAD_TET_CN_FEM_ORTHDEF_HPP
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

/** \file   Intrepid_HGRAD_TET_Cn_FEM_ORTHDef.hpp
    \brief  Definition file for FEM orthogonal basis functions of arbitrary degree 
    for H(grad) functions on TET.
    \author Created by R. Kirby
*/

namespace Intrepid {

  template<class Scalar, class ArrayScalar>
  Basis_HGRAD_TET_Cn_FEM_ORTH<Scalar,ArrayScalar>::Basis_HGRAD_TET_Cn_FEM_ORTH( int degree )
  {
    this -> basisCardinality_  = (degree+1)*(degree+2)*(degree+3)/6;
    this -> basisDegree_       = degree;
    this -> basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<4> >() );
    this -> basisType_         = BASIS_FEM_HIERARCHICAL;
    this -> basisCoordinates_  = COORDINATES_CARTESIAN;
    this -> basisTagsAreSet_   = false;
  }
  
  
  
  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_TET_Cn_FEM_ORTH<Scalar, ArrayScalar>::initializeTags() {
  
    // Basis-dependent initializations
    int tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
    int posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
    int posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
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
  void Basis_HGRAD_TET_Cn_FEM_ORTH<Scalar, ArrayScalar>::getValues(ArrayScalar &        outputValues,
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
    const int deg = this->getDegree();
  
    switch (operatorType) {
    case OPERATOR_VALUE:
      {
        TabulatorTet<Scalar,ArrayScalar,0>::tabulate( outputValues ,
                                                      deg ,
                                                      inputPoints );
      }
      break;
    case OPERATOR_GRAD:
    case OPERATOR_D1:
      {
        TabulatorTet<Scalar,ArrayScalar,1>::tabulate( outputValues ,
                                                      deg ,
                                                      inputPoints );
      }
      break;
    default:
      TEST_FOR_EXCEPTION( true , std::invalid_argument,
                          ">>> ERROR (Basis_HGRAD_TET_Cn_FEM_ORTH): invalid or unsupported operator" );
    }

    return;
  }
  
  template<class Scalar, class ArrayScalar>
  void Basis_HGRAD_TET_Cn_FEM_ORTH<Scalar, ArrayScalar>::getValues(ArrayScalar&           outputValues,
                                                                  const ArrayScalar &    inputPoints,
                                                                  const ArrayScalar &    cellVertices,
                                                                  const EOperator        operatorType) const {
    TEST_FOR_EXCEPTION( (true), std::logic_error,
                        ">>> ERROR (Basis_HGRAD_TET_Cn_FEM_ORTH): FEM Basis calling an FVD member function");
  }

  template<class Scalar, class ArrayScalar>
  void TabulatorTet<Scalar,ArrayScalar,0>::tabulate( ArrayScalar &outputValues ,
                                                    const int deg ,
                                                    const ArrayScalar &z )
  {
    const int np = z.dimension( 0 );
    int idxcur;
  
    // each point needs to be transformed from Pavel's element
    // z(i,0) --> (2.0 * z(i,0) - 1.0)
    // z(i,1) --> (2.0 * z(i,1) - 1.0)
    // z(i,2) --> (2.0 * z(i,2) - 1.0)
  
    Teuchos::Array<Scalar> f1(np),f2(np),f3(np),f4(np),f5(np);
  
    for (int i=0;i<np;i++) {
      f1[i] = 0.5 * ( 2.0 + 2.0*(2.0*z(i,0)-1.0) + (2.0*z(i,1)-1.0) + (2.0*z(i,2)-1.0) );
      Scalar foo =  0.5 * ( (2.0*z(i,1)-1.0) + (2.0*z(i,2)-1.0) );
      f2[i] = foo * foo;
      f3[i] = 0.5 * ( 1.0 + 2.0 * (2.0*z(i,1)-1.0) + (2.0*z(i,2)-1.0) );
      f4[i] = 0.5 * ( 1.0 - (2.0*z(i,2)-1.0) );
      f5[i] = f4[i] * f4[i];
    }

    // constant term
    idxcur = TabulatorTet<Scalar,ArrayScalar,0>::idx(0,0,0);
    for (int i=0;i<np;i++) {
      outputValues(idxcur,i) = 1.0 + z(i,0) - z(i,0) + z(i,1) - z(i,1) + z(i,2) - z(i,2);
    }
  
    if (deg > 0) {

      // D^{1,0,0}
      idxcur = TabulatorTet<Scalar,ArrayScalar,0>::idx(1,0,0);
      for (int i=0;i<np;i++) {
        outputValues(idxcur,i) = f1[i];
      }
  
      // p recurrence
      for (int p=1;p<deg;p++) {
        Scalar a1 = (2.0 * p + 1.0) / ( p + 1.0);
        Scalar a2 = p / ( p + 1.0 );
        int idxp = TabulatorTet<Scalar,ArrayScalar,0>::idx(p,0,0);
        int idxpp1 = TabulatorTet<Scalar,ArrayScalar,0>::idx(p+1,0,0);
        int idxpm1 = TabulatorTet<Scalar,ArrayScalar,0>::idx(p-1,0,0);
        for (int i=0;i<np;i++) {
          outputValues(idxpp1,i) = a1 * f1[i] * outputValues(idxp,i) - a2 * f2[i] * outputValues(idxpm1,i);
        }
      }
      // q = 1
      for (int p=0;p<deg;p++) {
        int idx0 = TabulatorTet<Scalar,ArrayScalar,0>::idx(p,0,0);
        int idx1 = TabulatorTet<Scalar,ArrayScalar,0>::idx(p,1,0);
        for (int i=0;i<np;i++) {
          outputValues(idx1,i) = outputValues(idx0,i) * ( p * ( 1.0 + (2.0*z(i,1)-1.0) ) +
                                                          0.5 * ( 2.0 + 3.0 * (2.0*z(i,1)-1.0) + (2.0*z(i,2)-1.0) ) );
        }
      }
  
      // q recurrence
      for (int p=0;p<deg-1;p++) {
        for (int q=1;q<deg-p;q++) {
          Scalar aq,bq,cq;

	  TabulatorTet<Scalar,ArrayScalar,0>::jrc(2.0*p+1.0 ,0 ,q, aq, bq, cq);
          int idxpqp1 = TabulatorTet<Scalar,ArrayScalar,0>::idx(p,q+1,0);
          int idxpq = TabulatorTet<Scalar,ArrayScalar,0>::idx(p,q,0);
          int idxpqm1 = TabulatorTet<Scalar,ArrayScalar,0>::idx(p,q-1,0);
          for (int i=0;i<np;i++) {
            outputValues(idxpqp1,i) = ( aq * f3[i] + bq * f4[i] ) * outputValues(idxpq,i) 
              - ( cq * f5[i] ) * outputValues(idxpqm1,i);
          }
        }
      }
  
      // r = 1
      for (int p=0;p<deg;p++) {
        for (int q=0;q<deg-p;q++) {
          int idxpq1 = TabulatorTet<Scalar,ArrayScalar,0>::idx(p,q,1);
          int idxpq0 = TabulatorTet<Scalar,ArrayScalar,0>::idx(p,q,0);
          for (int i=0;i<np;i++) {
            outputValues(idxpq1,i) = outputValues(idxpq0,i) * ( 1.0 + p + q + ( 2.0 + q + 
                                                                                p ) * (2.0*z(i,2)-1.0) );
          }
        }
      }
      // general r recurrence
      for (int p=0;p<deg-1;p++) {
        for (int q=0;q<deg-p-1;q++) {
          for (int r=1;r<deg-p-q;r++) {
            Scalar ar,br,cr;
            int idxpqrp1 = TabulatorTet<Scalar,ArrayScalar,0>::idx(p,q,r+1);
            int idxpqr = TabulatorTet<Scalar,ArrayScalar,0>::idx(p,q,r);
            int idxpqrm1 = TabulatorTet<Scalar,ArrayScalar,0>::idx(p,q,r-1);
            jrc(2.0*p+2.0*q+2.0, 0.0, r, ar, br, cr);
            for (int i=0;i<np;i++) {
              outputValues(idxpqrp1,i) = (ar * (2.0*z(i,2)-1.0) + br) * outputValues( idxpqr , i ) - cr * outputValues(idxpqrm1,i);
            }
          }
        }
      }

    }  
    // normalize
    for (int p=0;p<=deg;p++) {
      for (int q=0;q<=deg-p;q++) {
        for (int r=0;r<=deg-p-q;r++) {
          int idxcur = TabulatorTet<Scalar,ArrayScalar,0>::idx(p,q,r);
          Scalar scal = sqrt( (p+0.5)*(p+q+1.0)*(p+q+r+1.5) );
          for (int i=0;i<np;i++) {
            outputValues(idxcur,i) *= scal;
          }
        }
      }
    }
  
    return;
  
  }


  template<typename Scalar, typename ArrayScalar>
  void TabulatorTet<Scalar,ArrayScalar,1>::tabulate( ArrayScalar &outputValues ,
                                                    const int deg ,
                                                    const ArrayScalar &z ) 
  {
    const int np = z.dimension(0);
    const int card = outputValues.dimension(0);
    FieldContainer<Sacado::Fad::DFad<Scalar> > dZ( z.dimension(0) , z.dimension(1) );
    for (int i=0;i<np;i++) {
      for (int j=0;j<3;j++) {
        dZ(i,j) = Sacado::Fad::DFad<Scalar>( z(i,j) );
        dZ(i,j).diff(j,3);
      }
    }
    FieldContainer<Sacado::Fad::DFad<Scalar> > dResult(card,np);

    TabulatorTet<Sacado::Fad::DFad<Scalar>,FieldContainer<Sacado::Fad::DFad<Scalar> >,0>::tabulate( dResult ,
                                                                                                    deg ,
                                                                                                    dZ );

    for (int i=0;i<card;i++) {
      for (int j=0;j<np;j++) {
        for (int k=0;k<3;k++) {
          outputValues(i,j,k) = dResult(i,j).dx(k);
        }
      }
    }

    return;

  }


}// namespace Intrepid


#endif
