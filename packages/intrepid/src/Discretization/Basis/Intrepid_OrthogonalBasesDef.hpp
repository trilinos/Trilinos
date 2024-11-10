/*
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
*/

namespace Intrepid {
  
  template<class Scalar>
  void OrthogonalBases::jrc(const Scalar &alpha , const Scalar &beta , 
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
  
  
  template<class Scalar, class ScalarArray1, class ScalarArray2>
  void OrthogonalBases::tabulateTriangle( const ScalarArray1& z ,
                                          const int n ,
                                          ScalarArray2 & poly_val )
  {
    const int np = z.dimension( 0 );

    // each point needs to be transformed from Pavel's element
    // z(i,0) --> (2.0 * z(i,0) - 1.0)
    // z(i,1) --> (2.0 * z(i,1) - 1.0)

    // set up constant term
    int idx_cur = OrthogonalBases::idxtri(0,0);
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
    idx_cur = OrthogonalBases::idxtri(1,0);
    for (int i=0;i<np;i++) {
      poly_val(idx_cur,i) = f1[i];
    }

    // recurrence in p
    for (int p=1;p<n;p++) {
      idx_cur = OrthogonalBases::idxtri(p,0);
      idx_curp1 = OrthogonalBases::idxtri(p+1,0);
      idx_curm1 = OrthogonalBases::idxtri(p-1,0);
      Scalar a = (2.0*p+1.0)/(1.0+p);
      Scalar b = p / (p+1.0);

      for (int i=0;i<np;i++) {
        poly_val(idx_curp1,i) = a * f1[i] * poly_val(idx_cur,i)
          - b * f3[i] * poly_val(idx_curm1,i);
      }
    }
    
    // D^{p,1}
    for (int p=0;p<n;p++) {
      int idxp0 = OrthogonalBases::idxtri(p,0);
      int idxp1 = OrthogonalBases::idxtri(p,1);
      for (int i=0;i<np;i++) {
        poly_val(idxp1,i) = poly_val(idxp0,i)
          *0.5*(1.0+2.0*p+(3.0+2.0*p)*(2.0*z(i,1)-1.0));
      }
    }

    // recurrence in q
    for (int p=0;p<n-1;p++) {
      for (int q=1;q<n-p;q++) {
        int idxpqp1=OrthogonalBases::idxtri(p,q+1);
        int idxpq=OrthogonalBases::idxtri(p,q);
        int idxpqm1=OrthogonalBases::idxtri(p,q-1);
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

  template<class Scalar, class ScalarArray1, class ScalarArray2>
  void OrthogonalBases::tabulateTetrahedron(const ScalarArray1 &z , 
                                            const int n ,
                                            ScalarArray2 &poly_val )
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
      f2[i] = pow( 0.5 * ( (2.0*z(i,1)-1.0) + (2.0*z(i,2)-1.0) ) , 2 );
      f3[i] = 0.5 * ( 1.0 + 2.0 * (2.0*z(i,1)-1.0) + (2.0*z(i,2)-1.0) );
      f4[i] = 0.5 * ( 1.0 - (2.0*z(i,2)-1.0) );
      f5[i] = f4[i] * f4[i];
    }
    
    // constant term
    idxcur = idxtet(0,0,0);
    for (int i=0;i<np;i++) {
      poly_val(idxcur,i) = 1.0;
    }

    // D^{1,0,0}
    idxcur = idxtet(1,0,0);
    for (int i=0;i<np;i++) {
      poly_val(idxcur,i) = f1[i];
    }

    // p recurrence
    for (int p=1;p<n;p++) {
      Scalar a1 = (2.0 * p + 1.0) / ( p + 1.0);
      Scalar a2 = p / ( p + 1.0 );
      int idxp = idxtet(p,0,0);
      int idxpp1 = idxtet(p+1,0,0);
      int idxpm1 = idxtet(p-1,0,0);
      //cout << idxpm1 << " " << idxp << " " << idxpp1 << endl;
      for (int i=0;i<np;i++) {
        poly_val(idxpp1,i) = a1 * f1[i] * poly_val(idxp,i) - a2 * f2[i] * poly_val(idxpm1,i);
      }
    }
    // q = 1
    for (int p=0;p<n;p++) {
      int idx0 = idxtet(p,0,0);
      int idx1 = idxtet(p,1,0);
      for (int i=0;i<np;i++) {
        poly_val(idx1,i) = poly_val(idx0,i) * ( p * ( 1.0 + (2.0*z(i,1)-1.0) ) + 0.5 * ( 2.0 + 3.0 * (2.0*z(i,1)-1.0) + (2.0*z(i,2)-1.0) ) );
      }
    }

    // q recurrence
    for (int p=0;p<n-1;p++) {
      for (int q=1;q<n-p;q++) {
        Scalar aq,bq,cq;
        jrc((Scalar)(2.0*p+1.0),(Scalar)(0),q,aq,bq,cq);
        int idxpqp1 = idxtet(p,q+1,0);
        int idxpq = idxtet(p,q,0);
        int idxpqm1 = idxtet(p,q-1,0);
        for (int i=0;i<np;i++) {
          poly_val(idxpqp1,i) = ( aq * f3[i] + bq * f4[i] ) * poly_val(idxpq,i) 
            - ( cq * f5[i] ) * poly_val(idxpqm1,i);
        }
      }
    }
    
    // r = 1
    for (int p=0;p<n;p++) {
      for (int q=0;q<n-p;q++) {
        int idxpq1 = idxtet(p,q,1);
        int idxpq0 = idxtet(p,q,0);
        for (int i=0;i<np;i++) {
          poly_val(idxpq1,i) = poly_val(idxpq0,i) * ( 1.0 + p + q + ( 2.0 + q + p ) * (2.0*z(i,2)-1.0) );
        }
      }
    }
    
    // general r recurrence
    for (int p=0;p<n-1;p++) {
      for (int q=0;q<n-p-1;q++) {
        for (int r=1;r<n-p-q;r++) {
          Scalar ar,br,cr;
          int idxpqrp1 = idxtet(p,q,r+1);
          int idxpqr = idxtet(p,q,r);
          int idxpqrm1 = idxtet(p,q,r-1);
          jrc(2.0*p+2.0*q+2.0,0.0,r,ar,br,cr);
          for (int i=0;i<np;i++) {
            poly_val(idxpqrp1,i) = (ar * (2.0*z(i,2)-1.0) + br) * poly_val( idxpqr , i ) - cr * poly_val(idxpqrm1,i);
          }
        }
      }
    }
    
    return;
    
  }

} // namespace Intrepid;

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

