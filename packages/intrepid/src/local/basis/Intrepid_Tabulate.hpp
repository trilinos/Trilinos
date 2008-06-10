#ifndef INTREPID_RECURRENCE_BASIS_HPP
#define INTREPID_RECURRENCE_BASIS_HPP
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_Array.hpp"

using std::cout;
using std::endl;

namespace Intrepid {

  template<class Scalar>
  class OrthogonalExpansions {
  private:
    static void jrc( const Scalar &a , const Scalar &b , const int &n ,
	      Scalar &aj , Scalar &bj, Scalar &cj )
    {
      aj = (2.0 * n + 1.0 + a + b) * ( 2.0 * n + 2.0 + a + b ) 
	/ ( 2.0 * ( n + 1 ) * ( n + 1 + a + b ) );
      bj = (a*a-b*b)*(2.0*n+1.0+a+b) / ( 2.0*(n+1.0)*(2.0*n+a+b)*(n+1.0+a+b) );
      cj = (n+a)*(n+b)*(2.0*n+2.0+a+b) / ( (n+1.0)*(n+1.0+a+b)*(2.0*n+a+b) );
      return;
    }
    
    static inline int idx2d(int p, int q)
    {
      return (p+q)*(p+q+1)/2+q;
    }
    
    static inline int idx3d(const int &p,
			    const int &q,
			    const int &r)
    {
      return (p+q+r)*(p+q+r+1)*(p+q+r+2)/6+(q+r)*(q+r+1)/2+r;
    }
    
    static void tabulate_triangle( const int &n , 
			    const Teuchos::Array<Point<Scalar> > &xin , 
			    Teuchos::SerialDenseMatrix<int,Scalar> &results )
    {
      const int num_pts = xin.size();

      // transformation from Pavel's element
      Teuchos::Array<Point<Scalar> > x;
      Scalar xcur_arr[2];

      for (int i=0;i<num_pts;i++) {
	xcur_arr[0] = 2.0 * xin[i][0] - 1.0;
	xcur_arr[1] = 2.0 * xin[i][1] - 1.0;
	Point<Scalar> xcur( xcur_arr[0],xcur_arr[1] );
        x.append( xcur );
      }

      // set up constant term
      int idx_cur = idx2d(0,0);
      int idx_curp1,idx_curm1;

      
      // set D^{0,0} = 1.0
      for (int i=0;i<num_pts;i++) {
	results(idx_cur,i) = 1.0;
      }

      Teuchos::Array<Scalar> f1(num_pts),f2(num_pts),f3(num_pts);
      
      for (int i=0;i<num_pts;i++) {
	f1[i] = 0.5 * (1.0+2.0*x[i][0]+x[i][1]);
	f2[i] = 0.5 * (1.0-x[i][1]);
	f3[i] = f2[i] * f2[i];
      }

      // set D^{1,0} = f1
      idx_cur = idx2d(1,0);
      for (int i=0;i<num_pts;i++) {
	results(idx_cur,i) = f1[i];
      }

      // recurrence in p
      for (int p=1;p<n;p++) {
	idx_cur = idx2d(p,0);
	idx_curp1 = idx2d(p+1,0);
	idx_curm1 = idx2d(p-1,0);
	Scalar a = (2.0*p+1.0)/(1.0+p);
	Scalar b = p / (p+1.0);

	for (int i=0;i<num_pts;i++) {
	  results(idx_curp1,i) = a * f1[i] * results(idx_cur,i)
	    - b * f3[i] * results(idx_curm1,i);
	}
      }

      // D^{p,1}
      for (int p=0;p<n;p++) {
	int idxp0 = idx2d(p,0);
	int idxp1 = idx2d(p,1);
	for (int i=0;i<num_pts;i++) {
	  results(idxp1,i) = results(idxp0,i)
	    *0.5*(1.0+2.0*p+(3.0+2.0*p)*x[i][1]);
	}
      }

      // recurrence in q
      for (int p=0;p<n-1;p++) {
	for (int q=1;q<n-p;q++) {
	  int idxpqp1=idx2d(p,q+1);
	  int idxpq=idx2d(p,q);
	  int idxpqm1=idx2d(p,q-1);
	  Scalar a,b,c;
	  jrc(2*p+1,0,q,a,b,c);
	  for (int i=0;i<num_pts;i++) {
	    results(idxpqp1,i)
	      = (a*x[i][1]+b)*results(idxpq,i)
	      - c*results(idxpqm1,i);
	  }
	}
      }

      return;
    }

    static void tabulate_tetrahedron( const int &n , 
				      const Teuchos::Array<Point<Scalar> > &xin , 
				      Teuchos::SerialDenseMatrix<int,Scalar> &results )
    {
      const int num_pts = xin.size();

      // transformation from Pavel's element
      Teuchos::Array<Point<Scalar> > x;
      Scalar xcur_arr[3];
      int idxcur;

      for (int i=0;i<num_pts;i++) {
	xcur_arr[0] = 2.0 * xin[i][0] - 1.0;
	xcur_arr[1] = 2.0 * xin[i][1] - 1.0;
	xcur_arr[2] = 2.0 * xin[i][2] - 1.0;
	Point<Scalar> xcur( xcur_arr[0],xcur_arr[1],xcur_arr[2] );
        x.append( xcur );
      }

      Teuchos::Array<Scalar> f1(num_pts),f2(num_pts),f3(num_pts),f4(num_pts),f5(num_pts);

      for (int i=0;i<num_pts;i++) {
	f1[i] = 0.5 * ( 2.0 + 2.0*x[i][0] + x[i][1] + x[i][2] );
	f2[i] = pow( 0.5 * ( x[i][1] + x[i][2] ) , 2 );
	f3[i] = 0.5 * ( 1.0 + 2.0 * x[i][1] + x[i][2] );
	f4[i] = 0.5 * ( 1.0 - x[i][2] );
	f5[i] = f4[i] * f4[i];
      }

     
      // constant term
      idxcur = idx3d(0,0,0);
      for (int i=0;i<num_pts;i++) {
	results(idxcur,i) = 1.0;
      }

      // D^{1,0,0}
      idxcur = idx3d(1,0,0);
      for (int i=0;i<num_pts;i++) {
	results(idxcur,i) = f1[i];
      }

      // p recurrence
      for (int p=1;p<n;p++) {
	Scalar a1 = (2.0 * p + 1.0) / ( p + 1.0);
	Scalar a2 = p / ( p + 1.0 );
	int idxp = idx3d(p,0,0);
	int idxpp1 = idx3d(p+1,0,0);
	int idxpm1 = idx3d(p-1,0,0);
	//cout << idxpm1 << " " << idxp << " " << idxpp1 << endl;
	for (int i=0;i<num_pts;i++) {
	  results(idxpp1,i) = a1 * f1[i] * results(idxp,i) - a2 * f2[i] * results(idxpm1,i);
	}
      }

      // q = 1
      for (int p=0;p<n;p++) {
	int idx0 = idx3d(p,0,0);
	int idx1 = idx3d(p,1,0);
	for (int i=0;i<num_pts;i++) {
	  results(idx1,i) = results(idx0,i) * ( p * ( 1.0 + x[i][1] ) + 0.5 * ( 2.0 + 3.0 * x[i][1] + x[i][2] ) );
	}
      }

      // q recurrence
      for (int p=0;p<n-1;p++) {
	for (int q=1;q<n-p;q++) {
	  Scalar aq,bq,cq;
	  jrc(2.0*p+1.0,0,q,aq,bq,cq);
	  int idxpqp1 = idx3d(p,q+1,0);
	  int idxpq = idx3d(p,q,0);
	  int idxpqm1 = idx3d(p,q-1,0);
	  for (int i=0;i<num_pts;i++) {
	    results(idxpqp1,i) = ( aq * f3[i] + bq * f4[i] ) * results(idxpq,i) 
	      - ( cq * f5[i] ) * results(idxpqm1,i);
	  }
	}
      }

      // r = 1
      for (int p=0;p<n;p++) {
	for (int q=0;q<n-p;q++) {
	  int idxpq1 = idx3d(p,q,1);
	  int idxpq0 = idx3d(p,q,0);
	  for (int i=0;i<num_pts;i++) {
	    results(idxpq1,i) = results(idxpq0,i) * ( 1.0 + p + q + ( 2.0 + q + p ) * x[i][2] );
	  }
	}
      }

      // general r recurrence
      for (int p=0;p<n-1;p++) {
	for (int q=0;q<n-p-1;q++) {
	  for (int r=1;r<n-p-q;r++) {
	    Scalar ar,br,cr;
	    int idxpqrp1 = idx3d(p,q,r+1);
	    int idxpqr = idx3d(p,q,r);
	    int idxpqrm1 = idx3d(p,q,r-1);
	    jrc(2.0*p+2.0*q+2.0,0.0,r,ar,br,cr);
	    for (int i=0;i<num_pts;i++) {
	      results(idxpqrp1,i) = (ar * x[i][2] + br) * results( idxpqr , i ) - cr * results(idxpqrm1,i);
	    }
	  }
	}
      }
      
      return;

    }


  public:
    static void tabulate( ECell cell , int n , const Teuchos::Array<Point<Scalar> > &x ,
	      Teuchos::SerialDenseMatrix<int,Scalar> &results ) {
      if ( cell == CELL_TRI ) {
	tabulate_triangle( n , x , results );
      }
      else if ( cell == CELL_TET ) {
	tabulate_tetrahedron( n , x , results );
	//std::cout << results << std::endl;
      }
      
    }
  };
}

#endif


