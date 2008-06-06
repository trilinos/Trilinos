#ifndef INTREPID_RECURRENCE_BASIS_HPP
#define INTREPID_RECURRENCE_BASIS_HPP
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_Array.hpp"

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


  public:
    static void tabulate( int dim , int n , const Teuchos::Array<Point<Scalar> > &x ,
	      Teuchos::SerialDenseMatrix<int,Scalar> &results ) {
      if ( dim == 2 ) {
	tabulate_triangle( n , x , results );
      }
    }
  };
}

#endif


