#include "Teuchos_BLAS.hpp"
#include "Teuchos_Version.hpp"

int main(int argc, char* argv[])
{
  cout << Teuchos::Teuchos_Version() << endl << endl;

  // Creating an instance of the BLAS class for double-precision kernels looks like:
  Teuchos::BLAS<int, double> blas;

  // This instance provides the access to all the BLAS kernels listed in Figure \ref{blas_kernels}:
  const int n = 10;
  double alpha = 2.0;
  double x[ n ];
  for ( int i=0; i<n; i++ ) { x[i] = i; }
  blas.SCAL( n, alpha, x, 1 );
  int max_idx = blas.IAMAX( n, x, 1 );
  cout<< "The index of the maximum magnitude entry of x[] is the "
      <<  max_idx <<"-th and x[ " << max_idx-1 << " ] = "<< x[max_idx-1] 
      << endl;
  
  return 0;
}
