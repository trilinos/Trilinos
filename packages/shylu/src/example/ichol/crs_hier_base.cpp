#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef size_t size_type;

typedef Example::CrsMatrixBase<value_type,ordinal_type,size_type> CrsMatrixBase;
typedef Example::CrsMatrixView<CrsMatrixBase> CrsMatrixView;

typedef Example::CrsMatrixBase<CrsMatrixView,ordinal_type,size_type> CrsHierBase;

typedef Example::Uplo Uplo;

int main (int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " filename" << endl;
    return -1;
  }
  
  CrsMatrixBase A;

  ifstream in;
  in.open(argv[1]);
  if (!in.good()) {
    cout << "Error in open the file: " << argv[1] << endl;
    return -1;
  }
  A.importMatrixMarket(in);

  CrsMatrixBase L(A, Uplo::Lower);
  cout << "Lower Triangular L = " << endl
       << L << endl;

  // only 1x1 block matrix is implemented 
  // later scotch interface will be placed here 
  CrsHierBase H(L, 1, 1); 
  cout << "Hierarchical Block Matrix of L = " << endl
       << H << endl;

  return 0;
}
