#include "util.hpp"
#include "crs_matrix_base.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef size_t size_type;

typedef Example::CrsMatrixBase<value_type,ordinal_type,size_type> CrsMatrixBase;
typedef Example::Uplo Uplo;

int main (int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " filename" << endl;
    return -1;
  }
  
  CrsMatrixBase A;

  cout << "Empty Matrix A = " << endl;
  A.showMe(cout);

  ifstream in;
  in.open(argv[1]);
  if (!in.good()) {
    cout << "Error in open the file: " << argv[1] << endl;
    return -1;
  }
  A.importMatrixMarket(in);

  cout << "Imported Matrix A from " << argv[0] << endl
       << A << endl;

  CrsMatrixBase L(A, Uplo::Lower);
  cout << "Lower Triangular L = " << endl 
       << L << endl;

  return 0;
}
