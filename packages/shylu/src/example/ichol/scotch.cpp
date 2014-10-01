#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "graph_helper_scotch.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;
typedef Example::CrsMatrixBase<value_type,ordinal_type,size_type> CrsMatrix;
typedef Example::GraphHelper_Scotch<CrsMatrix> GraphHelper;

int main (int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " filename" << endl;
    return -1;
  }
  
  // --------------------------------------------------------------------
  CrsMatrix A;

  ifstream in;
  in.open(argv[1]);
  if (!in.good()) {
    cout << "Error in open the file: " << argv[1] << endl;
    return -1;
  }
  A.importMatrixMarket(in);
  A.showMe(cout);

  // --------------------------------------------------------------------
  GraphHelper S(A);

  S.computeOrdering();
  S.showMe(cout);
  // --------------------------------------------------------------------
  
  

  return 0;
}
