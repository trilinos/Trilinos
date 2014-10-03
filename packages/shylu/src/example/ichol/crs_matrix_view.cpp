#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_row_view.hpp"
#include "crs_matrix_view.hpp"


using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef size_t size_type;

typedef Example::CrsMatrixBase<value_type,ordinal_type,size_type> CrsMatrixBase;
typedef Example::CrsMatrixView<CrsMatrixBase> CrsMatrixView;

int main (int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " filename" << endl;
    return -1;
  }
  
  // --------------------------------------------------------------------
  CrsMatrixBase Abase;

  ifstream in;
  in.open(argv[1]);
  if (!in.good()) {
    cout << "Error in open the file: " << argv[1] << endl;
    return -1;
  }
  Abase.importMatrixMarket(in);
  //Abase.showMe(cout);

  // --------------------------------------------------------------------
  CrsMatrixView A(Abase, 
                  2, 6, 
                  3, 8);
  A.showMe(cout);


  return 0;
}
