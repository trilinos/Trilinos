#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_row_view.hpp"
#include "crs_matrix_view.hpp"

#include "dot.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef size_t size_type;

typedef Example::CrsMatrixBase<value_type,ordinal_type,size_type> CrsMatrixBase;
typedef Example::CrsMatrixView<CrsMatrixBase> CrsMatrixView;
typedef Example::CrsRowView<value_type,ordinal_type> CrsRowView;

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

  CrsMatrixView A2(A,   2, 6, 
                   /**/ 3, 8);
  
  CrsRowView r2 = A2.extractRow(2);
  
  CrsMatrixView A3(A,   3, 7, 
                   /**/ 3, 8);
  CrsRowView r3 = A3.extractRow(2);

  value_type result = Example::dot<CrsRowView>(r2,r3);

  cout << " r2 = " << endl << r2 << endl;
  cout << " r3 = " << endl << r3 << endl;
  cout << " Dot Result = " << result << endl;
  
  return 0;
}
