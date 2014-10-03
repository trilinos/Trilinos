#include "util.hpp"

#include "crs_matrix_base.hpp"

#include "crs_row_view.hpp"
#include "crs_matrix_view.hpp"

#include "partition.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

typedef Example::CrsMatrixBase<value_type,ordinal_type,size_type> CrsMatrixBase;
typedef Example::CrsMatrixView<CrsMatrixBase> CrsMatrixView;

typedef Example::Partition Partition;

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

  // --------------------------------------------------------------------
  CrsMatrixView A(Abase);

  CrsMatrixView 
    ATL, ATR,
    ABL, ABR;
  
  CrsMatrixView
    A00,  a01,     A02,
    a10t, alpha11, a12t,
    A20,  a21,     A22;
  
  Part_2x2(A,  ATL, ATR,
           /**/ABL, ABR, 
           0, 0, Partition::TopLeft);
  
  while (ATL.NumRows() < A.NumRows()) {
    Part_2x2_to_3x3(ATL, ATR, /**/  A00,  a01,     A02,
                    /*******/ /**/  a10t, alpha11, a12t,
                    ABL, ABR, /**/  A20,  a21,     A22,  
                    1, 1, Partition::BottomRight);
    // -----------------------------------------------------
    CrsMatrixView AB0, ab1;
    
    Merge_2x1(alpha11,
              a21,     ab1);
    
    Merge_2x1(a10t,
              A20,     AB0);

    cout << "alpha11 = " << endl;
    alpha11.showMe(cout);

    cout << "a21 = " << endl;
    a21.showMe(cout);

    cout << "A22 = " << endl;
    A22.showMe(cout);

    // -----------------------------------------------------
    Merge_3x3_to_2x2(A00,  a01,     A02,  /**/ ATL, ATR,
                     a10t, alpha11, a12t, /**/ /******/
                     A20,  a21,     A22,  /**/ ABL, ABR,
                     Partition::TopLeft);
  }
  
  return 0;
}
