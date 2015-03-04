#include <Kokkos_Core.hpp> 
#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "partition.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

typedef Kokkos::OpenMP host_type;

using namespace Example;

typedef CrsMatrixBase<value_type,ordinal_type,size_type,host_type> CrsMatrixBaseType;
typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;

int main (int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " filename" << endl;
    return -1;
  }

  Kokkos::initialize();
  cout << "Default execution space initialized = "
       << typeid(Kokkos::DefaultExecutionSpace).name()
       << endl;

  CrsMatrixBaseType AA("AA");

  ifstream in;
  in.open(argv[1]);
  if (!in.good()) {
    cout << "Error in open the file: " << argv[1] << endl;
    return -1;
  }
  AA.importMatrixMarket(in);

  {
    CrsMatrixViewType A(AA);
    
    CrsMatrixViewType ATL, ATR,     A00,  a01,     A02,
      /**/            ABL, ABR,     a10t, alpha11, a12t,
      /**/                          A20,  a21,     A22;  
  
    Part_2x2(A,  ATL, ATR,
             /**/ABL, ABR, 
             0, 0, Partition::TopLeft);
    
    while (ATL.NumRows() < A.NumRows()) {
      Part_2x2_to_3x3(ATL, ATR, /**/  A00,  a01,     A02,
                      /*******/ /**/  a10t, alpha11, a12t,
                      ABL, ABR, /**/  A20,  a21,     A22,  
                      1, 1, Partition::BottomRight);
      // -----------------------------------------------------
      CrsMatrixViewType AB0, ab1;
      
      Merge_2x1(alpha11,
                a21,     ab1);
      
      Merge_2x1(a10t,
                A20,     AB0);
      
      cout << "alpha11 = " << endl
           << alpha11 << endl;
      
      cout << "a21 = " << endl
           << a21 << endl;

      cout << "A22 = " << endl
           << A22 << endl;
      
      // -----------------------------------------------------
      Merge_3x3_to_2x2(A00,  a01,     A02,  /**/ ATL, ATR,
                       a10t, alpha11, a12t, /**/ /******/
                       A20,  a21,     A22,  /**/ ABL, ABR,
                       Partition::TopLeft);
    }
  }

  Kokkos::finalize();
  
  return 0;
}
