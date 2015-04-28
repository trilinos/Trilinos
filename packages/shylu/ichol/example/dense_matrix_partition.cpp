#include <Kokkos_Core.hpp> 
#include "util.hpp"

#include "dense_matrix_base.hpp"
#include "dense_matrix_view.hpp"

#include "partition.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

typedef Kokkos::OpenMP host_type;

using namespace Example;

typedef DenseMatrixBase<value_type,ordinal_type,size_type,host_type> DenseMatrixBaseType;
typedef DenseMatrixView<DenseMatrixBaseType> DenseMatrixViewType;

int main (int argc, char *argv[]) {
  Kokkos::initialize();
  cout << "Default execution space initialized = "
       << typeid(Kokkos::DefaultExecutionSpace).name()
       << endl;

  DenseMatrixBaseType AA("AA", 10, 10);
  {
    DenseMatrixViewType A(AA);
    
    DenseMatrixViewType ATL, ATR,     A00,  a01,     A02,
      /**/              ABL, ABR,     a10t, alpha11, a12t,
      /**/                            A20,  a21,     A22;  
  
    Part_2x2(A,  ATL, ATR,
             /**/ABL, ABR, 
             0, 0, Partition::TopLeft);
    
    while (ATL.NumRows() < A.NumRows()) {
      Part_2x2_to_3x3(ATL, ATR, /**/  A00,  a01,     A02,
                      /*******/ /**/  a10t, alpha11, a12t,
                      ABL, ABR, /**/  A20,  a21,     A22,  
                      1, 1, Partition::BottomRight);
      // -----------------------------------------------------
      DenseMatrixViewType AB0, ab1;
      
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
