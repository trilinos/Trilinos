#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "scale.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef size_t size_type;

typedef Kokkos::OpenMP host_type; 

typedef Example::CrsMatrixBase<value_type,ordinal_type,size_type,host_type> CrsMatrixBase;
typedef Example::CrsMatrixView<CrsMatrixBase> CrsMatrixView;
typedef Example::CrsRowView<CrsMatrixBase> CrsRowView;

int main (int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " filename" << endl;
    return -1;
  }

  Kokkos::initialize();
  cout << "Default execution space initialized = "
       << typeid(Kokkos::DefaultExecutionSpace).name()
       << endl;

  CrsMatrixBase AA("AA");

  ifstream in;
  in.open(argv[1]);
  if (!in.good()) {
    cout << "Error in open the file: " << argv[1] << endl;
    return -1;
  }
  AA.importMatrixMarket(in);

  {
    CrsMatrixView A(&AA, /**/ 2, 6, /**/ 3, 8);

    CrsRowView row;

    cout << " A = " << endl;
    for (ordinal_type i=0;i<A.NumRows();++i) {
      row.setView(A, i);
      cout << row << endl;
    }

    Example::scale(100, A); 

    cout << " scaled A by 100 = " << endl;
    for (ordinal_type i=0;i<A.NumRows();++i) {
      row.setView(A, i);
      cout << row << endl;
    }

  }

  Kokkos::finalize(); 

  return 0;
}
