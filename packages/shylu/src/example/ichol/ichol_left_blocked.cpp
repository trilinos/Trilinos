#include <Kokkos_Core.hpp>
#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "ichol_left_blocked.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

typedef Kokkos::OpenMP host_type; 

typedef Example::CrsMatrixBase<value_type,ordinal_type,size_type,host_type> CrsMatrixBase;
typedef Example::CrsMatrixView<CrsMatrixBase> CrsMatrixView;

typedef Example::Uplo Uplo;

int main (int argc, char *argv[]) {
  if (argc < 3) {
    cout << "Usage: " << argv[0] << " filename" << " blksize" << endl;
    return -1;
  }

  Kokkos::initialize();
  cout << "Default execution space initialized = "
       << typeid(Kokkos::DefaultExecutionSpace).name()
       << endl;
  
  CrsMatrixBase Abase("Abase");

  ifstream in;
  in.open(argv[1]);
  if (!in.good()) {
    cout << "Error in open the file: " << argv[1] << endl;
    return -1;
  }
  Abase.importMatrixMarket(in);
  cout << Abase << endl;

  CrsMatrixBase Lbase("Lower Triangular of Abase");
  Lbase.copy(Uplo::Lower, Abase);

  {
    CrsMatrixView L(Lbase);

    int r_val = Example::ichol_left_blocked_lower(L, atoi(argv[2]));
    if (r_val != 0) 
      cout << " Error = " << r_val << endl;
  }

  cout << Lbase << endl;

  Kokkos::finalize();

  return 0;
}
