#include <Kokkos_Core.hpp> 
#include <Kokkos_TaskPolicy.hpp>
#include <impl/Kokkos_Serial_TaskPolicy.hpp>
#include <Kokkos_Qthread.hpp>
#include <Qthread/Kokkos_Qthread_TaskPolicy.hpp>

#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "crs_matrix_helper.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef size_t size_type;

//typedef Kokkos::Qthread space_type;
typedef Kokkos::Serial space_type;

typedef Example::CrsMatrixBase<value_type,ordinal_type,size_type,space_type> CrsMatrixBase;
typedef Example::CrsMatrixView<CrsMatrixBase> CrsMatrixView;

typedef Example::CrsMatrixBase<CrsMatrixView,ordinal_type,size_type,space_type> CrsHierBase;
typedef Example::CrsMatrixHelper CrsMatrixHelper;

typedef Example::Uplo Uplo;

int main (int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " filename" << endl;
    return -1;
  }
  
  Kokkos::initialize();
  cout << "Default execution space initialized = "
       << typeid(Kokkos::DefaultExecutionSpace).name()
       << endl;

  // create a base matrix and import a sparse matrix from matrix market
  CrsMatrixBase A("A");

  ifstream in;
  in.open(argv[1]);
  if (!in.good()) {
    cout << "Error in open the file: " << argv[1] << endl;
    return -1;
  }
  A.importMatrixMarket(in);

  // copy and create for lower triangular 
  CrsMatrixBase L("L");
  L.copy(Uplo::Lower, A);
  cout << L << endl;

  // test hierarchical matrix
  CrsHierBase H("H");

  CrsMatrixHelper::flat2hier(L, H);
  cout << H << endl;

  Kokkos::finalize();

  return 0;
}
