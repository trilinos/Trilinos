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
#include "graph_helper_scotch.hpp" 

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

//typedef Kokkos::Qthread space_type;
typedef Kokkos::Serial space_type;

typedef Example::CrsMatrixBase<value_type,ordinal_type,size_type,space_type> CrsMatrixBase;
typedef Example::CrsMatrixView<CrsMatrixBase> CrsMatrixView;

typedef Example::CrsMatrixBase<CrsMatrixView,ordinal_type,size_type,space_type> CrsHierBase;
typedef Example::CrsMatrixHelper CrsMatrixHelper;

typedef Example::GraphHelper_Scotch<CrsMatrixBase> GraphHelper;    

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

  CrsMatrixBase AA("AA");

  ifstream in;
  in.open(argv[1]);
  if (!in.good()) {
    cout << "Error in open the file: " << argv[1] << endl;
    return -1;
  }
  AA.importMatrixMarket(in);

  GraphHelper S(AA);
  S.computeOrdering();
  cout << S << endl;

  CrsMatrixBase PA("Permuted AA");
  PA.copy(S.PermVector(), S.InvPermVector(), AA);

  CrsMatrixBase LL("LL");
  LL.copy(Uplo::Lower, PA);
  cout << LL << endl;

  CrsHierBase HH("HH");

  CrsMatrixHelper::flat2hier(LL, HH, 
                             S.NumBlocks(), 
                             S.RangeVector(),
                             S.TreeVector());
  cout << HH << endl;

  Kokkos::finalize();

  return 0;
}
