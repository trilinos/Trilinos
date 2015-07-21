#include <Kokkos_Core.hpp> 
#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "graph_helper_scotch.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

typedef Kokkos::Serial space_type;

using namespace Example;

typedef CrsMatrixBase<value_type,ordinal_type,size_type,space_type> CrsMatrixBaseType;
typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;

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
    GraphHelperType S(AA);
    
    S.computeOrdering();
    cout << S << endl;
    
    CrsMatrixBaseType PA("Permuted AA");
    PA.copy(S.PermVector(), S.InvPermVector(), AA);
  }

  Kokkos::finalize();

  return 0;
}
