#include <Kokkos_Core.hpp>
#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "ichol.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

typedef Kokkos::OpenMP space_type; 

typedef Example::CrsMatrixBase<value_type,ordinal_type,size_type,space_type> CrsMatrixBase;
typedef Example::CrsMatrixView<CrsMatrixBase> CrsMatrixView;

typedef Example::Uplo Uplo;
typedef Example::AlgoIChol AlgoIChol;

using Example::IChol;

int main (int argc, char *argv[]) {
  if (argc < 3) {
    cout << "Usage: " << argv[0] << " filename" << " blksize" << endl;
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
  cout << AA << endl;

  CrsMatrixBase LL("Lower Triangular of AA");
  LL.copy(Uplo::Lower, AA);

  {
    CrsMatrixView L(LL);

    //int r_val = IChol<Uplo::Lower,AlgoIChol::LeftBlocked>::invoke(L);
    int r_val = 0;
    IChol<Uplo::Lower,AlgoIChol::LeftBlocked>::blocksize = atoi(argv[2]);
    IChol<Uplo::Lower,AlgoIChol::LeftBlocked>::TaskFunctor<CrsMatrixView>(L).apply(r_val);
    if (r_val != 0) { 
      cout << " Error = " << r_val << endl;
      return r_val;
    }
  }

  cout << LL << endl;

  Kokkos::finalize();

  return 0;
}
