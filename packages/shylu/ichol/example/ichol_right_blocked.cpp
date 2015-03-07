#include <Kokkos_Core.hpp>

#include "util.hpp"
#include "graph_helper_scotch.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "crs_team_view.hpp"

#include "sequential_for.hpp"
#include "team_factory.hpp"

#include "ichol.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef int    size_type;

typedef Kokkos::OpenMP space_type;

using namespace Example;

typedef CrsMatrixBase<value_type,ordinal_type,size_type,space_type> CrsMatrixBaseType;
typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;

typedef TeamFactory<TeamPolicy,TeamThreadLoopRegion> TeamFactoryType;
typedef CrsTeamView<CrsMatrixBaseType,TeamFactoryType> CrsTeamViewType;

typedef GraphHelper_Scotch<CrsMatrixBaseType> GraphHelperType;

int main (int argc, char *argv[]) {
  if (argc < 3) {
    cout << "Usage: " << argv[0] << " filename" << " blksize" << endl;
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
  cout << AA << endl;

  GraphHelperType S(AA);
  S.computeOrdering();

  CrsMatrixBaseType PA("Permuted AA");
  PA.copy(S.PermVector(), S.InvPermVector(), AA);

  CrsMatrixBaseType UU("Upper Triangular of AA");
  UU.copy(Uplo::Upper, PA);

  CrsTeamViewType U(UU);

  {
    IChol<Uplo::Upper,AlgoIChol::RightBlocked>::blocksize = atoi(argv[2]);

    int r_val = 0;
    typedef typename CrsTeamViewType::policy_type::member_type member_type;

    IChol<Uplo::Upper,AlgoIChol::RightBlocked>
      ::TaskFunctor<CrsTeamViewType,SequentialFor>(U).apply(member_type(), r_val);

    if (r_val != 0) {
      cout << " Error = " << r_val << endl;
      return r_val;
    }
    cout << UU << endl;
  }

  Kokkos::finalize();

  return 0;
}
