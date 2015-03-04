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

#include "team_factory.hpp"
#include "crs_team_view.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef size_t size_type;

typedef Kokkos::Serial space_type; 

using namespace Example;

typedef CrsMatrixBase<value_type,ordinal_type,size_type,space_type> CrsMatrixBaseType;
typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;

typedef TeamFactory<Kokkos::TeamPolicy<space_type>, Kokkos::Impl::TeamThreadLoopBoundariesStruct> TeamFactoryType;
typedef CrsTeamView<CrsMatrixBaseType,TeamFactoryType> CrsTeamViewType;

typedef CrsMatrixBase<CrsTeamViewType,ordinal_type,size_type,space_type> CrsHierBaseType;
typedef CrsTeamView<CrsHierBaseType,TeamFactoryType> CrsHierViewType;

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

  CrsMatrixBaseType LL("LL");
  LL.copy(Uplo::Lower, AA);

  CrsHierBaseType HH("HH");
  CrsMatrixHelper::flat2hier(LL, HH);

  cout << "Hier Matrix HH = " << endl
       << HH << endl;

  CrsHierViewType H;
  H.setView(&HH, 2, 3, 2, 3);

  cout << "Block Partitioned Matrix H = " << endl
       << H << endl;

  Kokkos::finalize();

  return 0;
}

