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

#include "task_factory.hpp"
#include "crs_task_view.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef size_t size_type;

//typedef Kokkos::Qthread space_type;
typedef Kokkos::Serial space_type; 

typedef Example::CrsMatrixBase<value_type,ordinal_type,size_type,space_type> CrsMatrixBase;
typedef Example::CrsMatrixView<CrsMatrixBase> CrsMatrixView;

typedef Example::TaskFactory<Kokkos::TaskPolicy<space_type>,
                             Kokkos::Future<int,space_type> > TaskFactory;
typedef Example::CrsTaskView<CrsMatrixBase,TaskFactory> CrsTaskView;

typedef Example::CrsMatrixBase<CrsTaskView,ordinal_type,size_type,space_type> CrsHierBase;
typedef Example::CrsTaskView<CrsHierBase,TaskFactory> CrsHierView;

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
  
  CrsMatrixBase AA("AA");

  ifstream in;
  in.open(argv[1]);
  if (!in.good()) {
    cout << "Error in open the file: " << argv[1] << endl;
    return -1;
  }
  AA.importMatrixMarket(in);

  CrsMatrixBase LL("LL");
  LL.copy(Uplo::Lower, AA);

  CrsHierBase HH("HH");
  CrsMatrixHelper::flat2hier(LL, HH);

  CrsHierView H;
  H.setView(&HH, 2, 3, 2, 3);
  cout << "Block Matrix H = " << endl;
  H.showMeDetail(cout);

  Kokkos::finalize();

  return 0;
}

