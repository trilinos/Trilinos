#include <Kokkos_Core.hpp>
#include "util.hpp"

#include "dense_matrix_base.hpp"
#include "dense_matrix_view.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef size_t size_type;

typedef Kokkos::OpenMP space_type;

using namespace Example;

typedef DenseMatrixBase<value_type,ordinal_type,size_type,space_type> DenseMatrixBaseType;
typedef DenseMatrixView<DenseMatrixBaseType> DenseMatrixViewType;

int main (int argc, char *argv[]) {
  Kokkos::initialize();
  cout << "Default execution space initialized = "
       << typeid(Kokkos::DefaultExecutionSpace).name()
       << endl;

  DenseMatrixBaseType AA("AA", 10, 10);

  ordinal_type cnt = 0;  
  for (ordinal_type i=0;i<AA.NumRows();++i)
    for (ordinal_type j=0;j<AA.NumCols();++j)
      AA.Value(i, j) = cnt++;
  cout << AA << endl;

  DenseMatrixViewType A(&AA,   2, 6, 
                       /**/    3, 5);
  cout << A << endl;

  Kokkos::finalize();

  return 0;
}
