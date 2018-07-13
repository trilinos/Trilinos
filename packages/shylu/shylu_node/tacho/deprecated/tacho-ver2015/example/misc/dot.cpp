#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "dot.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef size_t size_type;

typedef Kokkos::OpenMP space_type; 

using namespace Example;

typedef CrsMatrixBase<value_type,ordinal_type,size_type,space_type> CrsMatrixBaseType;
typedef CrsMatrixView<CrsMatrixBaseType> CrsMatrixViewType;
typedef CrsRowView<CrsMatrixBaseType> CrsRowViewType;

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
    CrsRowViewType a, b;
    a.setView(CrsMatrixViewType(&AA, /**/ 2, 6, /**/ 3, 8), 2);
    b.setView(CrsMatrixViewType(&AA, /**/ 3, 7, /**/ 3, 8), 2);
    
    cout << " a = " << endl << a << endl;
    cout << " b = " << endl << b << endl;
    cout << " dot(a,b) = " << dot(a,b) << endl;
    cout << " dot(a,a) = " << dot(a) << endl;
  }

  Kokkos::finalize(); 

  return 0;
}
