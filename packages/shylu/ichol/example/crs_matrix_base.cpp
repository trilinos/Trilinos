#include <Kokkos_Core.hpp>
#include "util.hpp"

#include "crs_matrix_base.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef size_t size_type;

typedef Kokkos::OpenMP space_type;
typedef Kokkos::OpenMP device_type;

using namespace Example;

typedef CrsMatrixBase<value_type,ordinal_type,size_type,space_type> CrsMatrixBaseType;
typedef CrsMatrixBase<value_type,ordinal_type,size_type,device_type> CrsMatrixBaseDeviceType;

int main (int argc, char *argv[]) {
  if (argc < 2) {
    cout << "Usage: " << argv[0] << " filename" << endl;
    return -1;
  }

  Kokkos::initialize();
  cout << "Default execution space initialized = " 
       << typeid(Kokkos::DefaultExecutionSpace).name()
       << endl;

  { // Test on an empty matrix
    CrsMatrixBaseType A("Empty A");
    cout << A << endl;
  }

  { // Test on matrix allocation
    CrsMatrixBaseType A("A, 3x3 Allocated", 3, 3, 9);
    cout << A << endl;
  }

  { // Test on attaching buffers
    ordinal_type m = 3, n = 3, nnz = 9, cnt = 0;

    CrsMatrixBaseType::size_type_array    ap("External::RowPtrArray", m+1);
    CrsMatrixBaseType::ordinal_type_array aj("External::ColIdxArray", nnz);
    CrsMatrixBaseType::value_type_array   ax("External::ValueArray",  nnz);

    // column major 3x3 dense filling
    ap[0] = cnt;
    for (ordinal_type i=0;i<m;++i) {
      for (ordinal_type j=0;j<n;++j,++cnt) {
        aj[cnt] = j;
        ax[cnt] = (cnt + 10);
      }
      ap[i+1] = cnt;
    }

    CrsMatrixBaseType A("A, External buffer wrapped", 
                    m, n, nnz,
                    ap, aj, ax);
    
    cout << A << endl;

    // Test on copying operations
    CrsMatrixBaseType B(A);
    B.setLabel("B, shallow-copy of A");
    cout << B << endl;

    CrsMatrixBaseType C("C, deep-copy of A");
    C.copy(A);
    cout << C << endl;

    CrsMatrixBaseType Dl("D, deep-copy of A lower triangular");
    Dl.copy(Uplo::Lower, A);
    cout << Dl << endl;

    CrsMatrixBaseType Du("D, deep-copy of A upper triangular");
    Du.copy(Uplo::Upper, A);
    cout << Du << endl;

    // Test on permutation operators
    CrsMatrixBaseType::ordinal_type_array p ("External::PermVector", n);
    CrsMatrixBaseType::ordinal_type_array ip("External::InvPermVector", m);    

    for (ordinal_type i=0;i<m;++i)
      ip[i] = (m - i - 1);

    for (ordinal_type j=0;j<n;++j)
      p[j] = (n - j - 1);

    CrsMatrixBaseType E("E, permuted in A");
    E.copy(p, ip, A);
    cout << E << endl;

    // Heterogeneous copy
    CrsMatrixBaseDeviceType F("F, allocated in the device");
    F.copy(A);
    
    CrsMatrixBaseType G("G, allocated in the host");    
    G.copy(F);

    cout << G << endl;
  }

  { // Test on file io
    ifstream in;
    in.open(argv[1]);
    if (!in.good()) {
      cout << "Error in open the file: " << argv[1] << endl;
      return -1;
    }
    CrsMatrixBaseType A("A, imported from a file");
    A.importMatrixMarket(in);

    cout << A << endl;
  }

  Kokkos::finalize();

  return 0;
}
