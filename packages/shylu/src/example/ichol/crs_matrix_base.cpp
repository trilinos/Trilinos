#include <Kokkos_Core.hpp>
#include <iostream>
#include <typeinfo>

#include "util.hpp"
#include "crs_matrix_base.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef size_t size_type;

typedef Example::CrsMatrixBase<value_type,ordinal_type,size_type> CrsMatrixBase;
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

  { // Test on an empty matrix
    CrsMatrixBase A("Empty A");
    A.showMe(cout);
  }

  { // Test on matrix allocation
    CrsMatrixBase A("A, 3x3 Allocated", 3, 3, 9);
    A.showMe(cout);
  }

  { // Test on attaching buffers
    ordinal_type m = 3, n = 3, nnz = 9, cnt = 0;

    CrsMatrixBase::size_type_array    ap("External::RowPtrArray", m+1);
    CrsMatrixBase::ordinal_type_array aj("External::ColIdxArray", nnz);
    CrsMatrixBase::value_type_array   ax("External::ValueArray",  nnz);

    // column major 3x3 dense filling
    ap[0] = cnt;
    for (ordinal_type i=0;i<m;++i) {
      for (ordinal_type j=0;j<n;++j,++cnt) {
        aj[cnt] = j;
        ax[cnt] = (cnt + 10);
      }
      ap[i+1] = cnt;
    }

    CrsMatrixBase A("A, External buffer wrapped", 
                    m, n, nnz,
                    ap, aj, ax);
    
    A.showMe(cout);

    // Test on copying operations
    CrsMatrixBase B(A);
    B.setLabel("B, shallow-copy of A");
    B.showMe(cout);

    CrsMatrixBase C("C, deep-copy of A");
    C.copy(A);
    C.showMe(cout);

    CrsMatrixBase D("D, deep-copy of A lower triangular");
    D.copy(Uplo::Lower, A);
    D.showMe(cout);

    CrsMatrixBase::ordinal_type_array p ("External::PermVector", n);
    CrsMatrixBase::ordinal_type_array ip("External::InvPermVector", m);    

    for (ordinal_type i=0;i<m;++i)
      ip[i] = (m - i - 1);

    for (ordinal_type j=0;j<n;++j)
      p[j] = (n - j - 1);

    CrsMatrixBase E("E, permuted in A");
    E.copy(p, ip, A);
    E.showMe(cout);
  }

  { // File input
    ifstream in;
    in.open(argv[1]);
    if (!in.good()) {
      cout << "Error in open the file: " << argv[1] << endl;
      return -1;
    }
    CrsMatrixBase A("Imported A");
    A.importMatrixMarket(in);
    A.showMe(cout);
  }

  Kokkos::finalize();

  return 0;
}
