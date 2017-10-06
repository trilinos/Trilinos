#include <Kokkos_Core.hpp>
#include "util.hpp"

#include "dense_matrix_base.hpp"

using namespace std;

typedef double value_type;
typedef int    ordinal_type;
typedef size_t size_type;

typedef Kokkos::OpenMP space_type;
typedef Kokkos::OpenMP device_type;

using namespace Example;

typedef DenseMatrixBase<value_type,ordinal_type,size_type,space_type> DenseMatrixBaseType;
typedef DenseMatrixBase<value_type,ordinal_type,size_type,device_type> DenseMatrixBaseDeviceType;

int main (int argc, char *argv[]) {
  Kokkos::initialize();
  cout << "Default execution space initialized = " 
       << typeid(Kokkos::DefaultExecutionSpace).name()
       << endl;

  { // Test on an empty matrix
    DenseMatrixBaseType A("Empty A");
    cout << A << endl;
  }

  { // Test on matrix allocation
    DenseMatrixBaseType A("A, 3x3 Allocated", 3, 3);
    cout << A << endl;
  }

  { // Test on attaching buffers
    ordinal_type m = 3, n = 3;
    DenseMatrixBaseType::value_type_array a("External::ValueArray", m*n);

    // column major 3x3 dense filling
    ordinal_type cnt = 0;
    for (ordinal_type i=0;i<m;++i) 
      for (ordinal_type j=0;j<n;++j,++cnt) 
        a(i+j*m) = (cnt + 10);
    
    DenseMatrixBaseType A("A, External buffer wrapped", 
                          m, n, m, 1, a);
    
    cout << A << endl;

    // Test on copying operations
    DenseMatrixBaseType B(A);
    B.setLabel("B, shallow-copy of A");
    cout << B << endl;

    DenseMatrixBaseType C("C, deep-copy of A");
    C.copy(A);
    cout << C << endl;

    DenseMatrixBaseType Dl("D, deep-copy of A lower triangular");
    Dl.copy(Uplo::Lower, A);
    cout << Dl << endl;

    DenseMatrixBaseType Du("D, deep-copy of A upper triangular");
    Du.copy(Uplo::Upper, A);
    cout << Du << endl;

    // Test on permutation operators
    DenseMatrixBaseType::ordinal_type_array p ("External::PermVector", n);
    DenseMatrixBaseType::ordinal_type_array ip("External::InvPermVector", m);    

    for (ordinal_type i=0;i<m;++i)
      ip[i] = (m - i - 1);

    for (ordinal_type j=0;j<n;++j)
      p[j] = (n - j - 1);

    DenseMatrixBaseType E("E, permuted in A");
    E.copy(ip, A);
    cout << E << endl;

    // Heterogeneous copy
    DenseMatrixBaseDeviceType F("F, allocated in the device");
    F.copy(A);
    
    DenseMatrixBaseType G("G, allocated in the host");    
    G.copy(F);

    cout << G << endl;
  }

  Kokkos::finalize();

  return 0;
}
