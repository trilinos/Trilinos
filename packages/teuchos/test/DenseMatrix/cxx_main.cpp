#include <iostream>
#include <string>
#include "Teuchos_DenseMatrix.hpp"

#include "mp/mpint.h"
#include "mp/mpreal.h"

using namespace std;
using namespace Teuchos;

int main(int argc, char *argv[])
{
  mp::mp_init(2000);

  int i;
  mp_real a[9], c[9];
  for(i = 0; i < 9; i++)
    {
      a[i] = sqrt(i);
      c[i] = 0;
    }

  DenseMatrix<int, mp_real> A(Copy, a, 3, 3, 3);
  DenseMatrix<int, mp_real> B(Copy, a, 3, 3, 3);
  DenseMatrix<int, mp_real> C;
  C.shape(3, 3);

  mp_real one = 1.0;
  C.multiply('N', 'N', one, A, B, one);

  C.TempPrint();

  mp::mp_finalize();

  return 0;
}
