#include <iostream>
#include "Teuchos_BLAS.hpp"

#define cout std::cout
#define endl std::endl

int main()
{
  Teuchos::BLAS<int, double> D;
  Teuchos::BLAS<int, float> F;
  Teuchos::BLAS<int, int> I;

  double* dx = new double[3];
  double* dy = new double[3];
  float* fx = new float[3];
  float* fy = new float[3];
  int* ix = new int[3];
  int* iy = new int[3];

  int i;
  for(i = 0; i < 3; i++)
    {
      dx[i] = i + 1;
      fx[i] = i + 1;
      ix[i] = i + 1;
      dy[i] = i + 4;
      fy[i] = i + 4;
      iy[i] = i + 4;
    }

  D.AXPY(3, 1, dx, 1, dy, 1);
  F.AXPY(3, 1, fx, 1, fy, 1);
  I.AXPY(3, 1, ix, 1, iy, 1);

  for(i = 0; i < 3; i++)
    {
      cout << dy[i] << " ";
    }
  cout << endl;

  for(i = 0; i < 3; i++)
    {
      cout << fy[i] << " ";
    }
  cout << endl;

  for(i = 0; i < 3; i++)
    {
      cout << iy[i] << " ";
    }
  cout << endl;

  cout << I.DOT(3, ix, 1, iy, 1) << endl;

  return 0;
}
