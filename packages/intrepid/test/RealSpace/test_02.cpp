#include "Intrepid_RealSpace.hpp"

using namespace std;
using namespace Intrepid;


int main(int argc, char *argv[]) {

  cout << "\nTEST 1: class LinearMap in 3D\n\n";

  double mat3[3][3] = {{1,2,3},{4,5,7},{7,8,10}};

  cout << "Created linear map lmap3\n";
  LinearMap<double> lmap3(&mat3[0][0],3);
  cout << lmap3 << endl;

  cout << "Compute the transpose of lmap3\n";
  cout << lmap3.getTranspose() << endl;

  cout << "Transpose lmap3 (in place -- mutator)\n";
  lmap3.Transpose();
  cout << lmap3 << endl;

  cout << "Computing: lmap3 * lmap3\n";
  LinearMap<double> prod3 = lmap3 * lmap3;
  cout << prod3 << endl;

  cout << "Compute the inverse of lmap3\n";
  cout << lmap3.getInverse() << endl;

  cout << "Invert lmap3 (in place -- mutator)\n";
  lmap3.Invert();
  cout << lmap3 << endl;

  cout << "Invert lmap3 again (in place -- mutator)\n";
  lmap3.Invert();
  cout << lmap3 << endl;

  cout << "\nEND TEST 1: class LinearMap in 3D\n\n";


  cout << "\nTEST 2: class LinearMap in 2D\n\n";

  double mat2[] = {1,2,4,6};

  cout << "Created linear map lmap2\n";
  LinearMap<double> lmap2(mat2,2);
  cout << lmap2 << endl;

  cout << "Compute the transpose of lmap2\n";
  cout << lmap2.getTranspose() << endl;

  cout << "Transpose lmap2 (in place -- mutator)\n";
  lmap2.Transpose();
  cout << lmap2 << endl;

  cout << "Computing: lmap2 * lmap2\n";
  LinearMap<double> prod2 = lmap2 * lmap2;
  cout << prod2 << endl;

  cout << "Compute the inverse of lmap2\n";
  cout << lmap2.getInverse() << endl;

  cout << "Invert lmap2 (in place -- mutator)\n";
  lmap2.Invert();
  cout << lmap2 << endl;

  cout << "Invert lmap2 again (in place -- mutator)\n";
  lmap2.Invert();
  cout << lmap2 << endl;

  cout << "\nEND TEST 2: class LinearMap in 2D\n\n";


  cout << "\nTEST 3: class LinearMap in 1D\n\n";

  double mat1 = 4;

  cout << "Created linear map lmap1\n";
  LinearMap<double> lmap1(&mat1,1);
  cout << lmap1 << endl;

  cout << "Compute the transpose of lmap1\n";
  cout << lmap1.getTranspose() << endl;

  cout << "Transpose lmap1 (in place -- mutator)\n";
  lmap1.Transpose();
  cout << lmap1 << endl;

  cout << "Computing: lmap1 * lmap1\n";
  LinearMap<double> prod1 = lmap1 * lmap1;
  cout << prod1 << endl;

  cout << "Compute the inverse of lmap1\n";
  cout << lmap1.getInverse() << endl;

  cout << "Invert lmap1 (in place -- mutator)\n";
  lmap1.Invert();
  cout << lmap1 << endl;

  cout << "Invert lmap1 again (in place -- mutator)\n";
  lmap1.Invert();
  cout << lmap1 << endl;

  cout << "\nEND TEST 3: class LinearMap in 1D\n\n";

  return 0;
}
