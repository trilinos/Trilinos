#include "Intrepid_RealSpace.hpp"

using namespace std;
using namespace Intrepid;


int main(int argc, char *argv[]) {

  cout << "\nTEST 1: class Point in 3D\n\n";

  double vec[]  = {1.0, 2.0, 3.0};
  double vec2[] = {-4.0, -1.0, 10.0};

  Point<double> v1(vec, 3);
  cout << "Created vector v1:\n" << v1 << endl;
  Point<double> v2(vec2, 3);
  cout << "Created vector v2:\n" << v2 << endl;
  cout << "Computing: v1 + v2\n";
  cout << v1+v2 << endl;
  cout << "Computing: v1 - v2\n";
  cout << v1-v2 << endl;
  cout << "Computing: v1 dot v2\n    ";
  cout << v1*v2 << endl;
  cout << "Computing: (v1 cross v2)\n";
  cout << (v1^v2) << endl;
  cout << "Computing: 5.0 * v1\n";
  cout << 5.0*v1 << endl;
  cout << "Computing: v3 = (1.0 / (v1 dot v2)) * v1 cross v2\n";
  Point<double> v3 = (1.0 / (v1*v2)) * v1 ^ v2;
  cout << v3 << "\n";

  cout << "\nEND TEST 1: class Point in 3D\n\n";


  cout << "\nTEST 2: class Point in 2D\n\n";

  Point<double> v4(vec, 2);
  cout << "Created vector v4:\n" << v4 << endl;
  Point<double> v5(vec2, 2);
  cout << "Created vector v5:\n" << v5 << endl;
  cout << "Computing: v4 + v5\n";
  cout << v4+v5 << endl;
  cout << "Computing: v4 - v5\n";
  cout << v4-v5 << endl;
  cout << "Computing: v4 dot v5\n    ";
  cout << v4*v5 << endl;
  cout << "Computing: 5.0 * v4\n";
  cout << 5.0*v4 << endl;
  cout << "Computing: v6 = (1.0 / (v4 dot v5)) * (v4 + v5)\n";
  Point<double> v6 = (1.0 / (v4*v5)) * (v4 + v5);
  cout << v6 << "\n";

  cout << "\nEND TEST 2: class Point in 2D\n\n";


  cout << "\nTEST 3: class Point in 1D\n\n";

  Point<double> v7(vec, 1);
  cout << "Created vector v7:\n" << v7 << endl;
  Point<double> v8(vec2, 1);
  cout << "Created vector v8:\n" << v8 << endl;
  cout << "Computing: v7 + v8\n";
  cout << v7+v8 << endl;
  cout << "Computing: v7 - v8\n";
  cout << v7-v8 << endl;
  cout << "Computing: v7 dot v8\n    ";
  cout << v7*v8 << endl;
  cout << "Computing: 5.0 * v7\n";
  cout << 5.0*v7 << endl;
  cout << "Computing: v9 = (1.0 / (v7 dot v8)) * (v7 + v8)\n";
  Point<double> v9 = (1.0 / (v7*v8)) * (v7 + v8);
  cout << v9 << "\n";

  cout << "\nEND TEST 3: class Point in 1D\n\n";



  return 0;
}
