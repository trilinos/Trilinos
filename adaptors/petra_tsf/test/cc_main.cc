#include <cstdlib>
#include <cassert>
#include <iostream>
#include <strstream>
#include <cstring>
#include <cmath>
#include <iomanip>
using namespace std;

#include "Petra_Comm.h"
#include "Petra_Map.h"
#include "Petra_TSF_RDP_MultiVector.h"
#include "WhizBangVectorAdder.h"

int main(int argc, char *argv[]) {

  Petra_Comm Comm;
  Petra_Map Map(10, 0, Comm);
  Petra_TSF_RDP_MultiVector V1(Map, 1);
  Petra_TSF_RDP_MultiVector V2(Map, 1);

  V1.Random(); // Fill V1 with random numbers
  V2.Random(); // Same for V2

  cout << "\nContents of V1:\n" << V1 << endl;
  cout << "\nContents of V2:\n" << V2 << endl;

  
  V1.Update(1.0, V2, 1.0);

  cout << "\nContents of V1:\n" << V1 << endl;
  cout << "\nContents of V2:\n" << V2 << endl;

  WhizBangVectorAdder Adder(&V1, &V2);

  Adder.Add(); // Compute V1 <- V1 + V2

  cout << "\nContents of V1:\n" << V1 << endl;
  cout << "\nContents of V2:\n" << V2 << endl;

  return 0; // All done
}
