#include <iostream>
#include <string>
#include "Teuchos_ParameterList.hpp"

using namespace std;
using namespace Teuchos;

class Foo { };

int main(int argc, char *argv[])
{
  ParameterList P;

  char c = 'c';
  double d = 3.14159;
  float f = 3.14;
  int i = 10;
  std::string s = "!dlrow, olleH";

  Foo F;

  P.SetParameter("c", c);
  P.SetParameter("d", d);
  P.SetParameter("f", f);
  P.SetParameter("i", i);
  P.SetParameter("s", s);

  P.Print(0);

  cout << P.GetParameter("c", ' ') << " " 
       << P.GetParameter("d", (double)0) << " " << P.GetParameter("f", (float)0)  << " " << P.GetParameter("i", (int)0) <<  " "
       << P.GetParameter("s", s) << endl;

  P.SetParameter("foo", F);

  P.GetParameter("does not exist", 42);

  P.Print(0);

  return 0;
}
