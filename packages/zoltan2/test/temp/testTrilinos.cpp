// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
// ***********************************************************************
//
// Temporary tests (for example to verify how a Trilinos class works)
// can be created and built here.  Feel free to overwrite this test.

#include <Teuchos_ArrayRCP.hpp>

using Teuchos::ArrayRCP;
using Teuchos::arcp;

void func()
{
  ArrayRCP<int> *l = new ArrayRCP<int> [10];
  for (int i=0; i < 10; i++){
    l[i] = arcp(new int [100],0, 100, true);
  }
  ArrayRCP<ArrayRCP<int> > ll(l, 0, 10, true);
}

int main()
{
  func();
}
