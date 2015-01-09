#include <iostream>
#include <stk_percept/Percept.hpp>

#include "MeshDifference.hpp"


using namespace stk_classic;
using namespace percept;

int main(int argc,  char **argv)
{
  
  MeshDifference md;
  md.run(argc, argv);


  return 0;
}
