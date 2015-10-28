#include <iostream>

int main(int narg, char **arg)
{
  if (narg > 1)
    std::cout << arg[1] << std::endl;
  else
    std::cout << "BUMMER" << std::endl;
  return 0;
}
