#include "Phalanx_DataLayout.hpp"

//**********************************************************************
std::ostream& PHX::operator<<(std::ostream& os, const PHX::DataLayout& t)
{
  t.print(os);
  return os;
}

//**********************************************************************
