#ifndef PHX_DATA_CONTAINER_BASE_DEF_HPP
#define PHX_DATA_CONTAINER_BASE_DEF_HPP

#include "Teuchos_TestForException.hpp"

// **************************************************************************
template<typename Traits>
PHX::DataContainerBase<Traits>::DataContainerBase()
{

}

// **************************************************************************
template<typename Traits>
PHX::DataContainerBase<Traits>::~DataContainerBase()
{

}

// **************************************************************************
template<typename Traits>
std::ostream& 
PHX::operator<<(std::ostream& os, const PHX::DataContainerBase<Traits>& dc)
{ 
  dc.print(os);
  return os;
}

// **************************************************************************

#endif
