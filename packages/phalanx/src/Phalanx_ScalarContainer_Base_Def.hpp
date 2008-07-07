#ifndef PHX_SCALAR_CONTAINER_BASE_DEF_HPP
#define PHX_SCALAR_CONTAINER_BASE_DEF_HPP

#include "Teuchos_TestForException.hpp"
// **************************************************************************
template<typename Traits>
PHX::ScalarContainerBase<Traits>::ScalarContainerBase()
{

}

// **************************************************************************
template<typename Traits>
PHX::ScalarContainerBase<Traits>::~ScalarContainerBase()
{

}

// **************************************************************************
template<typename Traits>
void PHX::ScalarContainerBase<Traits>::
requireField(const PHX::FieldTag& f) 
{ 
  vp_manager_.requireField(f);
}

// **************************************************************************
template<typename Traits>
void PHX::ScalarContainerBase<Traits>::
registerEvaluator(const Teuchos::RCP<PHX::Evaluator<Traits> >& e) 
{ 
  vp_manager_.registerEvaluator(e);
}
    
// **************************************************************************
template<typename Traits>
std::ostream&
PHX::operator<<(std::ostream& os, const PHX::ScalarContainerBase<Traits>& sc)
{ 
  sc.print(os);
  return os;
}

// **************************************************************************

#endif
