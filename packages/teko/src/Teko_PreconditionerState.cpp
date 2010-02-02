#include "Teko_PreconditionerState.hpp"

#include "Thyra_DefaultPreconditioner.hpp"

using namespace Thyra;

namespace Teko {
/////////////////////////////////////////////////////

//! Set parameters from a parameter list and return with default values.
void PreconditionerState::setParameterList(const RCP<Teuchos::ParameterList> & paramList)
{
   paramList_ = paramList;
}

//! Get the parameter list that was set using setParameterList().
RCP<Teuchos::ParameterList> PreconditionerState::getNonconstParameterList()
{
   if(paramList_==Teuchos::null)
      paramList_ = Teuchos::rcp(new Teuchos::ParameterList());

   return paramList_;
}

//! Unset the parameter list that was set using setParameterList(). 
RCP<Teuchos::ParameterList> PreconditionerState::unsetParameterList()
{
   RCP<Teuchos::ParameterList> paramList = paramList_;
   paramList_ = Teuchos::null;
   return paramList;
}

/////////////////////////////////////////////////////

} // end namespace Teko
