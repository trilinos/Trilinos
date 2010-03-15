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

//! Merge internal storage of another PreconditionerState object into this one
void PreconditionerState::merge(const PreconditionerState & ps,int position)
{
   // merge the two linearOps lists
   linearOps_.insert(ps.linearOps_.begin(),ps.linearOps_.end());
 
   // merge two parameter lists
   Teuchos::ParameterList::ConstIterator itr;
   if(ps.paramList_!=Teuchos::null) {
      Teuchos::RCP<Teuchos::ParameterList> paramList = getNonconstParameterList();
      for(itr=ps.paramList_->begin();itr!=ps.paramList_->end();++itr)
         paramList->setEntry(itr->first,itr->second);
   }
}

/////////////////////////////////////////////////////

} // end namespace Teko
