#include "Teko_DiagonalPreconditionerFactory.hpp"
#include "Teko_DiagonalPreconditionerOp.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MpiComm.h"
#include "EpetraExt_PointToBlockDiagPermute.h"

using Teuchos::rcp;
using Teuchos::RCP;

namespace Teko {

DiagonalPrecondState::DiagonalPrecondState(){}

/*****************************************************/


DiagonalPreconditionerFactory::DiagonalPreconditionerFactory(){
}


RCP<PreconditionerState>  DiagonalPreconditionerFactory::buildPreconditionerState() const{
   DiagonalPrecondState*  mystate = new DiagonalPrecondState(); 
   return rcp(mystate);
}


LinearOp DiagonalPreconditionerFactory::buildPreconditionerOperator(LinearOp & lo,PreconditionerState & state) const{

  // Sanity check the state
  DiagonalPrecondState *MyState = dynamic_cast<DiagonalPrecondState *> (&state);
  TEUCHOS_ASSERT(MyState != 0);

  // Get the underlying Epetra_CrsMatrix, if we have one
  RCP<const Epetra_Operator> eo=Thyra::get_Epetra_Operator(*lo);
  TEUCHOS_ASSERT(eo!=Teuchos::null);
  TEUCHOS_ASSERT(dynamic_cast<const Epetra_MpiComm*>(&eo->Comm()));
  const Epetra_CrsMatrix *MAT=dynamic_cast<const Epetra_CrsMatrix*>(&*eo);
  TEUCHOS_ASSERT(MAT);

  // Create a new EpetraExt_PointToBlockDiagPermute for the state object, if we don't have one
  EpetraExt_PointToBlockDiagPermute *BDP;
  if(MyState->BDP_==Teuchos::null){
    BDP=new EpetraExt_PointToBlockDiagPermute(*MAT);
    BDP->SetParameters(List_);
    BDP->Compute();
    MyState->BDP_=rcp(BDP);
  }

  // Build the LinearOp object  (NTS: swapping the range and domain)
  DiagonalPreconditionerOp *MyOp=new DiagonalPreconditionerOp(MyState->BDP_,lo->domain(),lo->range());
  return rcp(MyOp);
}




void DiagonalPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList & pl){
  List_=pl;
}


} // end namespace Teko

