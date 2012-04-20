#include "Panzer_EpetraLinearObjContainer.hpp"

#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_EpetraLinearOp.hpp"

namespace panzer {

void EpetraLinearObjContainer::set_A_th(const Teuchos::RCP<Thyra::LinearOpBase<double> > & in)
{
   Teuchos::RCP<Epetra_Operator> eOp = Thyra::get_Epetra_Operator(*in);
   
   A = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(eOp,true); // throw on failure...must be CrsMatrix!
}

const Teuchos::RCP<Thyra::LinearOpBase<double> > EpetraLinearObjContainer::get_A_th() const
{
   return Thyra::nonconstEpetraLinearOp(A); // wrap in a thyra operator and return
}

}
