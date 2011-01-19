#include "Panzer_AssemblyEngine_InArgs.hpp"

#include "Teuchos_RCP.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"

#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include "Teuchos_DefaultSerialComm.hpp"

namespace panzer {

// try to get a comm object from a vector space object
static Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > getTeuchosComm(const Thyra::VectorSpaceBase<double> & vs)
{
   Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > comm;
   try {
      comm = Teuchos::dyn_cast<const Thyra::SpmdVectorSpaceBase<double> >(vs).getComm();
   }
   catch(std::bad_cast) {
      comm = Teuchos::rcp(new Teuchos::SerialComm<Thyra::Ordinal>);
   }

   return comm;
}

void AssemblyEngineInArgs::thyraToEpetra()
{
   // convert the linear operator
   Teuchos::RCP<const Epetra_Map> range, domain;
   if(th_j!=Teuchos::null) {
      Teuchos::RCP<Thyra::EpetraLinearOp> epJac = Teuchos::rcp_dynamic_cast<Thyra::EpetraLinearOp>(th_j);
      j = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(epJac->epetra_op());

      range = Teuchos::rcpFromRef(j->OperatorRangeMap());
      domain = Teuchos::rcpFromRef(j->OperatorDomainMap());
   }
   else {
      // extract a Teuchos Commm object
      Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > comm;
      if(th_x!=Teuchos::null) 
         comm = getTeuchosComm(*th_x->range());
      else if(th_dxdt!=Teuchos::null)
         comm = getTeuchosComm(*th_dxdt->range());
      else if(th_f!=Teuchos::null)
         comm = getTeuchosComm(*th_f->range());
      else
         TEST_FOR_EXCEPTION(true,std::runtime_error,"Could not construct a Teuchos::Comm object");

      Teuchos::RCP<const Epetra_Comm> eComm = Thyra::get_Epetra_Comm(*comm);

      // build teh maps
      if(th_x!=Teuchos::null)
         range = Thyra::get_Epetra_Map(*th_x->range(),eComm);
      else if(th_dxdt!=Teuchos::null)
         range = Thyra::get_Epetra_Map(*th_dxdt->range(),eComm);
      if(th_f!=Teuchos::null)
         domain = Thyra::get_Epetra_Map(*th_f->range(),eComm);
   }

   // convert the vectors
   if(th_x!=Teuchos::null)    x    = Thyra::get_Epetra_Vector(j->OperatorDomainMap(),th_x->col(0)); 
   if(th_dxdt!=Teuchos::null) dxdt = Thyra::get_Epetra_Vector(j->OperatorDomainMap(),th_dxdt->col(0)); 
   if(th_f!=Teuchos::null)    f    = Thyra::get_Epetra_Vector(j->OperatorRangeMap(),th_f->col(0)); 
}

void AssemblyEngineInArgs::epetraToThyra(const Teuchos::RCP<const Epetra_Map> & rangeMap, 
                                         const Teuchos::RCP<const Epetra_Map> & domainMap)
{
   using Teuchos::RCP;

   if(j!=Teuchos::null) th_j = Thyra::nonconstEpetraLinearOp(j); 

   RCP<const Thyra::VectorSpaceBase<double> > rangeVS = Thyra::create_VectorSpace(rangeMap);
   RCP<const Thyra::VectorSpaceBase<double> > domainVS = Thyra::create_VectorSpace(domainMap);

   // convert the vectors
   if(x!=Teuchos::null) {
      th_x = Thyra::create_Vector(x,domainVS);
      Teuchos::set_extra_data(domainMap,"epetra_map",Teuchos::inOutArg(th_x));
   }
   if(dxdt!=Teuchos::null) {
      th_dxdt = Thyra::create_Vector(dxdt,domainVS);
      Teuchos::set_extra_data(domainMap,"epetra_map",Teuchos::inOutArg(th_dxdt));
   }
   if(f!=Teuchos::null) { 
      th_f = Thyra::create_Vector(f,rangeVS);
      Teuchos::set_extra_data(rangeMap,"epetra_map",Teuchos::inOutArg(th_f));
   }
}

}

