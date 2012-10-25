#ifndef __Panzer_ResponseMESupport_Default_impl_hpp__
#define __Panzer_ResponseMESupport_Default_impl_hpp__

#include "Thyra_DefaultSpmdVectorSpace.hpp"

#include "Epetra_LocalMap.h"
#include "Epetra_Map.h"

namespace panzer {

template <typename EvalT>
Epetra_Vector & ResponseMESupport_Default<EvalT>::
getEpetraVector() const
{
   TEUCHOS_ASSERT(useEpetra());

   return *eVector_;
}

template <typename EvalT>
Thyra::ArrayRCP<double> ResponseMESupport_Default<EvalT>::
getThyraVector() const
{
   TEUCHOS_ASSERT(useThyra());

   Teuchos::ArrayRCP<double> data;
   Teuchos::rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(tVector_)->getNonconstLocalData(Teuchos::outArg(data));
 
   return data;
}

template <typename EvalT>
Teuchos::RCP<const Epetra_Map> ResponseMESupport_Default<EvalT>::
getMap() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(useThyra_,std::logic_error,
                             "Reponse field \"" << this->getName() << "\" has previously been initialized as a "
                             "Thyra object, now trying to initalize as a Epetra! Error!");

  // lazily construct the map only as needed
  if(map_==Teuchos::null) {
    if(this->vectorIsDistributed())
      map_ = Teuchos::rcp(new Epetra_Map(-1,(int) this->localSizeRequired(),0,eComm_));
    else
      map_ = Teuchos::rcp(new Epetra_LocalMap((int) this->localSizeRequired(),0,eComm_));
  }

  return map_;
}

template <typename EvalT>
void ResponseMESupport_Default<EvalT>::
setVector(const Teuchos::RCP<Epetra_Vector> & destVec)
{
  TEUCHOS_TEST_FOR_EXCEPTION(useThyra_,std::logic_error,
                             "Reponse field \"" << this->getName() << "\" has previously been initialized as a "
                             "Thyra object, now trying to initalize as a Epetra! Error!");

  eVector_ = destVec;

  useEpetra_ = true;
}

template <typename EvalT>
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > ResponseMESupport_Default<EvalT>::
getVectorSpace() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(useEpetra_,std::logic_error,
                             "Reponse field \"" << this->getName() << "\" has previously been initialized as an "
                             "Epetra object, now trying to initalize as a Thyra object! Error!");

  // lazily build the space and return it
  if(vSpace_==Teuchos::null)
    vSpace_ = Thyra::defaultSpmdVectorSpace<double>(1);
  // lazily construct the map only as needed
  if(map_==Teuchos::null) {
    if(this->vectorIsDistributed())
      vSpace_ = Thyra::defaultSpmdVectorSpace<double>(tComm_,this->localSizeRequired(),-1);
    else
      vSpace_ = Thyra::defaultSpmdVectorSpace<double>(this->localSizeRequired());
  }

  return vSpace_;
}

template <typename EvalT>
void ResponseMESupport_Default<EvalT>::
setVector(const Teuchos::RCP<Thyra::VectorBase<double> > & destVec)
{
  TEUCHOS_TEST_FOR_EXCEPTION(useEpetra_,std::logic_error,
                             "Reponse field \"" << this->getName() << "\" has previously been initialized as an "
                             "Epetra object, now trying to initalize as a Thyra object! Error!");

  tVector_ = destVec;

  useThyra_ = true;
}

}

#endif
