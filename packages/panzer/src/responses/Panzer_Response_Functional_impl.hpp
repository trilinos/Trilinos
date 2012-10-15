#ifndef __Panzer_Response_Functional_impl_hpp__
#define __Panzer_Response_Functional_impl_hpp__

#include "Thyra_DefaultSpmdVectorSpace.hpp"

#include "Epetra_LocalMap.h"

#include "Sacado_Traits.hpp"

namespace panzer {

template <typename ScalarT>
Teuchos::RCP<const Epetra_Map> Response_Functional<ScalarT>::
getMap() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(useThyra_,std::logic_error,
                             "Reponse field \"" << getName() << "\" has previously been initialized as a "
                             "Thyra object, now trying to initalize as a Epetra! Error!");

  // lazily construct the map only as needed
  if(map_==Teuchos::null) 
    map_ = Teuchos::rcp(new Epetra_LocalMap(1,0,eComm_));

  return map_;
}

template <typename ScalarT>
void Response_Functional<ScalarT>::
setVector(const Teuchos::RCP<Epetra_Vector> & destVec)
{
  TEUCHOS_TEST_FOR_EXCEPTION(useThyra_,std::logic_error,
                             "Reponse field \"" << getName() << "\" has previously been initialized as a "
                             "Thyra object, now trying to initalize as a Epetra! Error!");

  eVector_ = destVec;

  useEpetra_ = true;
}

template <typename ScalarT>
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > Response_Functional<ScalarT>::
getVectorSpace() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(useEpetra_,std::logic_error,
                             "Reponse field \"" << getName() << "\" has previously been initialized as an "
                             "Epetra object, now trying to initalize as a Thyra object! Error!");

  // lazily build the space and return it
  if(vSpace_==Teuchos::null)
    vSpace_ = Thyra::defaultSpmdVectorSpace<double>(1);

  return vSpace_;
}

template <typename ScalarT>
void Response_Functional<ScalarT>::
setVector(const Teuchos::RCP<Thyra::VectorBase<double> > & destVec)
{
  TEUCHOS_TEST_FOR_EXCEPTION(useEpetra_,std::logic_error,
                             "Reponse field \"" << getName() << "\" has previously been initialized as an "
                             "Epetra object, now trying to initalize as a Thyra object! Error!");

  tVector_ = destVec;

  useThyra_ = true;
}

template <typename ScalarT>
void Response_Functional<ScalarT>::
scatterResponse() 
{
  double locValue = Sacado::ScalarValue<ScalarT>::eval(value);
  double glbValue = 0.0;

  // do global summation
  eComm_.SumAll(&locValue,&glbValue,1);

  // built data in vectors
  if(useEpetra_)
    (*eVector_)[0] = glbValue;
  else {
    Teuchos::ArrayRCP<double> data;
    Teuchos::rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(tVector_)->getNonconstLocalData(Teuchos::outArg(data));

    data[0] = glbValue;
  }
}

}

#endif
