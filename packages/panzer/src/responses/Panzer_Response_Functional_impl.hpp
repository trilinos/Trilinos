#ifndef __Panzer_Response_Functional_impl_hpp__
#define __Panzer_Response_Functional_impl_hpp__

#include "Teuchos_Comm.hpp"

#include "Epetra_LocalMap.h"

#include "Sacado_Traits.hpp"

namespace panzer {

template <typename EvalT>
void Response_Functional<EvalT>::
scatterResponse() 
{
  double locValue = Sacado::ScalarValue<ScalarT>::eval(value);
  double glbValue = 0.0;

  // do global summation
  Teuchos::reduceAll(*this->getComm(), Teuchos::REDUCE_SUM, static_cast<Thyra::Ordinal>(1), &locValue,&glbValue);

  value = glbValue;

  // built data in vectors
  if(this->useEpetra()) {
    // use epetra 
    this->getEpetraVector()[0] = glbValue;
  }
  else {
    // use thyra
    TEUCHOS_ASSERT(this->useThyra());

    this->getThyraVector()[0] = glbValue;
  }
}

}

#endif
