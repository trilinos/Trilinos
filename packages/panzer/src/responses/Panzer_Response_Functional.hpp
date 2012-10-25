#ifndef __Panzer_Response_Functional_hpp__
#define __Panzer_Response_Functional_hpp__

#include <string>


#include <mpi.h> // need for comm

#include "Teuchos_RCP.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MpiComm.h"

#include "Panzer_ResponseMESupport_Default.hpp"
#include "Panzer_GlobalEvaluationData.hpp"


namespace panzer {

/** This class provides a response storage for
  * simple functionals of the solution (i.e. scalar
  * values).
  */
template <typename EvalT>
class Response_Functional : public ResponseMESupport_Default<EvalT> {
public:
   typedef typename EvalT::ScalarT ScalarT;

   Response_Functional(const std::string & responseName,MPI_Comm comm)
     : ResponseMESupport_Default<EvalT>(responseName,comm) {}

   //! provide direct access, this thing is pretty simple
   ScalarT value;

   //! This simply does global summation, then shoves the result into a vector
   virtual void scatterResponse();

   virtual void initializeResponse()  
   { value = 0.0; }

   // from ResponseMESupport_Default

   //! What is the number of values you need locally
   virtual std::size_t localSizeRequired() const { return 1; }

   //! Is the vector distributed (or replicated)
   virtual bool vectorIsDistributed() const { return false; }
    
private:
   // hide these methods
   Response_Functional();
   Response_Functional(const Response_Functional &);
};

}

#include "Panzer_Response_Functional_impl.hpp"

#endif
