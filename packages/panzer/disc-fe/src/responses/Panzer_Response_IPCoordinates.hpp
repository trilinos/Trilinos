#ifndef __Panzer_Response_IPCoordinates_hpp__
#define __Panzer_Response_IPCoordinates_hpp__

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
class Response_IPCoordinates : public ResponseBase {
public:
   typedef typename EvalT::ScalarT ScalarT;

   Response_IPCoordinates(const std::string & responseName)
     : ResponseBase(responseName) {}

   virtual void scatterResponse() {}

   virtual void initializeResponse()  
   { 
     if(coords==Teuchos::null)
       coords = Teuchos::rcp(new std::vector<panzer::Traits::Residual::ScalarT>);

     coords->clear();
   }

   Teuchos::RCP<const std::vector<panzer::Traits::Residual::ScalarT> > getCoords() const
   { return coords; }

   Teuchos::RCP<std::vector<panzer::Traits::Residual::ScalarT> > getNonconstCoords()
   { return coords; }

private:
   // hide these methods
   Response_IPCoordinates();
   Response_IPCoordinates(const Response_IPCoordinates &);

   Teuchos::RCP<std::vector<panzer::Traits::Residual::ScalarT> > coords;
};

}

#endif
