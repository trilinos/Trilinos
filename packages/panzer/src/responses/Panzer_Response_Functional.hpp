#ifndef __Panzer_Response_Functional_hpp__
#define __Panzer_Response_Functional_hpp__

#include <string>

#include "Panzer_ResponseBase.hpp"

#include <mpi.h> // need for comm

#include "Teuchos_RCP.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MpiComm.h"

#include "Panzer_GlobalEvaluationData.hpp"


namespace panzer {

/** This class provides a response storage for
  * simple functionals of the solution (i.e. scalar
  * values).
  */
template <typename ScalarT>
class Response_Functional : public ResponseBase {
public:

   Response_Functional(const std::string & responseName,MPI_Comm comm)
     : ResponseBase(responseName), useEpetra_(false), eComm_(comm), useThyra_(false) {}

   //! provide direct access, this thing is pretty simple
   ScalarT value;

   //! This simply does global summation, then shoves the result into a vector
   virtual void scatterResponse();

   virtual void initializeResponse()  
   { value = 0.0; }

   // This is the epetra view of the world
   ///////////////////////////////////////////////////////////

   //! Get the <code>Epetra_Map</code> for this response, map is constructed lazily.
   Teuchos::RCP<const Epetra_Map> getMap() const;

   /** Set the vector (to be filled) for this response. This must be 
     * constructed from the vector space returned by <code>getMap</code>.
     */
   void setVector(const Teuchos::RCP<Epetra_Vector> & destVec);

   // This is the Thyra view of the world
   ///////////////////////////////////////////////////////////

   //! Get the vector space for this response, vector space is constructed lazily.
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > getVectorSpace() const;

   /** Set the vector (to be filled) for this response. This must be 
     * constructed from the vector space returned by <code>getVectorSpace</code>.
     */
   void setVector(const Teuchos::RCP<Thyra::VectorBase<double> > & destVec);

private:
   // hide these methods
   Response_Functional();
   Response_Functional(const Response_Functional &);

   bool useEpetra_;
   mutable Teuchos::RCP<const Epetra_Map> map_;
   Teuchos::RCP<Epetra_Vector> eVector_;
   Epetra_MpiComm eComm_;

   bool useThyra_;
   mutable Teuchos::RCP<const Thyra::VectorSpaceBase<double> > vSpace_;
   Teuchos::RCP<Thyra::VectorBase<double> > tVector_;
};

}

#include "Panzer_Response_Functional_impl.hpp"

#endif
