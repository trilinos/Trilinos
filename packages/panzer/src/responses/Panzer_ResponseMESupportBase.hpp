#ifndef __Panzer_ResponseMESupportBase_hpp__
#define __Panzer_ResponseMESupportBase_hpp__

#include <string>

#include "Teuchos_RCP.hpp"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"

#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"

#include "Panzer_ResponseBase.hpp"

namespace panzer {

template <typename EvalT>
class ResponseMESupportBase : public ResponseBase {
public:
   ResponseMESupportBase(const std::string & responseName)
     : ResponseBase(responseName) {}

   virtual ~ResponseMESupportBase() {}

   // This is the epetra view of the world
   ///////////////////////////////////////////////////////////

   //! Get the <code>Epetra_Map</code> for this response, map is constructed lazily.
   virtual Teuchos::RCP<const Epetra_Map> getMap() const = 0;

   /** Set the vector (to be filled) for this response. This must be 
     * constructed from the vector space returned by <code>getMap</code>.
     */
   virtual void setVector(const Teuchos::RCP<Epetra_Vector> & destVec) = 0;

   // This is the Thyra view of the world
   ///////////////////////////////////////////////////////////

   //! Get the vector space for this response, vector space is constructed lazily.
   virtual Teuchos::RCP<const Thyra::VectorSpaceBase<double> > getVectorSpace() const = 0;

   /** Set the vector (to be filled) for this response. This must be 
     * constructed from the vector space returned by <code>getVectorSpace</code>.
     */
   virtual void setVector(const Teuchos::RCP<Thyra::VectorBase<double> > & destVec) = 0;

private:
   // hide these methods
   ResponseMESupportBase();
   ResponseMESupportBase(const ResponseMESupportBase<EvalT> &);
};

template < >
class ResponseMESupportBase<panzer::Traits::Jacobian> : public ResponseBase {
public:
   ResponseMESupportBase(const std::string & responseName)
     : ResponseBase(responseName) {}

   virtual ~ResponseMESupportBase() {}
 
   //! Does this response support derivative evaluation?
   virtual bool supportsDerivative() const = 0;

   // This is the epetra view of the world
   ///////////////////////////////////////////////////////////

   //! Get the <code>Epetra_Vector</code> for this response, map is constructed lazily.
   virtual Teuchos::RCP<Epetra_Vector> buildEpetraDerivative() const = 0;

   /** Set the derivative (to be filled) for this response. 
     */
   virtual void setDerivative(const Teuchos::RCP<Epetra_Vector> & derivative) = 0;

   // This is the Thyra view of the world
   ///////////////////////////////////////////////////////////

   //! Get an <code>Epetra_Operator</code> for this response, map is constructed lazily.
   virtual Teuchos::RCP<Thyra::VectorBase<double> > buildDerivative() const = 0;

   /** Set the derivative (to be filled) for this response. This must be 
     * constructed from the vector space returned by <code>getMap</code>.
     */
   virtual void setDerivative(const Teuchos::RCP<Thyra::VectorBase<double> > & derivative) = 0;

private:
   // hide these methods
   ResponseMESupportBase();
   ResponseMESupportBase(const ResponseMESupportBase<panzer::Traits::Jacobian> &);
};

}

#endif
