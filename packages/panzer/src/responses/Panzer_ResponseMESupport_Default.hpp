#ifndef __Panzer_ResponseMESupport_Default_hpp__
#define __Panzer_ResponseMESupport_Default_hpp__

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include "Panzer_ResponseMESupportBase.hpp"

namespace panzer {

template <typename EvalT>
class ResponseMESupport_Default : public ResponseMESupportBase<EvalT> {
public:
   ResponseMESupport_Default(const std::string & responseName,MPI_Comm comm)
     : ResponseMESupportBase<EvalT>(responseName), useEpetra_(false), eComm_(comm), useThyra_(false) 
   {
     tComm_ = Teuchos::rcp(new Teuchos::MpiComm<Thyra::Ordinal>(Teuchos::opaqueWrapper(comm)));
   }

   virtual ~ResponseMESupport_Default() {}

   //! What is the number of values you need locally
   virtual std::size_t localSizeRequired() const = 0;

   //! Is the vector distributed (or replicated)
   virtual bool vectorIsDistributed() const = 0;

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

protected:
   //! Get the teuchos comm object
   Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > getComm() const { return tComm_; }

   //! Is Epetra the right vector
   bool useEpetra() const { return useEpetra_; }

   //! Is Thyra the right vector
   bool useThyra() const { return useThyra_; }

   //! Access the epetra vector
   Epetra_Vector & getEpetraVector() const;

   //! Access the thyra vector
   Thyra::ArrayRCP<double> getThyraVector() const;

private:
   // hide these methods
   ResponseMESupport_Default();
   ResponseMESupport_Default(const ResponseMESupport_Default<EvalT> &);

   bool useEpetra_;
   mutable Teuchos::RCP<const Epetra_Map> map_;
   Teuchos::RCP<Epetra_Vector> eVector_;
   Epetra_MpiComm eComm_;

   bool useThyra_;
   mutable Teuchos::RCP<const Thyra::VectorSpaceBase<double> > vSpace_;
   Teuchos::RCP<Thyra::VectorBase<double> > tVector_;
   Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > tComm_;
};

}

#include "Panzer_ResponseMESupport_Default_impl.hpp"

#endif
