#ifndef __Panzer_ResponseMESupport_Default_hpp__
#define __Panzer_ResponseMESupport_Default_hpp__

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include "Panzer_ResponseMESupportBase.hpp"

#include "Thyra_EpetraThyraWrappers.hpp"

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
   Epetra_MpiComm eComm_;
   mutable Teuchos::RCP<const Epetra_Map> map_;
   Teuchos::RCP<Epetra_Vector> eVector_;

   bool useThyra_;
   mutable Teuchos::RCP<const Thyra::VectorSpaceBase<double> > vSpace_;
   Teuchos::RCP<Thyra::VectorBase<double> > tVector_;
   Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > tComm_;
};

template < >
class ResponseMESupport_Default<panzer::Traits::Jacobian> : public ResponseMESupportBase<panzer::Traits::Jacobian> {
public:

   ResponseMESupport_Default(const std::string & responseName,MPI_Comm comm,
                             const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & derivVecSpace=Teuchos::null)
     : ResponseMESupportBase<panzer::Traits::Jacobian>(responseName), derivVecSpace_(derivVecSpace)
   { tComm_ = Teuchos::rcp(new Teuchos::MpiComm<Thyra::Ordinal>(Teuchos::opaqueWrapper(comm))); }

   virtual ~ResponseMESupport_Default() {}

   //! What is the number of values you need locally
   virtual std::size_t localSizeRequired() const = 0;

   //! Is the vector distributed (or replicated). For derivative assembly this must be false!
   virtual bool vectorIsDistributed() const = 0;

   //! Does this response support derivative evaluation?
   bool supportsDerivative() const { return getDerivativeVectorSpace()!=Teuchos::null; }

   /** Get the vector for this response. This must be 
     * constructed from the vector space returned by <code>getMap</code>.
     */
   Teuchos::RCP<Thyra::VectorBase<double> > getDerivative() const
   { return derivative_; }

   // This is the epetra view of the world
   ///////////////////////////////////////////////////////////

   //! Get the <code>Epetra_Map</code> for this response, map is constructed lazily.
   virtual Teuchos::RCP<Epetra_Vector> buildEpetraDerivative() const 
   {
     TEUCHOS_ASSERT(!vectorIsDistributed());
     TEUCHOS_ASSERT(localSizeRequired()==1);
     TEUCHOS_ASSERT(supportsDerivative());

     if(eMap_==Teuchos::null)
       eMap_ = Thyra::get_Epetra_Map(*getDerivativeVectorSpace(),Thyra::get_Epetra_Comm(*tComm_));

     return Teuchos::rcp(new Epetra_Vector(*eMap_));
   }

   /** Set the vector (to be filled) for this response. This must be 
     * constructed from the vector space returned by <code>getMap</code>.
     */
   virtual void setDerivative(const Teuchos::RCP<Epetra_Vector> & derivative) 
   { 
     TEUCHOS_ASSERT(!vectorIsDistributed());
     TEUCHOS_ASSERT(localSizeRequired()==1);
     TEUCHOS_ASSERT(supportsDerivative());
     TEUCHOS_ASSERT(eMap_!=Teuchos::null);

     derivative_ = Thyra::create_Vector(derivative,getDerivativeVectorSpace());
   }

   // This is the Thyra view of the world
   ///////////////////////////////////////////////////////////

   //! Get the <code>Epetra_Map</code> for this response, map is constructed lazily.
   virtual Teuchos::RCP<Thyra::VectorBase<double> > buildDerivative() const 
   { 
     TEUCHOS_ASSERT(!vectorIsDistributed());
     TEUCHOS_ASSERT(localSizeRequired()==1);
     TEUCHOS_ASSERT(supportsDerivative());
     return Thyra::createMember(*getDerivativeVectorSpace()); 
   }

   /** Set the vector (to be filled) for this response. This must be 
     * constructed from the vector space returned by <code>getMap</code>.
     */
   virtual void setDerivative(const Teuchos::RCP<Thyra::VectorBase<double> > & derivative) 
   { 
     TEUCHOS_ASSERT(!vectorIsDistributed());
     TEUCHOS_ASSERT(localSizeRequired()==1);
     TEUCHOS_ASSERT(supportsDerivative());
     derivative_ = derivative; 
   }

protected:
   //! Get the teuchos comm object
   Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > getComm() const { return tComm_; }

   //! Get the derivative vector space
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > getDerivativeVectorSpace() const
   { return derivVecSpace_; }

   //! Set the derivative vector space
   void setDerivativeVectorSpace(const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & vs)
   { derivVecSpace_ = vs; }

private:
   // hide these methods
   ResponseMESupport_Default();
   ResponseMESupport_Default(const ResponseMESupport_Default<panzer::Traits::Jacobian> &);

   Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > tComm_;
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > derivVecSpace_;
   mutable Teuchos::RCP<const Epetra_Map> eMap_;

   Teuchos::RCP<Thyra::VectorBase<double> > derivative_;
};

}

#include "Panzer_ResponseMESupport_Default_impl.hpp"

#endif
