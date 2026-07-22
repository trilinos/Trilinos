// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_ResponseMESupport_Default_hpp__
#define __Panzer_ResponseMESupport_Default_hpp__

#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

#include "Panzer_ResponseMESupportBase.hpp"

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_SpmdVectorBase.hpp"

#include "PanzerDiscFE_config.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Epetra_LocalMap.h"
#include "Epetra_Map.h"
#include "Thyra_EpetraThyraWrappers.hpp"
#endif

namespace panzer {

template <typename EvalT>
class ResponseMESupport_Default : public ResponseMESupportBase<EvalT> {
public:
   ResponseMESupport_Default(const std::string & responseName,MPI_Comm comm)
     : ResponseMESupportBase<EvalT>(responseName), useEpetra_(false),
#ifdef PANZER_HAVE_EPETRA_STACK
     eComm_(comm),
#endif
     useThyra_(false)
   {
     tComm_ = Teuchos::rcp(new Teuchos::MpiComm<Thyra::Ordinal>(Teuchos::opaqueWrapper(comm)));
   }

   virtual ~ResponseMESupport_Default() {}

   //! What is the number of values you need locally
   virtual std::size_t localSizeRequired() const = 0;

   //! Is the vector distributed (or replicated)
   virtual bool vectorIsDistributed() const = 0;

#ifdef PANZER_HAVE_EPETRA_STACK
   // This is the epetra view of the world
   ///////////////////////////////////////////////////////////

   //! Get the <code>Epetra_Map</code> for this response, map is constructed lazily.
   Teuchos::RCP<const Epetra_Map> getMap() const;

   /** Set the vector (to be filled) for this response. This must be
     * constructed from the vector space returned by <code>getMap</code>.
     */
   void setVector(const Teuchos::RCP<Epetra_Vector> & destVec);
#endif

   // This is the Thyra view of the world
   ///////////////////////////////////////////////////////////

   //! Get the vector space for this response, vector space is constructed lazily.
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > getVectorSpace() const;

   /** Set the vector (to be filled) for this response. This must be
     * constructed from the vector space returned by <code>getVectorSpace</code>.
     */
   void setVector(const Teuchos::RCP<Thyra::VectorBase<double> > & destVec);

   //! set the vector space for this response
   void setVectorSpace(Teuchos::RCP<const Thyra::VectorSpaceBase<double> > vs)
   { vSpace_ = vs; }

   //! Access the response vector
   Teuchos::RCP<Thyra::VectorBase<double> > getVector() const
   { return tVector_; }

protected:
   //! Get the teuchos comm object
   Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > getComm() const { return tComm_; }

   //! Is Epetra the right vector
   bool useEpetra() const { return useEpetra_; }

   //! Is Thyra the right vector
   bool useThyra() const { return useThyra_; }

#ifdef PANZER_HAVE_EPETRA_STACK
   //! Access the epetra vector
   Epetra_Vector & getEpetraVector() const;
#endif

   //! Access the thyra vector
   Thyra::ArrayRCP<double> getThyraVector() const;

   //! Access the thyra MultiVector
    Teuchos::RCP<Thyra::MultiVectorBase<double> > getThyraMultiVector() const
    { return tVector_;}


private:
   // hide these methods
   ResponseMESupport_Default();
   ResponseMESupport_Default(const ResponseMESupport_Default<EvalT> &);

   bool useEpetra_;
#ifdef PANZER_HAVE_EPETRA_STACK
   Epetra_MpiComm eComm_;
   mutable Teuchos::RCP<const Epetra_Map> map_;
   Teuchos::RCP<Epetra_Vector> eVector_;
#endif

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
   Teuchos::RCP<Thyra::MultiVectorBase<double> > getDerivative() const
   { return derivative_; }

#ifdef PANZER_HAVE_EPETRA_STACK
   // This is the epetra view of the world
   ///////////////////////////////////////////////////////////

   //! Get the <code>Epetra_Map</code> for this response, map is constructed lazily.
   virtual Teuchos::RCP<Epetra_MultiVector> buildEpetraDerivative() const
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
   virtual void setDerivative(const Teuchos::RCP<Epetra_MultiVector> & derivative)
   {
     TEUCHOS_ASSERT(!vectorIsDistributed());
     TEUCHOS_ASSERT(localSizeRequired()==1);
     TEUCHOS_ASSERT(supportsDerivative());
     TEUCHOS_ASSERT(eMap_!=Teuchos::null);

     derivative_ = Thyra::create_MultiVector(derivative,getDerivativeVectorSpace());
   }
#endif

   // This is the Thyra view of the world
   ///////////////////////////////////////////////////////////

   //! Get the <code>Epetra_Map</code> for this response, map is constructed lazily.
   virtual Teuchos::RCP<Thyra::MultiVectorBase<double> > buildDerivative() const
   {
     TEUCHOS_ASSERT(!vectorIsDistributed());
     TEUCHOS_ASSERT(localSizeRequired()==1);
     TEUCHOS_ASSERT(supportsDerivative());
     return Thyra::createMember(*getDerivativeVectorSpace());
   }

   /** Set the vector (to be filled) for this response. This must be
     * constructed from the vector space returned by <code>getMap</code>.
     */
   virtual void setDerivative(const Teuchos::RCP<Thyra::MultiVectorBase<double> > & derivative)
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
#ifdef PANZER_HAVE_EPETRA_STACK
   mutable Teuchos::RCP<const Epetra_Map> eMap_;
#endif

   Teuchos::RCP<Thyra::MultiVectorBase<double> > derivative_;
};

template < >
class ResponseMESupport_Default<panzer::Traits::Tangent> : public ResponseMESupportBase<panzer::Traits::Tangent> {
public:
  typedef panzer::Traits::Tangent EvalT;

    ResponseMESupport_Default(const std::string & responseName,MPI_Comm comm)
     : ResponseMESupportBase<EvalT>(responseName), useEpetra_(false),
#ifdef PANZER_HAVE_EPETRA_STACK
      eComm_(comm),
#endif
     useThyra_(false)
   {
     tComm_ = Teuchos::rcp(new Teuchos::MpiComm<Thyra::Ordinal>(Teuchos::opaqueWrapper(comm)));
   }

   virtual ~ResponseMESupport_Default() {}

   //! What is the number of values you need locally
   virtual std::size_t localSizeRequired() const = 0;

   //! Is the vector distributed (or replicated)
   virtual bool vectorIsDistributed() const = 0;

#ifdef PANZER_HAVE_EPETRA_STACK
   // This is the epetra view of the world
   ///////////////////////////////////////////////////////////

   //! Get the <code>Epetra_Map</code> for this response, map is constructed lazily.
   Teuchos::RCP<const Epetra_Map> getMap() const {
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

   /** Set the vector (to be filled) for this response. This must be
     * constructed from the vector space returned by <code>getMap</code>.
     */
   void setVector(const Teuchos::RCP<Epetra_MultiVector> & destVec) {
     TEUCHOS_TEST_FOR_EXCEPTION(useThyra_,std::logic_error,
                                "Reponse field \"" << this->getName() << "\" has previously been initialized as a "
                                "Thyra object, now trying to initalize as a Epetra! Error!");
     eVector_ = destVec;
     useEpetra_ = true;
   }
#endif

   // This is the Thyra view of the world
   ///////////////////////////////////////////////////////////

   //! Get the vector space for this response, vector space is constructed lazily.
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > getVectorSpace() const {
     TEUCHOS_TEST_FOR_EXCEPTION(useEpetra_,std::logic_error,
                                "Reponse field \"" << this->getName() << "\" has previously been initialized as an "
                                "Epetra object, now trying to initalize as a Thyra object! Error!");
     // lazily build the space and return it
     if(vSpace_==Teuchos::null) {
       if(this->vectorIsDistributed())
         vSpace_ = Thyra::defaultSpmdVectorSpace<double>(tComm_,this->localSizeRequired(),-1);
       else
         vSpace_ = Thyra::locallyReplicatedDefaultSpmdVectorSpace<double>(tComm_,this->localSizeRequired());
     }
     return vSpace_;
   }

   /** Set the vector (to be filled) for this response. This must be
     * constructed from the vector space returned by <code>getVectorSpace</code>.
     */
   void setVector(const Teuchos::RCP<Thyra::MultiVectorBase<double> > & destVec) {
     TEUCHOS_TEST_FOR_EXCEPTION(useEpetra_,std::logic_error,
                                "Reponse field \"" << this->getName() << "\" has previously been initialized as an "
                                "Epetra object, now trying to initalize as a Thyra object! Error!");
     tVector_ = destVec;
     useThyra_ = true;
   }

protected:
   //! Get the teuchos comm object
   Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > getComm() const { return tComm_; }

   //! Is Epetra the right vector
   bool useEpetra() const { return useEpetra_; }

   //! Is Thyra the right vector
   bool useThyra() const { return useThyra_; }

#ifdef PANZER_HAVE_EPETRA_STACK
   //! Access the epetra vector
   Epetra_MultiVector & getEpetraMultiVector() const {
     TEUCHOS_ASSERT(useEpetra());
     return *eVector_;
   }
#endif

   //! Access the thyra vector
   Thyra::ArrayRCP< Thyra::ArrayRCP<double> > getThyraMultiVector() const {
     TEUCHOS_ASSERT(useThyra());
     const int num_col = tVector_->domain()->dim();
     Thyra::ArrayRCP< Thyra::ArrayRCP<double> > data(num_col);
     for (int i=0; i<num_col; ++i)
       Teuchos::rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(tVector_->col(i),true)->getNonconstLocalData(Teuchos::outArg(data[i]));
     return data;
  }

  //! Return the number of columns in the multivector
  int numDeriv() const {
#ifdef PANZER_HAVE_EPETRA_STACK
    if (useEpetra())
      return eVector_->NumVectors();
    else
#endif
      return tVector_->domain()->dim();
  }


private:
   // hide these methods
   ResponseMESupport_Default();
   ResponseMESupport_Default(const ResponseMESupport_Default<EvalT> &);

   bool useEpetra_;
#ifdef PANZER_HAVE_EPETRA_STACK
   Epetra_MpiComm eComm_;
   mutable Teuchos::RCP<const Epetra_Map> map_;
   Teuchos::RCP<Epetra_MultiVector> eVector_;
#endif

   bool useThyra_;
   mutable Teuchos::RCP<const Thyra::VectorSpaceBase<double> > vSpace_;
   Teuchos::RCP<Thyra::MultiVectorBase<double> > tVector_;
   Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > tComm_;
};

#ifdef Panzer_BUILD_HESSIAN_SUPPORT

template < >
class ResponseMESupport_Default<panzer::Traits::Hessian> : public ResponseMESupportBase<panzer::Traits::Hessian> {
public:

   ResponseMESupport_Default(const std::string & responseName,MPI_Comm comm,
                             const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > & derivVecSpace=Teuchos::null)
     : ResponseMESupportBase<panzer::Traits::Hessian>(responseName), derivVecSpace_(derivVecSpace)
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
   Teuchos::RCP<Thyra::MultiVectorBase<double> > getDerivative() const
   { return derivative_; }

   //! Get the <code>Epetra_Map</code> for this response, map is constructed lazily.
   virtual Teuchos::RCP<Thyra::MultiVectorBase<double> > buildDerivative() const
   {
     TEUCHOS_ASSERT(!vectorIsDistributed());
     TEUCHOS_ASSERT(localSizeRequired()==1);
     TEUCHOS_ASSERT(supportsDerivative());
     return Thyra::createMember(*getDerivativeVectorSpace());
   }

   /** Set the vector (to be filled) for this response. This must be
     * constructed from the vector space returned by <code>getMap</code>.
     */
   virtual void setDerivative(const Teuchos::RCP<Thyra::MultiVectorBase<double> > & derivative)
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
   ResponseMESupport_Default(const ResponseMESupport_Default<panzer::Traits::Hessian> &);

   Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > tComm_;
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > derivVecSpace_;

   Teuchos::RCP<Thyra::MultiVectorBase<double> > derivative_;
};

#endif

}

#include "Panzer_ResponseMESupport_Default_impl.hpp"

#endif
