#ifndef __PB_BlockPreconditionerFactory_hpp__
#define __PB_BlockPreconditionerFactory_hpp__

#include "Teuchos_ParameterListAcceptor.hpp"

// Thyra includes
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_LinearOpSourceBase.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"

// Teko includes
#include "PB_Utilities.hpp"
#include "PB_InverseLibrary.hpp"
#include "PB_CloneFactory.hpp"

namespace Teko {

using Thyra::LinearOpBase;
using Thyra::DefaultPreconditioner;
using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::ParameterList;

/** \brief Skeleton implementation of a state object for block
  *        preconditioners.
  *
  * A skeleton implementation of a state object for block
  * preconditioners. This implementation takes a parameter list.
  * However, it is easy to override this class and fill it
  * with what ever is neccessary to build the preconditioner.
  * This is essentially an emtpy "bag" for filling.
  */
class BlockPreconditionerState : public Teuchos::ParameterListAcceptor {
public:
   //! \name Default and copy constructors
   //@{ 
   BlockPreconditionerState() : isInitialized_(false) {}
   BlockPreconditionerState(const BlockPreconditionerState & src) 
      : paramList_(rcp(new ParameterList(*src.paramList_))), 
        srcVector_(src.srcVector_),
        inverses_(src.inverses_),
        isInitialized_(src.isInitialized_) {}
   //@}

   //! \name for ParameterListAcceptor
   //@{

   //! Set parameters from a parameter list and return with default values.
   void setParameterList(const RCP<ParameterList> & paramList); 

   //! Get the parameter list that was set using setParameterList().
   RCP<ParameterList> getNonconstParameterList();

   //! Unset the parameter list that was set using setParameterList(). 
   RCP<ParameterList> unsetParameterList();
 
   //@}

   /** Has this state object been initialized for the particular operator
     * being used?
     */
   virtual bool isInitialized() const
   { return isInitialized_; }

   /** Set that this state object has been initialized for this operator.
     */
   virtual void setInitialized(bool init=true) 
   { isInitialized_ = init; }

   //! Set the vector associated with this operator (think nonlinear system)
   virtual void setSourceVector(const Teko::BlockedMultiVector & srcVec)
   { srcVector_ = srcVec; }

   //! Set the vector associated with this operator (think nonlinear system)
   virtual const Teko::BlockedMultiVector getSourceVector() const
   { return srcVector_; }

   //! Add a named inverse to the state object
   virtual void addInverse(const std::string & name,const Teko::InverseLinearOp & ilo)
   { inverses_[name] = ilo; }

   //! Get a named inverse from the state object
   virtual Teko::InverseLinearOp getInverse(const std::string & name) const
   { std::map<std::string,Teko::InverseLinearOp>::const_iterator itr;
     itr =  inverses_.find(name);
     if(itr==inverses_.end()) return Teuchos::null; 
     return itr->second; }

protected:
   //! for ParameterListAcceptor
   RCP<ParameterList>          paramList_;

   //! Store a source vector
   Teko::BlockedMultiVector srcVector_;

   //! Store a map of inverse linear operators
   std::map<std::string,Teko::InverseLinearOp> inverses_;

   //! Stores the initialization state 
   bool isInitialized_;
};

/** \brief An extension of the <code>Thyra::DefaultPreconditioner</code>
  *        class with some specializations useful for use within Teko.
  *
  * An extension of the <code>Thyra::DefaultPreconditioner</code>
  * class with some specializations useful for use within Teko. This includes
  * having facilities to store the source vector and the state object.
  */
class BlockPreconditioner : public DefaultPreconditioner<double> {
public:
   //! \name Constructors based from Thyra::DefaultPreconditioner
   //@{
   BlockPreconditioner()
      : DefaultPreconditioner<double>() {}
   BlockPreconditioner(const RCP<LinearOpBase<double> > & leftPrecOp,
                       const RCP<LinearOpBase<double> > & rightPrecOp)
      : DefaultPreconditioner<double>(leftPrecOp,rightPrecOp) {}
   BlockPreconditioner(const RCP<const LinearOpBase<double> > & leftPrecOp,
                       const RCP<const LinearOpBase<double> > & rightPrecOp)
      : DefaultPreconditioner<double>(leftPrecOp,rightPrecOp) {}
   BlockPreconditioner(const RCP<LinearOpBase<double> > & unspecifiedPrecOp)
      : DefaultPreconditioner<double>(unspecifiedPrecOp) {}
   BlockPreconditioner(const RCP<const LinearOpBase<double> > & unspecifiedPrecOp)
      : DefaultPreconditioner<double>(unspecifiedPrecOp) {}
   //@}

   /** Set the vector associated with this operator (think nonlinear system)
     *
     * \param[in] srcVec The source vector associated with this preconditioner.
     */
   virtual void setSourceVector(const RCP<Thyra::MultiVectorBase<double> > & srcVec)
   { if(srcVec!=Teuchos::null) state_->setSourceVector(Teko::toBlockedMultiVector(srcVec));
     else                      state_->setSourceVector(Teuchos::null); }

   /** Set the state object associated with this preconditioner
     *
     * \param[in] state The state object to use
     */
   virtual void setStateObject(const RCP<BlockPreconditionerState> & state)
   { state_ = state; }

   /** Get the state object associated with this preconditioner
     *
     * \returns The state object used by this preconditioner
     */
   virtual const RCP<BlockPreconditionerState> getStateObject()
   { return state_; }

   /** Get the state object associated with this preconditioner
     *
     * \returns The state object used by this preconditioner
     */
   virtual const RCP<const BlockPreconditionerState> getStateObject() const
   { return state_; }

protected:
   //! User defined state for this preconditioner
   RCP<BlockPreconditionerState> state_;
};


/** \brief Abstract class which block preconditioner factories in Teko
  *        should be based on.
  *
  * Abstract class which block preconditioner factories in Teko should
  * be based on. All that is needed is the implementation of 
  * "buildPreconditionerOperator".
  */
class BlockPreconditionerFactory 
   : public virtual Thyra::PreconditionerFactoryBase<double> {
public:

   /** \brief Function that is called to build the preconditioner
     *        for the linear operator that is passed in.
     *
     * This function builds a preconditioner based on the passed
     * in BlockedLinearOp. 
     *
     * \param[in] blo   Source linear operator that is to be preconditioned.
     * \param[in] state An object associated with this operator to store
     *                  the preconditioner state.
     * 
     * \returns The preconditioner as a linear operator (i.e. to perform
    *           a matrix-vector operation simply call "apply").
     */
   virtual LinearOp buildPreconditionerOperator(BlockedLinearOp & blo,BlockPreconditionerState & state) const = 0;

   /** \brief Function that permits the construction of an arbitrary
     *        BlockPreconditionerState object.
     *
     * Function that permits the construction of an arbitrary
     * BlockPreconditionerState object. If the basic state object,
     * which takes a parameter list, is sufficient the default behavior
     * does precisely what is needed. Otherwise, an author of a
     * PreconditionerFactory would need to reimplement this method to 
     * return a new state object.
     *
     * \returns A state object associated with this factory.
     */
   virtual RCP<BlockPreconditionerState> buildPreconditionerState() const
   { return rcp(new BlockPreconditionerState()); }

   /** \brief This function builds the internals of the preconditioner factory
     *        from a parameter list.
     *        
     * This function builds the internals of the preconditioner factory
     * from a parameter list. Furthermore, it allows a preconditioner factory
     * developer to easily add a factory to the build system. This function
     * is required for building a preconditioner from a parameter list.
     *
     * \param[in] settings Parmaeter list to use as the internal settings
     *
     * \note The default implementation does nothing.
     */
   virtual void initializeFromParameterList(const Teuchos::ParameterList & settings)
   { }

   /** \brief Request the additional parameters this preconditioner factory
     *        needs. 
     *
     * Request the additonal parameters needed by this preconditioner factory.
     * The parameter list will have a set of fields that can be filled with 
     * the requested values. These fields include all requirements, even those
     * of the sub-solvers if there are any.  Once correctly filled the object
     * can be updated by calling the updateRequestedParameters with the filled
     * parameter list.
     *
     * \returns A parameter list with the requested parameters.
     *
     * \note The default implementation returns Teuchos::null.
     */
   virtual Teuchos::RCP<Teuchos::ParameterList> getRequestedParameters() const
   { return Teuchos::null; }
   
   /** \brief Update this object with the fields from a parameter list.
     *
     * Update the requested fields using a parameter list. This method is
     * expected to pair with the getRequestedParameters method (i.e. the fields
     * requested are going to be update using this method).
     *
     * \param[in] pl Parameter list containing the requested parameters.
     *
     * \returns If the method succeeded (found all its required parameters) this
     *          method returns true, otherwise it returns false.
     *
     * \note The default implementation returns true (it does nothing!).
     */
   virtual bool updateRequestedParameters(const Teuchos::ParameterList & pl)
   { return true; }
 
   //! @name Ignore me methods.
   //@{
 
   virtual LinearOp buildPreconditionerOperator(LinearOp & blo,BlockPreconditionerState & state) const;

   //! is this operator compatiable with the preconditioner factory?
   bool isCompatible(const Thyra::LinearOpSourceBase<double> &fwdOpSrc) const;

   //! create an instance of the preconditioner
   RCP<Thyra::PreconditionerBase<double> > createPrec() const;

   /** \brief initialize a newly created preconditioner object
     *
     * Initialize a newly created preconditioner object. For use with
     * nonlinear solvers.
     *
     * \param[in] fwdOpSrc Forward operator to be preconditioned
     * \param[in] solnVec Vector associated with this linear operator.
     * \param[in,out] precOp Return location for the preconditioner
     * \param[in] supportSolveUse Thyra information (?)
     */
   void initializePrec(const RCP<const Thyra::LinearOpSourceBase<double> > & fwdOpSrc,
                       const RCP<const Thyra::MultiVectorBase<double> > & solnVec,
                       Thyra::PreconditionerBase<double> * precOp,
                       const Thyra::ESupportSolveUse supportSolveUse) const;

   //! initialize a newly created preconditioner object
   void initializePrec(const RCP<const Thyra::LinearOpSourceBase<double> > & fwdOpSrc,
                       Thyra::PreconditionerBase<double> * precOp,
                       const Thyra::ESupportSolveUse supportSolveUse) const;

   //! wipe clean a already initialized preconditioner object
   void uninitializePrec(Thyra::PreconditionerBase<double> * prec, 
                         RCP<const Thyra::LinearOpSourceBase<double> > * fwdOpSrc,
                         Thyra::ESupportSolveUse *supportSolveUse) const;

   // for ParameterListAcceptor

   //! Set parameters from a parameter list and return with default values.
   void setParameterList(const RCP<ParameterList> & paramList); 

   //! Get the parameter list that was set using setParameterList().
   RCP< ParameterList > getNonconstParameterList();

   //! Unset the parameter list that was set using setParameterList(). 
   RCP< ParameterList > unsetParameterList();
   //@}

   //! Set the inverse library used by this preconditioner factory
   void setInverseLibrary(const RCP<const InverseLibrary> & il);

   //! Get the inverse library used by this preconditioner factory
   RCP<const InverseLibrary> getInverseLibrary() const;

protected:
   //! for ParameterListAcceptor
   RCP<ParameterList>          paramList_;

private:
   //! Inverse library to be used by this factory
   RCP<const InverseLibrary> inverseLibrary_;

public:

   /** \brief Builder function for creating preconditioner factories (yes
     *        this is a factory factory).
     *
     * Builder function for creating preconditioner factories (yes
     * this is a factory factory).
     * 
     * \param[in] name     String name of factory to build
     * \param[in] settings Parameter list describing the parameters for the
     *                     factory to build
     * \param[in] invLib   Inverse library for the factory to use.
     *
     * \returns If the name is associated with a preconditioner
     *          a pointer is returned, otherwise Teuchos::null is returned.
     */
   static RCP<BlockPreconditionerFactory> 
   buildPreconditionerFactory(const std::string & name, 
                              const Teuchos::ParameterList & settings,
                              const RCP<const InverseLibrary> & invLib=Teuchos::null);

   /** \brief Add a preconditioner factory to the builder. This is done using the
     *        clone pattern. 
     *
     * Add a preconditioner factory to the builder. This is done using the
     * clone pattern. If your class does not support the Cloneable interface then
     * you can use the AutoClone class to construct your object.
     *
     * \note If this method is called twice with the same string, the latter clone pointer
     *       will be used.
     *
     * \param[in] name String to associate with this object
     * \param[in] clone Pointer to Cloneable object
     */
   static void addPreconditionerFactory(const std::string & name,const RCP<Cloneable> & clone);

private:

   //! for creating the preconditioner factories objects
   static CloneFactory<BlockPreconditionerFactory> precFactoryBuilder_;

   //! This is where the default objects are put into the precFactoryBuilder_
   static void initializePrecFactoryBuilder();
};

} // end namespace Teko

#endif
