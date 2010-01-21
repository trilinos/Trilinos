#ifndef __Teko_InverseFactory_hpp__
#define __Teko_InverseFactory_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"

#include "Teko_Utilities.hpp"

namespace Teko {

/** \brief Abstract class for building an inverse operator
  *
  * Abstract class for building an inverse operator. It pairs
  * with a linear operator and gives you a new operator that
  * behaves like its inverse.
  */
class InverseFactory {
public:
   virtual ~InverseFactory() {}

   /** \brief Build an inverse operator
     *
     * Build the inverse operator using this factory.
     *
     * \param[in] linearOp Linear operator needing to be inverted.
     *
     * \returns New linear operator that functions as the inverse
     *          of <code>linearOp</code>.
     */
   virtual InverseLinearOp buildInverse(const LinearOp & linearOp) const = 0;

   /** \brief Pass in an already constructed inverse operator. Update
     *        the inverse operator based on the new source operator.
     *
     * Pass in an already constructed inverse operator. Update
     * the inverse operator based on the new source operator.
     *
     * \param[in]     source Source operator to be inverted.
     * \param[in,out] dest   Pre constructed inverse operator to be
     *                        rebuilt using the <code>source</code>
     *                        object.
     */
   virtual void rebuildInverse(const LinearOp & source,InverseLinearOp & dest) const = 0;

   /** \brief A function that permits inspection of the parameters used to create
     *        this object.
     *
     * A function that permits inspection of the parameters used to create this
     * object. Useful for determining defaults and settings used.
     *
     * \returns A list used to parameterize this object.
     */
   virtual Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const = 0;

   /** Return a string that describes this factory */
   virtual std::string toString() const = 0;

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
};

class SolveInverseFactory : public InverseFactory {
public:
   //! \name Constructors
   //@{ 
   
   /** \brief Constructor that takes a Thyra solve factory and 
     *        makes it look like an InverseFactory
     *
     * Constructor that takes a Thyra solve factory and 
     * makes it look like an InverseFactory.
     * 
     * \param[in] lowsFactory Thyra LineaerOpWithSolveFactoryBase used for building 
     *                        the inverse.
     */
   SolveInverseFactory(const Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > & lowsFactory);

   //! Copy constructor
   SolveInverseFactory(const SolveInverseFactory & siFactory);
   //@}

   virtual ~SolveInverseFactory() {}

   /** \brief Build an inverse operator
     *
     * Build the inverse operator using this factory.
     *
     * \param[in] linearOp Linear operator needing to be inverted.
     *
     * \returns New linear operator that functions as the inverse
     *          of <code>linearOp</code>.
     */
   virtual InverseLinearOp buildInverse(const LinearOp & linearOp) const;

   /** \brief Pass in an already constructed inverse operator. Update
     *        the inverse operator based on the new source operator.
     *
     * Pass in an already constructed inverse operator. Update
     * the inverse operator based on the new source operator.
     *
     * \param[in]     source Source operator to be inverted.
     * \param[in,out] dest   Pre constructed inverse operator to be
     *                        rebuilt using the <code>source</code>
     *                        object.
     */
   virtual void rebuildInverse(const LinearOp & source,InverseLinearOp & dest) const;

   /** \brief A function that permits inspection of the parameters used to create
     *        this object.
     *
     * A function that permits inspection of the parameters used to create this
     * object. Useful for determining defaults and settings used.
     *
     * \returns A list used to parameterize this object.
     */
   virtual Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;

   /** Return a string that describes this factory */
   virtual std::string toString() const { return lowsFactory_->description(); }

protected:
   Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory_;

private:
   // hide me!
   SolveInverseFactory();
};

class PreconditionerInverseFactory : public InverseFactory {
public:
   //! \name Constructors
   //@{ 
   
   /** \brief Constructor that takes a Thyra solve factory and 
     *        makes it look like an InverseFactory
     *
     * Constructor that takes a Thyra solve factory and 
     * makes it look like an InverseFactory.
     * 
     * \param[in] precFactory Thyra PreconditionerFactoryBase used for building 
     *                        the inverse.
     */
   PreconditionerInverseFactory(const Teuchos::RCP<Thyra::PreconditionerFactoryBase<double> > & precFactory);

   /** \brief Constructor that takes a Thyra solve factory and 
     *        makes it look like an InverseFactory. This constructor
     *        also permits the passing of an "Extra Parameters" parameter
     *        list.
     *
     * Constructor that takes a Thyra solve factory and 
     * makes it look like an InverseFactory.  This constructor
     * also permits the passing of an "Extra Parameters" parameter
     * list to be used and updated through the "RequestedParameters" function.
     * 
     * \param[in] precFactory Thyra PreconditionerFactoryBase used for building 
     *                        the inverse.
     * \param[in] xtraParam Parameter list containing extra parameters.
     */
   PreconditionerInverseFactory(const Teuchos::RCP<Thyra::PreconditionerFactoryBase<double> > & precFactory,
                                const Teuchos::RCP<const Teuchos::ParameterList> & xtraParam);

   //! Copy constructor
   PreconditionerInverseFactory(const PreconditionerInverseFactory & pFactory);
   //@}

   virtual ~PreconditionerInverseFactory() {}

   /** \brief Build an inverse operator
     *
     * Build the inverse operator using this factory. This also tacks
     * on extra data to the RCP called "prec". This is the
     * PreconditionerBase object, and it is need for <code>rebuildInverse</code>
     *
     * \param[in] linearOp Linear operator needing to be inverted.
     *
     * \returns New linear operator that functions as the inverse
     *          of <code>linearOp</code>.
     */
   virtual InverseLinearOp buildInverse(const LinearOp & linearOp) const;

   /** \brief Pass in an already constructed inverse operator. Update
     *        the inverse operator based on the new source operator.
     *
     * Pass in an already constructed inverse operator. Update
     * the inverse operator based on the new source operator. This
     * method assumes the <code>dest</code> object also contains
     * the associated PreconditionerBase object as "prec" as extra
     * data in the RCP.
     *
     * \param[in]     source Source operator to be inverted.
     * \param[in,out] dest   Pre constructed inverse operator to be
     *                        rebuilt using the <code>source</code>
     *                        object.
     */
   virtual void rebuildInverse(const LinearOp & source,InverseLinearOp & dest) const;

   /** \brief A function that permits inspection of the parameters used to create
     *        this object.
     *
     * A function that permits inspection of the parameters used to create this
     * object. Useful for determining defaults and settings used.
     *
     * \returns A list used to parameterize this object.
     */
   virtual Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;

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
   virtual Teuchos::RCP<Teuchos::ParameterList> getRequestedParameters() const;
   
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
   virtual bool updateRequestedParameters(const Teuchos::ParameterList & pl);

   /** Return a string that describes this factory */
   virtual std::string toString() const { return precFactory_->description(); }

   /** Get the preconditioner factroy */
   Teuchos::RCP<const Thyra::PreconditionerFactoryBase<double> > getPrecFactory() const
   { return precFactory_; }

   /** Get the preconditioner factroy */
   Teuchos::RCP<Thyra::PreconditionerFactoryBase<double> > getPrecFactory()
   { return precFactory_; }

protected:
   Teuchos::RCP<Thyra::PreconditionerFactoryBase<double> > precFactory_;
   Teuchos::RCP<Teuchos::ParameterList> extraParams_;

private:
   // hide me!
   PreconditionerInverseFactory();
};

class StaticOpInverseFactory : public InverseFactory {
public:
   //! \name Constructors
   //@{ 
   
   /** \brief Constructor that takes a linear operator and
     *        uses it as a static inverse
     *
     * Constructor that takes a linear operator and
     * uses it as a static inverse
     * 
     * \param[in] inv Linear operator to use as the inverse.
     */
   StaticOpInverseFactory(const LinearOp inv) 
      : inverse_(inv) {}

   //! Copy constructor
   StaticOpInverseFactory(const StaticOpInverseFactory & saFactory)
      : inverse_(saFactory.inverse_) {}
   //@}

   virtual ~StaticOpInverseFactory() {}

   /** \brief Build an inverse operator
     *
     * Build the inverse operator using this factory. This also tacks
     * on extra data to the RCP called "prec". This is the
     * PreconditionerBase object, and it is need for <code>rebuildInverse</code>
     *
     * \param[in] linearOp Linear operator needing to be inverted.
     *
     * \returns New linear operator that functions as the inverse
     *          of <code>linearOp</code>.
     */
   virtual InverseLinearOp buildInverse(const LinearOp & linearOp) const
   { return Teuchos::rcp_const_cast<Thyra::LinearOpBase<double> >(inverse_); }

   /** \brief Pass in an already constructed inverse operator. Update
     *        the inverse operator based on the new source operator.
     *
     * Pass in an already constructed inverse operator. Update
     * the inverse operator based on the new source operator. This
     * method assumes the <code>dest</code> object also contains
     * the associated PreconditionerBase object as "prec" as extra
     * data in the RCP.
     *
     * \param[in]     source Source operator to be inverted.
     * \param[in,out] dest   Pre constructed inverse operator to be
     *                        rebuilt using the <code>source</code>
     *                        object.
     */
   virtual void rebuildInverse(const LinearOp & source,InverseLinearOp & dest) const
   { }

   /** \brief A function that permits inspection of the parameters used to create
     *        this object.
     *
     * A function that permits inspection of the parameters used to create this
     * object. Useful for determining defaults and settings used.
     *
     * \returns A list used to parameterize this object.
     */
   virtual Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const
   { return Teuchos::null; }

   /** Return a string that describes this factory */
   virtual std::string toString() const { return inverse_->description(); }

protected:
   Teko::LinearOp inverse_;

private:
   // hide me!
   StaticOpInverseFactory();
};

//! @name Functions for constructing and initializing solvers
//@{

/** Build an inverse operator using a factory and a linear operator
  *
  * \param[in] factory The inverse factory used to construct the inverse
  *                    operator
  * \param[in] A       Linear operator whose inverse is required
  *
  * \returns An (approximate) inverse operator is returned for the operator <code>A</code>.
  *
  * \relates InverseFactory
  */
InverseLinearOp buildInverse(const InverseFactory & factory,const LinearOp & A);

/** Using a prebuilt linear operator, use factory to build an inverse operator
  * given a new forward operator.
  *
  * \note This function sometimes fails depending on the underlying type
  *       of the inverse factory.  Use with caution.
  *
  * \param[in] factory The inverse factory used to construct the inverse
  *                    operator
  * \param[in] A       Linear operator whose inverse is required
  * \param[in] invA    The inverse operator that is to be rebuilt using
  *                    the <code>A</code> operator.
  *
  * \relates InverseFactory
  */
void rebuildInverse(const InverseFactory & factory, const LinearOp & A, InverseLinearOp & invA);

/** \brief Build an InverseFactory object from a ParameterList, as specified in Stratimikos.
  *
  * Build an InverseFactory object from a ParameterList, as specified in Stratimikos.
  * The specific inverse routine (either solver or preconditioner) to be chosen is specified
  * by a string.
  *
  * \note It is preferred that the <code>InverseLibrary</code> is used to construct an
  *       <code>InverseFactory</code> instead.
  *
  * \param[in] list ParameterList that describes the available solvers/preconditioners.
  * \param[in] type String saying which solver/preconditioner to use.
  *
  * \returns An inverse factory using the specified inverse operation.
  *
  * \relates InverseFactory
  */
Teuchos::RCP<InverseFactory> invFactoryFromParamList(const Teuchos::ParameterList & list,const std::string & type);

/** \brief Get a valid parameter list for the inverse factory class.
  *
  * Get a valid parameter list for the inverse factory class. This will
  * specify the set of parameters for each possible "inverse".
  *
  * \note It is preferred that the <code>InverseLibrary</code> is used 
  *       to get paramter lists for <code>InverseFactory</code> construction.
  *
  * \returns A parameter list is returned that is suitable to be passed
  *          to <code>invFactoryFromParamList</code>.
  *
  * \relates InverseFactory
  */
Teuchos::RCP<const Teuchos::ParameterList> invFactoryValidParameters();

//@}

} // end namespace Teko

#endif
