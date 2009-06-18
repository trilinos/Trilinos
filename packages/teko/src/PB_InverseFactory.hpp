#ifndef __PB_InverseFactory_hpp__
#define __PB_InverseFactory_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"

#include "PB_Utilities.hpp"

namespace PB {

/** \brief Abstract class for building an inverse operator
  *
  * Abstract class for building an inverse operator. It pairs
  * with a linear operator and gives you a new operator that
  * behaves like its inverse.
  */
class InverseFactory {
public:
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

   //! Copy constructor
   PreconditionerInverseFactory(const PreconditionerInverseFactory & pFactory);
   //@}

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

   /** Return a string that describes this factory */
   virtual std::string toString() const { return precFactory_->description(); }

protected:
   Teuchos::RCP<Thyra::PreconditionerFactoryBase<double> > precFactory_;

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
   PB::LinearOp inverse_;

private:
   // hide me!
   StaticOpInverseFactory();
};

//! @name Functions for constructing and initializing solvers
//@{
// typedef Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > InverseFactory;

//! Build an inverse operator using a factory and a linear operator
InverseLinearOp buildInverse(const InverseFactory & factory,const LinearOp & A);

/** Using a prebuilt linear operator, use factory to build an inverse operator
  * given a new forward operator.
  */
void rebuildInverse(const InverseFactory & factory, const LinearOp & A, InverseLinearOp & invA);

/** \brief Build an InverseFactory object from a ParameterList, as specified in Stratimikos.
  *
  * Build an InverseFactory object from a ParameterList, as specified in Stratimikos.
  * The specific inverse routine (either solver or preconditioner) to be chosen is specified
  * by a string.
  *
  * \param[in] list ParameterList that describes the available solvers/preconditioners.
  * \param[in] type String saying which solver/preconditioner to use.
  *
  * \returns An inverse factory using the specified inverse operation.
  */
Teuchos::RCP<const InverseFactory> invFactoryFromParamList(const Teuchos::ParameterList & list,const std::string & type);

/** \brief Get a valid parameter list for the inverse factory class.
  *
  * Get a valid parameter list for the inverse factory class. This will
  * specify the set of parameters for each possible "inverse".
  *
  * \returns A parameter list is returned that is suitable to be passed
  *          to <code>invFactoryFromParamList</code>.
  */
Teuchos::RCP<const Teuchos::ParameterList> invFactoryValidParameters();

//@}

} // end namespace PB

#endif
