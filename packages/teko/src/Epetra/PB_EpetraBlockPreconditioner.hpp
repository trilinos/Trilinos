#ifndef __PB_EpetraBlockPreconditioner_hpp__
#define __PB_EpetraBlockPreconditioner_hpp__

#include "PB_BlockPreconditionerFactory.hpp"
#include "Epetra/PB_EpetraInverseOpWrapper.hpp"

namespace PB {
namespace Epetra {

/** \brief A single Epetra wrapper for all the BlockPreconditioners.
  *
  * This class uses the Thyra based preconditioner factories to
  * build an Epetra_Operator that behaves like a preconditioner.
  * This is done by using the BlockPreconditionerFactory, and letting
  * it build whatever preconditioner is neccessary. Thus the Epetra
  * "layer" is just a single class that handles any generic
  * BlockPreconditionerFactory.
  */
class EpetraBlockPreconditioner : public EpetraInverseOpWrapper {
public:
   /** \brief Constructor that takes the BlockPreconditionerFactory that will
     *        build the preconditioner.
     *
     * Constructor that takes the BlockPreconditionerFactory that will
     * build the preconditioner.
     */
   EpetraBlockPreconditioner(const Teuchos::RCP<const BlockPreconditionerFactory> & bfp); 

   /** \brief Build this preconditioner from an Epetra_Operator 
     * passed in to this object.
     *
     * Build this preconditioner from an Epetra_Operator 
     * passed in to this object. It is assumed that this Epetra_Operator
     * will be a EpetraOperatorWrapper object, so the block Thyra components
     * can be easily extracted.
     *
     * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
     */
   virtual void buildPreconditioner(const Epetra_Operator & A);

   /** \brief Build this preconditioner from an Epetra_Operator 
     * passed in to this object. It is assumed that this Epetra_Operator
     *
     * Build this preconditioner from an Epetra_Operator 
     * passed in to this object. It is assumed that this Epetra_Operator
     * will be a EpetraOperatorWrapper object, so the block Thyra components
     * can be easily extracted.
     *
     * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
     * \param[in] mv A vector that was used to build the source operator.
     */
   virtual void buildPreconditioner(const Epetra_Operator & A,const Epetra_MultiVector & mv);

   /** \brief Rebuild this preconditioner from an Epetra_Operator passed
     * in this to object. 
     *
     * Rebuild this preconditioner from an Epetra_Operator passed
     * in this to object.  If <code>buildPreconditioner</code> has not been called
     * the preconditioner will be built instead. Otherwise efforts are taken
     * to only rebuild what is neccessary. Also, it is assumed that this Epetra_Operator
     * will be an EpetraOperatorWrapper object, so the block Thyra components
     * can be easily extracted.
     *
     * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
     */
   virtual void rebuildPreconditioner(const Epetra_Operator & A);

   /** \brief Rebuild this preconditioner from an Epetra_Operator passed
     * in this to object. 
     *
     * Rebuild this preconditioner from an Epetra_Operator passed
     * in this to object.  If <code>buildPreconditioner</code> has not been called
     * the preconditioner will be built instead. Otherwise efforts are taken
     * to only rebuild what is neccessary. Also, it is assumed that this Epetra_Operator
     * will be an EpetraOperatorWrapper object, so the block Thyra components
     * can be easily extracted.
     *
     * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
     * \param[in] mv A vector that was used to build the source operator.
     */
   virtual void rebuildPreconditioner(const Epetra_Operator & A,const Epetra_MultiVector & mv);

   /** Try to get a <code>Teko::BlockPreconditionerState</code> object. This method
     * attempts to cast its internal representation of a preconditioner 
     * object to a <code>Teko::BlockPreconditioner</code> object.  If it suceeds a 
     * state object is returned.  Otherwise, <code>Teuchos::null</code> is returned.
     *
     * \returns Get the state object associated with this preconditioner.
     *          If it doesn't exist for this type of preconditioner factory
     *          this method returns null.
     */
   virtual Teuchos::RCP<BlockPreconditionerState> getPreconditionerState();

   /** Try to get a <code>Teko::BlockPreconditionerState</code> object. This method
     * attempts to cast its internal representation of a preconditioner 
     * object to a <code>Teko::BlockPreconditioner</code> object.  If it suceeds a 
     * state object is returned.  Otherwise, <code>Teuchos::null</code> is returned.
     *
     * \returns Get the state object associated with this preconditioner.
     *          If it doesn't exist for this type of preconditioner factory
     *          this method returns null.
     */
   virtual Teuchos::RCP<const BlockPreconditionerState> getPreconditionerState() const;

protected:
   EpetraBlockPreconditioner(); 
   EpetraBlockPreconditioner(const EpetraBlockPreconditioner &); 

   Teuchos::RCP<const BlockPreconditionerFactory> preconFactory_;
   Teuchos::RCP<Thyra::PreconditionerBase<double> > preconObj_;
};

} // end namespace Epetra
} // end namespace PB

#endif
