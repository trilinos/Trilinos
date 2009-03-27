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
     * passed in to this object. It is assume that this Epetra_Operator
     *
     * Build this preconditioner from an Epetra_Operator 
     * passed in to this object. It is assume that this Epetra_Operator
     * will be a EpetraOperatorWrapper object, so the block Thyra components
     * can be easily extracted.
     *
     * \param[in] A The Epetra source operator. (Should be a EpetraOperatorWrapper!)
     */
   virtual void buildPreconditioner(const Epetra_Operator & A);

protected:
   EpetraBlockPreconditioner(); 
   EpetraBlockPreconditioner(const EpetraBlockPreconditioner &); 

   Teuchos::RCP<const BlockPreconditionerFactory> preconFactory_;
};

} // end namespace Epetra
} // end namespace PB

#endif
