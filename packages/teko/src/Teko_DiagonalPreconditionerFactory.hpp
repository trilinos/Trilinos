#ifndef __Teko_DiagonalPreconditionerFactory_hpp__
#define __Teko_DiagonalPreconditionerFactory_hpp__

// Teko includes
#include "Teko_PreconditionerState.hpp"
#include "Teko_PreconditionerFactory.hpp"

class EpetraExt_PointToBlockDiagPermute; 


namespace Teko {

using Thyra::LinearOpBase;
using Thyra::DefaultPreconditioner;


/** Preconditioning factories must supply a 'State' class which
  * is where data specific to the preconditioner construction 
  * is stored. The constructor will be invoked elsewhere.
  */
class DiagonalPrecondState : public Teko::PreconditionerState {
public:
  DiagonalPrecondState();

  Teuchos::RCP<EpetraExt_PointToBlockDiagPermute> BDP_;
};



/** \brief Preconditioner factory that for (block) diagonals of explicit operators.
  *
  * Preconditioner factory that for (block) diagonals of explicit operators.
  * These operators need to be Epetra_CrsMatrices under the hood or this will bomb.
  * Uses EpetraExt_PointToBlockDiagPermute.
  */
class DiagonalPreconditionerFactory 
   : public virtual Teko::PreconditionerFactory {
public:

  //! @name Constructors.
  //@{
  
  /** Build an empty Gauss-Seidel preconditioner factory
   */
  DiagonalPreconditionerFactory();


  //@}
  
  //! Builds a preconditioner state object
  Teuchos::RCP<PreconditionerState> buildPreconditionerState() const;


  /** Create the diagonal preconditioner operator.
   */
  LinearOp buildPreconditionerOperator(LinearOp & lo,PreconditionerState & state) const;

  //! Initialize from a parameter list
  virtual void initializeFromParameterList(const Teuchos::ParameterList & pl);

protected: 
  //! some members
  mutable Teuchos::ParameterList List_;
  

};

} // end namespace Teko

#endif
