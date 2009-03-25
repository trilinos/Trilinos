#ifndef __PB_BlockJacobiPreconditionerFactory_hpp__
#define __PB_BlockJacobiPreconditionerFactory_hpp__

#include "Teuchos_ParameterListAcceptor.hpp"

#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_LinearOpSourceBase.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_PreconditionerBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"

#include "PB_BlockInvDiagonalStrategy.hpp"

namespace PB {

using namespace Teuchos;
using namespace Thyra;

//! BlockJacobiPreconditionerFactory: Factory for building a block Jacobi preconditioner.
/*!
 Builds matrix 
\f[
    M^{-1} = \left( \begin{array}{cccc}
                   D_0^{-1} & 0        & 0         & 0      \\
                    0       & D_1^{-1} & 0         & 0      \\
                    0       & 0        & D_2^{-1}  & 0      \\
                    0       & 0        & 0         & \cdots
                    \end{array} \right)
\f]
              
as an approximate inverse of
\f[
    A = \left( \begin{array}{cccc}
          D_0     & U_{01}  & U_{02} & \dots \\
          L_{10}  & D_1     & U_{12} & \dots \\
          L_{20}  & L_{21}  & D_2    & \dots \\
           .      &  .      & .      & \dots
               \end{array} \right)
\f]
*/

class BlockJacobiPreconditionerFactory 
   : public virtual PreconditionerFactoryBase<double> {
   public:
      //! @name Constructors.

      /*! Construct a PreconditionerFactory assuming a specific block
          \f$2\times2\f$ matrix. This case is a common one.
      */ 
      BlockJacobiPreconditionerFactory(const RCP<const LinearOpBase<double> > & invD0,
                                       const RCP<const LinearOpBase<double> > & invD1);

      /*! The most flexible BlockJacobiPreconditionerFactory constructor.
          Pass in a generally defined BlockInvDiagonalStrategy to use the
          full generality of this class.
      */
      BlockJacobiPreconditionerFactory(const RCP<const BlockInvDiagonalStrategy> & strategy);
      //@}

      // for PreconditionerFactoryBase
      ///////////////////////////////////////////////////////////////////////
 
      //! @name Public methods.
      //@{
      //! is this operator compatiable with the preconditioner factory?
      bool isCompatible(const LinearOpSourceBase<double> &fwdOpSrc) const;

      //! create an instance of the preconditioner
      RCP<PreconditionerBase<double> > createPrec() const;

      //! initialize a newly created preconditioner object
      void initializePrec(const RCP<const LinearOpSourceBase<double> > & fwdOpSrc,
                          PreconditionerBase<double> * precOp,
                          const ESupportSolveUse supportSolveUse=SUPPORT_SOLVE_UNSPECIFIED) const;

      //! wipe clean a already initialized preconditioner object
      void uninitializePrec(PreconditionerBase<double> * prec, 
                            RCP<const LinearOpSourceBase<double> > * fwdOpSrc=NULL,
                            ESupportSolveUse *supportSolveUse=NULL) const;

      // for ParameterListAcceptor
      ///////////////////////////////////////////////////////////////////////

      //! Set parameters from a parameter list and return with default values.
      void setParameterList(const RCP<ParameterList> & paramList); 

      //! Get the parameter list that was set using setParameterList().
      RCP< ParameterList > getNonconstParameterList();

      //! Unset the parameter list that was set using setParameterList(). 
      RCP< ParameterList > unsetParameterList();
      //@}

   protected:
      //! for use by inherited classes only
      BlockJacobiPreconditionerFactory() {}

      //! some members
      Teuchos::RCP<const BlockInvDiagonalStrategy> invOpsStrategy_;

      //! for ParameterListAcceptor
      mutable RCP<ParameterList>  validPL_;
      RCP<ParameterList>          paramList_;
};

} // end namespace PB

#endif
