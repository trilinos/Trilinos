#ifndef __PB_Block2x2PreconditionerFactory_hpp__
#define __PB_Block2x2PreconditionerFactory_hpp__

#include "Teuchos_ParameterListAcceptor.hpp"

#include "Thyra_LinearOpBaseDecl.hpp"
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_LinearOpSourceBase.hpp"
#include "Thyra_PreconditionerBase.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"

#include "PB_Block2x2Strategy.hpp"

namespace PB {

/** \brief Construct a preconditioner using a LDU dcomposition of a block
 *  2x2 matrix.
 *
 * This produces a preconditioner using the block-LDU decomposition of
 * the matrix. The general assumption made is that the matrix is 2x2 
 * and the block factorization can be constructed (i.e. assumptions about
 * the invertability of some blocks). The pattern used, and the one you
 * should follow if you want to use this software is
 *
 * \f$
 * A = \left[ 
 * \begin{array}{cc}
 * F & U \\
 * L & D
 * \end{array}
 * \right]
 * = \left[ 
 * \begin{array}{cc}
 * I & 0  \\
 * L F^{-1} & I
 * \end{array}
 * \right]
 * \left[ 
 * \begin{array}{cc}
 * F & 0  \\
 * 0 & -S
 * \end{array}
 * \right]
 * \left[ 
 * \begin{array}{cc}
 * I &  F^{-1} U \\
 * 0 & I
 * \end{array}
 * \right]
 * \f$
 *
 * where the Schur complement is \f$ S=-D+L F^{-1} U \f$ .
 *
 * In order to facilate using this class in a nonlinear solve (or for a 
 * time-dependent problem) the additional abstraction of a ``Strategy''
 * has been added. This strategy, abstractly represented as the Block2x2Strategy,
 * provides the inv(F) and inv(S) operators. Typical usage for this class
 * is to build a Block2x2Strategy and pass it into the primary constructor. 
 * Additional constructors are provided both for convenience and to ease
 * adoption. Underneath the hood all these constructors do is invoke the
 * corresponding strategy object.
 *
 * For example, assume that you have the particularly nice case that
 * your approximations of inv(F) and inv(S) are indedent of the source
 * operator. Then, one way to instantiate a Block2x2PreconditionerFactory
 * is

   <code>
      RCP<LinearOpBase<double> > invF = buildInvF(...);\n
      RCP<LinearOpBase<double> > invS = buildInvF(...);\n
      RCP<Block2x2PreconditionerFactory> precFact = rcp(new Block2x2PreconditionerFactory(invF,invS));
   </code>

   Now using the strategy constructor, an entirely equivalent factory
   object can be constructed by

   <code>
      RCP<LinearOpBase<double> > invF = buildInvF(...);\n
      RCP<LinearOpBase<double> > invS = buildInvF(...);\n
      RCP<Block2x2Strateghy> precStrat = rcp(new StaticBlock2x2Strategy(invF,invS));\n
      RCP<Block2x2PreconditionerFactory> precFact = rcp(new Block2x2PreconditionerFactory(precStrat));
   </code>
 
   Notice that the StaticBlock2x2Strategy takes the same objects 
   as the original constructor, it acts as an intermediary to tell the 
   Block2x2PreconditionerFactory what those operators are.
 **/
class Block2x2PreconditionerFactory 
   : public Thyra::PreconditionerFactoryBase<double> {
   public:
      //@{
 
      /** @brief Build a simple static Block2x2 preconditioner */
      Block2x2PreconditionerFactory(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invF,
                                    const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invPschur);

      /** @brief Constructor that permits the most generality in build inv(F) and
        *        inv(Schur).
        *
        * Constructor that permits the most generality in build inv(F) and inv(Schur).
        *
        * @param[in] strategy  Strategy object that takes a 2x2 block matrix and
        *                      and constructs the inv(F) and inv(Schur) objects.
        */
      Block2x2PreconditionerFactory(const Teuchos::RCP<const Block2x2Strategy> & strategy);

      /** @brief Gets the underlying block strategy object.
        * 
        * Gets the underlying block strategy object.
        */
      const Teuchos::RCP<const Block2x2Strategy> GetBlock2x2Strategy() const
      { return invOpsStrategy_; }
 
      // @}

      /** @name Methods inherited from Thyra::PreconditionerFactoryBase */
      //@{

      /** @brief Is this operator compatiable with the preconditioner factory? */
      bool isCompatible(const Thyra::LinearOpSourceBase<double> &fwdOpSrc) const;

      /** @brief Create an instance of the preconditioner. */
      Teuchos::RCP<Thyra::PreconditionerBase<double> > createPrec() const;

      /** @brief Initialize a newly created preconditioner object */
      void initializePrec(const Teuchos::RCP<const Thyra::LinearOpSourceBase<double> > & fwdOpSrc,
                          Thyra::PreconditionerBase<double> * precOp,
                          const Thyra::ESupportSolveUse supportSolveUse=Thyra::SUPPORT_SOLVE_UNSPECIFIED) const;

      /** @brief wipe clean a already initialized preconditioner object */
      void uninitializePrec(Thyra::PreconditionerBase<double> * prec, 
                            Teuchos::RCP<const Thyra::LinearOpSourceBase<double> > * fwdOpSrc=NULL,
                            Thyra::ESupportSolveUse *supportSolveUse=NULL) const;
 
      //@}

      /** @name Methods inherited from Teuchos::ParameterListAcceptor (unused) */
      //@{

      // Set parameters from a parameter list and return with default values.
      void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList); 

      // Get the parameter list that was set using setParameterList().
      Teuchos::RCP< Teuchos::ParameterList > getNonconstParameterList();

      // Unset the parameter list that was set using setParameterList(). 
      Teuchos::RCP< Teuchos::ParameterList > unsetParameterList();

      //@}

   protected:
      /** @brief An empty constructor restricted for use by child classes
        *        only.
        * 
        * An empty constructor restricted for use by child classes only.
        */
      Block2x2PreconditionerFactory() {}
                                   
      /** @brief Build a 2x2 inverse operator based on the LDU decomposition 
       *
       * Build a 2x2 inverse operator based on the LDU decomposition. This is the
       * core functionality of this class. Its provided in the protected section to
       * offer easy and convenient access to its children.
       * 
       * @param[in] A    Source operator to build LDU matrix from.
       * @param[in] invF An approximate inverse of inv(F).
       * @param[in] invS An approximate inverse of inv(S).
       *
       * @return An operator that will produce the approximate inverse
       *         of A using an LDU decomposition and inv(F) and inv(S). The
       *         operator will apply inv(F) twice, and inv(S) once.
       */
      const Teuchos::RCP<const Thyra::LinearOpBase<double> > 
      build2x2InverseOperator(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A, 
                              const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invF, 
                              const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invS) const;

      // for ParameterListAcceptor
      /** Member to support the parameter list: unused */
      mutable Teuchos::RCP<Teuchos::ParameterList>  validPL_; 

      /** Member to support the parameter list: unused */
      Teuchos::RCP<Teuchos::ParameterList>          paramList_; 

   private:
  
      /** @brief Strategy object that constructs inv(F) and inv(Schur) 
       *
       * Strategy object that constructs inv(F) and inv(Schur). 
       */
      Teuchos::RCP<const Block2x2Strategy> invOpsStrategy_; 
};

} // end namespace PB

#endif
