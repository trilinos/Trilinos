/*
 * Author: Zhen Wang
 * Email: wangz@ornl.gov
 *        zhen.wang@alum.emory.edu
 */

#ifndef __Teko_ModALPreconditionerFactory_hpp__
#define __Teko_ModALPreconditionerFactory_hpp__

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_InvModALStrategy.hpp"

namespace Teko
{

namespace NS
{

/** \brief Class for saving state variables for
 * ModALPreconditionerFactory.
 */
class ModALPrecondState : public BlockPreconditionerState
{
public:
   ModALPrecondState();

   LinearOp pressureMassMatrix_;
   LinearOp invPressureMassMatrix_;
   ModifiableLinearOp B1tMpB1_, B2tMpB2_, B3tMpB3_;
   ModifiableLinearOp A11p_, A22p_, A33p_;
   ModifiableLinearOp invA11p_, invA22p_, invA33p_, invS_;
   ModifiableLinearOp S_;

   double gamma_;
   bool isStabilized_;
};

/** \brief Modified augmented Lagrangian-based preconditioner
 *  for incompressible Navier-Stokes equations.
 *
 * This class implements sparse matrix vector multiplication
 * for augmented Lagrangian-based preconditioners.
 * Details can be found in the following papers:
 *
 * [1] Benzi, M. A. Olshanskii and Z. Wang,
 * Modified Augmented Lagrangian Preconditioners for the Incompressible Navier-Stokes Equations,
 * International Journal for Numerical Methods in Fluids, 66 (2011), pp. 486-508.
 *
 * [2] M. Benzi and Z. Wang,
 * Analysis of Augmented Lagrangian-Based Preconditioners
 * for the Steady Incompressible Navier-Stokes Equations,
 * SIAM J. Scientific Computing, 33 (2011), pp. 2761-2784.
 *
 * Suppose we are solving the following linear system:
 *
 * \f$
 * \left[ \begin{array}{cc}
 * A & B^T \\
 * B & -C
 * \end{array} \right]
 * \left[ \begin{array}{c}
 * u \\
 * p
 * \end{array} \right]
 * =
 * \left[ \begin{array}{c}
 * f \\
 * g
 * \end{array} \right].
 * \f$
 *
 * The equivalent augmented Lagrangian (AL) formulation is:
 *
 * \f$
 * \left[ \begin{array}{cc}
 * A + \gamma B^T W^{-1} B & B^T - \gamma B^T W^{-1} C \\
 * B & -C
 * \end{array} \right]
 * \left[ \begin{array}{c}
 * u \\
  * p
 * \end{array} \right]
 * =
 * \left[ \begin{array}{c}
 * f + \gamma B^T W^{-1} g \\
 * g
 * \end{array} \right].
 * \f$
 *
 * Here \f$ W \f$ can be take as the diagonal of the pressure
 * mass matrix and \f$ \gamma \f$ is a positive number.
 *
 * If
 * \f$ A = (A_{11}, A_{12}; A_{21}, A_{22})
 * \quad \mathrm{and} \quad
 * B = (B_1, B_2)
 * \f$,
 *
 * then the modified AL preconditioner is defined as
 *
 * \f$
 * \left[ \begin{array}{ccc}
 * A_{11} + \gamma B_1^T W^{-1} B_1^T & A_{12} + \gamma B_1^T W^{-1} B_2^T & B_1^T - \gamma B_1^T W^{-1} C \\
 * 0 & A_{22} + \gamma B_2^T W^{-1} B_2^T & B_1^T - \gamma B_2^T W^{-1} C \\
 * 0 & 0 & S
 * \end{array} \right].
 * \f$
 *
 * More details (analysis, numerical results) can be found in [2].
 */

class ModALPreconditionerFactory : public BlockPreconditionerFactory
{
public:

   ModALPreconditionerFactory();

   ModALPreconditionerFactory(const Teuchos::RCP<InverseFactory> & factory);

   ModALPreconditionerFactory(const Teuchos::RCP<InverseFactory> & invFactoryA,
         const Teuchos::RCP<InverseFactory> & invFactoryS);

   ModALPreconditionerFactory(const Teuchos::RCP<InverseFactory> & factory,
         LinearOp & pressureMassMatrix);

   ModALPreconditionerFactory(const Teuchos::RCP<InverseFactory> & invFactoryA,
         const Teuchos::RCP<InverseFactory> & invFactoryS,
         LinearOp & pressureMassMatrix);

   ModALPreconditionerFactory(const Teuchos::RCP<InvModALStrategy> & strategy);

   /** Set pressure mass matrix.
    *
    * \param[in] pressureMassMatrix
    *            Pressure mass matrix.
    */
   void
   setPressureMassMatrix(const LinearOp & pressureMassMatrix)
   {
      invOpsStrategy_->setPressureMassMatrix(pressureMassMatrix);
   }

   /** Set the augmentation parameter gamma.
    *
    * \param[in] gamma
    *            Augmentation paramter.
    */
   void
   setGamma(double gamma)
   {
      invOpsStrategy_->setGamma(gamma);
   }

   /**
    * Build modified AL preconditioner.
    *
    * \param[in] alOp
    *            The AL operator.
    * \param[in] state
    *            State object for storying reusable information about the AL operator.
    */
   virtual LinearOp
   buildPreconditionerOperator(BlockedLinearOp & alOp,
         BlockPreconditionerState & state) const;

   /**
    * Build the ModALPrecondState object.
    */
   virtual Teuchos::RCP<PreconditionerState>
   buildPreconditionerState() const
   {
      return Teuchos::rcp(new ModALPrecondState());
   }

protected:

   Teuchos::RCP<InvModALStrategy> invOpsStrategy_;

   bool isSymmetric_;
};

} // end namespace NS

} // end namespace Teko

#endif /* __Teko_ModALPreconditionerFactory_hpp__ */
