/*
 * Author: Zhen Wang
 * Email: wangz@ornl.gov
 *        zhen.wang@alum.emory.edu
 */

#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp_def.hpp"
#include "Thyra_DefaultIdentityLinearOp_decl.hpp"

#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

#include "Teuchos_Time.hpp"

#include "Teko_Utilities.hpp"

#include "Teko_InvModALStrategy.hpp"
#include "Teko_ModALPreconditionerFactory.hpp"

using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcp_const_cast;

namespace Teko
{

namespace NS
{

// Empty constructor.
InvModALStrategy::InvModALStrategy() :
      invFactoryA_(Teuchos::null), invFactoryS_(Teuchos::null),
      pressureMassMatrix_(Teuchos::null), gamma_(0.05),
      scaleType_(Diagonal), isSymmetric_(true)
{
}

// If there is only one InverseFactory, use it for all solves.
InvModALStrategy::InvModALStrategy(const Teuchos::RCP<InverseFactory> & factory) :
      invFactoryA_(factory), invFactoryS_(factory),
      pressureMassMatrix_(Teuchos::null), gamma_(0.05),
      scaleType_(Diagonal), isSymmetric_(true)
{
}

// If there are two InverseFactory's...
InvModALStrategy::InvModALStrategy(const Teuchos::RCP<InverseFactory> & invFactoryA,
      const Teuchos::RCP<InverseFactory> & invFactoryS) :
      invFactoryA_(invFactoryA), invFactoryS_(invFactoryS),
      pressureMassMatrix_(Teuchos::null), gamma_(0.05),
      scaleType_(Diagonal), isSymmetric_(true)
{
}

// If there are two InverseFactory's...
InvModALStrategy::InvModALStrategy(const Teuchos::RCP<InverseFactory> & invFactory,
      LinearOp & pressureMassMatrix) :
      invFactoryA_(invFactory), invFactoryS_(invFactory),
      pressureMassMatrix_(pressureMassMatrix), gamma_(0.05),
      scaleType_(Diagonal), isSymmetric_(true)
{
}

// If there are two InverseFactory's...
InvModALStrategy::InvModALStrategy(const Teuchos::RCP<InverseFactory> & invFactoryA,
      const Teuchos::RCP<InverseFactory> & invFactoryS, LinearOp & pressureMassMatrix) :
      invFactoryA_(invFactoryA), invFactoryS_(invFactoryS),
      pressureMassMatrix_(pressureMassMatrix), gamma_(0.05),
      scaleType_(Diagonal), isSymmetric_(true)
{
}

// Return "inverses".
LinearOp InvModALStrategy::getInvA11p(BlockPreconditionerState & state) const
{
   return state.getInverse("invA11p");
}

LinearOp InvModALStrategy::getInvA22p(BlockPreconditionerState & state) const
{
   return state.getInverse("invA22p");
}

LinearOp InvModALStrategy::getInvA33p(BlockPreconditionerState & state) const
{
   return state.getInverse("invA33p");
}

LinearOp InvModALStrategy::getInvS(BlockPreconditionerState & state) const
{
   return state.getInverse("invS");
}

// Set pressure mass matrix.
void InvModALStrategy::setPressureMassMatrix(const LinearOp & pressureMassMatrix)
{
   pressureMassMatrix_ = pressureMassMatrix;
}

// Set gamma.
void InvModALStrategy::setGamma(double gamma)
{
   TEUCHOS_ASSERT(gamma > 0.0);
   gamma_ = gamma;
}

void InvModALStrategy::buildState(const BlockedLinearOp & alOp,
      BlockPreconditionerState & state) const
{
   Teko_DEBUG_SCOPE("InvModALStrategy::buildState", 10);

   ModALPrecondState * modALState = dynamic_cast<ModALPrecondState*>(&state);
   TEUCHOS_ASSERT(modALState != NULL);

   // if necessary save state information
   if(not modALState->isInitialized())
   {
      Teko_DEBUG_EXPR(Teuchos::Time timer(""));

      {
         // construct operators
         Teko_DEBUG_SCOPE("ModAL::buildState: Initializing state object", 1);
         Teko_DEBUG_EXPR(timer.start(true));

         initializeState(alOp, modALState);

         Teko_DEBUG_EXPR(timer.stop());
         Teko_DEBUG_MSG("ModAL::buildState: BuildOpsTime = " << timer.totalElapsedTime(), 1);
      }

      {
         // Build the inverses
         Teko_DEBUG_SCOPE("ModAL::buildState: Computing inverses", 1);
         Teko_DEBUG_EXPR(timer.start(true));

         computeInverses(alOp, modALState);

         Teko_DEBUG_EXPR(timer.stop());
         Teko_DEBUG_MSG("ModAL::buildState: BuildInvTime = " << timer.totalElapsedTime(), 1);
      }
   }
}

// Initialize the state object using the ALOperator.
void InvModALStrategy::initializeState(const BlockedLinearOp & alOp,
      ModALPrecondState *state) const
{
   Teko_DEBUG_SCOPE("InvModALStrategy::initializeState", 10);

   // Extract sub-matrices from blocked linear operator.
   int dim = blockRowCount(alOp) - 1;
   TEUCHOS_ASSERT(dim == 2 || dim == 3);

   LinearOp lpA11 = getBlock(0, 0, alOp);
   LinearOp lpA22 = getBlock(1, 1, alOp);
   LinearOp lpA33, lpB1, lpB2, lpB3, lpB1t, lpB2t, lpB3t, lpC;

   // 2D problem.
   if(dim == 2)
   {
      lpB1  = getBlock(2, 0, alOp);
      lpB2  = getBlock(2, 1, alOp);
      lpB1t = getBlock(0, 2, alOp);
      lpB2t = getBlock(1, 2, alOp);
      lpC = getBlock(2, 2, alOp);
   }
   // 3D problem.
   else if(dim == 3)
   {
      lpA33 = getBlock(2, 2, alOp);
      lpB1  = getBlock(3, 0, alOp);
      lpB2  = getBlock(3, 1, alOp);
      lpB3  = getBlock(3, 2, alOp);
      lpB1t = getBlock(0, 3, alOp);
      lpB2t = getBlock(1, 3, alOp);
      lpB3t = getBlock(2, 3, alOp);
      lpC   = getBlock(3, 3, alOp);
   }

   // For problems using stabilized finite elements,
   // lpB1t, lpB2t and lpB3t are added linear operators. Extract original operators.
   LinearOp B1t = (rcp_dynamic_cast<const Thyra::DefaultAddedLinearOp<double> >(lpB1t))->getOp(0);
   LinearOp B2t = (rcp_dynamic_cast<const Thyra::DefaultAddedLinearOp<double> >(lpB2t))->getOp(0);
   LinearOp B3t;
   if(dim == 3)
   {
      B3t = (rcp_dynamic_cast<const Thyra::DefaultAddedLinearOp<double> >(lpB3t))->getOp(0);
   }

   //std::cout << Teuchos::describe(*lpC, Teuchos::VERB_EXTREME) << std::endl;
   // Check whether the finite elements are stable or not.
   state->isStabilized_ =(not isZeroOp(lpC));
   //state->isStabilized_ = false;
   //std::cout << state->isStabilized_ << std::endl;

   state->pressureMassMatrix_ = pressureMassMatrix_;
   // If pressure mass matrix is not set, use identity.
   if(state->pressureMassMatrix_ == Teuchos::null)
   {
      Teko_DEBUG_MSG("InvModALStrategy::initializeState(): Build identity type \""
            << getDiagonalName(scaleType_) << "\"", 1);
      state->invPressureMassMatrix_ = Thyra::identity<double>(lpB1->range());
   }
   // If the inverse of the pressure mass matrix is not set,
   // build it from the pressure mass matrix.
   else if(state->invPressureMassMatrix_ == Teuchos::null)
   {
      Teko_DEBUG_MSG("ModAL::initializeState(): Build Scaling <mass> type \""
            << getDiagonalName(scaleType_) << "\"", 1);
      state->invPressureMassMatrix_ = getInvDiagonalOp(pressureMassMatrix_, scaleType_);
   }
   // Else "invPressureMassMatrix_" should be set and there is no reason to rebuild it
   state->gamma_ = gamma_;
   //S_ = scale(1.0/gamma_, pressureMassMatrix_);
   std::cout << Teuchos::describe(*(state->invPressureMassMatrix_), Teuchos::VERB_EXTREME) << std::endl;

   // Build state variables: B_1^T*W^{-1}*B_1, A11p, etc.
   // B_1^T*W^{-1}*B_1 may not change so save it in the state.
   if(state->B1tMpB1_ == Teuchos::null)
      state->B1tMpB1_ = explicitMultiply(B1t, state->invPressureMassMatrix_, lpB1, state->B1tMpB1_);
   // Get the(1,1) block of the non-augmented matrix.
   // Recall alOp is augmented. So lpA11 = A11 + gamma B_1^T W^{-1} B_1.
   // Cast lpA11 as an added linear operator and get the first item.
   LinearOp A11 = (rcp_dynamic_cast<const Thyra::DefaultAddedLinearOp<double> >(lpA11))->getOp(0);
   state->A11p_ = explicitAdd(A11, scale(state->gamma_, state->B1tMpB1_), state->A11p_);
   //std::cout << Teuchos::describe(*(state->B1tMpB1_), Teuchos::VERB_EXTREME) << std::endl;
   Teko_DEBUG_MSG("Computed A11p", 10);

   if(state->B2tMpB2_ == Teuchos::null)
      state->B2tMpB2_ = explicitMultiply(B2t, state->invPressureMassMatrix_, lpB2, state->B2tMpB2_);
   LinearOp A22 = (rcp_dynamic_cast<const Thyra::DefaultAddedLinearOp<double> >(lpA22))->getOp(0);
   state->A22p_ = explicitAdd(A22, scale(state->gamma_, state->B2tMpB2_), state->A22p_);
   Teko_DEBUG_MSG("Computed A22p", 10);

   if(dim == 3)
   {
      if(state->B3tMpB3_ == Teuchos::null)
         state->B3tMpB3_ = explicitMultiply(B3t, state->invPressureMassMatrix_, lpB3, state->B3tMpB3_);
      LinearOp A33 = (rcp_dynamic_cast<const Thyra::DefaultAddedLinearOp<double> >(lpA33))->getOp(0);
      state->A33p_ = explicitAdd(A33, scale(state->gamma_, state->B3tMpB3_), state->A33p_);
      Teko_DEBUG_MSG("Computed A33p", 10);
   }

   // Inverse the Schur complement.
   if(state->isStabilized_)
   {
      if(state->S_ == Teuchos::null)
      {
         state->S_ = explicitAdd(scale(-1.0, lpC), scale(1.0/state->gamma_, pressureMassMatrix_), state->S_);
      }
      Teko_DEBUG_MSG("Computed S", 10);
   }

   state->setInitialized(true);
}

// Compute inverses.
void InvModALStrategy::computeInverses(const BlockedLinearOp & alOp,
      ModALPrecondState *state) const
{
   int dim = blockRowCount(alOp) - 1;
   TEUCHOS_ASSERT(dim == 2 || dim == 3);

   Teko_DEBUG_SCOPE("InvModALStrategy::computeInverses", 10);
   Teko_DEBUG_EXPR(Teuchos::Time invTimer(""));

   //(re)build the inverse of A11
   Teko_DEBUG_MSG("ModAL::computeInverses(): Building inv(A11)", 1);
   Teko_DEBUG_EXPR(invTimer.start(true));

   InverseLinearOp invA11p = state->getInverse("invA11p");
   if(invA11p == Teuchos::null)
   {
      invA11p = buildInverse(*invFactoryA_, state->A11p_);
      state->addInverse("invA11p", invA11p);
   }
   else
   {
      rebuildInverse(*invFactoryA_, state->A11p_, invA11p);
   }

   Teko_DEBUG_EXPR(invTimer.stop());
   Teko_DEBUG_MSG("ModAL::computeInverses GetInvA11 = " << invTimer.totalElapsedTime(), 1);

   //(re)build the inverse of A22
   Teko_DEBUG_MSG("ModAL::computeInverses(): Building inv(A22)", 2);
   Teko_DEBUG_EXPR(invTimer.start(true));

   InverseLinearOp invA22p = state->getInverse("invA22p");
   if(invA22p == Teuchos::null)
   {
      invA22p = buildInverse(*invFactoryA_, state->A22p_);
      state->addInverse("invA22p", invA22p);
   }
   else
   {
      rebuildInverse(*invFactoryA_, state->A22p_, invA22p);
   }

   Teko_DEBUG_EXPR(invTimer.stop());
   Teko_DEBUG_MSG("ModAL::computeInverses(): GetInvA22 = " << invTimer.totalElapsedTime(), 2);

   //(re)build the inverse of A33
   if(dim == 3)
   {
      Teko_DEBUG_MSG("ModAL::computeInverses Building inv(A33)", 3);
      Teko_DEBUG_EXPR(invTimer.start(true));

      InverseLinearOp invA33p = state->getInverse("invA33p");
      if(invA33p == Teuchos::null)
      {
         invA33p = buildInverse(*invFactoryA_, state->A33p_);
         state->addInverse("invA33p", invA33p);
      }
      else
      {
         rebuildInverse(*invFactoryA_, state->A33p_, invA33p);
      }

      Teko_DEBUG_EXPR(invTimer.stop());
      Teko_DEBUG_MSG("ModAL::computeInverses GetInvA33 = " << invTimer.totalElapsedTime(), 3);
   }

   //(re)build the inverse of S
   Teko_DEBUG_MSG("ModAL::computeInverses Building inv(S)", 4);
   Teko_DEBUG_EXPR(invTimer.start(true));

   // There are two ways to "invert" S.
   // The following method construct invS by InverseFactory invFactoryS_.
   // The other way is to use diagonal approximation,
   // which is done in ModALPreconditionerFactory.cpp.
   if(state->isStabilized_)
   {
      InverseLinearOp invS = state->getInverse("invS");
      if(invS == Teuchos::null)
      {
         invS = buildInverse(*invFactoryS_, state->S_);
         state->addInverse("invS", invS);
      }
      else
      {
         rebuildInverse(*invFactoryS_, state->S_, invS);
      }
   }

   Teko_DEBUG_EXPR(invTimer.stop());
   Teko_DEBUG_MSG("ModAL::computeInverses GetInvS = " << invTimer.totalElapsedTime(), 4);

}

} // end namespace NS

} // end namespace Teko
