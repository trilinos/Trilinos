// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MiniEM_FullMaxwellPreconditionerFactory_Augmentation.hpp"

#include "Teko_BlockLowerTriInverseOp.hpp"
#include "Teko_BlockUpperTriInverseOp.hpp"

#include "Teko_SolveInverseFactory.hpp"

#include "Thyra_DiagonalLinearOpBase.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"

#include "Teuchos_as.hpp"
#include "Teuchos_Time.hpp"

#include "Teko_TpetraHelpers.hpp"
#include "Teko_TpetraOperatorWrapper.hpp"

#include "Thyra_TpetraLinearOp.hpp"
#include "Panzer_NodeType.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_LinearObjContainer.hpp"
#include "Panzer_ThyraObjContainer.hpp"

#include "Thyra_DefaultDiagonalLinearOp.hpp"

#include "MiniEM_Utils.hpp"

using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

namespace mini_em {


////////////////////////////////////////////////////
// FullMaxwellPreconditionerFactory_Augmentation  //
////////////////////////////////////////////////////

Teko::LinearOp FullMaxwellPreconditionerFactory_Augmentation::buildPreconditionerOperator(Teko::BlockedLinearOp & blo, Teko::BlockPreconditionerState & /* state */) const
{
   Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("MaxwellPreconditioner::build")));

   // Output stream for debug information
   Teuchos::RCP<Teuchos::FancyOStream> debug = Teuchos::null;
   if (doDebug)
     debug = Teko::getOutputStream();

   // Check that system is right size
   int rows = Teko::blockRowCount(blo);
   int cols = Teko::blockColCount(blo);
   TEUCHOS_ASSERT(rows==cols);
   TEUCHOS_ASSERT(rows==2);

   // Notation:
   // 0 - Hgrad
   // 1 - Hcurl
   // 2 - HDiv

   // M_k(a) - mass matrix on space k with weight a
   // D_k - derivative from space k to k+1

   // The block matrix is
   //
   // | Q_B  K   |
   // | Kt   Q_E |
   //
   // where
   // Q_B = 1/dt * M_2(1)
   // K   = M_2(1) * D_1
   // Kt  = - D_1^T * M_2(1/mu)
   // Q_E = 1/dt * M_1(eps)

   // Extract the blocks
   Teko::LinearOp Q_B   = Teko::getBlock(0,0,blo);
   Teko::LinearOp K     = Teko::getBlock(0,1,blo);
   Teko::LinearOp Kt    = Teko::getBlock(1,0,blo);
   Teko::LinearOp Q_E   = Teko::getBlock(1,1,blo);

   // nodal mass matrix
   Teko::LinearOp Q_rho = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix AUXILIARY_NODE"));

   // Set up the Schur complement
   Teko::LinearOp S_E;
   {
     Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("MaxwellPreconditioner: Schur complement"));
     S_E = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("SchurComplement AUXILIARY_EDGE"));
   }

   /////////////////////////////////////////////////
   // Debug and matrix dumps                      //
   /////////////////////////////////////////////////

   if (dump) {
     writeOut("Q_B.mm",*Q_B);
     writeOut("K.mm",*K);
     writeOut("Kt.mm",*Kt);
     writeOut("Q_E.mm",*Q_E);
     writeOut("Q_rho.mm",*Q_rho);
   }
   describeMatrix("Q_B",*Q_B,debug);
   describeMatrix("K",*K,debug);
   describeMatrix("Kt",*Kt,debug);
   describeMatrix("Q_E",*Q_E,debug);
   describeMatrix("Q_rho",*Q_rho,debug);


   /////////////////////////////////////////////////
   // Set up inverses for sub-blocks              //
   /////////////////////////////////////////////////

   // Inverse of B mass matrix
   Teko::LinearOp invQ_B;
   {
     Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("MaxwellPreconditioner: Inverse Q_B"));
     // Are we building a solver or a preconditioner?
     if (useAsPreconditioner) {
       invQ_B = Teko::buildInverse(*invLib.getInverseFactory("Q_B Preconditioner"),Q_B);
     } else {
       Teko::LinearOp invDiagQ_B = Teko::buildInverse(*invLib.getInverseFactory("Q_B Preconditioner"),Q_B);
       describeMatrix("invDiagQ_B",*invDiagQ_B,debug);
       invQ_B = Teko::buildInverse(*invLib.getInverseFactory("Q_B Solve"),Q_B, invDiagQ_B);
     }
   }

   // Solver for S_E
   Teko::LinearOp invS_E;
   {
     Teuchos::TimeMonitor tm1(*Teuchos::TimeMonitor::getNewTimer("MaxwellPreconditioner: Solver S_E"));
     // Get auxiliary operators for gradient and nodal mass matrix
     Teko::LinearOp G;
     if (use_discrete_gradient_) {
       Teko::LinearOp T = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Gradient"));
       if (dump)
         writeOut("DiscreteGradient.mm",*T);
       describeMatrix("DiscreteGradient",*T,debug);
       G = Teko::explicitMultiply(Q_E,Teko::scale(-1.0,T));
     } else {
       G = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Weak Gradient"));
     }
     if (dump)
       writeOut("WeakGradient.mm",*G);
     describeMatrix("WeakGradient",*G,debug);

     Teko::LinearOp Gt = Teko::explicitTranspose(G);
     // Compute grad-div term
     Teko::LinearOp invDiagQ_rho = Teko::getInvDiagonalOp(Q_rho,Teko::AbsRowSum);
     Teko::LinearOp GGt = Teko::explicitMultiply(G,invDiagQ_rho,Gt);

     // Rescale such that grad-div is large enough to fix curl-curl null-space while not dominating
     double scaling = Teko::infNorm(Q_E)/Teko::infNorm(GGt);

     // Augmented Schur complement and its inverse
     if (doDebug)
       *debug << "Adding up T_E" << std::endl;
     Teko::LinearOp T_E = Teko::explicitAdd(S_E, Thyra::scale(scaling,GGt));
     if (debug != Teuchos::null)
       *debug << "Added up T_E" << std::endl;
     *Teko::getOutputStream() << "Building T_E inverse operator" << std::endl;
     Teko::LinearOp invT_E = Teko::buildInverse(*invLib.getInverseFactory("T_E Solve"),T_E);

     // Correction term
     Teko::LinearOp Z_E = Thyra::add(Q_E, Thyra::scale(scaling,GGt));

     // if (dump)
     //   writeOut("Z_E.mm",*Z_E);
     describeMatrix("Z_E",*Z_E,debug);

     // Mass inverse - diagonal approximation
     Teko::LinearOp invQ_E = Teko::getInvDiagonalOp(Q_E,Teko::AbsRowSum);

     // Approximate inverse of S_E
     invS_E = Teko::multiply(invQ_E,Z_E,invT_E);
   }


   /////////////////////////////////////////////////
   // Build block  inverse matrices               //
   /////////////////////////////////////////////////

   {
     Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("MaxwellPreconditioner: Block preconditioner"));

     // Inverse blocks
     std::vector<Teko::LinearOp> diag(2);
     diag[0] = invQ_B;
     diag[1] = invS_E;

     // Upper tri blocks
     Teko::BlockedLinearOp U = Teko::createBlockedOp();
     Teko::beginBlockFill(U,rows,rows);
     Teko::setBlock(0,0,U,Q_B);
     Teko::setBlock(1,1,U,S_E);
     Teko::setBlock(0,1,U,K);
     Teko::endBlockFill(U);

     Teko::LinearOp invU = Teko::createBlockUpperTriInverseOp(U,diag);

     if (!useAsPreconditioner) {
       Teko::BlockedLinearOp invL = Teko::createBlockedOp();
       Teko::LinearOp id_B = Teko::identity(Teko::rangeSpace(Q_B));
       Teko::LinearOp id_E = Teko::identity(Teko::rangeSpace(Q_E));
       Teko::beginBlockFill(invL,rows,rows);
       Teko::setBlock(0,0,invL,id_B);
       Teko::setBlock(1,0,invL,Teko::multiply(Thyra::scale(-1.0, Kt), invQ_B));
       Teko::setBlock(1,1,invL,id_E);
       Teko::endBlockFill(invL);

       return Teko::multiply(invU, Teko::toLinearOp(invL));
     } else
       return invU;
   }
}

//! Initialize from a parameter list
void FullMaxwellPreconditionerFactory_Augmentation::initializeFromParameterList(const Teuchos::ParameterList & pl)
{
   /////////////////////
   // Solver options  //
   // //////////////////

   params = pl;

   use_discrete_gradient_ = params.get("Use discrete gradient",false);
   dump                   = params.get("Dump",false);
   doDebug                = params.get("Debug",false);
   useAsPreconditioner    = params.get("Use as preconditioner",false);

   //////////////////////////////////
   // Set up sub-solve factories   //
   //////////////////////////////////

   // New inverse lib to add inverse factories to
   invLib = *getInverseLibrary();

   // Q_B solve
   if (pl.isParameter("Q_B Solve")) {
     Teuchos::ParameterList cg_pl = pl.sublist("Q_B Solve");
     invLib.addInverse("Q_B Solve",cg_pl);
   }

   // Q_B preconditioner
   Teuchos::ParameterList Q_B_prec_pl = pl.sublist("Q_B Preconditioner");
   invLib.addStratPrecond("Q_B Preconditioner",
                          Q_B_prec_pl.get<std::string>("Prec Type"),
                          Q_B_prec_pl.sublist("Prec Types"));

   // T_E solve
   Teuchos::ParameterList T_E_pl = pl.sublist("T_E Solve");
   invLib.addInverse("T_E Solve",T_E_pl);

}

} // namespace mini_em
