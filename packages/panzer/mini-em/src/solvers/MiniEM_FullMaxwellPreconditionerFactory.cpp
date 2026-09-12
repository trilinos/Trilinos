// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MiniEM_FullMaxwellPreconditionerFactory.hpp"

#include "Teko_BlockLowerTriInverseOp.hpp"
#include "Teko_BlockUpperTriInverseOp.hpp"

#include "Teko_SolveInverseFactory.hpp"

#include "Teko_Utilities.hpp"
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
#include "MueLu_Maxwell1.hpp"
#include "MueLu_RefMaxwell.hpp"

#include "Thyra_DefaultDiagonalLinearOp.hpp"

#include "MiniEM_Utils.hpp"

using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

namespace mini_em {


///////////////////////////////////////
// FullMaxwellPreconditionerFactory  //
///////////////////////////////////////

Teko::LinearOp FullMaxwellPreconditionerFactory::buildPreconditionerOperator(Teko::BlockedLinearOp & blo, Teko::BlockPreconditionerState & /* state */) const
{
   typedef double Scalar;
   typedef int LocalOrdinal;
   typedef panzer::GlobalOrdinal GlobalOrdinal;
   typedef panzer::TpetraNodeType Node;

   Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("MaxwellPreconditioner::build")));

   // Check that system is right size
   int rows = Teko::blockRowCount(blo);
   int cols = Teko::blockColCount(blo);
   TEUCHOS_ASSERT(rows==cols);
   TEUCHOS_ASSERT(rows==2);

  auto rh = getRequestHandler();

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

   // S_E = Q_E - Kt * Q_B^-1 * K
   //     = 1/dt * M_1(eps) + dt * D_1^T * M_2(1/mu) * D_1
   //
   // -> curl-curl term scales like dt / mu
   //
   // for refmaxwell: Q_rho = M_0(mu / dt) so that the addon is:
   // M_1(1) * D_0 * M_0(mu / dt)^-1 * D_0^T * M_1(1)

   // Modify the system
   if (simplifyFaraday_) {
     RCP<Teuchos::FancyOStream> out = Teko::getOutputStream();
     *out << std::endl;
     *out << "*** WARNING ***" << std::endl;
     *out << "We are modifying the linear system. That's not a friendly thing to do." << std::endl;
     *out << std::endl;

     Teko::LinearOp Q_B  = Teko::getBlock(0,0,blo);
     Teko::LinearOp id_B = getIdentityMatrix(Q_B, 1/dt);
     Teko::LinearOp hoC  = rh->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Curl"));
     Teko::LinearOp Kt    = Teko::getBlock(1,0,blo);
     Teko::LinearOp Q_E   = Teko::getBlock(1,1,blo);
     blo->beginBlockFill(2,2);
     Teko::setBlock(0,0,blo,id_B);
     Teko::setBlock(0,1,blo,hoC);
     Teko::setBlock(1,0,blo,Kt);
     Teko::setBlock(1,1,blo,Q_E);
     blo->endBlockFill();
   }

   // Extract the blocks
   Teko::LinearOp Q_B   = Teko::getBlock(0,0,blo);
   Teko::LinearOp K     = Teko::getBlock(0,1,blo);
   Teko::LinearOp Kt    = Teko::getBlock(1,0,blo);
   Teko::LinearOp Q_E   = Teko::getBlock(1,1,blo);

   // discrete curl and its transpose
   Teko::LinearOp C;
   if (use_discrete_curl_) {
     C = rh->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Curl"));
   }

   // Set up the Schur complement
   Teko::LinearOp S_E;
   {
     Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("MaxwellPreconditioner: Schur complement"));
     if (simplifyFaraday_) {
       // Q_B is a scaled identy matrix
       // K is the discrete curl
       // Kt is the weak curl
       auto inv_Q_B = Teko::getInvDiagonalOp(Q_B, Teko::Diagonal);
       auto curl_curl = Teko::explicitMultiply(Kt, inv_Q_B, K);
       S_E = Teko::explicitAdd(Q_E, Teko::scale(-1.0, curl_curl), Teuchos::null);
     } else
       S_E = rh->request<Teko::LinearOp>(Teko::RequestMesg("SchurComplement AUXILIARY_EDGE"));
   }

   /////////////////////////////////////////////////
   // Debug and matrix dumps                      //
   /////////////////////////////////////////////////

   debugInfo("Q_B", Q_B);
   debugInfo("K", K);
   debugInfo("Kt", Kt);
   debugInfo("Q_E", Q_E);
   debugInfo("S_E", S_E);
   if (!C.is_null()) {
     debugInfo("DiscreteCurl", C);

     if (dump || doDebug) {
       Teko::LinearOp K2 = Teko::explicitScale(dt, Teko::explicitMultiply(Q_B,C));
       debugInfo("K2", K2);
       Teko::LinearOp diffK = Teko::explicitAdd(K, Teko::scale(-1.0,K2));
       debugInfo("diff", diffK);
       TEUCHOS_ASSERT(Teko::infNorm(diffK) < 1.0e-8 * Teko::infNorm(K));
     }
   }

   /////////////////////////////////////////////////
   // Set up inverses for sub-blocks              //
   /////////////////////////////////////////////////

   // Inverse of B mass matrix
   Teko::LinearOp invQ_B;
   if (!simplifyFaraday_) {
     Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("MaxwellPreconditioner: Inverse Q_B"));
     // Are we building a solver or a preconditioner?
     if (useAsPreconditioner) {
       invQ_B = Teko::buildInverse(*invLib.getInverseFactory("Q_B Preconditioner"),Q_B);
     } else {
       Teko::LinearOp invDiagQ_B = Teko::buildInverse(*invLib.getInverseFactory("Q_B Preconditioner"),Q_B);
       debugInfo("invDiagQ_B", invDiagQ_B);
       invQ_B = Teko::buildInverse(*invLib.getInverseFactory("Q_B Solve"),Q_B, invDiagQ_B);
     }
   }

   // Solver for S_E
   Teko::LinearOp invS_E;
   {
     Teuchos::TimeMonitor tm1(*Teuchos::TimeMonitor::getNewTimer("MaxwellPreconditioner: Solver S_E"));

     RCP<Teko::InverseFactory> S_E_prec_factory;
     if (S_E_prec_type_ != "none") {
       Teuchos::ParameterList S_E_prec_pl;
       S_E_prec_factory = invLib.getInverseFactory("S_E Preconditioner");
       S_E_prec_pl = *S_E_prec_factory->getParameterList();
       Teuchos::ParameterList& muelu_pl = S_E_prec_pl.sublist("Preconditioner Types").sublist(S_E_prec_type_);

       std::set<std::string> requiredUserData;
       std::set<std::string> optionalUserData;

       if (S_E_prec_type_ == "MueLuRefMaxwell")
         std::tie(requiredUserData, optionalUserData) = MueLu::RefMaxwell<Scalar, LocalOrdinal, GlobalOrdinal, Node>::requiredAndOptionalUserData(muelu_pl);
       else if (S_E_prec_type_ == "MueLuMaxwell1")
         std::tie(requiredUserData, optionalUserData) = MueLu::Maxwell1<Scalar, LocalOrdinal, GlobalOrdinal, Node>::requiredAndOptionalUserData(muelu_pl);

       // add nodal coordinates
       if (requiredUserData.contains("Coordinates")) {
         auto Coordinates = S_E_prec_pl.get<RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >("Coordinates");
         muelu_pl.sublist("user data").set("Coordinates",Coordinates);
       }

       // add discrete gradient
       if (requiredUserData.contains("Dk_1") || requiredUserData.contains("D0")) {
         auto T = rh->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Gradient"));
         debugInfo("DiscreteGradient", T);
         if (requiredUserData.contains("Dk_1"))
           muelu_pl.sublist("user data").set("Dk_1",T);
         if (requiredUserData.contains("D0"))
           muelu_pl.sublist("user data").set("D0",T);
       }

       // add edge mass matrix
       Teko::LinearOp Q_E_aux;
       if (requiredUserData.contains("Mk_one")) {
         Q_E_aux = rh->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix AUXILIARY_EDGE"));
         muelu_pl.sublist("user data").set("Mk_one",Q_E_aux);
       }

       // add weighted edge matrix
       if (requiredUserData.contains("M1_beta")) {
         Teko::LinearOp Q_E_aux_weighted;
         try {
           Q_E_aux_weighted = rh->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix weighted AUXILIARY_EDGE"));
         } catch (std::runtime_error&) {
           if (Q_E_aux.is_null())
             Q_E_aux = rh->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix AUXILIARY_EDGE"));
           Q_E_aux_weighted = Q_E_aux;
         }
         muelu_pl.sublist("user data").set("M1_beta",Q_E_aux_weighted);
       }

       // inverse lumped diagonal of weighted nodal mass matrix
       if (requiredUserData.contains("invMk_1_invBeta")) {
         // nodal mass matrix
         auto Q_rho = rh->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix AUXILIARY_NODE"));
         debugInfo("Q_rho", Q_rho);
         auto invDiagQ_rho = Teuchos::rcp_dynamic_cast<const Thyra::DiagonalLinearOpBase<Scalar> >(Teko::getInvLumpedMatrix(Q_rho), true);
         muelu_pl.sublist("user data").set("invMk_1_invBeta",invDiagQ_rho);
       }

       // Form the curl-curl matrix
       if (requiredUserData.contains("CurlCurl")) {
         auto curl_curl = Teko::explicitAdd(S_E, Teko::scale(-1.0, Q_E), Teuchos::null);
         muelu_pl.sublist("user data").set("CurlCurl",curl_curl);
       }

       if ((S_E_prec_type_ == "MueLuRefMaxwell") || (S_E_prec_type_ == "MueLuMaxwell1")) {
         Teko::InverseLibrary myInvLib = invLib;
         muelu_pl.set("Type",S_E_prec_type_);
         myInvLib.addInverse("S_E Preconditioner",muelu_pl);
         S_E_prec_factory = myInvLib.getInverseFactory("S_E Preconditioner");
       }

       if (S_E_prec_factory.is_null())
         S_E_prec_factory = invLib.getInverseFactory("S_E Preconditioner");
     }

     if (useAsPreconditioner) {
       invS_E = Teko::buildInverse(*S_E_prec_factory,S_E);
     } else {
       if (S_E_prec_.is_null()) {
         if (S_E_prec_type_ != "none") {
           S_E_prec_ = Teko::buildInverse(*S_E_prec_factory,S_E);
         }
       } else
         Teko::rebuildInverse(*S_E_prec_factory,S_E, S_E_prec_);
       if (!S_E_prec_.is_null())
         invS_E = Teko::buildInverse(*invLib.getInverseFactory("S_E Solve"),S_E,S_E_prec_);
       else
         invS_E = Teko::buildInverse(*invLib.getInverseFactory("S_E Solve"),S_E);
     }
   }


   /////////////////////////////////////////////////
   // Build block inverse matrices                //
   /////////////////////////////////////////////////

   if (!use_discrete_curl_) {
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
   } else {
     Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("MaxwellPreconditioner: Block preconditioner"));

     Teko::LinearOp id_B = Teko::identity(Teko::rangeSpace(Q_B));

     // Inverse blocks
     std::vector<Teko::LinearOp> diag(2);
     diag[0] = Teko::scale(dt,id_B);
     diag[1] = invS_E;

     // Upper tri blocks
     Teko::BlockedLinearOp U = Teko::createBlockedOp();
     Teko::beginBlockFill(U,rows,rows);
     Teko::setBlock(0,0,U,Teko::scale(1/dt,id_B));
     Teko::setBlock(1,1,U,S_E);
     Teko::setBlock(0,1,U,C);
     Teko::endBlockFill(U);

     Teko::LinearOp invU = Teko::createBlockUpperTriInverseOp(U,diag);

     if (!useAsPreconditioner) {
       Teko::BlockedLinearOp invL = Teko::createBlockedOp();
       Teko::LinearOp id_E = Teko::identity(Teko::rangeSpace(Q_E));
       Teko::beginBlockFill(invL,rows,rows);
       Teko::setBlock(0,0,invL,id_B);
       Teko::setBlock(1,0,invL,Thyra::scale(-dt, Kt));
       Teko::setBlock(1,1,invL,id_E);
       Teko::endBlockFill(invL);

       if (!simplifyFaraday_) {
         Teko::BlockedLinearOp invDiag = Teko::createBlockedOp();
         Teko::beginBlockFill(invDiag,rows,rows);
         Teko::setBlock(0,0,invDiag,Teko::scale(1/dt,invQ_B));
         Teko::setBlock(1,1,invDiag,id_E);
         Teko::endBlockFill(invDiag);

         return Teko::multiply(invU, Teko::multiply(Teko::toLinearOp(invL), Teko::toLinearOp(invDiag)));
       } else
         return Teko::multiply(invU, Teko::toLinearOp(invL));
     } else
       return invU;
   }
}

//! Initialize from a parameter list
void FullMaxwellPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList & pl)
{
  /////////////////////
  // Solver options  //
  // //////////////////

  params = pl;

  use_discrete_curl_     = params.get("Use discrete curl",false);
  dump                   = params.get("Dump",false);
  doDebug                = params.get("Debug",false);
  useAsPreconditioner    = params.get("Use as preconditioner",false);
  simplifyFaraday_       = params.get("Simplify Faraday",false) && use_discrete_curl_;
  dt                     = params.get<double>("dt");

  if(pl.isSublist("S_E Preconditioner") && pl.sublist("S_E Preconditioner").isParameter("Type"))
    S_E_prec_type_ = pl.sublist("S_E Preconditioner").get<std::string>("Type");
  else
    S_E_prec_type_ = "";

  //////////////////////////////////
  // Set up sub-solve factories   //
  //////////////////////////////////

  // New inverse lib to add inverse factories to
  invLib = *getInverseLibrary();

  // Q_B solve
  if (pl.isParameter("Q_B Solve")) {
    Teuchos::ParameterList Q_B_solve_pl = pl.sublist("Q_B Solve");
    invLib.addInverse("Q_B Solve",Q_B_solve_pl);
  }

  // Q_B preconditioner
  Teuchos::ParameterList Q_B_prec_pl = pl.sublist("Q_B Preconditioner");
  invLib.addStratPrecond("Q_B Preconditioner",
                         Q_B_prec_pl.get<std::string>("Prec Type"),
                         Q_B_prec_pl.sublist("Prec Types"));

  // S_E solve
  if (pl.isSublist("S_E Solve")) {
    Teuchos::ParameterList S_E_solve_pl = pl.sublist("S_E Solve");
    invLib.addInverse("S_E Solve",S_E_solve_pl);
  }

  // S_E preconditioner
  Teuchos::ParameterList S_E_prec_pl = pl.sublist("S_E Preconditioner");
  if (S_E_prec_type_ == "MueLuRefMaxwell" || S_E_prec_type_ == "MueLuMaxwell1") {
    invLib.addInverse("S_E Preconditioner",S_E_prec_pl);
  } else {
    if (S_E_prec_pl.isType<std::string>("Prec Type") && S_E_prec_pl.isSublist("Prec Types"))
      invLib.addStratPrecond("S_E Preconditioner",
                             S_E_prec_pl.get<std::string>("Prec Type"),
                             S_E_prec_pl.sublist("Prec Types"));
    else
      S_E_prec_type_ = "none";
  }
}

void FullMaxwellPreconditionerFactory::debugInfo(const std::string &name, const Teko::LinearOp &mat) const {
  if (dump)
    writeOut(name+".mm",*mat);
  if (doDebug) {
     auto debug = Teko::getOutputStream();
     describeMatrix(name,*mat,debug);
  }
}

} // namespace mini_em
