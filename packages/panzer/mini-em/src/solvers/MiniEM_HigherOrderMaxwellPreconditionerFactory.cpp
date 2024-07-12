// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MiniEM_HigherOrderMaxwellPreconditionerFactory.hpp"

#include "Teko_BlockLowerTriInverseOp.hpp"
#include "Teko_BlockUpperTriInverseOp.hpp"

#include "Teko_SolveInverseFactory.hpp"

#include "Teuchos_Time.hpp"

#include "Teko_TpetraHelpers.hpp"

#include "Thyra_TpetraLinearOp.hpp"
#include "Panzer_NodeType.hpp"
#include "Panzer_String_Utilities.hpp"

#include "MiniEM_Utils.hpp"
#include <stdexcept>

using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

namespace mini_em {


///////////////////////////////////////
// HigherOrderMaxwellPreconditionerFactory  //
///////////////////////////////////////

Teko::LinearOp HigherOrderMaxwellPreconditionerFactory::buildPreconditionerOperator(Teko::BlockedLinearOp & blo, Teko::BlockPreconditionerState & /* state */) const
{
   using Scalar        = double;
   using LocalOrdinal  = int;
   using GlobalOrdinal = panzer::GlobalOrdinal;
   using Node          = panzer::TpetraNodeType;
   using TpMV          = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;

   Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("MaxwellPreconditioner::build")));

   // Output stream for debug information
   RCP<Teuchos::FancyOStream> debug = Teuchos::null;
   if (doDebug)
     debug = Teko::getOutputStream();

   // Check that system is right size
   int rows = Teko::blockRowCount(blo);
   int cols = Teko::blockColCount(blo);
   TEUCHOS_ASSERT(rows==cols);
   TEUCHOS_ASSERT(rows==2);
   // Interpolation between FE spaces not implemented for Epetra
   TEUCHOS_ASSERT(Teko::TpetraHelpers::isTpetraLinearOp(Teko::getBlock(0,0,blo)));

   // Notation:
   // 0 - Hgrad
   // 1 - Hcurl
   // 2 - HDiv

   // M_k(a) - mass matrix on space k with weight a
   // I_k - identity matrix on space k
   // D_k - derivative from space k to k+1

   // The input block matrix is
   //
   // | Q_B  K   |
   // | Kt   Q_E |
   //
   // where
   //  Q_B = 1/dt * M_2(1)
   //  K   = M_2(1) * D_1
   //  Kt  = - D_1^T * M_2(1/mu)
   //  Q_E = 1/dt * M_1(eps)

   // We form a block factorization resulting in the Schur complement:
   //
   // | Q_B  K   |  = | dt*Q_B      | * | I_2        | * | 1/dt*I_2  D_1 |
   // | Kt   Q_E |    |         I_1 |   | dt*Kt  I_1 |   |           S_E |
   //
   // S_E = Q_E - Kt * Q_B^-1 * K
   //     = 1/dt * M_1(eps) + dt * D_1^T * M_2(1/mu) * D_1
   //
   // -> curl-curl term scales like dt / mu
   //
   // We assemble the curl-curl term to preserve its sparsity.
   // We also note that we need the strong form of the curl operator, D_1 = hoC.

   // The inverse is given by
   //
   // | 1/dt*I_2  D_1 |^-1 * | I_2         | * |1/dt*Q_B^-1     |
   // |           S_E |      | -dt*Kt  I_1 |   |             I_1|
   //
   //
   // The inversion of S_E is the only difficult one, the rest are mass and identity solves.

   // For "Simplify Faraday"" = true, we apply the inverse of the diagonal
   // matrix to the linear system directly and obtain the modified system:
   //
   // | Q_B  K   |
   // | Kt   Q_E |
   //
   // where
   //  Q_B = 1/dt * I_2
   //  K   = D_1
   //  Kt  = - D_1^T * M_2(1/mu)
   //  Q_E = 1/dt * M_1(eps)

   // For first order finite elements, we apply RefMaxwell directly to S_E. (See FullMaxwellPreconditionerFactory.)
   // For higher order elements, we p-coarsen down to first order elements, and then apply RefMaxwell.

   // We need to assemble:
   // - interpolation operators between high-order and low-order Hgrad and Hcurl spaces (interpHgrad and interpHcurl),
   // - strong gradients in high- and low-order spaces (hoT and loT)
   // - Hcurl mass matrices with weights 1 and dt/mu in the low-order space (M1 and Ms)
   // - inverse lumped diagonal of Hgrad mass matrix with weight mu / dt in the low-order space (M0inv)
   //
   // We directly go from high-order to low-order space, but we could allow taking several p-coarsening steps.

   // for refmaxwell: Q_rho = M_0(mu / dt) so that the addon is:
   // M_1(1) * D_0 * M_0(mu / dt)^-1 * D_0^T * M_1(1)


   // Modify the system
   if (simplifyFaraday_) {
     RCP<Teuchos::FancyOStream> out = Teko::getOutputStream();
     *out << std::endl;
     *out << "*** WARNING ***" << std::endl;
     *out << "We are modifying the linear system. That's not a friendly thing to do." << std::endl;
     *out << std::endl;

     Teko::LinearOp Q_B  = Teko::getBlock(0, 0, blo);
     Teko::LinearOp id_B = getIdentityMatrix(Q_B, 1/dt);
     Teko::LinearOp hoC  = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Curl"));
     Teko::LinearOp Kt   = Teko::getBlock(1, 0, blo);
     Teko::LinearOp Q_E  = Teko::getBlock(1, 1, blo);
     blo->beginBlockFill(2,2);
     Teko::setBlock(0, 0, blo, id_B);
     Teko::setBlock(0, 1, blo, hoC);
     Teko::setBlock(1, 0, blo, Kt);
     Teko::setBlock(1, 1, blo, Q_E);
     blo->endBlockFill();
   }

   // Extract the blocks
   Teko::LinearOp Q_B = Teko::getBlock(0, 0, blo);
   Teko::LinearOp K   = Teko::getBlock(0, 1, blo);
   Teko::LinearOp Kt  = Teko::getBlock(1, 0, blo);
   Teko::LinearOp Q_E = Teko::getBlock(1, 1, blo);

   // discrete curl between high-order spaces
   Teko::LinearOp hoC = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Curl"));

   describeAndWriteMatrix("Q_B", *Q_B, debug, dump);
   describeAndWriteMatrix("K",   *K,   debug, dump);
   describeAndWriteMatrix("Kt",  *Kt,  debug, dump);
   describeAndWriteMatrix("Q_E", *Q_E, debug, dump);
   describeAndWriteMatrix("hoC", *hoC, debug, dump);

   // interpolations between HCurl spaces of different orders
   std::vector<Teko::LinearOp> interpolationsHCurl;
   // discrete gradients HGrad -> HCurl
   std::vector<Teko::LinearOp> discreteGradients;
   // Schur complements
   std::vector<Teko::LinearOp> schurComplements;
   // projection of Schur complements into HGrad using discrete gradients
   std::vector<Teko::LinearOp> projectedSchurComplements;
   {
     using assembleType = std::pair<std::string, std::vector<Teko::LinearOp>& >;

     std::vector<assembleType> assemble = {{"Discrete Gradient",                       discreteGradients},
                                           {"SchurComplement AUXILIARY_EDGE",          schurComplements},
                                           {"ProjectedSchurComplement AUXILIARY_NODE", projectedSchurComplements}};

     for (auto ait = assemble.begin(); ait != assemble.end(); ++ait) {
       std::string operatorName = ait->first;
       std::vector<Teko::LinearOp>& operatorVector = ait->second;
       operatorVector.push_back(getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg(operatorName)));
       describeAndWriteMatrix(operatorName, *operatorVector.back(), debug, dump);
     }

     auto it = pCoarsenSchedule_.begin();
     int p = *it;
     ++it;
     while (it != pCoarsenSchedule_.end()) {
       int q = *it;

       assemble = {{"Discrete Gradient "+std::to_string(q),                              discreteGradients},
                   {"SchurComplement AUXILIARY_EDGE_"+std::to_string(q),                 schurComplements},
                   {"ProjectedSchurComplement AUXILIARY_NODE_"+std::to_string(q),        projectedSchurComplements},
                   {"Interpolation Hcurl "+std::to_string(q) + "->" + std::to_string(p), interpolationsHCurl}};

       for (auto ait = assemble.begin(); ait != assemble.end(); ++ait) {
         std::string operatorName = ait->first;
         std::vector<Teko::LinearOp>& operatorVector = ait->second;
         operatorVector.push_back(getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg(operatorName)));
         describeAndWriteMatrix(operatorName, *operatorVector.back(), debug, dump);
       }

       p = q;
       ++it;
     }
   }

   // Schur complement
   Teko::LinearOp S_E = schurComplements.front();
   // Hgrad mass matrix, mu / dt weight
   Teko::LinearOp M0 = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix AUXILIARY_NODE_1"));
   // Get inverse lumped diagonal of M0 -> M0inv
   Teko::LinearOp M0inv = getLumpedInverseDiagonal(M0);
   // Hcurl mass matrix, unit weight
   Teko::LinearOp M1 = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix AUXILIARY_EDGE_1"));
   // Hcurl mass matrix, dt/mu weight
   Teko::LinearOp Ms;
   try {
     Ms = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix weighted AUXILIARY_EDGE 1"));
   } catch (std::runtime_error&) {
     Ms = M1;
   }

   describeAndWriteMatrix("S_E",*S_E,debug,dump);
   describeAndWriteMatrix("M0",*M0,debug,dump);
   describeAndWriteMatrix("M1",*M1,debug,dump);
   describeAndWriteMatrix("Ms",*Ms,debug,dump);

   /////////////////////////////////////////////////
   // Set up inverses for sub-blocks              //
   /////////////////////////////////////////////////

   // Inverse of B mass matrix
   Teko::LinearOp invQ_B;
   if (!simplifyFaraday_) {
     Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("MaxwellPreconditioner: Inverse Q_B"));
     // Are we building a solver or a preconditioner?
     if (useAsPreconditioner) {
       invQ_B = invLib.getInverseFactory("Q_B Preconditioner")->buildInverse(Q_B);
     } else {
       Teko::LinearOp invDiagQ_B = invLib.getInverseFactory("Q_B Preconditioner")->buildInverse(Q_B);
       describeAndWriteMatrix("invDiagQ_B",*invDiagQ_B,debug,dump);
       invQ_B = invLib.getInverseFactory("Q_B Solve")->buildInverse(Q_B, invDiagQ_B);
     }
   }

   // Solver for S_E
   Teko::LinearOp invS_E;
   {
     Teuchos::TimeMonitor tm1(*Teuchos::TimeMonitor::getNewTimer("MaxwellPreconditioner: Inverse S_E"));

     RCP<Teko::InverseFactory> S_E_prec_factory = invLib.getInverseFactory("S_E Preconditioner");
     Teuchos::ParameterList  S_E_prec_pl  = *S_E_prec_factory->getParameterList();
     Teuchos::ParameterList& muelulist = S_E_prec_pl.sublist("Preconditioner Types").sublist(S_E_prec_type_);

     int maxLevels = interpolationsHCurl.size()+1;

     // Make sure MueLu only creates levels for the operators that we pass in
     muelulist.set("max levels", maxLevels);


     bool implicitTranspose = muelulist.get("transpose: use implicit", false);

     // Cannot explicitly transpose matrix-free operators
     if (!implicitTranspose && isMatrixFreeOperator(interpolationsHCurl[0])) {
       implicitTranspose = true;
       muelulist.set("transpose: use implicit", true);
     }
     std::string smootherType = muelulist.get("smoother: type", "HIPTMAIR");
     std::transform(smootherType.begin(), smootherType.end(), smootherType.begin(), ::toupper);
     if (smootherType == "HIPTMAIR") {
       if (isMatrixFreeOperator(interpolationsHCurl[0]))
         muelulist.sublist("smoother: params").set("hiptmair: implicit transpose", true);
       TEUCHOS_ASSERT_EQUALITY(muelulist.get<bool>("transpose: use implicit"),
                               muelulist.sublist("smoother: params").get<bool>("hiptmair: implicit transpose"));
     }

     std::vector<Teko::LinearOp> interpolationsHCurlT;
     if (!implicitTranspose)
       // Get restriction operators
       for (int lvl = 1; lvl < maxLevels; ++lvl) {
         Teko::LinearOp interpT = Teko::explicitTranspose(interpolationsHCurl[lvl-1]);
         interpolationsHCurlT.push_back(interpT);
       }

     for (int lvl = 0; lvl < maxLevels; ++lvl) {
       Teuchos::ParameterList& lvlList = muelulist.sublist("level " + std::to_string(lvl) + " user data");
       if (lvl > 0) {
         lvlList.set("A",schurComplements[lvl]);
         lvlList.set("P",interpolationsHCurl[lvl-1]);
         if (!implicitTranspose)
           lvlList.set("R",interpolationsHCurlT[lvl-1]);
       }
       // Operators for Hiptmair smoothing
       lvlList.set("NodeMatrix",projectedSchurComplements[lvl]);
       lvlList.set("D0",discreteGradients[lvl]);
     }

     // Operators for RefMaxwell coarse grid solve
     // ("A" and "D0" are already set above.)
     Teuchos::ParameterList& lvlList = muelulist.sublist("level " + std::to_string(maxLevels-1) + " user data");
     lvlList.set("Dk_1",discreteGradients.back());
     lvlList.set("Mk_one",M1);
     lvlList.set("M1_beta",Ms);
     lvlList.set("Coordinates", S_E_prec_pl.get<RCP<TpMV> >("Coordinates"));
     lvlList.set("invMk_1_invBeta",M0inv);

     // make sure we build all levels
     muelulist.set("coarse: max size", Teuchos::as<int>(S_E_prec_pl.get<RCP<TpMV> >("Coordinates")->getMap()->getGlobalNumElements()-1));

     Teko::InverseLibrary myInvLib = invLib;

     muelulist.set("Type",S_E_prec_type_);
     myInvLib.addInverse("S_E Preconditioner", muelulist);
     S_E_prec_factory = myInvLib.getInverseFactory("S_E Preconditioner");

     // Are we building a solver or a preconditioner?
     if (useAsPreconditioner)
       invS_E = S_E_prec_factory->buildInverse(S_E);
     else {
       if (S_E_prec_.is_null())
         S_E_prec_ = S_E_prec_factory->buildInverse(S_E);
       else
         S_E_prec_factory->rebuildInverse(S_E, S_E_prec_);
       invS_E = invLib.getInverseFactory("S_E Solve")->buildInverse(S_E,S_E_prec_);
     }
   }


   /////////////////////////////////////////////////
   // Build block inverse matrices                //
   /////////////////////////////////////////////////

   {
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
     Teko::setBlock(0,1,U,hoC);
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
void HigherOrderMaxwellPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList & pl)
{
   /////////////////////
   // Solver options  //
   // //////////////////

   params = pl;

   dump                   = params.get("Dump",false);
   doDebug                = params.get("Debug",false);
   useAsPreconditioner    = params.get("Use as preconditioner",false);
   simplifyFaraday_       = params.get("Simplify Faraday",true);
   dt                     = params.get<double>("dt");

   std::string pCoarsenScheduleStr = params.get<std::string>("p coarsen schedule");
   std::vector<std::string> pCoarsenScheduleVecStr;
   panzer::StringTokenizer(pCoarsenScheduleVecStr, pCoarsenScheduleStr, ",");
   panzer::TokensToInts(pCoarsenSchedule_, pCoarsenScheduleVecStr);

   S_E_prec_type_ = pl.sublist("S_E Preconditioner").get<std::string>("Type");
   TEUCHOS_ASSERT(S_E_prec_type_ == "MueLu");

   // Output stream for debug information
   RCP<Teuchos::FancyOStream> debug = Teuchos::null;
   if (doDebug)
     debug = Teko::getOutputStream();

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

   // S_E solve
   if (pl.isParameter("S_E Solve")) {
     Teuchos::ParameterList ml_pl = pl.sublist("S_E Solve");
     invLib.addInverse("S_E Solve",ml_pl);
   }

   // S_E preconditioner
   Teuchos::ParameterList S_E_prec_pl = pl.sublist("S_E Preconditioner");
   invLib.addInverse("S_E Preconditioner",S_E_prec_pl);
}

} // namespace mini_em
