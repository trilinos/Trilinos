#include "MiniEM_HigherOrderMaxwellPreconditionerFactory.hpp"

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
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_LinearObjContainer.hpp"
#include "Panzer_ThyraObjContainer.hpp"
#include "Panzer_String_Utilities.hpp"

#include "Thyra_DefaultDiagonalLinearOp.hpp"

#include "MiniEM_Utils.hpp"
#include "Xpetra_ThyraUtils.hpp"
#include "Xpetra_MatrixFactory.hpp"

using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

namespace mini_em {


  // template<class SC, class LO, class GO,class NT>
  // RCP<Tpetra::CrsMatrix<SC, LO, GO, NT> >
  // getIdentityMatrixTpetra (Teuchos::RCP<const Tpetra::Map<LO,GO,NT> >& identityRowMap)
  // {
  //   using Teuchos::RCP;
  //   typedef Tpetra::CrsMatrix<SC, LO, GO, NT> Matrix_t;

  //   RCP<Matrix_t> identityMatrix = Teuchos::rcp(new Matrix_t(identityRowMap, 1));
  //   Teuchos::ArrayView<const GO> gblRows = identityRowMap->getLocalElementList ();
  //   Teuchos::Array<SC> val (1, Teuchos::ScalarTraits<SC>::one ());
  //   for (auto it = gblRows.begin (); it != gblRows.end (); ++it) {
  //     Teuchos::Array<GO> col (1, *it);
  //     identityMatrix->insertGlobalValues (*it, col (), val ());
  //   }
  //   identityMatrix->fillComplete ();
  //   return identityMatrix;
  // }


///////////////////////////////////////
// HigherOrderMaxwellPreconditionerFactory  //
///////////////////////////////////////

Teko::LinearOp HigherOrderMaxwellPreconditionerFactory::buildPreconditionerOperator(Teko::BlockedLinearOp & blo, Teko::BlockPreconditionerState & /* state */) const
{
   typedef double Scalar;
   typedef int LocalOrdinal;
   typedef panzer::GlobalOrdinal GlobalOrdinal;
   typedef panzer::TpetraNodeType Node;

   using XpThyUtils   = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
   using ThyLinOpBase = Thyra::LinearOpBase<Scalar>;
   using ThyDiagLinOpBase = Thyra::DiagonalLinearOpBase<Scalar>;
   using XpMat = Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
   using TpMV  = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;
   using XpMV  = Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>;

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

   // Check whether we are using Tpetra or Epetra
   bool useTpetra;
   {
     Teko::LinearOp blo00 = Teko::getBlock(0,0,blo);
     RCP<const Thyra::EpetraLinearOp> EOp = rcp_dynamic_cast<const Thyra::EpetraLinearOp>(blo00);
     useTpetra = (EOp == Teuchos::null);
   }
   TEUCHOS_ASSERT(useTpetra);  // Interpolation between FE spaces not implemented for Epetra

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

     Teko::LinearOp Q_B  = Teko::getBlock(0,0,blo);
     Teko::LinearOp id_B = getIdentityMatrix(Q_B, 1/dt);
     Teko::LinearOp hoC  = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Curl"));
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

   // discrete curl between high-order spaces
   Teko::LinearOp hoC = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Curl"));

   // Schur complement
   Teko::LinearOp S_E;
   {
     Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("MaxwellPreconditioner: Schur complement"));
     Teko::LinearOp CurlCurl = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Curl Curl AUXILIARY_EDGE"));
     S_E = Teko::explicitAdd(Q_E, CurlCurl);
   }

   // interpolations between FE spaces of different orders
   std::vector<Teko::LinearOp> interpolationsHGrad;
   std::vector<Teko::LinearOp> interpolationsHCurl;
   // discrete gradients
   std::vector<Teko::LinearOp> discreteGradients;
   {
     std::string discGradName, interpHcurlName, interpHgradName;
     auto it = pCoarsenSchedule_.begin();
     int p = *it;
     discGradName = "Discrete Gradient";
     discreteGradients.push_back(getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg(discGradName)));
     describeMatrix(discGradName, *discreteGradients.back(), debug);
     if (dump) {
       writeOut(discGradName+".mm", *discreteGradients.back());
     }
     ++it;
     while (it != pCoarsenSchedule_.end()) {
       int q = *it;
       discGradName = "Discrete Gradient "+std::to_string(*it);
       interpHcurlName = "Interpolation Hcurl "+std::to_string(q) + "->" + std::to_string(p);
       interpHgradName = "Interpolation Hgrad "+std::to_string(q) + "->" + std::to_string(p);
       discreteGradients.push_back(getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg(discGradName)));
       interpolationsHGrad.push_back(getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg(interpHgradName)));
       interpolationsHCurl.push_back(getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg(interpHcurlName)));

       describeMatrix(discGradName, *discreteGradients.back(), debug);
       describeMatrix(interpHcurlName, *interpolationsHCurl.back(), debug);
       describeMatrix(interpHgradName, *interpolationsHGrad.back(), debug);
       if (dump) {
         writeOut(discGradName+".mm", *discreteGradients.back());
         writeOut(interpHcurlName+".mm", *interpolationsHCurl.back());
         writeOut(interpHgradName+".mm", *interpolationsHGrad.back());
       }
       p = q;
       ++it;
     }
   }

   // Hgrad mass matrix, mu / dt weight
   Teko::LinearOp M0 = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix AUXILIARY_NODE_1"));
   // Hcurl mass matrix, unit weight
   Teko::LinearOp M1 = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix AUXILIARY_EDGE_1"));

   // Get inverse lumped diagonal of M0 -> M0inv
   RCP<const ThyDiagLinOpBase> M0inv;
   {
     // Get inverse of lumped M0
     RCP<Thyra::VectorBase<Scalar> > ones = Thyra::createMember(M0->domain());
     RCP<Thyra::VectorBase<Scalar> > diagonal = Thyra::createMember(M0->range());
     Thyra::assign(ones.ptr(),1.0);
     // compute lumped diagonal
     Thyra::apply(*M0,Thyra::NOTRANS,*ones,diagonal.ptr());
     Thyra::reciprocal(*diagonal,diagonal.ptr());
     M0inv = rcp(new Thyra::DefaultDiagonalLinearOp<Scalar>(diagonal));
   }

   /////////////////////////////////////////////////
   // Debug and matrix dumps                      //
   /////////////////////////////////////////////////

   describeMatrix("Q_B",*Q_B,debug);
   describeMatrix("K",*K,debug);
   describeMatrix("Kt",*Kt,debug);
   describeMatrix("Q_E",*Q_E,debug);

   describeMatrix("hoC",*hoC,debug);
   describeMatrix("S_E",*S_E,debug);

   describeMatrix("M0",*M0,debug);
   describeMatrix("M1",*M1,debug);

   if (dump) {
     writeOut("Q_B.mm",*Q_B);
     writeOut("K.mm",*K);
     writeOut("Kt.mm",*Kt);
     writeOut("Q_E.mm",*Q_E);

     writeOut("hoC.mm",*hoC);
     writeOut("S_E.mm",*S_E);

     writeOut("M0.mm",*M0);
     writeOut("M1.mm",*M1);

     Teko::LinearOp K2 = Teko::explicitMultiply(Q_B,hoC);
     Teko::LinearOp diffK;

     {
       RCP<const Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tOp = rcp_dynamic_cast<const Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(K2,true);
       RCP<Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tOp2 = Teuchos::rcp_const_cast<Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>>(tOp);
       RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsOp = rcp_dynamic_cast<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(tOp2->getTpetraOperator(),true);
       crsOp->scale(dt);
       diffK = Teko::explicitAdd(K, Teko::scale(-1.0,K2));

       writeOut("K2.mm",*K2);
       writeOut("diff.mm",*diffK);

     }

     TEUCHOS_ASSERT(Teko::infNorm(diffK) < 1.0e-8 * Teko::infNorm(K));
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
       describeMatrix("invDiagQ_B",*invDiagQ_B,debug);
       invQ_B = Teko::buildInverse(*invLib.getInverseFactory("Q_B Solve"),Q_B, invDiagQ_B);
     }
   }

   // Solver for S_E
   Teko::LinearOp invS_E;
   {
     Teuchos::TimeMonitor tm1(*Teuchos::TimeMonitor::getNewTimer("MaxwellPreconditioner: Inverse S_E"));

     RCP<Teko::InverseFactory> S_E_prec_factory = invLib.getInverseFactory("S_E Preconditioner");
     Teuchos::ParameterList  S_E_prec_pl  = *S_E_prec_factory->getParameterList();
     Teuchos::ParameterList& maxwell1list = S_E_prec_pl.sublist("Preconditioner Types").sublist(S_E_prec_type_);
     Teuchos::ParameterList& list11       = maxwell1list.sublist("maxwell1: 11list");
     Teuchos::ParameterList& list22       = maxwell1list.sublist("maxwell1: 22list");

     // Maxwell1 list
     maxwell1list.set("D0", discreteGradients.front());

     int maxLevels = interpolationsHGrad.size()+1;

     list22.set("max levels", maxLevels);

     for (int lvl = 1; lvl < maxLevels; ++lvl) {
       // (2,2) list
       auto xpInterpHgrad = XpThyUtils::toXpetra(Teuchos::rcp_const_cast<ThyLinOpBase>(interpolationsHGrad[lvl-1]));
       list22.sublist("level " + std::to_string(lvl) + " user data").set("P",xpInterpHgrad);

       // (1,1) list
       auto xpInterpHcurl = XpThyUtils::toXpetra(Teuchos::rcp_const_cast<ThyLinOpBase>(interpolationsHCurl[lvl-1]));
       list11.sublist("level " + std::to_string(lvl) + " user data").set("P",xpInterpHcurl);

       auto xpLoT = XpThyUtils::toXpetra(Teuchos::rcp_const_cast<ThyLinOpBase>(discreteGradients[lvl]));
       list11.sublist("level " + std::to_string(lvl) + " user data").set("D0",xpLoT);
     }


     auto xpM1 = XpThyUtils::toXpetra(Teuchos::rcp_const_cast<ThyLinOpBase>(M1));
     list11.sublist("level " + std::to_string(maxLevels-1) + " user data").set("M1",xpM1);
     list11.sublist("level " + std::to_string(maxLevels-1) + " user data").set("Ms",xpM1);

     RCP<XpMV> CoordinatesLO = Xpetra::toXpetra(S_E_prec_pl.get<RCP<TpMV> >("Coordinates"));
     list11.sublist("level " + std::to_string(maxLevels-1) + " user data").set("Coordinates", CoordinatesLO);

     RCP<XpMat> xpM0inv = XpThyUtils::toXpetra(Teuchos::rcp_const_cast<ThyDiagLinOpBase>(M0inv));;
     list11.sublist("level " + std::to_string(maxLevels-1) + " user data").set("M0inv",xpM0inv);

     Teko::InverseLibrary myInvLib = invLib;

     maxwell1list.set("Type",S_E_prec_type_);
     myInvLib.addInverse("S_E Preconditioner", maxwell1list);
     S_E_prec_factory = myInvLib.getInverseFactory("S_E Preconditioner");

     // Are we building a solver or a preconditioner?
     if (useAsPreconditioner)
       invS_E = Teko::buildInverse(*S_E_prec_factory,S_E);
     else {
       if (S_E_prec_.is_null())
         S_E_prec_ = Teko::buildInverse(*S_E_prec_factory,S_E);
       else
         Teko::rebuildInverse(*S_E_prec_factory,S_E, S_E_prec_);
       invS_E = Teko::buildInverse(*invLib.getInverseFactory("S_E Solve"),S_E,S_E_prec_);
     }
   }


   /////////////////////////////////////////////////
   // Build block  inverse matrices               //
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

   std::string pCoarsenScheduleStr = params.get<std::string>("p coarsen schedule");
   std::vector<std::string> pCoarsenScheduleVecStr;
   panzer::StringTokenizer(pCoarsenScheduleVecStr, pCoarsenScheduleStr, ",");
   panzer::TokensToInts(pCoarsenSchedule_, pCoarsenScheduleVecStr);

   S_E_prec_type_ = pl.sublist("S_E Preconditioner").get<std::string>("Type");
   TEUCHOS_ASSERT(S_E_prec_type_ == "MueLuMaxwell1");

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

   dt = params.get<double>("dt");

   // S_E solve
   Teuchos::ParameterList ml_pl = pl.sublist("S_E Solve");
   invLib.addInverse("S_E Solve",ml_pl);

   // S_E preconditioner
   Teuchos::ParameterList S_E_prec_pl = pl.sublist("S_E Preconditioner");

   invLib.addInverse("S_E Preconditioner",S_E_prec_pl);
}

} // namespace mini_em
