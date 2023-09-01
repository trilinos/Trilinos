#include "MiniEM_FullDarcyPreconditionerFactory.hpp"

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
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Thyra_EpetraThyraWrappers.hpp"
#endif
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_LinearObjContainer.hpp"
#include "Panzer_ThyraObjContainer.hpp"

#include "Thyra_DefaultDiagonalLinearOp.hpp"

#include "MiniEM_Utils.hpp"

using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

namespace mini_em {


///////////////////////////////////////
// FullDarcyPreconditionerFactory  //
///////////////////////////////////////

Teko::LinearOp FullDarcyPreconditionerFactory::buildPreconditionerOperator(Teko::BlockedLinearOp & blo, Teko::BlockPreconditionerState & /* state */) const
{
   typedef double Scalar;
   typedef int LocalOrdinal;
   typedef panzer::GlobalOrdinal GlobalOrdinal;
   typedef panzer::TpetraNodeType Node;

   Teuchos::TimeMonitor tM(*Teuchos::TimeMonitor::getNewTimer(std::string("DarcyPreconditioner::build")));

   // Output stream for debug information
   RCP<Teuchos::FancyOStream> debug = Teuchos::null;
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
   // 3 - HVol

   // M_k(a) - mass matrix on space k with weight a
   // D_k - derivative from space k to k+1

   // The block matrix is
   //
   // | Q_u  K   |
   // | Kt   Q_sigma |
   //
   // where
   // Q_u = 1/dt * M_3(1)
   // K   = -M_3(1) * D_2
   // Kt  = D_2^T * M_3(1)
   // Q_sigma = M_2(1/kappa)

   // S_sigma = Q_sigma - Kt * Q_u^-1 * K
   //     = M_2(1/kappa) + dt * D_2^T * M_3(1) * D_2
   //
   // -> grad-div term scales like dt
   //
   // for refmaxwell: Q_rho = M_1(1/dt) so that the addon is:
   // M_2(1) * D_1 * M_1(1/dt)^-1 * D_1^T * M_2(1)

   // Modify the system
   if (simplifyFaraday_) {
     RCP<Teuchos::FancyOStream> out = Teko::getOutputStream();
     *out << std::endl;
     *out << "*** WARNING ***" << std::endl;
     *out << "We are modifying the linear system. That's not a friendly thing to do." << std::endl;
     *out << std::endl;

     Teko::LinearOp Q_u  = Teko::getBlock(0,0,blo);
     Teko::LinearOp id_u = getIdentityMatrix(Q_u, 1/dt);
     Teko::LinearOp hoDiv  = Teko::scale(-1.0, getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Div")));
     Teko::LinearOp Kt    = Teko::getBlock(1,0,blo);
     Teko::LinearOp Q_sigma   = Teko::getBlock(1,1,blo);
     blo->beginBlockFill(2,2);
     Teko::setBlock(0,0,blo,id_u);
     Teko::setBlock(0,1,blo,hoDiv);
     Teko::setBlock(1,0,blo,Kt);
     Teko::setBlock(1,1,blo,Q_sigma);
     blo->endBlockFill();
   }

   // Extract the blocks
   Teko::LinearOp Q_u   = Teko::getBlock(0,0,blo);
   Teko::LinearOp K     = Teko::getBlock(0,1,blo);
   Teko::LinearOp Kt    = Teko::getBlock(1,0,blo);
   Teko::LinearOp Q_sigma   = Teko::getBlock(1,1,blo);

   // discrete curl and its transpose
   Teko::LinearOp Div, DivT;
   if (use_discrete_div_) {
     Div = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Div"));
     DivT = Teko::explicitTranspose(Div);
   }

   // Set up the Schur complement
   Teko::LinearOp S_sigma;
   {
     Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("DarcyPreconditioner: Schur complement"));
     S_sigma = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("DarcySchurComplement AUXILIARY_FACE"));
   }

   // Check whether we are using Tpetra or Epetra
   RCP<const Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> > checkTpetra = Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>>(Q_sigma);
   bool useTpetra = nonnull(checkTpetra);

   /////////////////////////////////////////////////
   // Debug and matrix dumps                      //
   /////////////////////////////////////////////////

   if (dump) {
     writeOut("Q_u.mm",*Q_u);
     writeOut("K.mm",*K);
     writeOut("Kt.mm",*Kt);
     writeOut("Q_sigma.mm",*Q_sigma);
     writeOut("S_sigma.mm",*S_sigma);

     if (Div != Teuchos::null) {
       Teko::LinearOp K2 = Teko::explicitMultiply(Q_u, Div);
       Teko::LinearOp diffK;

       if (useTpetra) {
         typedef panzer::TpetraNodeType Node;
         typedef int LocalOrdinal;
         typedef panzer::GlobalOrdinal GlobalOrdinal;

         RCP<const Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tOp = rcp_dynamic_cast<const Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(K2,true);
         RCP<Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tOp2 = Teuchos::rcp_const_cast<Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>>(tOp);
         RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > crsOp = rcp_dynamic_cast<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(tOp2->getTpetraOperator(),true);
         crsOp->scale(dt);
         diffK = Teko::explicitAdd(K, Teko::scale(-1.0,K2));

         writeOut("DiscreteDiv.mm",*Div);
         writeOut("K2.mm",*K2);
         writeOut("diff.mm",*diffK);

       } else {
         diffK = Teko::explicitAdd(K, Teko::scale(-dt,K2));

         writeOut("DiscreteDiv.mm",*Div);
         writeOut("K2.mm",*K2);
         writeOut("diff.mm",*diffK);
       }

       TEUCHOS_ASSERT(Teko::infNorm(diffK) < 1.0e-8 * Teko::infNorm(K));
     }
   }
   describeMatrix("Q_u",*Q_u,debug);
   describeMatrix("K",*K,debug);
   describeMatrix("Kt",*Kt,debug);
   describeMatrix("Q_sigma",*Q_sigma,debug);
   if (Div != Teuchos::null)
     describeMatrix("Div",*Div,debug);
   describeMatrix("S_sigma",*S_sigma,debug);


   /////////////////////////////////////////////////
   // Set up inverses for sub-blocks              //
   /////////////////////////////////////////////////

   // Inverse of B mass matrix
   Teko::LinearOp invQ_u;
   if (!simplifyFaraday_) {
     Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("DarcyPreconditioner: Inverse Q_u"));
     // Are we building a solver or a preconditioner?
     if (useAsPreconditioner) {
       invQ_u = Teko::buildInverse(*invLib.getInverseFactory("Q_u Preconditioner"),Q_u);
     } else {
       Teko::LinearOp invDiagQ_u = Teko::buildInverse(*invLib.getInverseFactory("Q_u Preconditioner"),Q_u);
       describeMatrix("invDiagQ_u",*invDiagQ_u,debug);
       invQ_u = Teko::buildInverse(*invLib.getInverseFactory("Q_u Solve"),Q_u, invDiagQ_u);
     }
   }

   // Solver for S_sigma
   Teko::LinearOp invS_sigma;
   {
     Teuchos::TimeMonitor tm1(*Teuchos::TimeMonitor::getNewTimer("DarcyPreconditioner: Solver S_sigma"));

     if (S_sigma_prec_type_ == "MueLuRefDarcy") {// refDarcy

       // edge mass matrix
       Teko::LinearOp Q_rho = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix AUXILIARY_EDGE"));
       describeMatrix("Q_rho",*Q_rho,debug);
       if (dump)
         writeOut("Q_rho.mm",*Q_rho);

       // Teko::LinearOp T = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Gradient"));
       // Teko::LinearOp KT = Teko::explicitMultiply(K,T);
       // TEUCHOS_ASSERT(Teko::infNorm(KT) < 1.0e-14 * Teko::infNorm(T) * Teko::infNorm(K));

       RCP<Teko::InverseFactory> S_sigma_prec_factory;
       Teuchos::ParameterList S_sigma_prec_pl;
       S_sigma_prec_factory = invLib.getInverseFactory("S_sigma Preconditioner");
       S_sigma_prec_pl = *S_sigma_prec_factory->getParameterList();

       // Get coordinates
       {
         if (useTpetra) {
           RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Coordinates = S_sigma_prec_pl.get<RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > >("Coordinates");
           S_sigma_prec_pl.sublist("Preconditioner Types").sublist(S_sigma_prec_type_).set("Coordinates",Coordinates);
         }
#ifdef PANZER_HAVE_EPETRA_STACK
         else {
           RCP<Epetra_MultiVector> Coordinates = S_sigma_prec_pl.get<RCP<Epetra_MultiVector> >("Coordinates");
           S_sigma_prec_pl.sublist("Preconditioner Types").sublist(S_sigma_prec_type_).set("Coordinates",Coordinates);
         }
#else
         else
           TEUCHOS_ASSERT(false);
#endif
       }

       RCP<const Thyra::DiagonalLinearOpBase<Scalar> > invDiagQ_rho;
       {
         Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("MaxwellPreconditioner: Lumped diagonal Q_rho"));

         // Get inverse of lumped Q_rho
         RCP<Thyra::VectorBase<Scalar> > ones = Thyra::createMember(Q_rho->domain());
         RCP<Thyra::VectorBase<Scalar> > diagonal = Thyra::createMember(Q_rho->range());
         Thyra::assign(ones.ptr(),1.0);
         // compute lumped diagonal
         Thyra::apply(*Q_rho,Thyra::NOTRANS,*ones,diagonal.ptr());
         Thyra::reciprocal(*diagonal,diagonal.ptr());
         invDiagQ_rho = rcp(new Thyra::DefaultDiagonalLinearOp<Scalar>(diagonal));
       }

       {
         Teko::InverseLibrary myInvLib = invLib;
         S_sigma_prec_pl.sublist("Preconditioner Types").sublist(S_sigma_prec_type_).set("M1inv",invDiagQ_rho);
         S_sigma_prec_pl.sublist("Preconditioner Types").sublist(S_sigma_prec_type_).set("Type",S_sigma_prec_type_);
         myInvLib.addInverse("S_sigma Preconditioner",S_sigma_prec_pl.sublist("Preconditioner Types").sublist(S_sigma_prec_type_));
         S_sigma_prec_factory = myInvLib.getInverseFactory("S_sigma Preconditioner");
       }

       // Are we building a solver or a preconditioner?
       {
         Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("MaxwellPreconditioner: Build S_sigma preconditioner"));

         if (useAsPreconditioner)
           invS_sigma = Teko::buildInverse(*S_sigma_prec_factory,S_sigma);
         else {
           if (S_sigma_prec_.is_null())
             S_sigma_prec_ = Teko::buildInverse(*S_sigma_prec_factory,S_sigma);
           else
             Teko::rebuildInverse(*S_sigma_prec_factory,S_sigma, S_sigma_prec_);
           invS_sigma = Teko::buildInverse(*invLib.getInverseFactory("S_sigma Solve"),S_sigma,S_sigma_prec_);
         }
       }
     } else {
       if (useAsPreconditioner)
         invS_sigma = Teko::buildInverse(*invLib.getInverseFactory("S_sigma Preconditioner"),S_sigma);
       else {
         if (S_sigma_prec_.is_null())
           S_sigma_prec_ = Teko::buildInverse(*invLib.getInverseFactory("S_sigma Preconditioner"),S_sigma);
         else
           Teko::rebuildInverse(*invLib.getInverseFactory("S_sigma Preconditioner"),S_sigma, S_sigma_prec_);
         invS_sigma = Teko::buildInverse(*invLib.getInverseFactory("S_sigma Solve"),S_sigma,S_sigma_prec_);
       }
     }
   }


   /////////////////////////////////////////////////
   // Build block  inverse matrices               //
   /////////////////////////////////////////////////

   if (!use_discrete_div_) {
     Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("MaxwellPreconditioner: Block preconditioner"));

     // Inverse blocks
     std::vector<Teko::LinearOp> diag(2);
     diag[0] = invQ_u;
     diag[1] = invS_sigma;

     // Upper tri blocks
     Teko::BlockedLinearOp U = Teko::createBlockedOp();
     Teko::beginBlockFill(U,rows,rows);
     Teko::setBlock(0,0,U,Q_u);
     Teko::setBlock(1,1,U,S_sigma);
     Teko::setBlock(0,1,U,K);
     Teko::endBlockFill(U);

     Teko::LinearOp invU = Teko::createBlockUpperTriInverseOp(U,diag);

     if (!useAsPreconditioner) {
       Teko::BlockedLinearOp invL = Teko::createBlockedOp();
       Teko::LinearOp id_u = Teko::identity(Teko::rangeSpace(Q_u));
       Teko::LinearOp id_sigma = Teko::identity(Teko::rangeSpace(Q_sigma));
       Teko::beginBlockFill(invL,rows,rows);
       Teko::setBlock(0,0,invL,id_u);
       Teko::setBlock(1,0,invL,Teko::multiply(Thyra::scale(-1.0, Kt), invQ_u));
       Teko::setBlock(1,1,invL,id_sigma);
       Teko::endBlockFill(invL);

       return Teko::multiply(invU, Teko::toLinearOp(invL));
     } else
       return invU;
   } else {
     Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("MaxwellPreconditioner: Block preconditioner"));

     Teko::LinearOp id_u = Teko::identity(Teko::rangeSpace(Q_u));

     // Inverse blocks
     std::vector<Teko::LinearOp> diag(2);
     diag[0] = Teko::scale(dt,id_u);
     diag[1] = invS_sigma;

     // Upper tri blocks
     Teko::BlockedLinearOp U = Teko::createBlockedOp();
     Teko::beginBlockFill(U,rows,rows);
     Teko::setBlock(0,0,U,Teko::scale(1/dt,id_u));
     Teko::setBlock(1,1,U,S_sigma);
     Teko::setBlock(0,1,U,Div);
     Teko::endBlockFill(U);

     Teko::LinearOp invU = Teko::createBlockUpperTriInverseOp(U,diag);

     if (!useAsPreconditioner) {
       Teko::BlockedLinearOp invL = Teko::createBlockedOp();
       Teko::LinearOp id_sigma = Teko::identity(Teko::rangeSpace(Q_sigma));
       Teko::beginBlockFill(invL,rows,rows);
       Teko::setBlock(0,0,invL,id_u);
       Teko::setBlock(1,0,invL,Thyra::scale(-dt, Kt));
       Teko::setBlock(1,1,invL,id_sigma);
       Teko::endBlockFill(invL);

       if (!simplifyFaraday_) {
         Teko::BlockedLinearOp invDiag = Teko::createBlockedOp();
         Teko::beginBlockFill(invDiag,rows,rows);
         Teko::setBlock(0,0,invDiag,Teko::scale(1/dt,invQ_u));
         Teko::setBlock(1,1,invDiag,id_sigma);
         Teko::endBlockFill(invDiag);

         return Teko::multiply(invU, Teko::multiply(Teko::toLinearOp(invL), Teko::toLinearOp(invDiag)));
       } else
         return Teko::multiply(invU, Teko::toLinearOp(invL));
     } else
       return invU;
   }
}

//! Initialize from a parameter list
void FullDarcyPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList & pl)
{
   /////////////////////
   // Solver options  //
   // //////////////////

   params = pl;

   use_discrete_div_     = params.get("Use discrete div",false);
   dump                   = params.get("Dump",false);
   doDebug                = params.get("Debug",false);
   useAsPreconditioner    = params.get("Use as preconditioner",false);
   simplifyFaraday_       = params.get("Simplify Faraday",false) && use_discrete_div_;

   if(pl.isSublist("S_sigma Preconditioner") && pl.sublist("S_sigma Preconditioner").isParameter("Type"))
     S_sigma_prec_type_ = pl.sublist("S_sigma Preconditioner").get<std::string>("Type");
   else
     S_sigma_prec_type_ = "";

   // Output stream for debug information
   RCP<Teuchos::FancyOStream> debug = Teuchos::null;
   if (doDebug)
     debug = Teko::getOutputStream();

   //////////////////////////////////
   // Set up sub-solve factories   //
   //////////////////////////////////

   // New inverse lib to add inverse factories to
   invLib = *getInverseLibrary();

   // Q_u solve
   if (pl.isParameter("Q_u Solve")) {
     Teuchos::ParameterList cg_pl = pl.sublist("Q_u Solve");
     invLib.addInverse("Q_u Solve",cg_pl);
   }

   // Q_u preconditioner
   Teuchos::ParameterList Q_u_prec_pl = pl.sublist("Q_u Preconditioner");
   invLib.addStratPrecond("Q_u Preconditioner",
                          Q_u_prec_pl.get<std::string>("Prec Type"),
                          Q_u_prec_pl.sublist("Prec Types"));

   dt = params.get<double>("dt");

   if (S_sigma_prec_type_ == "MueLuRefDarcy" || S_sigma_prec_type_ == "ML") { // RefDarcy based solve

     // S_sigma solve
     Teuchos::ParameterList ml_pl = pl.sublist("S_sigma Solve");
     invLib.addInverse("S_sigma Solve",ml_pl);

     // S_sigma preconditioner
     Teuchos::ParameterList S_sigma_prec_pl = pl.sublist("S_sigma Preconditioner");

     // add discrete curl and face mass matrix
     Teko::LinearOp Q_sigma_aux = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix AUXILIARY_FACE"));
     Teko::LinearOp Q_sigma_aux_weighted = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix weighted AUXILIARY_FACE"));
     Teko::LinearOp Curl = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Curl"));
     if (S_sigma_prec_type_ == "ML") {
#ifdef PANZER_HAVE_EPETRA_STACK
       RCP<const Epetra_CrsMatrix> eCurl = get_Epetra_CrsMatrix(*Curl);
       RCP<const Epetra_CrsMatrix> eQ_sigma_aux = get_Epetra_CrsMatrix(*Q_sigma_aux);
       RCP<const Epetra_CrsMatrix> eQ_sigma_aux_weighted = get_Epetra_CrsMatrix(*Q_sigma_aux_weighted);
       S_sigma_prec_pl.sublist("ML Settings").set("D1",eCurl);
       S_sigma_prec_pl.sublist("ML Settings").set("M2",eQ_sigma_aux);
       S_sigma_prec_pl.sublist("ML Settings").set("Ms",eQ_sigma_aux_weighted);
#else
       TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"ERROR: MiniEM_FullDarcyPreconditionerFactory: ML is not supported in this build! Requires Epetra support be enabled!");
#endif
     } else {
       S_sigma_prec_pl.sublist("Preconditioner Types").sublist(S_sigma_prec_type_).set("D1",Curl);
       S_sigma_prec_pl.sublist("Preconditioner Types").sublist(S_sigma_prec_type_).set("M2",Q_sigma_aux);
       S_sigma_prec_pl.sublist("Preconditioner Types").sublist(S_sigma_prec_type_).set("Ms",Q_sigma_aux_weighted);
     }

     if (dump) {
       writeOut("DiscreteCurl.mm",*Curl);
     }
     describeMatrix("DiscreteCurl",*Curl,debug);

     invLib.addInverse("S_sigma Preconditioner",S_sigma_prec_pl);
   } else {
     // S_sigma solve
     if (pl.isParameter("S_sigma Solve")) {
       Teuchos::ParameterList cg_pl = pl.sublist("S_sigma Solve");
       invLib.addInverse("S_sigma Solve",cg_pl);
     }

     // S_sigma preconditioner
     Teuchos::ParameterList S_sigma_prec_pl = pl.sublist("S_sigma Preconditioner");
     invLib.addStratPrecond("S_sigma Preconditioner",
                            S_sigma_prec_pl.get<std::string>("Prec Type"),
                            S_sigma_prec_pl.sublist("Prec Types"));
   }

}

} // namespace mini_em
