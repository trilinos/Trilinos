#include "MiniEM_FullMaxwellPreconditionerFactory.hpp"

#include "Teko_BlockLowerTriInverseOp.hpp"
#include "Teko_BlockUpperTriInverseOp.hpp"

#include "Teko_SolveInverseFactory.hpp"

#include "Thyra_DiagonalLinearOpBase.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"

#include "Teuchos_as.hpp"
#include "Teuchos_Time.hpp"

#include "Teko_TpetraHelpers.hpp"

#include "Thyra_TpetraLinearOp.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_LinearObjContainer.hpp"
#include "Panzer_ThyraObjContainer.hpp"

#include "Thyra_DefaultDiagonalLinearOp.hpp"

using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

namespace mini_em {

void writeOut(const std::string & s,const Thyra::LinearOpBase<double> & op)
{
  using Teuchos::RCP;

  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType NT;
  const RCP<const Thyra::TpetraLinearOp<double,int,panzer::Ordinal64,NT> > tOp = rcp_dynamic_cast<const Thyra::TpetraLinearOp<double,int,panzer::Ordinal64,NT> >(Teuchos::rcpFromRef(op));
  if(tOp != Teuchos::null) {
    *Teko::getOutputStream() << "Dumping matrix \'" << s << "\'" << std::endl;
    const RCP<const Tpetra::CrsMatrix<double,int,panzer::Ordinal64,NT> > crsOp = rcp_dynamic_cast<const Tpetra::CrsMatrix<double,int,panzer::Ordinal64,NT> >(tOp->getConstTpetraOperator(),true);
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::Ordinal64,NT> >::writeMapFile(("rowmap_"+s).c_str(),*(crsOp->getRowMap()));
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::Ordinal64,NT> >::writeMapFile(("colmap_"+s).c_str(),*(crsOp->getColMap()));
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::Ordinal64,NT> >::writeMapFile(("domainmap_"+s).c_str(),*(crsOp->getDomainMap()));
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::Ordinal64,NT> >::writeMapFile(("rangemap_"+s).c_str(),*(crsOp->getRangeMap()));
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::Ordinal64,NT> >::writeSparseFile(s.c_str(),crsOp);
  }
}


void describeMatrix(const std::string & s,const Thyra::LinearOpBase<double> & op,Teuchos::RCP<Teuchos::FancyOStream> out)
{
  using Teuchos::RCP;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType NT;
  if (out!=Teuchos::null) {
    const RCP<const Thyra::TpetraLinearOp<double,int,panzer::Ordinal64,NT> > tOp = rcp_dynamic_cast<const Thyra::TpetraLinearOp<double,int,panzer::Ordinal64,NT> >(Teuchos::rcpFromRef(op));
    if(tOp != Teuchos::null) {
      const RCP<const Tpetra::CrsMatrix<double,int,panzer::Ordinal64,NT> > crsOp = rcp_dynamic_cast<const Tpetra::CrsMatrix<double,int,panzer::Ordinal64,NT> >(tOp->getConstTpetraOperator(),true);
      *out << "\nDebug: " << s << std::endl;
      crsOp->describe(*out,Teuchos::VERB_MEDIUM);
    }
  }
}


///////////////////////////////////////
// FullMaxwellPreconditionerFactory  //
///////////////////////////////////////

Teko::LinearOp FullMaxwellPreconditionerFactory::buildPreconditionerOperator(Teko::BlockedLinearOp & blo, Teko::BlockPreconditionerState & /* state */) const
{
   Teuchos::RCP<Teuchos::TimeMonitor> tM = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("MaxwellPreconditioner::build"))));

   // Output stream for debug information
   Teuchos::RCP<Teuchos::FancyOStream> debug;
   // print debug information
   if (params.isParameter("Debug") && params.get<bool>("Debug"))
     debug = Teko::getOutputStream();
   else
     debug = Teuchos::null;

   // Check that system is right size
   int rows = Teko::blockRowCount(blo);
   int cols = Teko::blockColCount(blo);
   TEUCHOS_ASSERT(rows==cols);
   TEUCHOS_ASSERT(rows==2);

   // Extract the blocks
   Teko::LinearOp Q_B   = Teko::getBlock(0,0,blo);  // actually 1/dt * Q_B = mu/dt * M_2(1/mu)
   Teko::LinearOp K     = Teko::getBlock(0,1,blo);  // actually K = Q_B * D_1 = mu * M_2(1/mu) * D_1
   Teko::LinearOp Kt    = Teko::getBlock(1,0,blo);  // actually -Kt  = - D_1^T * M_2(1/mu)
   Teko::LinearOp Q_E   = Teko::getBlock(1,1,blo);  // actually 1/(c^2*dt) * Q_E = 1/dt * M_1(eps)

   // pull in nodal mass matrix
   Teko::LinearOp Q_rho = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix AUXILIARY_NODE"));

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

   //for refmaxwell: Q_rho = M_0(epsilon / dt / cfl^2 / min_dx^2)
   // S_E = Q_E - Kt * Q_B^-1 * K
   //     = 1/dt * M_1(eps) + dt * D_1^T * M_2(1/mu) * D_1
   // addon in paper: dt * M_1(1) * D_0 * M_0(mu)^-1 * D_0^T * M_1(1)
   // addon here:     1/dt * M_1(eps) * D_0 * M_0(epsilon / dt / cfl^2 / min_dx^2)^-1 * D_0^T * 1/dt * M_1(eps)
   //               = Q_E * D_0 * Q_rho^-1 * D_0^T * Q_E


   bool useDiscreteGradient = false;
   if (params.isParameter("Use discrete gradient"))
     useDiscreteGradient = params.get<bool>("Use discrete gradient");
   bool solveLowerTriangular = false;
   if (params.isParameter("Solve lower triangular"))
     solveLowerTriangular = params.get<bool>("Solve lower triangular");
   bool useAsPreconditioner = true;
   if (params.isParameter("Use as preconditioner"))
     useAsPreconditioner = params.get<bool>("Use as preconditioner");

   // Compute the approximate Schur complement
   Teko::LinearOp idQ_B = Teko::getInvDiagonalOp(Q_B,Teko::AbsRowSum);
   Teko::LinearOp KtK   = Teko::explicitMultiply(Kt,idQ_B,K);
   if (debug != Teuchos::null) 
     *debug << "Adding up S_E" << std::endl;
   Teko::LinearOp S_E   = Teko::explicitAdd(Q_E, Thyra::scale(-1.0,KtK));
   if (debug != Teuchos::null) 
     *debug << "Added up S_E" << std::endl;
   if (dump)
     writeOut("S_E.mm",*S_E);
   describeMatrix("S_E",*S_E,debug);

   Teko::LinearOp invQ_B,invS_E;

   // Inverse of B mass matrix
   *Teko::getOutputStream() << "Building Q_B inverse operator" << std::endl;
   // Are we building a solver or a preconditioner?
   if (useAsPreconditioner) {
     invQ_B = Teko::buildInverse(*invLib.getInverseFactory("Q_B Preconditioner"),Q_B);
   } else {
     Teko::LinearOp invDiagQ_B = Teko::buildInverse(*invLib.getInverseFactory("Q_B Preconditioner"),Q_B);
     describeMatrix("invDiagQ_B",*invDiagQ_B,debug);
     invQ_B = Teko::buildInverse(*invLib.getInverseFactory("Q_B Solve"),Q_B, invDiagQ_B);
   }

   if(!use_refmaxwell) // Augmentation based solver
   {
     // Get auxiliary operators for gradient and nodal mass matrix
     Teko::LinearOp G;
     if (useDiscreteGradient) {
       Teko::LinearOp T = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Gradient"));
       if (dump)
         writeOut("DiscreteGradient.mm",*T);
       describeMatrix("DiscreteGradient",*T,debug);
       double dt = params.get<double>("dt");
       double epsilon = params.get<double>("epsilon");
       G = Teko::explicitMultiply(Q_E,Teko::scale(-dt/epsilon,T));
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
     if (debug != Teuchos::null)
       *debug << "Adding up T_E" << std::endl;
     Teko::LinearOp T_E = Teko::explicitAdd(S_E, Thyra::scale(scaling,GGt));
     if (debug != Teuchos::null)
       *debug << "Added up T_E" << std::endl;
     *Teko::getOutputStream() << "Building T_E inverse operator" << std::endl;
     Teko::LinearOp invT_E = Teko::buildInverse(*invLib.getInverseFactory("T_E Solve"),T_E);

     // Correction term
     Teko::LinearOp Z_E = Thyra::add(Q_E, Thyra::scale(scaling,GGt));

     if (dump)
       writeOut("Z_E.mm",*Z_E);
     describeMatrix("Z_E",*Z_E,debug);

     // Mass inverse - diagonal approximation
     Teko::LinearOp invQ_E = Teko::getInvDiagonalOp(Q_E,Teko::AbsRowSum);

     // Approximate inverse of S_E
     invS_E = Teko::multiply(invQ_E,Z_E,invT_E);
   }
   else {// refMaxwell

     // Teko::LinearOp T = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Gradient"));
     // Teko::LinearOp KT = Teko::explicitMultiply(K,T);
     // TEUCHOS_ASSERT(Teko::infNorm(KT) < 1.0e-14 * Teko::infNorm(T) * Teko::infNorm(K));

     // Inverse of Schur complement
     *Teko::getOutputStream() << "Building S_E inverse operator" << std::endl;
     Teuchos::RCP<Teko::InverseFactory> S_E_prec_factory = invLib.getInverseFactory("S_E Preconditioner");
     Teuchos::ParameterList S_E_prec_pl = *S_E_prec_factory->getParameterList();

     // Get coordinates
     Teuchos::RCP<Tpetra::MultiVector<double, int, panzer::Ordinal64> > Coordinates = S_E_prec_pl.get<Teuchos::RCP<Tpetra::MultiVector<double, int, panzer::Ordinal64> > >("Coordinates");
     S_E_prec_pl.sublist("Preconditioner Types").sublist("MueLuRefMaxwell-Tpetra").set("Coordinates",Coordinates);
     S_E_prec_pl.remove("Coordinates");

     // Set M1 = Q_E.
     // We do this here, since we cannot get it from the request handler.
     S_E_prec_pl.sublist("Preconditioner Types").sublist("MueLuRefMaxwell-Tpetra").set("M1",Q_E);

     // Get inverse of lumped Q_rho
     RCP<Thyra::VectorBase<double> > ones = Thyra::createMember(Q_rho->domain());
     RCP<Thyra::VectorBase<double> > diagonal = Thyra::createMember(Q_rho->range());
     Thyra::assign(ones.ptr(),1.0);
     // compute lumped diagonal
     Thyra::apply(*Q_rho,Thyra::NOTRANS,*ones,diagonal.ptr());
     Thyra::reciprocal(*diagonal,diagonal.ptr());
     RCP<const Thyra::DiagonalLinearOpBase<double> > invDiagQ_rho = rcp(new Thyra::DefaultDiagonalLinearOp<double>(diagonal));
     S_E_prec_pl.sublist("Preconditioner Types").sublist("MueLuRefMaxwell-Tpetra").set("M0inv",invDiagQ_rho);

     Teko::InverseLibrary myInvLib = invLib;
     S_E_prec_pl.sublist("Preconditioner Types").sublist("MueLuRefMaxwell-Tpetra").set("Type","MueLuRefMaxwell-Tpetra");
     myInvLib.addInverse("S_E Preconditioner",S_E_prec_pl.sublist("Preconditioner Types").sublist("MueLuRefMaxwell-Tpetra"));
     S_E_prec_factory = myInvLib.getInverseFactory("S_E Preconditioner");

     // Are we building a solver or a preconditioner?
     if (useAsPreconditioner)
       invS_E = Teko::buildInverse(*S_E_prec_factory,S_E);
     else {
       Teko::LinearOp S_E_prec = Teko::buildInverse(*S_E_prec_factory,S_E);
       invS_E = Teko::buildInverse(*invLib.getInverseFactory("S_E Solve"),S_E,S_E_prec);
     }
   }

   /////////////////////////////////////////////////
   // Build block upper triangular inverse matrix //
   /////////////////////////////////////////////////

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

   if (solveLowerTriangular) {
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

//! Initialize from a parameter list
void FullMaxwellPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList & pl)
{
   /////////////////////
   // Solver options  //
   // //////////////////

   params = pl;

   // Don't augment and use refMaxwell for S_E solve
   use_refmaxwell = false;
   if(pl.isParameter("Use refMaxwell"))
     use_refmaxwell = pl.get<bool>("Use refMaxwell");

   // Output stream for debug information
   Teuchos::RCP<Teuchos::FancyOStream> debug;
   // print debug information
   bool doDebug = false;
   if (pl.isParameter("Debug"))
     doDebug = pl.get<bool>("Debug");
   if (doDebug)
     debug = Teko::getOutputStream();
   else
     debug = Teuchos::null;

   // dump matrices
   dump = false;
   if (pl.isParameter("Dump"))
     dump = pl.get<bool>("Dump");

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

   if (!use_refmaxwell){
     // T_E solve
     Teuchos::ParameterList T_E_pl = pl.sublist("T_E Solve");
     invLib.addInverse("T_E Solve",T_E_pl);

   } else { // RefMaxwell based solve
     // S_E solve
     Teuchos::ParameterList ml_pl = pl.sublist("S_E Solve");
     invLib.addInverse("S_E Solve",ml_pl);

     // S_E preconditioner
     Teuchos::ParameterList S_E_prec_pl = pl.sublist("S_E Preconditioner");

     // add discrete gradient
     Teko::LinearOp T = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Gradient"));
     S_E_prec_pl.sublist("Preconditioner Types").sublist("MueLuRefMaxwell-Tpetra").set("D0",T);

     if (dump) {
       writeOut("DiscreteGradient.mm",*T);
     }
     describeMatrix("DiscreteGradient",*T,debug);

     // add edge mass matrix
     // commented out, since the edge mass matrix isn't registered in the request handler
     // Teko::LinearOp Q_E = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix E_edge"));
     // S_E_prec_pl.sublist("Preconditioner Types").sublist("MueLuRefMaxwell-Tpetra-Tpetra").set("M1",Q_E);

     invLib.addInverse("S_E Preconditioner",S_E_prec_pl);
   }
}

} // namespace mini_em
