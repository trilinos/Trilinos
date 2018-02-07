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

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include <Stratimikos_MueLuHelpers.hpp>
#include "Thyra_DefaultDiagonalLinearOp.hpp"

#include <BelosTypes.hpp>

using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

namespace mini_em {

void writeOut(const std::string & s,const Thyra::LinearOpBase<double> & op)
{
  using Teuchos::RCP;

  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType NT;
  const RCP<const Thyra::TpetraLinearOp<double,int,panzer::Ordinal64,NT> > tOp = rcp_dynamic_cast<const Thyra::TpetraLinearOp<double,int,panzer::Ordinal64,NT> >(Teuchos::rcpFromRef(op));
  if(tOp != Teuchos::null) {
    const RCP<const Tpetra::CrsMatrix<double,int,panzer::Ordinal64,NT> > crsOp = rcp_dynamic_cast<const Tpetra::CrsMatrix<double,int,panzer::Ordinal64,NT> >(tOp->getConstTpetraOperator(),true);
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<double,int,panzer::Ordinal64,NT> >::writeSparseFile(s.c_str(),crsOp);
  }
}


///////////////////////////////////////
// FullMaxwellPreconditionerFactory  //
///////////////////////////////////////

Teko::LinearOp FullMaxwellPreconditionerFactory::buildPreconditionerOperator(Teko::BlockedLinearOp & blo, Teko::BlockPreconditionerState & /* state */) const
{
   Teuchos::RCP<Teuchos::TimeMonitor> tM = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("MaxwellPreconditioner::build"))));

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

   //for refmaxwell: Q_rho = M_0(epsilon / dt / cfl^2 / min_dx^2)
   // S_E = Q_E - Kt * Q_B^-1 * K = 1/dt * M_1(eps) + dt * D_1^T * M_2(1/mu) * D_1
   // addon: dt * M_1(1) * D_0 * M_0(mu)^-1 * D_0^T * M_1(1) 
   
   if(!use_refmaxwell) // Augmentation based solver
   {
     // Inverse of B mass matrix
     *Teko::getOutputStream() << "Building B inverse operator" << std::endl;
     Teko::LinearOp invQ_B = Teko::buildInverse(*invLib.getInverseFactory("Q_B Solve"),Q_B);

     // Compute the approximate Schur complement
     Teko::LinearOp idQ_B = Teko::getInvDiagonalOp(Q_B,Teko::AbsRowSum);
     Teko::LinearOp KtK   = Teko::explicitMultiply(Kt,idQ_B,K);
     Teko::LinearOp S_E   = Teko::explicitAdd(Q_E, Thyra::scale(-1.0,KtK));

     // Get auxiliary operators for gradient and nodal mass matrix
     Teko::LinearOp G     = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Weak Gradient"));
     Teko::LinearOp Gt    = Teko::explicitTranspose(G);
     Teko::LinearOp Q_rho = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix AUXILIARY_NODE"));

     // Compute grad-div term
     Teko::LinearOp idQ_rho = Teko::getInvDiagonalOp(Q_rho,Teko::AbsRowSum);
     Teko::LinearOp GGt     = Teko::explicitMultiply(G,idQ_rho,Gt);

     // Rescale such that grad-div is large enough to fix curl-curl null-space while not dominating
     double scaling = Teko::infNorm(Q_E)/Teko::infNorm(GGt);

     // Augmented Schur complement and its inverse
     Teko::LinearOp T_E = Teko::explicitAdd(S_E, Thyra::scale(scaling,GGt));
     *Teko::getOutputStream() << "Building T_E inverse operator" << std::endl;
     Teko::LinearOp invT_E = Teko::buildInverse(*invLib.getInverseFactory("T_E Solve"),T_E);

     // Correction term
     Teko::LinearOp Z_E = Thyra::add(Q_E, Thyra::scale(scaling,GGt));

     // Mass inverse - diagonal approximation
     Teko::LinearOp invQ_E = Teko::getInvDiagonalOp(Q_E,Teko::AbsRowSum);
 
     /////////////////////////////////////////////////
     // Build block upper triangular inverse matrix //
     /////////////////////////////////////////////////
     Teko::LinearOp invU;

     // Inverse blocks
     std::vector<Teko::LinearOp> diag(2);
     diag[0] = invQ_B;
     diag[1] = Teko::multiply(invQ_E,Z_E,invT_E);

     // Upper tri blocks
     Teko::BlockedLinearOp U = Teko::createBlockedOp();
     Teko::beginBlockFill(U,rows,rows);
        Teko::setBlock(0,0,U,Q_B);
        Teko::setBlock(1,1,U,S_E);
        Teko::setBlock(0,1,U,K);
     Teko::endBlockFill(U);

     // return upper tri preconditioner
     return(Teko::createBlockUpperTriInverseOp(U,diag));
   }
   else {// refMaxwell

     // Inverse of B mass matrix
     *Teko::getOutputStream() << "Building Q_B inverse operator" << std::endl;
     Teko::LinearOp invDiagQ_B = Teko::getInvDiagonalOp(Q_B,Teko::Diagonal);
     Teko::LinearOp invQ_B = Teko::buildInverse(*invLib.getInverseFactory("Q_B Solve"),Q_B, invDiagQ_B);

     // Compute the approximate Schur complement
     Teko::LinearOp idQ_B = Teko::getInvDiagonalOp(Q_B,Teko::AbsRowSum);
     Teko::LinearOp KtK   = Teko::explicitMultiply(Kt,idQ_B,K);
     Teko::LinearOp S_E   = Teko::explicitAdd(Q_E, Thyra::scale(-1.0,KtK));

     // Get nodal mass matrix and discrete gradient
     // Q_rho = M_0(mu)
     Teko::LinearOp Q_rho = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Mass Matrix AUXILIARY_NODE"));
     Teko::LinearOp T     = getRequestHandler()->request<Teko::LinearOp>(Teko::RequestMesg("Discrete Gradient"));
     // Teko::LinearOp KT = Teko::explicitMultiply(K,T);
     // TEUCHOS_ASSERT(Teko::infNorm(KT) < 1.0e-14 * Teko::infNorm(T) * Teko::infNorm(K));

     // Get inverse of lumped Q_rho
     RCP<Thyra::VectorBase<double> > ones = Thyra::createMember(Q_rho->domain());
     RCP<Thyra::VectorBase<double> > diagonal = Thyra::createMember(Q_rho->range());
     // set to all ones
     Thyra::assign(ones.ptr(),1.0);
     // compute lumped diagonal
     Thyra::apply(*Q_rho,Thyra::NOTRANS,*ones,diagonal.ptr());
     Thyra::reciprocal(*diagonal,diagonal.ptr());
     RCP<const Thyra::DiagonalLinearOpBase<double> > invDiagQ_rho = rcp(new Thyra::DefaultDiagonalLinearOp<double>(diagonal));

     // Get coordinates
     Teuchos::ParameterList SList2 = *invLib.getInverseFactory("S_E Solve")->getParameterList();
     Teuchos::RCP<Tpetra::MultiVector<double, int, panzer::Ordinal64> > Coordinates = SList2.get<Teuchos::RCP<Tpetra::MultiVector<double, int, panzer::Ordinal64> > >("Coordinates");

     *Teko::getOutputStream() << "Building S_E inverse operator" << std::endl;

     // Build the rest of the Stratimikos list
     Teuchos::ParameterList SList;
     SList.set("Linear Solver Type","Belos");
     SList.sublist("Linear Solver Types").sublist("Belos").set("Solver Type", "Pseudo Block CG");
     // SList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist("Pseudo Block CG").set("Output Frequency",1);
     SList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist("Pseudo Block CG").set("Maximum Iterations",500);
     SList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist("Pseudo Block CG").set("Convergence Tolerance",1e-5);
     // SList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist("Pseudo Block CG").set("Output Style",1);
     // SList.sublist("Linear Solver Types").sublist("Belos").sublist("Solver Types").sublist("Pseudo Block CG").set("Verbosity",33);
     // SList.sublist("Linear Solver Types").sublist("Belos").sublist("VerboseObject").set("Verbosity Level", "medium");

     SList.set("Preconditioner Type","MueLuRefMaxwell");
     // Teuchos::ParameterList refMaxwellPL = SList.sublist("Preconditioner Types").sublist("MueLuRefMaxwell");
     Teuchos::ParameterList refMaxwellPL;
     refMaxwellPL.set("parameterlist: syntax","muelu");
     refMaxwellPL.set("refmaxwell: mode","additive");
     refMaxwellPL.set("refmaxwell: disable addon",false);
     refMaxwellPL.set("refmaxwell: dump matrices",true);
     refMaxwellPL.set("refmaxwell: max coarse size",25);
     refMaxwellPL.set("refmaxwell: max levels",4);
     refMaxwellPL.set("smoother: type","CHEBYSHEV");

     Teuchos::ParameterList params11 = refMaxwellPL.sublist("refmaxwelll: 11list");
     params11.set("coarse: max size", 128);
     params11.set("number of equations",3);

     Teuchos::ParameterList params22 = refMaxwellPL.sublist("refmaxwelll: 22list");
     params22.set("coarse: max size", 128);

     refMaxwellPL.set("D0",T);
     refMaxwellPL.set("M0inv",invDiagQ_rho);
     refMaxwellPL.set("M1",Q_E);
     refMaxwellPL.set("Coordinates",Coordinates);

     SList.sublist("Preconditioner Types").set("MueLuRefMaxwell",refMaxwellPL);

     /* Stratimikos setup */
     Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
     Stratimikos::enableMueLuRefMaxwell<int,panzer::Ordinal64>(linearSolverBuilder);                // Register MueLu RefMaxwell as a Stratimikos preconditioner strategy.
     linearSolverBuilder.setParameterList(rcp(&SList,false));
     RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory = createLinearSolveStrategy(linearSolverBuilder);
     Teko::SolveInverseFactory siFactory = Teko::SolveInverseFactory(lowsFactory);
     Teko::LinearOp invS_E = Teko::buildInverse(siFactory, S_E);

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

     Teko::BlockedLinearOp invL = Teko::createBlockedOp();
     Teko::LinearOp id_B = Teko::identity(Teko::rangeSpace(Q_B));
     Teko::LinearOp id_E = Teko::identity(Teko::rangeSpace(Q_E));
     Teko::beginBlockFill(invL,rows,rows);
        Teko::setBlock(0,0,invL,id_B);
        Teko::setBlock(1,0,invL,Teko::multiply(Thyra::scale(-1.0, Kt), invQ_B));
        Teko::setBlock(1,1,invL,id_E);
     Teko::endBlockFill(invL);

     Teko::LinearOp prec = Teko::multiply(invU, Teko::toLinearOp(invL));
     return(prec);
   }

}

//! Initialize from a parameter list
void FullMaxwellPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList & pl)
{
   /////////////////////
   // Solver options  //
   // //////////////////            

   // Print residual for each sub-solve
   bool print_diagnostics = false;
   if(pl.isParameter("Print Diagnostics"))
     print_diagnostics = pl.get<bool>("Print Diagnostics");
   std::string name_append = "";
   if(print_diagnostics)
     name_append = " Base";

   // Use ILU smoother for Schur complement solve
   bool use_ilu = false;
   if(pl.isParameter("Use ILU"))
     use_ilu = pl.get<bool>("Use ILU");

   // Don't augment and use refMaxwell for S_E solve
   use_refmaxwell = false;
   if(pl.isParameter("Use refMaxwell"))
     use_refmaxwell = pl.get<bool>("Use refMaxwell");


   //////////////////////////////////
   // Set up sub-solve factories   //
   //////////////////////////////////

   // New inverse lib to add inverse factories to
   invLib = *getInverseLibrary();


   if (!use_refmaxwell){
     { // MueLu with Gauss-Seidel smoother
       Teuchos::ParameterList ml_pl("MueLu GS");
       ml_pl.set("Type", "MueLu-Tpetra");
       ml_pl.set("verbosity", "high");
       ml_pl.set("multigrid algorithm",      "unsmoothed");
       ml_pl.set("coarse: type",             "KLU2");
       ml_pl.set("coarse: max size",         2500);
       ml_pl.set("aggregation: type",        "uncoupled");
       ml_pl.set("aggregation: drop scheme", "classical");
       ml_pl.set("aggregation: drop tol",    0.0);
       ml_pl.set("smoother: pre or post",    "both");
       ml_pl.set("smoother: type",           "RELAXATION");
       {
         Teuchos::ParameterList& smoother = ml_pl.sublist("smoother: params");
         smoother.set("relaxation: type",           "MT Gauss-Seidel");
         smoother.set("relaxation: symmetric matrix structure",         true);
         smoother.set("relaxation: sweeps",         4);
         smoother.set("relaxation: damping factor", 1.0);
       }
       ml_pl.set("repartition: enable",true);
       ml_pl.set("repartition: partitioner","zoltan2");
       ml_pl.set("repartition: start level",2);
       ml_pl.set("repartition: min rows per proc",1024);
       ml_pl.set("repartition: max imbalance",1.327);
       ml_pl.set("repartition: remap parts",true);
       ml_pl.set("repartition: rebalance P and R",true);
       {
         Teuchos::ParameterList& repartition = ml_pl.sublist("repartition: params");
         repartition.set("algorithm","multijagged");
       }
       // add coordinates to parameter list
       {
         Teuchos::ParameterList& required = ml_pl.sublist("Required Parameters");
         required.set("Coordinates","B_face");
       }
       invLib.addInverse("Q_B Solve"+name_append,ml_pl);
     }

     if(!use_ilu)
       { // MueLu with Chebyshev smoother
         Teuchos::ParameterList ml_pl("MueLu Cheb");
         ml_pl.set("Type", "MueLu-Tpetra");
         ml_pl.set("verbosity", "high");
         ml_pl.set("multigrid algorithm",      "unsmoothed");
         ml_pl.set("coarse: type",             "KLU2");
         ml_pl.set("coarse: max size",         2500);
         ml_pl.set("aggregation: type",        "uncoupled");
         ml_pl.set("aggregation: drop scheme", "classical");
         ml_pl.set("aggregation: drop tol",    0.0);
         ml_pl.set("smoother: pre or post",    "both");
         ml_pl.set("smoother: type",           "CHEBYSHEV");
         {
           Teuchos::ParameterList& smoother = ml_pl.sublist("smoother: params");
           smoother.set("chebyshev: degree",2);
           smoother.set("chebyshev: ratio eigenvalue",20.0);
           smoother.set("chebyshev: min eigenvalue",1.0);
           smoother.set("chebyshev: eigenvalue max iterations",15);
         }
         ml_pl.set("repartition: enable",true);
         ml_pl.set("repartition: partitioner","zoltan2");
         ml_pl.set("repartition: start level",2);
         ml_pl.set("repartition: min rows per proc",2500);
         ml_pl.set("repartition: max imbalance",1.327);
         ml_pl.set("repartition: remap parts",true);
         ml_pl.set("repartition: rebalance P and R",true);
         {
           Teuchos::ParameterList& repartition = ml_pl.sublist("repartition: params");
           repartition.set("algorithm","multijagged");
         }
         // add coordinates to parameter list
         {
           Teuchos::ParameterList& required = ml_pl.sublist("Required Parameters");
           required.set("Coordinates","E_edge");
         }
         invLib.addInverse("T_E Solve"+name_append,ml_pl);
       }
     else
       { // MueLu with ILU smoother
         Teuchos::ParameterList ml_pl("MueLu ILU");
         ml_pl.set("Type", "MueLu-Tpetra");
         ml_pl.set("verbosity", "high");
         // ml_pl.set("number of equations", 3);
         ml_pl.set("multigrid algorithm",      "unsmoothed");
         ml_pl.set("coarse: type",             "KLU2");
         ml_pl.set("coarse: max size",         2500);
         ml_pl.set("aggregation: type",        "uncoupled");
         ml_pl.set("aggregation: drop scheme", "classical");
         ml_pl.set("aggregation: drop tol",    0.0);
         ml_pl.set("smoother: pre or post",    "both");
         ml_pl.set("smoother: type",           "SCHWARZ");
         {
           Teuchos::ParameterList& smoother = ml_pl.sublist("smoother: params");
           smoother.set("schwarz: overlap level", 1);
           smoother.set("schwarz: combine mode", "Zero");
           smoother.set("subdomain solver name", "RILUK");
           {
             Teuchos::ParameterList& subdomain = smoother.sublist("subdomain solver parameters");
             subdomain.set("fact: iluk level-of-fill", 1);
           }
         }
         ml_pl.set("repartition: enable",true);
         ml_pl.set("repartition: partitioner","zoltan2");
         ml_pl.set("repartition: start level",2);
         ml_pl.set("repartition: min rows per proc",2500);
         ml_pl.set("repartition: max imbalance",1.327);
         ml_pl.set("repartition: remap parts",true);
         ml_pl.set("repartition: rebalance P and R",true);
         {
           Teuchos::ParameterList& repartition = ml_pl.sublist("repartition: params");
           repartition.set("algorithm","multijagged");
         }
         // add coordinates to parameter list
         {
           Teuchos::ParameterList& required = ml_pl.sublist("Required Parameters");
           required.set("Coordinates","E_edge");
         }
         invLib.addInverse("T_E Solve"+name_append,ml_pl);
       }

     if(print_diagnostics){
       { // Diagnostic Q_B solve
         Teuchos::ParameterList diag_pl("Q_B Solve");
         diag_pl.set("Type","Diagnostic Inverse");
         diag_pl.set("Inverse Factory","Q_B Solve Base");
         diag_pl.set("Descriptive Label","Q_B");
         diag_pl.set("Print Residual",true);
         invLib.addInverse("Q_B Solve",diag_pl);
       }
       { // Diagnostic T_E solve
         Teuchos::ParameterList diag_pl("T_E Solve");
         diag_pl.set("Type","Diagnostic Inverse");
         diag_pl.set("Inverse Factory","T_E Solve Base");
         diag_pl.set("Descriptive Label","T_E");
         diag_pl.set("Print Residual",true);
         diag_pl.set("Convergence Tolerance", 1.0e-5);
         invLib.addInverse("T_E Solve",diag_pl);
       }
     }
   } else {
     { // Q_B solve
       Teuchos::ParameterList cg_pl("Belos CG");
       cg_pl.set("Type", "Belos");
       cg_pl.set("Solver Type", "Block CG");
       {
         Teuchos::ParameterList& st_pl = cg_pl.sublist("Solver Types");
         {
           Teuchos::ParameterList& bcg_pl = st_pl.sublist("Block CG");
           // bcg_pl.set("Verbosity", Belos::StatusTestDetails+Belos::FinalSummary+Belos::Warnings);
           // bcg_pl.set("Output Frequency", 1);
           // bcg_pl.set("Output Style", Belos::Brief);
           bcg_pl.set("Convergence Tolerance", 1.0e-5);
           bcg_pl.set("Maximum Iterations", 100);
         }
       }
       {
         Teuchos::ParameterList& verb = cg_pl.sublist("VerboseObject");
         verb.set("Output File", "none");
         // verb.set("Verbosity Level", "medium");
       }
       invLib.addInverse("Q_B Solve"+name_append,cg_pl);
     }

     { // MueLu RefMaxwell
       Teuchos::ParameterList ml_pl("Belos Block CG");
       ml_pl.set("Type", "MueLu-Tpetra");

       // add coordinates to parameter list
       {
         Teuchos::ParameterList& required = ml_pl.sublist("Required Parameters");
         required.set("Coordinates","AUXILIARY_NODE");
       }
       invLib.addInverse("S_E Solve"+name_append,ml_pl);
     }
     if(print_diagnostics){
       { // Diagnostic Q_B solve
         Teuchos::ParameterList diag_pl("Q_B Solve");
         diag_pl.set("Type","Diagnostic Inverse");
         diag_pl.set("Inverse Factory","Q_B Solve Base");
         diag_pl.set("Descriptive Label","Q_B");
         diag_pl.set("Print Residual",true);
         invLib.addInverse("Q_B Solve",diag_pl);
       }
       { // Diagnostic S_E solve
         Teuchos::ParameterList diag_pl("S_E Solve");
         diag_pl.set("Type","Diagnostic Inverse");
         diag_pl.set("Inverse Factory","S_E Solve Base");
         diag_pl.set("Descriptive Label","S_E");
         diag_pl.set("Print Residual",true);
         diag_pl.set("Convergence Tolerance", 1.0e-5);
         invLib.addInverse("S_E Solve",diag_pl);
       }
     }
   }
}
 
}
