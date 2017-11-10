#include "MiniEM_FullMaxwellPreconditionerFactory.hpp"

#include "Teko_BlockLowerTriInverseOp.hpp"
#include "Teko_BlockUpperTriInverseOp.hpp"

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

using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;

namespace mini_em {

///////////////////////////////////////
// FullMaxwellPreconditionerFactory  //
///////////////////////////////////////

Teko::LinearOp FullMaxwellPreconditionerFactory::buildPreconditionerOperator(Teko::BlockedLinearOp & blo, Teko::BlockPreconditionerState & state) const
{
   Teuchos::RCP<Teuchos::TimeMonitor> tM = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer(std::string("MaxwellPreconditioner::build"))));

   // Check that system is right size
   int rows = Teko::blockRowCount(blo);
   int cols = Teko::blockColCount(blo);
   TEUCHOS_ASSERT(rows==cols);
   TEUCHOS_ASSERT(rows==2);

   // Extract the blocks
   Teko::LinearOp Q_B   = Teko::getBlock(0,0,blo);
   Teko::LinearOp K     = Teko::getBlock(0,1,blo);
   Teko::LinearOp Kt    = Teko::getBlock(1,0,blo);
   Teko::LinearOp Q_E   = Teko::getBlock(1,1,blo);

   // Inverse of B mass matrix
   *Teko::getOutputStream() << "Building B inverse operator" << std::endl;
   Teko::LinearOp invQ_B = Teko::buildInverse(*invLib.getInverseFactory("Q_B Solve"),Q_B);

   // Compute the approximate Schur complement
   Teko::LinearOp idQ_B = Teko::getInvDiagonalOp(Q_B,Teko::AbsRowSum);
   Teko::LinearOp KtK   = Teko::explicitMultiply(Kt,idQ_B,K);
   Teko::LinearOp S_E   = Teko::explicitAdd(Q_E, Thyra::scale(-1.0,KtK));

   if(!use_refmaxwell) // Augmentation based solver
   {
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
   else // refMaxwell
     TEUCHOS_ASSERT(false);

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
     // TODO: implement refMaxwell version


   //////////////////////////////////
   // Set up sub-solve factories   //
   //////////////////////////////////

   // New inverse lib to add inverse factories to
   invLib = *getInverseLibrary();

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
       invLib.addInverse("T_E Solve",diag_pl);
     }
   }
}
 
}
