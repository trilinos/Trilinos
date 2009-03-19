#include "PB_Helpers.hpp"

// Thyra includes
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"

#include <cmath>

// Anasazi includes
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziThyraAdapter.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziStatusTestMaxIters.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace PB {

Teuchos::RCP<Thyra::LinearOpBase<double> > inverse(const Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > & lowsFact,
                                                   const Teuchos::RCP<const Thyra::PreconditionerFactoryBase<double> > & precFact,
                                                   const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A)
{
   using Teuchos::RCP;

   // construct the preconditioner
   RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec(*precFact,A);

   // initialize the solver with a preconditioner
   RCP<Thyra::LinearOpWithSolveBase<double> > A_lows = lowsFact->createOp(); 
   Thyra::initializePreconditionedOp<double>(*lowsFact,A,prec,&*A_lows);

   // construct a Handle level linear operator for direct use and manipulation
   return Thyra::inverse<double>(A_lows);
}

double computeSpectralRad(const RCP<const Thyra::LinearOpBase<double> > & A, double tol,
                          bool isHermitian,int numBlocks,int restart,int verbosity)
{
   typedef Thyra::LinearOpBase<double> OP;
   typedef Thyra::MultiVectorBase<double> MV;

   int startVectors = 1;

   // construct an initial guess
   const RCP<MV> ivec = Thyra::createMember(A->domain());
   Thyra::randomize(-1.0,1.0,&*ivec);
   
   RCP<Anasazi::BasicEigenproblem<double,MV,OP> > eigProb
         = rcp(new Anasazi::BasicEigenproblem<double,MV,OP>(A,ivec));
   eigProb->setNEV(1);
   eigProb->setHermitian(isHermitian);

   // set the problem up
   if(not eigProb->setProblem()) {
      // big time failure!
      return Teuchos::ScalarTraits<double>::nan();
   }

   // we want largert magnitude eigenvalue
   std::string which("LM"); // largest magnitude

   // Create the parameter list for the eigensolver
   Teuchos::ParameterList MyPL;
   MyPL.set( "Verbosity", verbosity );
   MyPL.set( "Which", which );
   MyPL.set( "Block Size", startVectors );
   MyPL.set( "Num Blocks", numBlocks );
   MyPL.set( "Maximum Restarts", restart ); 
   MyPL.set( "Convergence Tolerance", tol );

   // build status test manager
   RCP<Anasazi::StatusTestMaxIters<double,MV,OP> > statTest
         = rcp(new Anasazi::StatusTestMaxIters<double,MV,OP>(10));

   // Create the Block Krylov Schur solver
   // This takes as inputs the eigenvalue problem and the solver parameters
   Anasazi::BlockKrylovSchurSolMgr<double,MV,OP> MyBlockKrylovSchur(eigProb, MyPL );
 
   // Solve the eigenvalue problem, and save the return code
   Anasazi::ReturnType solverreturn = MyBlockKrylovSchur.solve();

   if(solverreturn==Anasazi::Unconverged) {
      // cout << "Anasazi::BlockKrylovSchur::solve() did not converge!" << endl;
      return -std::abs(MyBlockKrylovSchur.getRitzValues().begin()->realpart);
   }
   else { // solverreturn==Anasazi::Converged
      // cout << "Anasazi::BlockKrylovSchur::solve() converged!" << endl;
      return std::abs(eigProb->getSolution().Evals.begin()->realpart);
   }
}

}
