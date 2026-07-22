// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** 
    \class BelosKrylov
    \brief Provides interface for using ROL::Vector with Belos solvers
                
    \author Created by Greg von Winckel
*/


#ifndef ROL_BELOS_KRYLOV_HPP
#define ROL_BELOS_KRYLOV_HPP

// The GenericSolverFactory will register managers for any types.
// To be used when the types are not part of the standard set.
// You must explicitly use GenericSolverFactory instead of SolverFactory.
#include "BelosSolverFactory_Generic.hpp"

#include "ROL_Krylov.hpp"
#include "ROL_BelosMultiVector.hpp"
#include "ROL_BelosOperator.hpp"
#include "ROL_MultiVectorDefault.hpp"

namespace ROL {

    template<class Real>
    class BelosKrylov : public Krylov<Real> {
        typedef Real                      ST;
        typedef LinearOperator<ST>        OP;
        typedef Vector<Real>              V; 
        typedef MultiVector<ST>           MV;
        typedef MultiVectorDefault<ST>    MVD;

        // For testing
	typedef Belos::MultiVecTraits<ST,MV>    MVT;
        typedef Belos::OperatorTraits<ST,MV,OP> OPT;

        private:

            Belos::GenericSolverFactory<ST,MV,OP> factory_;
            Teuchos::RCP<Belos::SolverManager<ST,MV,OP> > solver_;
            Teuchos::RCP<Belos::LinearProblem<ST,MV,OP> > problem_;  

        public:
           
            /// \brief Create a Belos solver 
            BelosKrylov(Teuchos::ParameterList &parlist) : 
                problem_(ROL::makePtr<Belos::LinearProblem<ST,MV,OP>>()) {

                Teuchos::RCP<Teuchos::ParameterList> solverParams = Teuchos::rcp( new Teuchos::ParameterList() );

                // Options likely to be of interest include CG, MINRES, GMRES, and RCG
                int blockSize          = 1; // Only support single solution & single RHS for now 
                std::string solverName = parlist.get("Krylov Method","MINRES");  
                int maxit              = parlist.get("Maximum Number of Krylov Iterations",50);
                Real abstol            = parlist.get("Absolute Krylov Tolerance",1.e-4);
                int numVectors         = parlist.get("Number of Stored Vectors",3);
 
                solverParams->setName("Belos input parameters"); 
                solverParams->set("Block Size",blockSize);
                solverParams->set("Maximum Iterations",maxit);
                solverParams->set("Convergence Tolerance",abstol);  
                solverParams->set("Num Blocks",numVectors);

                solver_ = factory_.create(solverName,solverParams);                 
            }


            /// \brief Compute solution vector
            Real run( V &x, OP& A, const V &b, OP &M, int &iter, int &flag )  {

                
                
                ;


                // Get pointers to ROL::Vectors
                Teuchos::RCP<V>        xp = Teuchos::rcpFromRef(x);

                // Wasteful, but have not yet implemented const case for MV
                Teuchos::RCP<V>        bp = b.clone();
                bp->set(b);

                // Make ROL::MultiVectors from the pointers to ROL::Vectors
                Teuchos::RCP<MV> xmvp = Teuchos::rcp( new MultiVectorDefault<Real>(xp) );
                Teuchos::RCP<MV> bmvp = Teuchos::rcp( new MultiVectorDefault<Real>(bp) );

                Teuchos::RCP<OP> Ap = Teuchos::rcpFromRef(A);
                Teuchos::RCP<OP> Mp = Teuchos::rcpFromRef(M);

                // Wrap x and b in ROL::MultiVector objects 
                MVD xmv(xp);
                MVD bmv(bp);
 
                problem_->setOperator(Ap);
                problem_->setLeftPrec(Mp);
                problem_->setProblem(xmvp,bmvp);

                solver_->setProblem(problem_);

                flag = static_cast<int>(solver_->solve());

                iter = solver_->getNumIters();

                return solver_->achievedTol();
            }
    };
}

#endif
