// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef ANASAZI_GENERALIZED_DAVIDSON_SOLMGR_HPP
#define ANASAZI_GENERALIZED_DAVIDSON_SOLMGR_HPP

/*! \file AnasaziGeneralizedDavidsonSolMgr.hpp
 *  \brief The Anasazi::GeneralizedDavidsonSolMgr provides a solver manager for the GeneralizedDavidson eigensolver.
*/

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCPDecl.hpp"

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziEigenproblem.hpp"
#include "AnasaziSolverManager.hpp"
#include "AnasaziBasicOrthoManager.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziICGSOrthoManager.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziGeneralizedDavidson.hpp"
#include "AnasaziStatusTestResNorm.hpp"
#include "AnasaziStatusTestWithOrdering.hpp"

using Teuchos::RCP;

namespace Anasazi {

/*!
 * \class GeneralizedDavidsonSolMgr
 * \brief Solver Manager for GeneralizedDavidson
 *
 * This class provides a simple interface to the GeneralizedDavidson
 * eigensolver.  This manager creates
 * appropriate orthogonalization/sort/output managers based on user
 * specified ParameterList entries (or selects suitable defaults),
 * provides access to solver functionality, and manages the restarting
 * process.
 *
 * This class is currently only implemented for real scalar types
 * (i.e. float, double).
 */
template <class ScalarType, class MV, class OP>
class GeneralizedDavidsonSolMgr : public SolverManager<ScalarType,MV,OP>
{
    public:

        /*!
         * \brief Basic constructor for GeneralizedDavidsonSolMgr
         *
         * This constructor accepts the Eigenproblem to be solved and a parameter list of options
         * for the solver.
         * The following options control the behavior
         * of the solver:
         * - "Which" -- a string specifying the desired eigenvalues: SM, LM, SR, LR, SI, or LI. Default: "LM."
         * - "Block Size" -- block size used by algorithm.  Default: 1.
         * - "Maximum Subspace Dimension" -- maximum number of basis vectors for subspace.  Two
         *  (for standard eigenvalue problems) or three (for generalized eigenvalue problems) sets of basis
         *  vectors of this size will be required. Default: 3*problem->getNEV()*"Block Size"
         * - "Restart Dimension" -- Number of vectors retained after a restart.  Default: NEV
         * - "Maximum Restarts" -- an int specifying the maximum number of restarts the underlying solver
         *  is allowed to perform.  Default: 20
         * - "Orthogonalization" -- a string specifying the desired orthogonalization: DGKS, SVQB, ICGS.
         *   Default: "SVQB"
         * - "Verbosity" -- a sum of MsgType specifying the verbosity.  Default: AnasaziErrors
         * - "Convergence Tolerance" -- a MagnitudeType specifying the level that residual norms must
         *  reach to decide convergence.  Default: machine precision
         * - "Relative Convergence Tolerance" -- a bool specifying whether residual norms should be
         *  scaled by the magnitude of the corresponding Ritz value.  Care should be taken when performing
         *  scaling for problems where the eigenvalue can be very large or very small.  Default: "false".
         * - "Initial Guess" -- how should initial vector be selected: "Random" or "User".
         *   If "User," the value in problem->getInitVec() will be used.  Default: "Random".
         * - "Print Number of Ritz Values" -- an int specifying how many Ritz values should be printed
         *   at each iteration.  Default: "NEV".
         */
        GeneralizedDavidsonSolMgr( const RCP< Eigenproblem<ScalarType,MV,OP> > &problem,
                                   Teuchos::ParameterList &pl );

        /*!
         * \brief Return the eigenvalue problem.
         */
        const Eigenproblem<ScalarType,MV,OP> & getProblem() const { return *d_problem; }

        /*!
         * \brief Get the iteration count for the most recent call to solve()
         */
        int getNumIters() const { return d_solver->getNumIters(); }

        /*!
         * \brief This method performs possibly repeated calls to the underlying eigensolver's iterate()
         *  routine until the problem has been solved (as decided by the StatusTest) or the solver manager decides to quit.
         */
        ReturnType solve();

    private:

        void getRestartState( GeneralizedDavidsonState<ScalarType,MV> &state );

        typedef MultiVecTraits<ScalarType,MV>        MVT;
        typedef Teuchos::ScalarTraits<ScalarType>    ST;
        typedef typename ST::magnitudeType           MagnitudeType;
        typedef Teuchos::ScalarTraits<MagnitudeType> MT;

        RCP< Eigenproblem<ScalarType,MV,OP> >           d_problem;
        RCP< GeneralizedDavidson<ScalarType,MV,OP> >    d_solver;
        RCP< OutputManager<ScalarType> >                d_outputMan;
        RCP< OrthoManager<ScalarType,MV> >              d_orthoMan;
        RCP< SortManager<MagnitudeType> >               d_sortMan;
        RCP< StatusTest<ScalarType,MV,OP> >             d_tester;
        int d_maxRestarts;
        int d_restartDim;

}; // class GeneralizedDavidsonSolMgr

//---------------------------------------------------------------------------//
// Prevent instantiation on complex scalar type
//---------------------------------------------------------------------------//
template <class MagnitudeType, class MV, class OP>
class GeneralizedDavidsonSolMgr<std::complex<MagnitudeType>,MV,OP>
{
  public:

    typedef std::complex<MagnitudeType> ScalarType;
    GeneralizedDavidsonSolMgr(
            const RCP<Eigenproblem<ScalarType,MV,OP> > &problem,
            Teuchos::ParameterList &pl )
    {
        // Provide a compile error when attempting to instantiate on complex type
        MagnitudeType::this_class_is_missing_a_specialization();
    }
};

//---------------------------------------------------------------------------//
// Start member definitions
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
GeneralizedDavidsonSolMgr<ScalarType,MV,OP>::GeneralizedDavidsonSolMgr(
        const RCP<Eigenproblem<ScalarType,MV,OP> > &problem,
        Teuchos::ParameterList &pl )
   : d_problem(problem)
{
    TEUCHOS_TEST_FOR_EXCEPTION( d_problem == Teuchos::null,                std::invalid_argument, "Problem not given to solver manager." );
    TEUCHOS_TEST_FOR_EXCEPTION( !d_problem->isProblemSet(),                std::invalid_argument, "Problem not set." );
    TEUCHOS_TEST_FOR_EXCEPTION( d_problem->getA() == Teuchos::null &&
                                d_problem->getOperator() == Teuchos::null, std::invalid_argument, "A operator not supplied on Eigenproblem." );
    TEUCHOS_TEST_FOR_EXCEPTION( d_problem->getInitVec() == Teuchos::null,  std::invalid_argument, "No vector to clone from on Eigenproblem." );
    TEUCHOS_TEST_FOR_EXCEPTION( d_problem->getNEV() <= 0,                  std::invalid_argument, "Number of requested eigenvalues must be positive.");

    if( !pl.isType<int>("Block Size") )
    {
        pl.set<int>("Block Size",1);
    }

    if( !pl.isType<int>("Maximum Subspace Dimension") )
    {
        pl.set<int>("Maximum Subspace Dimension",3*problem->getNEV()*pl.get<int>("Block Size"));
    }

    if( !pl.isType<int>("Print Number of Ritz Values") )
    {
        int numToPrint = std::max( pl.get<int>("Block Size"), d_problem->getNEV() );
        pl.set<int>("Print Number of Ritz Values",numToPrint);
    }

    // Get convergence info
    MagnitudeType tol = pl.get<MagnitudeType>("Convergence Tolerance", MT::eps() );
    TEUCHOS_TEST_FOR_EXCEPTION( pl.get<MagnitudeType>("Convergence Tolerance") <= MT::zero(),
                                std::invalid_argument, "Convergence Tolerance must be greater than zero." );

    // Get maximum restarts
    if( pl.isType<int>("Maximum Restarts") )
    {
        d_maxRestarts = pl.get<int>("Maximum Restarts");
        TEUCHOS_TEST_FOR_EXCEPTION( d_maxRestarts < 0, std::invalid_argument, "Maximum Restarts must be non-negative" );
    }
    else
    {
        d_maxRestarts = 20;
    }

    // Get maximum restarts
    d_restartDim = pl.get<int>("Restart Dimension",d_problem->getNEV());
    TEUCHOS_TEST_FOR_EXCEPTION( d_restartDim < d_problem->getNEV(),
            std::invalid_argument, "Restart Dimension must be at least NEV" );

    // Get initial guess type
    std::string initType;
    if( pl.isType<std::string>("Initial Guess") )
    {
        initType = pl.get<std::string>("Initial Guess");
        TEUCHOS_TEST_FOR_EXCEPTION( initType!="User" && initType!="Random", std::invalid_argument,
                                    "Initial Guess type must be 'User' or 'Random'." );
    }
    else
    {
        initType = "User";
    }

    // Get sort type
    std::string which;
    if( pl.isType<std::string>("Which") )
    {
        which = pl.get<std::string>("Which");
        TEUCHOS_TEST_FOR_EXCEPTION( which!="LM" && which!="SM" && which!="LR" && which!="SR" && which!="LI" && which!="SI",
                                    std::invalid_argument,
                                    "Which must be one of LM,SM,LR,SR,LI,SI." );
    }
    else
    {
        which = "LM";
    }

    // Build sort manager (currently must be stored as pointer to derived class)
    d_sortMan = Teuchos::rcp( new BasicSort<MagnitudeType>(which) );

    // Build orthogonalization manager
    std::string ortho = pl.get<std::string>("Orthogonalization","SVQB");
    TEUCHOS_TEST_FOR_EXCEPTION( ortho!="DGKS" && ortho!= "SVQB" && ortho!="ICGS", std::invalid_argument,
                                "Anasazi::GeneralizedDavidsonSolMgr::constructor: Invalid orthogonalization type" );

    if( ortho=="DGKS" )
    {
        d_orthoMan = Teuchos::rcp( new BasicOrthoManager<ScalarType,MV,OP>() );
    }
    else if( ortho=="SVQB" )
    {
        d_orthoMan = Teuchos::rcp( new SVQBOrthoManager<ScalarType,MV,OP>() );
    }
    else if( ortho=="ICGS" )
    {
        d_orthoMan = Teuchos::rcp( new ICGSOrthoManager<ScalarType,MV,OP>() );
    }

    // Build StatusTest
    bool scaleRes  = false; // Always false, scaling the residual is handled by the solver
    bool failOnNaN = false;
    RCP<StatusTest<ScalarType,MV,OP> > resNormTest = Teuchos::rcp(
            new StatusTestResNorm<ScalarType,MV,OP>(tol,d_problem->getNEV(),
                                    StatusTestResNorm<ScalarType,MV,OP>::RES_2NORM,scaleRes,failOnNaN) );
    d_tester = Teuchos::rcp( new StatusTestWithOrdering<ScalarType,MV,OP>(resNormTest,d_sortMan,d_problem->getNEV()) );

    // Build output manager
    int verbosity = pl.get<int>("Verbosity",Errors);
    d_outputMan = Teuchos::rcp( new BasicOutputManager<ScalarType>() );
    d_outputMan->setVerbosity( verbosity );

    // Build solver
    d_outputMan->stream(Debug) << " >> Anasazi::GeneralizedDavidsonSolMgr: Building solver" << std::endl;
    d_solver = Teuchos::rcp( new GeneralizedDavidson<ScalarType,MV,OP>( problem, d_sortMan, d_outputMan, d_tester, d_orthoMan, pl ) );
}

//---------------------------------------------------------------------------//
// Solve
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
ReturnType GeneralizedDavidsonSolMgr<ScalarType,MV,OP>::solve()
{
    Eigensolution<ScalarType,MV> sol;
    sol.numVecs = 0;
    d_problem->setSolution(sol);

    d_solver->initialize();
    int restarts = 0;
    while( 1 )
    {
        // Call iterate on the solver
        d_solver->iterate();

        // If the solver converged, we're done
        if( d_tester->getStatus() == Passed )
            break;

        // If we're already at maximum number of restarts, wrap it up
        if( restarts == d_maxRestarts )
            break;

        // We need to restart
        d_solver->sortProblem( d_restartDim );
        GeneralizedDavidsonState<ScalarType,MV> state = d_solver->getState();
        getRestartState( state );
        d_solver->initialize( state );
        restarts++;
    }

    // Output final state
    if( d_outputMan->isVerbosity(FinalSummary) )
        d_solver->currentStatus(d_outputMan->stream(FinalSummary));

    // Fill solution struct
    sol.numVecs = d_tester->howMany();
    if( sol.numVecs > 0 )
    {
        std::vector<int> whichVecs = d_tester->whichVecs();
        std::vector<int> origIndex = d_solver->getRitzIndex();

        // Make sure no conjugate pairs are split
        // Because these are not sorted we have to check all values
        for( int i=0; i<sol.numVecs; ++i )
        {
            if( origIndex[ whichVecs[i] ] == 1 )
            {
                if( std::find( whichVecs.begin(), whichVecs.end(), whichVecs[i]+1 ) == whichVecs.end() )
                {
                    whichVecs.push_back( whichVecs[i]+1 );
                    sol.numVecs++;
                }
            }
            else if( origIndex[ whichVecs[i] ] == -1 )
            {
                if( std::find( whichVecs.begin(), whichVecs.end(), whichVecs[i]-1 ) == whichVecs.end() )
                {
                    whichVecs.push_back( whichVecs[i]-1 );
                    sol.numVecs++;
                }
            }
        }

        if( d_outputMan->isVerbosity(Debug) )
        {
            d_outputMan->stream(Debug) << " >> Anasazi::GeneralizedDavidsonSolMgr: "
                << sol.numVecs << " eigenpairs converged" << std::endl;
        }

        // Sort converged values
        std::vector< Value<ScalarType> > origVals = d_solver->getRitzValues();
        std::vector<MagnitudeType> realParts;
        std::vector<MagnitudeType> imagParts;
        for( int i=0; i<sol.numVecs; ++i )
        {
            realParts.push_back( origVals[whichVecs[i]].realpart );
            imagParts.push_back( origVals[whichVecs[i]].imagpart );
        }

        std::vector<int> permVec(sol.numVecs);
        d_sortMan->sort( realParts, imagParts, Teuchos::rcpFromRef(permVec), sol.numVecs );

        // Create new which vector
        std::vector<int> newWhich;
        for( int i=0; i<sol.numVecs; ++i )
            newWhich.push_back( whichVecs[permVec[i]] );

        // Check if converged vectors are ordered
        bool ordered = true;
        for( int i=0; i<sol.numVecs; ++i )
        {
            if( newWhich[i]!=i )
            {
                ordered = false;
                break;
            }
        }

        if( ordered  )
        {
            // Everything is ordered, pull directly from solver and resize
            sol.index = origIndex;
            sol.index.resize(sol.numVecs);
            sol.Evals = d_solver->getRitzValues();
            sol.Evals.resize(sol.numVecs);
        }
        else
        {
            // Manually copy values into sol

            sol.index.resize(sol.numVecs);
            sol.Evals.resize(sol.numVecs);

            for( int i=0; i<sol.numVecs; ++i )
            {
                sol.index[i] = origIndex[ newWhich[i] ];
                sol.Evals[i] = origVals[  newWhich[i] ];
            }
        }
        sol.Evecs = MVT::CloneCopy( *(d_solver->getRitzVectors()), newWhich );
    }
    d_problem->setSolution(sol);

    // Return convergence status
    if( sol.numVecs < d_problem->getNEV() )
        return Unconverged;

    return Converged;
}

//---------------------------------------------------------------------------//
// Update GeneralizedDavidson state for restarting
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
void GeneralizedDavidsonSolMgr<ScalarType,MV,OP>::getRestartState(
        GeneralizedDavidsonState<ScalarType,MV> &state )
{
    TEUCHOS_TEST_FOR_EXCEPTION( state.curDim <= d_restartDim, std::runtime_error,
            "Anasazi::GeneralizedDavidsonSolMgr: State dimension at restart is smaller than Restart Dimension" );

    std::vector<int> ritzIndex = d_solver->getRitzIndex();

    // Don't split conjugate pair when restarting
    int restartDim = d_restartDim;
    if( ritzIndex[d_restartDim-1]==1 )
        restartDim++;

    d_outputMan->stream(Debug) << " >> Anasazi::GeneralizedDavidsonSolMgr: Restarting with "
        << restartDim << " vectors" << std::endl;

    // We have already sorted the problem with d_restartDim "best" values
    //  in the leading position.  If we partition the Schur vectors (Z)
    //  of the projected problem as Z = [Z_wanted Z_unwanted], then the
    //  search subspace after the restart is V_restart = V*Z_wanted
    //  (same for AV,BV)

    // Get view of wanted portion of Z
    const Teuchos::SerialDenseMatrix<int,ScalarType> Z_wanted =
        Teuchos::SerialDenseMatrix<int,ScalarType>(Teuchos::View,*state.Z,state.curDim,restartDim);

    // Get indices for restart
    std::vector<int> allIndices(state.curDim);
    for( int i=0; i<state.curDim; ++i )
        allIndices[i] = i;

    RCP<const MV>  V_orig = MVT::CloneView( *state.V,  allIndices );

    // Get indices for restart
    std::vector<int> restartIndices(restartDim);
    for( int i=0; i<restartDim; ++i )
        restartIndices[i] = i;

    // Views of subspace vectors to be updated
    RCP<MV>  V_restart  = MVT::CloneViewNonConst( *state.V,  restartIndices );

    // Temp storage
    RCP<MV> restartVecs = MVT::Clone(*state.V,restartDim);

    // Reset V
    MVT::MvTimesMatAddMv(ST::one(),*V_orig,Z_wanted,ST::zero(),*restartVecs);
    MVT::SetBlock(*restartVecs,restartIndices,*V_restart);

    // V, Z each have orthonormal columns, therefore V*Z should as well
    if( d_outputMan->isVerbosity(Debug) )
    {
        MagnitudeType orthErr = d_orthoMan->orthonormError(*V_restart);
        std::stringstream os;
        os << " >> Anasazi::GeneralizedDavidsonSolMgr: Error in V^T V == I after restart : " << orthErr << std::endl;
        d_outputMan->print(Debug,os.str());
    }

    // Reset AV
    RCP<MV> AV_restart  = MVT::CloneViewNonConst( *state.AV, restartIndices );
    RCP<const MV> AV_orig = MVT::CloneView( *state.AV, allIndices );

    MVT::MvTimesMatAddMv(ST::one(),*AV_orig,Z_wanted,ST::zero(),*restartVecs);
    MVT::SetBlock(*restartVecs,restartIndices,*AV_restart);

    int err;

    // Update matrix projection as Z^{*}(V^{*}AV)Z
    const Teuchos::SerialDenseMatrix<int,ScalarType> VAV_orig( Teuchos::View, *state.VAV, state.curDim, state.curDim );
    Teuchos::SerialDenseMatrix<int,ScalarType> tmpMat(state.curDim, restartDim);
    err = tmpMat.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, ST::one(), VAV_orig, Z_wanted, ST::zero() );
    TEUCHOS_TEST_FOR_EXCEPTION( err!=0, std::runtime_error, "GeneralizedDavidsonSolMgr::getRestartState: multiply returned nonzero error code" );

    Teuchos::SerialDenseMatrix<int,ScalarType> VAV_restart( Teuchos::View, *state.VAV, restartDim, restartDim );
    err = VAV_restart.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, ST::one(), Z_wanted, tmpMat, ST::zero() );
    TEUCHOS_TEST_FOR_EXCEPTION( err!=0, std::runtime_error, "GeneralizedDavidsonSolMgr::getRestartState: multiply returned nonzero error code" );

    if( d_problem->getM() != Teuchos::null )
    {
        // Reset BV
        RCP<const MV> BV_orig     = MVT::CloneView( *state.BV, allIndices );
        RCP<MV>       BV_restart  = MVT::CloneViewNonConst( *state.BV, restartIndices );

        MVT::MvTimesMatAddMv(ST::one(),*BV_orig,Z_wanted,ST::zero(),*restartVecs);
        MVT::SetBlock(*restartVecs,restartIndices,*BV_restart);


        // Update matrix projection as Z^{*}(V^{*}BV)Z
        const Teuchos::SerialDenseMatrix<int,ScalarType> VBV_orig( Teuchos::View, *state.VBV, state.curDim, state.curDim );
        err = tmpMat.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, ST::one(), VBV_orig, Z_wanted, ST::zero() );
        TEUCHOS_TEST_FOR_EXCEPTION( err!=0, std::runtime_error, "GeneralizedDavidsonSolMgr::getRestartState: multiply returned nonzero error code" );

        Teuchos::SerialDenseMatrix<int,ScalarType> VBV_restart( Teuchos::View, *state.VBV, restartDim, restartDim );
        VBV_restart.multiply( Teuchos::TRANS, Teuchos::NO_TRANS, ST::one(), Z_wanted, tmpMat, ST::zero() );
        TEUCHOS_TEST_FOR_EXCEPTION( err!=0, std::runtime_error, "GeneralizedDavidsonSolMgr::getRestartState: multiply returned nonzero error code" );
    }

    // Set Q,Z to identity
    state.Q->putScalar( ST::zero() );
    state.Z->putScalar( ST::zero() );
    for( int ii=0; ii<restartDim; ii++ )
    {
       (*state.Q)(ii,ii)= ST::one();
       (*state.Z)(ii,ii)= ST::one();
    }

    // Update current dimension
    state.curDim = restartDim;
}

} // namespace Anasazi

#endif // ANASAZI_GENERALIZED_DAVIDSON_SOLMGR_HPP

