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

#ifndef ANASAZI_GENERALIZED_DAVIDSON_HPP
#define ANASAZI_GENERALIZED_DAVIDSON_HPP

/*! \file AnasaziGeneralizedDavidson.hpp
    \brief Implementation of a block Generalized Davidson eigensolver.

    \author Steven Hamilton
*/

#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziEigenproblem.hpp"
#include "AnasaziEigensolver.hpp"
#include "AnasaziOrthoManager.hpp"
#include "AnasaziOutputManager.hpp"
#include "AnasaziSortManager.hpp"
#include "AnasaziStatusTest.hpp"

using Teuchos::RCP;

namespace Anasazi {

/*!
 * \brief Structure to contain pointers to GeneralizedDavidson state variables.
 */
template <class ScalarType, class MV>
struct GeneralizedDavidsonState {
    /*! \brief The current subspace dimension. */
    int curDim;

    /*! \brief Orthonormal basis for search subspace. */
    RCP<MV> V;

    /*! \brief Image of V under A. */
    RCP<MV> AV;

    /*! \brief Image of V under B. */
    RCP<MV> BV;

    /*! \brief Projection of A onto V. */
    RCP< Teuchos::SerialDenseMatrix<int,ScalarType> > VAV;

    /*! \brief Projection of B onto V. */
    RCP< Teuchos::SerialDenseMatrix<int,ScalarType> > VBV;

    /*! \brief Left quasi upper triangular matrix from QZ decomposition of (VAV,VBV) */
    RCP< Teuchos::SerialDenseMatrix<int,ScalarType> > S;

    /*! \brief Right quasi upper triangular matrix from QZ decomposition of (VAV,VBV) */
    RCP< Teuchos::SerialDenseMatrix<int,ScalarType> > T;

    /*! \brief Left generalized Schur vectors from QZ decomposition of (VAV,VBV) */
    RCP< Teuchos::SerialDenseMatrix<int,ScalarType> > Q;

    /*! \brief Right generalized Schur vectors from QZ decomposition of (VAV,VBV) */
    RCP< Teuchos::SerialDenseMatrix<int,ScalarType> > Z;

    /*! \brief Vector of generalized eigenvalues */
    std::vector< Value<ScalarType> > eVals;

    GeneralizedDavidsonState() : curDim(0), V(Teuchos::null), AV(Teuchos::null),
                                 BV(Teuchos::null), VAV(Teuchos::null),
                                 VBV(Teuchos::null), S(Teuchos::null),
                                 T(Teuchos::null), Q(Teuchos::null),
                                 Z(Teuchos::null), eVals(0) {}

};


/*!
 * \class GeneralizedDavidson
 * \brief Solves eigenvalue problem using generalized Davidson method
 *
 * This class searches for a few eigenvalues and corresponding eigenvectors
 * for either a standard eigenvalue problem \f$Ax=\lambda x\f$
 * or a generalized eigenvalue problem \f$Ax=\lambda B x\f$
 * Note that unlike some other solvers, the generalized Davidson method places
 * no restrictions on either matrix in a generalized eigenvalue problem.
 *
 * Tips for preconditioning:  A good preconditioner usually approximates
 * \f$(A-\sigma I)^{-1}\f$ or \f$(A-\sigma B)^{-1}\f$, where \f$\sigma\f$
 * is close to the target eigenvalue.  When searching for largest magnitude
 * eigenvalues, selecting a preconditioner \f$P^{-1} \approx B^{-1}\f$
 * usually works well and when searching for smallest magnitude eigenvalues
 * selecting \f$P^{-1} \approx A^{-1}\f$ is usually appropriate.
 *
 * This class is currently only implemented for real scalar types
 * (i.e. float, double).
 */
template <class ScalarType, class MV, class OP>
class GeneralizedDavidson : public Eigensolver<ScalarType,MV,OP>
{
  private:
    // Convenience Typedefs
    typedef MultiVecTraits<ScalarType,MV>            MVT;
    typedef OperatorTraits<ScalarType,MV,OP>         OPT;
    typedef Teuchos::ScalarTraits<ScalarType>        ST;
    typedef typename ST::magnitudeType               MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType>     MT;

  public:

    /*!
     *  \brief Constructor.
     *
     * GeneralizedDavidson constructor with eigenproblem, parameters, and
     * solver utilities.
     *
     * Behavior of the solver is controlled by the following ParameterList
     * entries:
     * - "Block Size" -- block size used by algorithm.  Default: 1.
     * - "Maximum Subspace Dimension" -- maximum number of basis vectors for subspace.  Two
     *  for standard eigenvalue problem) or three (for generalized eigenvalue problem) sets of basis
     *  vectors of this size will be required. Default: 3*problem->getNEV()*"Block Size"
     * - "Initial Guess" -- how should initial vector be selected: "Random" or "User".
     *   If "User," the value in problem->getInitVec() will be used.  Default: "Random".
     * - "Print Number of Ritz Values" -- an int specifying how many Ritz values should be printed
     *   at each iteration.  Default: "NEV".
     * - "Relative Convergence Tolerance" -- should residual be scaled by corresponding Ritz value
     *   to measure convergence.  Default: "false"
     *
     */
    GeneralizedDavidson(const RCP<Eigenproblem<ScalarType,MV,OP> >  &problem,
                        const RCP<SortManager<MagnitudeType> >      &sortman,
                        const RCP<OutputManager<ScalarType> >       &outputman,
                        const RCP<StatusTest<ScalarType,MV,OP> >    &tester,
                        const RCP<OrthoManager<ScalarType,MV> >     &orthoman,
                        Teuchos::ParameterList                      &pl);

    /*!
     * \brief Solves the eigenvalue problem.
     */
    void iterate();

    /*!
     * \brief Initialize the eigenvalue problem
     *
     * Anything on the state that is not null is assumed to be valid.
     * Anything not present on the state will be generated.
     * Very limited error checking can be performed to ensure the validity of
     * state components (e.g. we cannot verify that <tt> state.AV </tt>actually corresponds
     * to <tt>A*state.V</tt>), so this function should be used carefully.
     */
    void initialize();

    /*!
     * \brief Initialize solver from state
     */
    void initialize( GeneralizedDavidsonState<ScalarType,MV> state );

    /*!
     * \brief Get number of iterations
     */
    int getNumIters() const { return d_iteration; }

    /*!
     * \brief Reset the number of iterations
     */
    void resetNumIters() { d_iteration=0; d_opApplies=0; }

    /*!
     * \brief Get the current Ritz vectors
     */
    RCP<const MV> getRitzVectors()
    {
        if( !d_ritzVectorsValid )
            computeRitzVectors();
        return d_ritzVecs;
    }

    /*!
     * \brief Get the current Ritz values
     */
    std::vector< Value<ScalarType> > getRitzValues();

    /*!
     * \brief Get the current Ritz index vector
     */
    std::vector<int> getRitzIndex()
    {
        if( !d_ritzIndexValid )
            computeRitzIndex();
        return d_ritzIndex;
    }

    /*!
     * \brief Get indices of current block
     *
     * Number of entries is equal to getBlockSize()
     */
    std::vector<int> getBlockIndex() const
    {
        return d_expansionIndices;
    }

    /*!
     * \brief Get the current residual norms (w.r.t. norm defined by OrthoManager)
     */
    std::vector<MagnitudeType> getResNorms();

    /*!
     * \brief Get the current residual norms (w.r.t. norm defined by OrthoManager)
     */
    std::vector<MagnitudeType> getResNorms(int numWanted);

    /*!
     * \brief Get the current residual norms (2-norm)
     */
    std::vector<MagnitudeType> getRes2Norms() { return d_resNorms; }

    /*!
     * \brief Get the current Ritz residual norms (2-norm)
     *
     * GeneralizedDavidson doesn't compute Ritz residual norms
     * so this is equivalent to calling getRes2Norms()
     */
    std::vector<MagnitudeType> getRitzRes2Norms() { return d_resNorms; }

    /*!
     * \brief Get current subspace dimension
     */
    int getCurSubspaceDim() const { return d_curDim; }

    /*!
     * \brief Get maximum subspace dimension
     */
    int getMaxSubspaceDim() const { return d_maxSubspaceDim; }

    /*!
     * \brief Set status test
     */
    void setStatusTest( RCP<StatusTest<ScalarType,MV,OP> > tester) { d_tester = tester; }

    /*!
     * \brief Get status test
     */
    RCP<StatusTest<ScalarType,MV,OP> > getStatusTest() const { return d_tester; }

    /*!
     * \brief Get eigenproblem
     */
    const Eigenproblem<ScalarType,MV,OP> & getProblem() const { return *d_problem; }

    /*!
     * \brief Get block size
     */
    int getBlockSize() const { return d_expansionSize; }

    /*!
     * \brief Set block size
     */
    void setBlockSize(int blockSize);

    /*!
     * \brief Set problem size.
     */
    void setSize(int blockSize, int maxSubDim);

    /*!
     * \brief Get the auxilliary vectors
     */
    Teuchos::Array< RCP<const MV> > getAuxVecs() const { return d_auxVecs; }

    /*!
     * \brief Set auxilliary vectors
     *
     * Manually setting the auxilliary vectors invalidates the current state
     * of the solver.  Reuse of any components of the solver requires extracting
     * the state, orthogonalizing V against the aux vecs and reinitializing.
     */
    void setAuxVecs( const Teuchos::Array< RCP<const MV> > &auxVecs );

    /*!
     * \brief Query if the solver is in an initialized state
     */
    bool isInitialized() const { return d_initialized; }

    /*!
     * \brief Print current status of solver
     */
    void currentStatus( std::ostream &myout );

    /*!
     * \brief Get the current state of the eigensolver.
     */
    GeneralizedDavidsonState<ScalarType,MV> getState();

    /*!
     * Reorder Schur form, bringing wanted values to front
     */
    void sortProblem( int numWanted );

  private:

    // Expand subspace
    void expandSearchSpace();

    // Apply Operators
    void applyOperators();

    // Update projections
    void updateProjections();

    // Solve projected eigenproblem
    void solveProjectedEigenproblem();

    // Compute eigenvectors of matrix pair
    void computeProjectedEigenvectors( Teuchos::SerialDenseMatrix<int,ScalarType> &X );

    // Scale projected eigenvectors by alpha/beta
    void scaleEigenvectors( const Teuchos::SerialDenseMatrix<int,ScalarType> &X,
                                  Teuchos::SerialDenseMatrix<int,ScalarType> &X_alpha,
                                  Teuchos::SerialDenseMatrix<int,ScalarType> &X_beta );

    // Sort vectors of pairs
    void sortValues( std::vector<MagnitudeType> &realParts,
                     std::vector<MagnitudeType> &imagParts,
                     std::vector<int>    &permVec,
                     int N);

    // Compute Residual
    void computeResidual();

    // Update the current Ritz index vector
    void computeRitzIndex();

    // Compute the current Ritz vectors
    void computeRitzVectors();

    // Operators
    RCP<Eigenproblem<ScalarType,MV,OP> > d_problem;
    Teuchos::ParameterList d_pl;
    RCP<const OP> d_A;
    RCP<const OP> d_B;
    RCP<const OP> d_P;
    bool d_haveB;
    bool d_haveP;

    // Parameters
    int d_blockSize;
    int d_maxSubspaceDim;
    int d_NEV;
    int d_numToPrint;
    std::string d_initType;
    int d_verbosity;
    bool d_relativeConvergence;

    // Managers
    RCP<OutputManager<ScalarType> >     d_outputMan;
    RCP<OrthoManager<ScalarType,MV> >   d_orthoMan;
    RCP<SortManager<MagnitudeType> >    d_sortMan;
    RCP<StatusTest<ScalarType,MV,OP> >  d_tester;

    // Eigenvalues
    std::vector< Value<ScalarType> > d_eigenvalues;

    // Residual Vector
    RCP<MV> d_R;
    std::vector<MagnitudeType> d_resNorms;

    // Subspace Vectors
    RCP<MV> d_V;
    RCP<MV> d_AV;
    RCP<MV> d_BV;
    RCP<MV> d_ritzVecSpace;
    RCP<MV> d_ritzVecs;
    Teuchos::Array< RCP<const MV> > d_auxVecs;

    // Serial Matrices
    RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > d_VAV;
    RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > d_VBV;
    RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > d_S;
    RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > d_T;
    RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > d_Q;
    RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > d_Z;

    // Arrays for holding Ritz values
    std::vector<MagnitudeType> d_alphar;
    std::vector<MagnitudeType> d_alphai;
    std::vector<MagnitudeType> d_betar;
    std::vector<int>    d_ritzIndex;
    std::vector<int>    d_convergedIndices;
    std::vector<int>    d_expansionIndices;

    // Current subspace dimension
    int d_curDim;

    // How many vectors are to be added to the subspace
    int d_expansionSize;

    // Should subspace expansion use leading vectors
    //  (if false, will use leading unconverged vectors)
    bool d_useLeading;

    // What should be used for test subspace (V, AV, or BV)
    std::string d_testSpace;

    // How many residual vectors are valid
    int d_residualSize;

    int  d_iteration;
    int  d_opApplies;
    bool d_initialized;
    bool d_ritzIndexValid;
    bool d_ritzVectorsValid;

};

//---------------------------------------------------------------------------//
// Prevent instantiation on complex scalar type
//---------------------------------------------------------------------------//
template <class MagnitudeType, class MV, class OP>
class GeneralizedDavidson<std::complex<MagnitudeType>,MV,OP>
{
  public:

    typedef std::complex<MagnitudeType> ScalarType;
    GeneralizedDavidson(
        const RCP<Eigenproblem<ScalarType,MV,OP> > &problem,
        const RCP<SortManager<MagnitudeType> >     &sortman,
        const RCP<OutputManager<ScalarType> >      &outputman,
        const RCP<StatusTest<ScalarType,MV,OP> >   &tester,
        const RCP<OrthoManager<ScalarType,MV> >    &orthoman,
        Teuchos::ParameterList                     &pl)
    {
        // Provide a compile error when attempting to instantiate on complex type
        MagnitudeType::this_class_is_missing_a_specialization();
    }
};

//---------------------------------------------------------------------------//
// PUBLIC METHODS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Constructor
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
GeneralizedDavidson<ScalarType,MV,OP>::GeneralizedDavidson(
        const RCP<Eigenproblem<ScalarType,MV,OP> > &problem,
        const RCP<SortManager<MagnitudeType> >     &sortman,
        const RCP<OutputManager<ScalarType> >      &outputman,
        const RCP<StatusTest<ScalarType,MV,OP> >   &tester,
        const RCP<OrthoManager<ScalarType,MV> >    &orthoman,
        Teuchos::ParameterList                     &pl )
{
    TEUCHOS_TEST_FOR_EXCEPTION(   problem == Teuchos::null, std::invalid_argument, "No Eigenproblem given to solver." );
    TEUCHOS_TEST_FOR_EXCEPTION( outputman == Teuchos::null, std::invalid_argument, "No OutputManager given to solver." );
    TEUCHOS_TEST_FOR_EXCEPTION(  orthoman == Teuchos::null, std::invalid_argument, "No OrthoManager given to solver." );
    TEUCHOS_TEST_FOR_EXCEPTION(   sortman == Teuchos::null, std::invalid_argument, "No SortManager given to solver." );
    TEUCHOS_TEST_FOR_EXCEPTION(    tester == Teuchos::null, std::invalid_argument, "No StatusTest given to solver." );
    TEUCHOS_TEST_FOR_EXCEPTION(   !problem->isProblemSet(), std::invalid_argument, "Problem has not been set." );

    d_problem = problem;
    d_pl = pl;
    TEUCHOS_TEST_FOR_EXCEPTION( problem->getA()==Teuchos::null && problem->getOperator()==Teuchos::null,
                                std::invalid_argument, "Either A or Operator must be non-null on Eigenproblem");
    d_A = problem->getA()!=Teuchos::null ? problem->getA() : problem->getOperator();
    d_B = problem->getM();
    d_P = problem->getPrec();
    d_sortMan = sortman;
    d_outputMan = outputman;
    d_tester = tester;
    d_orthoMan = orthoman;

    // Pull entries from the ParameterList and Eigenproblem
    d_NEV        = d_problem->getNEV();
    d_initType   = d_pl.get<std::string>("Initial Guess","Random");
    d_numToPrint = d_pl.get<int>("Print Number of Ritz Values",-1);
    d_useLeading = d_pl.get<bool>("Use Leading Vectors",false);

    if( d_B != Teuchos::null )
        d_haveB = true;
    else
        d_haveB = false;

    if( d_P != Teuchos::null )
        d_haveP = true;
    else
        d_haveP = false;

    d_testSpace = d_pl.get<std::string>("Test Space","V");
    TEUCHOS_TEST_FOR_EXCEPTION( d_testSpace!="V" && d_testSpace!="AV" && d_testSpace!="BV", std::invalid_argument,
        "Anasazi::GeneralizedDavidson: Test Space must be V, AV, or BV" );
    TEUCHOS_TEST_FOR_EXCEPTION( d_testSpace=="V" ? false : !d_haveB, std::invalid_argument,
        "Anasazi::GeneralizedDavidson: Test Space must be V for standard eigenvalue problem" );

    // Allocate space for subspace vectors, projected matrices
    int blockSize  = d_pl.get<int>("Block Size",1);
    int maxSubDim  = d_pl.get<int>("Maximum Subspace Dimension",3*d_NEV*blockSize);
    d_blockSize      = -1;
    d_maxSubspaceDim = -1;
    setSize( blockSize, maxSubDim );
    d_relativeConvergence = d_pl.get<bool>("Relative Convergence Tolerance",false);

    // Make sure subspace size is consistent with requested eigenvalues
    TEUCHOS_TEST_FOR_EXCEPTION( d_blockSize <= 0, std::invalid_argument, "Block size must be positive");
    TEUCHOS_TEST_FOR_EXCEPTION( d_maxSubspaceDim <= 0, std::invalid_argument, "Maximum Subspace Dimension must be positive" );
    TEUCHOS_TEST_FOR_EXCEPTION( d_problem->getNEV()+2 > pl.get<int>("Maximum Subspace Dimension"),
                                std::invalid_argument, "Maximum Subspace Dimension must be strictly greater than NEV");
    TEUCHOS_TEST_FOR_EXCEPTION( d_maxSubspaceDim > MVT::GetVecLength(*problem->getInitVec()),
                                std::invalid_argument, "Maximum Subspace Dimension cannot exceed problem size");


    d_curDim = 0;
    d_iteration = 0;
    d_opApplies = 0;
    d_ritzIndexValid = false;
    d_ritzVectorsValid = false;
}


//---------------------------------------------------------------------------//
// Iterate
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
void GeneralizedDavidson<ScalarType,MV,OP>::iterate()
{
    // Initialize Problem
    if( !d_initialized )
    {
        d_outputMan->stream(Warnings) << "WARNING: GeneralizedDavidson::iterate called without first calling initialize" << std::endl;
        d_outputMan->stream(Warnings) << "         Default initialization will be performed" << std::endl;
        initialize();
    }

    // Print current status
    if( d_outputMan->isVerbosity(Debug) )
    {
        currentStatus( d_outputMan->stream(Debug) );
    }
    else if( d_outputMan->isVerbosity(IterationDetails) )
    {
        currentStatus( d_outputMan->stream(IterationDetails) );
    }

    while( d_tester->getStatus() != Passed && d_curDim+d_expansionSize <= d_maxSubspaceDim )
    {
        d_iteration++;

        expandSearchSpace();

        applyOperators();

        updateProjections();

        solveProjectedEigenproblem();

        // Make sure the most significant Ritz values are in front
        // We want the greater of the block size and the number of
        //  requested values, but can't exceed the current dimension
        int numToSort = std::max(d_blockSize,d_NEV);
        numToSort = std::min(numToSort,d_curDim);
        sortProblem( numToSort );

        computeResidual();

        // Print current status
        if( d_outputMan->isVerbosity(Debug) )
        {
            currentStatus( d_outputMan->stream(Debug) );
        }
        else if( d_outputMan->isVerbosity(IterationDetails) )
        {
            currentStatus( d_outputMan->stream(IterationDetails) );
        }
    }
}

//---------------------------------------------------------------------------//
// Return the current state struct
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
GeneralizedDavidsonState<ScalarType,MV> GeneralizedDavidson<ScalarType,MV,OP>::getState()
{
    GeneralizedDavidsonState<ScalarType,MV> state;
    state.curDim = d_curDim;
    state.V      = d_V;
    state.AV     = d_AV;
    state.BV     = d_BV;
    state.VAV    = d_VAV;
    state.VBV    = d_VBV;
    state.S      = d_S;
    state.T      = d_T;
    state.Q      = d_Q;
    state.Z      = d_Z;
    state.eVals  = getRitzValues();
    return state;
}

//---------------------------------------------------------------------------//
// Set block size
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
void GeneralizedDavidson<ScalarType,MV,OP>::setBlockSize(int blockSize)
{
    setSize(blockSize,d_maxSubspaceDim);
}

//---------------------------------------------------------------------------//
// Set block size and maximum subspace dimension.
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
void GeneralizedDavidson<ScalarType,MV,OP>::setSize(int blockSize, int maxSubDim )
{
    if( blockSize != d_blockSize || maxSubDim != d_maxSubspaceDim )
    {
        d_blockSize = blockSize;
        d_maxSubspaceDim = maxSubDim;
        d_initialized = false;

        d_outputMan->stream(Debug) << " >> Anasazi::GeneralizedDavidson: Allocating eigenproblem"
            << " state with block size of " << d_blockSize
            << " and maximum subspace dimension of " << d_maxSubspaceDim << std::endl;

        // Resize arrays for Ritz values
        d_alphar.resize(d_maxSubspaceDim);
        d_alphai.resize(d_maxSubspaceDim);
        d_betar.resize(d_maxSubspaceDim);

        // Shorten for convenience here
        int msd = d_maxSubspaceDim;

        // Temporarily save initialization vector to clone needed vectors
        RCP<const MV> initVec = d_problem->getInitVec();

        // Allocate subspace vectors
        d_V            = MVT::Clone(*initVec, msd);
        d_AV           = MVT::Clone(*initVec, msd);

        // Allocate serial matrices
        d_VAV = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(msd,msd) );
        d_S   = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(msd,msd) );
        d_Q   = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(msd,msd) );
        d_Z   = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(msd,msd) );

        // If this is generalized eigenproblem, allocate B components
        if( d_haveB )
        {
            d_BV  = MVT::Clone(*initVec, msd);
            d_VBV = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(msd,msd) );
            d_T   = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(msd,msd) );
        }

        /* Allocate space for residual and Ritz vectors
         * The residual serves two purposes in the Davidson algorithm --
         *  subspace expansion (via the preconditioner) and convergence checking.
         * We need "Block Size" vectors for subspace expantion and NEV vectors
         *  for convergence checking.  Allocate space for max of these, one
         *  extra to avoid splitting conjugate pairs
         * Allocate one more than "Block Size" to avoid splitting a conjugate pair
         */
        d_R = MVT::Clone(*initVec,std::max(d_blockSize,d_NEV)+1);
        d_ritzVecSpace = MVT::Clone(*initVec,std::max(d_blockSize,d_NEV)+1);
    }
}

//---------------------------------------------------------------------------//
/*
 * Initialize the eigenvalue problem
 *
 * Anything on the state that is not null is assumed to be valid.
 * Anything not present on the state will be generated
 * Very limited error checking can be performed to ensure the validity of
 * state components (e.g. we cannot verify that state.AV actually corresponds
 * to A*state.V), so this function should be used carefully.
 */
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
void GeneralizedDavidson<ScalarType,MV,OP>::initialize( GeneralizedDavidsonState<ScalarType,MV> state )
{
    // If state has nonzero dimension, we initialize from that, otherwise
    //  we'll pick d_blockSize vectors to start with
    d_curDim = (state.curDim > 0 ? state.curDim : d_blockSize );

    d_outputMan->stream(Debug) << " >> Anasazi::GeneralizedDavidson: Initializing state"
        << " with subspace dimension " << d_curDim << std::endl;

    // Index for 1st block_size vectors
    std::vector<int> initInds(d_curDim);
    for( int i=0; i<d_curDim; ++i )
        initInds[i] = i;

    // View of vectors that need to be initialized
    RCP<MV>  V1 = MVT::CloneViewNonConst(*d_V,initInds);

    // If state's dimension is large enough, use state.V to initialize
    bool reset_V = false;;
    if( state.curDim > 0 && state.V != Teuchos::null && MVT::GetNumberVecs(*state.V) >= d_curDim )
    {
        if( state.V != d_V )
            MVT::SetBlock(*state.V,initInds,*V1);
    }
    // If there aren't enough vectors in problem->getInitVec() or the user specifically
    //  wants to use random data, set V to random
    else if( MVT::GetNumberVecs(*d_problem->getInitVec()) < d_blockSize || d_initType == "Random" )
    {
        MVT::MvRandom(*V1);
        reset_V = true;
    }
    // Use vectors in problem->getInitVec()
    else
    {
        RCP<const MV> initVec = MVT::CloneView(*d_problem->getInitVec(),initInds);
        MVT::SetBlock(*initVec,initInds,*V1);
        reset_V = true;
    }

    // If we reset V, it needs to be orthonormalized
    if( reset_V )
    {
        int rank = d_orthoMan->projectAndNormalize( *V1, d_auxVecs );
        TEUCHOS_TEST_FOR_EXCEPTION( rank < d_blockSize, std::logic_error,
            "Anasazi::GeneralizedDavidson::initialize(): Error generating initial orthonormal basis" );
    }

    if( d_outputMan->isVerbosity(Debug) )
    {
        d_outputMan->stream(Debug) << " >> Anasazi::GeneralizedDavidson: Error in V^T V == I: "
            << d_orthoMan->orthonormError( *V1 ) << std::endl;
    }

    // Now process AV
    RCP<MV> AV1 = MVT::CloneViewNonConst(*d_AV,initInds);

    // If AV in the state is valid and of appropriate size, use it
    // We have no way to check that AV is actually A*V
    if( !reset_V && state.AV != Teuchos::null && MVT::GetNumberVecs(*state.AV) >= d_curDim )
    {
        if( state.AV != d_AV )
            MVT::SetBlock(*state.AV,initInds,*AV1);
    }
    // Otherwise apply A to V
    else
    {
        OPT::Apply( *d_A, *V1, *AV1 );
        d_opApplies += MVT::GetNumberVecs( *V1 );
    }

    // Views of matrix to be updated
    Teuchos::SerialDenseMatrix<int,ScalarType> VAV1( Teuchos::View, *d_VAV, d_curDim, d_curDim );

    // If the state has a valid VAV, use it
    if( !reset_V && state.VAV != Teuchos::null && state.VAV->numRows() >= d_curDim && state.VAV->numCols() >= d_curDim )
    {
        if( state.VAV != d_VAV )
        {
            Teuchos::SerialDenseMatrix<int,ScalarType> state_VAV( Teuchos::View, *state.VAV, d_curDim, d_curDim );
            VAV1.assign( state_VAV );
        }
    }
    // Otherwise compute VAV from V,AV
    else
    {
        if( d_testSpace == "V" )
        {
            MVT::MvTransMv( ST::one(),  *V1, *AV1, VAV1 );
        }
        else if( d_testSpace == "AV" )
        {
            MVT::MvTransMv( ST::one(), *AV1, *AV1, VAV1 );
        }
        else if( d_testSpace == "BV" )
        {
            RCP<MV> BV1 = MVT::CloneViewNonConst(*d_BV,initInds);
            MVT::MvTransMv( ST::one(), *BV1, *AV1, VAV1 );
        }
    }

    // Process BV if we have it
    if( d_haveB )
    {
        RCP<MV> BV1 = MVT::CloneViewNonConst(*d_BV,initInds);

        // If BV in the state is valid and of appropriate size, use it
        // We have no way to check that BV is actually B*V
        if( !reset_V && state.BV != Teuchos::null && MVT::GetNumberVecs(*state.BV) >= d_curDim )
        {
            if( state.BV != d_BV )
                MVT::SetBlock(*state.BV,initInds,*BV1);
        }
        // Otherwise apply B to V
        else
        {
            OPT::Apply( *d_B, *V1, *BV1 );
        }

        // Views of matrix to be updated
        Teuchos::SerialDenseMatrix<int,ScalarType> VBV1( Teuchos::View, *d_VBV, d_curDim, d_curDim );

        // If the state has a valid VBV, use it
        if( !reset_V && state.VBV != Teuchos::null && state.VBV->numRows() >= d_curDim && state.VBV->numCols() >= d_curDim )
        {
            if( state.VBV != d_VBV )
            {
                Teuchos::SerialDenseMatrix<int,ScalarType> state_VBV( Teuchos::View, *state.VBV, d_curDim, d_curDim );
                VBV1.assign( state_VBV );
            }
        }
        // Otherwise compute VBV from V,BV
        else
        {
            if( d_testSpace == "V" )
            {
                MVT::MvTransMv( ST::one(),  *V1, *BV1, VBV1 );
            }
            else if( d_testSpace == "AV" )
            {
                MVT::MvTransMv( ST::one(), *AV1, *BV1, VBV1 );
            }
            else if( d_testSpace == "BV" )
            {
                MVT::MvTransMv( ST::one(), *BV1, *BV1, VBV1 );
            }
        }
    }

    // Update Ritz values
    solveProjectedEigenproblem();

    // Sort
    int numToSort = std::max(d_blockSize,d_NEV);
    numToSort = std::min(numToSort,d_curDim);
    sortProblem( numToSort );

    // Get valid residual
    computeResidual();

    // Set solver to initialized
    d_initialized = true;
}

//---------------------------------------------------------------------------//
// Initialize the eigenvalue problem with empty state
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
void GeneralizedDavidson<ScalarType,MV,OP>::initialize()
{
    GeneralizedDavidsonState<ScalarType,MV> empty;
    initialize( empty );
}

//---------------------------------------------------------------------------//
// Get current residual norms
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>
    GeneralizedDavidson<ScalarType,MV,OP>::getResNorms()
{
    return getResNorms(d_residualSize);
}

//---------------------------------------------------------------------------//
// Get current residual norms
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>
    GeneralizedDavidson<ScalarType,MV,OP>::getResNorms(int numWanted)
{
    std::vector<int> resIndices(numWanted);
    for( int i=0; i<numWanted; ++i )
        resIndices[i]=i;

    RCP<const MV> R_checked = MVT::CloneView( *d_R, resIndices );

    std::vector<MagnitudeType> resNorms;
    d_orthoMan->norm( *R_checked, resNorms );

    return resNorms;
}

//---------------------------------------------------------------------------//
// Get current Ritz values
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
std::vector< Value<ScalarType> > GeneralizedDavidson<ScalarType,MV,OP>::getRitzValues()
{
    std::vector< Value<ScalarType> > ritzValues;
    for( int ival=0; ival<d_curDim; ++ival )
    {
        Value<ScalarType> thisVal;
        thisVal.realpart = d_alphar[ival] / d_betar[ival];
        if( d_betar[ival] != MT::zero() )
            thisVal.imagpart = d_alphai[ival] / d_betar[ival];
        else
            thisVal.imagpart = MT::zero();

        ritzValues.push_back( thisVal );
    }

    return ritzValues;
}

//---------------------------------------------------------------------------//
/*
 * Set auxilliary vectors
 *
 * Manually setting the auxilliary vectors invalidates the current state
 * of the solver.  Reuse of any components of the solver requires extracting
 * the state, orthogonalizing V against the aux vecs and reinitializing.
 */
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
void GeneralizedDavidson<ScalarType,MV,OP>::setAuxVecs(
        const Teuchos::Array< RCP<const MV> > &auxVecs )
{
    d_auxVecs = auxVecs;

    // Set state to uninitialized if any vectors were set here
    typename Teuchos::Array< RCP<const MV> >::const_iterator arrItr;
    int numAuxVecs=0;
    for( arrItr=auxVecs.begin(); arrItr!=auxVecs.end(); ++arrItr )
    {
        numAuxVecs += MVT::GetNumberVecs( *(*arrItr) );
    }
    if( numAuxVecs > 0 )
        d_initialized = false;
}

//---------------------------------------------------------------------------//
// Reorder Schur form, bringing wanted values to front
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
void GeneralizedDavidson<ScalarType,MV,OP>::sortProblem( int numWanted )
{
    // Get permutation vector
    std::vector<MagnitudeType> realRitz(d_curDim), imagRitz(d_curDim);
    std::vector< Value<ScalarType> > ritzVals = getRitzValues();
    for( int i=0; i<d_curDim; ++i )
    {
        realRitz[i] = ritzVals[i].realpart;
        imagRitz[i] = ritzVals[i].imagpart;
    }

    std::vector<int> permVec;
    sortValues( realRitz, imagRitz, permVec, d_curDim );

    std::vector<int> sel(d_curDim,0);
    for( int ii=0; ii<numWanted; ++ii )
        sel[ permVec[ii] ]=1;

    if( d_haveB )
    {
        int ijob  = 0; // reorder only, no condition number estimates
        int wantq = 1; // keep left Schur vectors
        int wantz = 1; // keep right Schur vectors
        int work_size=10*d_maxSubspaceDim+16;
        std::vector<ScalarType> work(work_size);
        int sdim   = 0;
        int iwork_size = 1;
        int iwork;
        int info   = 0;

        Teuchos::LAPACK<int,ScalarType> lapack;
        lapack.TGSEN( ijob, wantq, wantz, &sel[0], d_curDim, d_S->values(), d_S->stride(), d_T->values(), d_T->stride(),
                      &d_alphar[0], &d_alphai[0], &d_betar[0], d_Q->values(), d_Q->stride(), d_Z->values(), d_Z->stride(),
                      &sdim, 0, 0, 0, &work[0], work_size, &iwork, iwork_size, &info );

        d_ritzIndexValid   = false;
        d_ritzVectorsValid = false;

        std::stringstream ss;
        ss << "Anasazi::GeneralizedDavidson: TGSEN returned error code " << info << std::endl;
        TEUCHOS_TEST_FOR_EXCEPTION( info<0, std::runtime_error, ss.str() );
        if( info > 0 )
        {
            // Only issue a warning for positive error code, this usually indicates
            //  that the system has not been fully reordered, presumably due to ill-conditioning.
            // This is usually not detrimental to the calculation.
            d_outputMan->stream(Warnings) << "WARNING: " << ss.str() << std::endl;
            d_outputMan->stream(Warnings) << "  Problem may not be correctly sorted" << std::endl;
        }
    }
    else
    {
        char getCondNum = 'N'; // no condition number estimates
        char getQ = 'V';       // keep Schur vectors
        int subDim;
        int work_size = d_curDim;
        std::vector<ScalarType> work(work_size);
        int iwork_size = 1;
        int iwork;
        int info;

        Teuchos::LAPACK<int,ScalarType> lapack;
        lapack.TRSEN( getCondNum, getQ, &sel[0], d_curDim, d_S->values(), d_S->stride(), d_Z->values(), d_Z->stride(),
                      &d_alphar[0], &d_alphai[0], &subDim, 0, 0, &work[0], work_size, &iwork, iwork_size, &info );


        std::fill( d_betar.begin(), d_betar.end(), ST::one() );

        d_ritzIndexValid = false;
        d_ritzVectorsValid = false;

        std::stringstream ss;
        ss << "Anasazi::GeneralizedDavidson: TRSEN returned error code " << info << std::endl;
        TEUCHOS_TEST_FOR_EXCEPTION( info!=0, std::runtime_error, ss.str() );
    }
}


//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Expand subspace using preconditioner and orthogonalize
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
void GeneralizedDavidson<ScalarType,MV,OP>::expandSearchSpace()
{
    // Get indices into relevant portion of residual and
    //  location to be added to search space
    std::vector<int> newIndices(d_expansionSize);
    for( int i=0; i<d_expansionSize; ++i )
    {
        newIndices[i] = d_curDim+i;
    }

    // Get indices into pre-existing search space
    std::vector<int> curIndices(d_curDim);
    for( int i=0; i<d_curDim; ++i )
        curIndices[i] = i;

    // Get View of vectors
    RCP<MV>       V_new    = MVT::CloneViewNonConst( *d_V, newIndices);
    RCP<const MV> V_cur    = MVT::CloneView(         *d_V, curIndices);
    RCP<const MV> R_active = MVT::CloneView(         *d_R, d_expansionIndices);

    if( d_haveP )
    {
        // Apply Preconditioner to Residual
        OPT::Apply( *d_P, *R_active, *V_new );
    }
    else
    {
        // Just copy the residual
        MVT::SetBlock( *R_active, newIndices, *d_V );
    }

    // Normalize new vector against existing vectors in V plus auxVecs
    Teuchos::Array< RCP<const MV> > against = d_auxVecs;
    against.push_back( V_cur );
    int rank = d_orthoMan->projectAndNormalize(*V_new,against);

    if( d_outputMan->isVerbosity(Debug) )
    {
        std::vector<int> allIndices(d_curDim+d_expansionSize);
        for( int i=0; i<d_curDim+d_expansionSize; ++i )
            allIndices[i]=i;

        RCP<const MV> V_all = MVT::CloneView( *d_V, allIndices );

        d_outputMan->stream(Debug) << " >> Anasazi::GeneralizedDavidson: Error in V^T V == I: "
            << d_orthoMan->orthonormError( *V_all ) << std::endl;
    }

    TEUCHOS_TEST_FOR_EXCEPTION( rank != d_expansionSize, std::runtime_error,
        "Anasazi::GeneralizedDavidson::ExpandSearchSpace(): Orthonormalization of new vectors failed" );

}

//---------------------------------------------------------------------------//
// Apply operators
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
void GeneralizedDavidson<ScalarType,MV,OP>::applyOperators()
{
    // Get indices for different components
    std::vector<int> newIndices(d_expansionSize);
    for( int i=0; i<d_expansionSize; ++i )
        newIndices[i] = d_curDim+i;

    // Get Views
    RCP<const MV>  V_new = MVT::CloneView(         *d_V,  newIndices );
    RCP<MV>       AV_new = MVT::CloneViewNonConst( *d_AV, newIndices );

    // Multiply by A
    OPT::Apply( *d_A, *V_new, *AV_new );
    d_opApplies += MVT::GetNumberVecs( *V_new );

    // Multiply by B
    if( d_haveB )
    {
        RCP<MV>       BV_new = MVT::CloneViewNonConst( *d_BV, newIndices );
        OPT::Apply( *d_B, *V_new, *BV_new );
    }
}

//---------------------------------------------------------------------------//
// Update projected matrices.
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
void GeneralizedDavidson<ScalarType,MV,OP>::updateProjections()
{
    // Get indices for different components
    std::vector<int> newIndices(d_expansionSize);
    for( int i=0; i<d_expansionSize; ++i )
        newIndices[i] = d_curDim+i;

    std::vector<int> curIndices(d_curDim);
    for( int i=0; i<d_curDim; ++i )
        curIndices[i] = i;

    std::vector<int> allIndices(d_curDim+d_expansionSize);
    for( int i=0; i<d_curDim+d_expansionSize; ++i )
        allIndices[i] = i;

    // Test subspace can be V, AV, or BV
    RCP<const MV> W_new, W_all;
    if( d_testSpace == "V" )
    {
        W_new = MVT::CloneView(*d_V, newIndices );
        W_all = MVT::CloneView(*d_V, allIndices );
    }
    else if( d_testSpace == "AV" )
    {
        W_new = MVT::CloneView(*d_AV, newIndices );
        W_all = MVT::CloneView(*d_AV, allIndices );
    }
    else if( d_testSpace == "BV" )
    {
        W_new = MVT::CloneView(*d_BV, newIndices );
        W_all = MVT::CloneView(*d_BV, allIndices );
    }

    // Get views of AV
    RCP<const MV>     AV_new = MVT::CloneView(*d_AV, newIndices);
    RCP<const MV> AV_current = MVT::CloneView(*d_AV, curIndices);

    // Last block_size rows of VAV (minus final entry)
    Teuchos::SerialDenseMatrix<int,ScalarType> VAV_lastrow( Teuchos::View, *d_VAV, d_expansionSize, d_curDim, d_curDim, 0 );
    MVT::MvTransMv( ST::one(), *W_new, *AV_current, VAV_lastrow );

    // Last block_size columns of VAV
    Teuchos::SerialDenseMatrix<int,ScalarType> VAV_lastcol( Teuchos::View, *d_VAV, d_curDim+d_expansionSize, d_expansionSize, 0, d_curDim );
    MVT::MvTransMv( ST::one(), *W_all, *AV_new, VAV_lastcol );

    if( d_haveB )
    {
        // Get views of BV
        RCP<const MV>     BV_new = MVT::CloneView(*d_BV, newIndices);
        RCP<const MV> BV_current = MVT::CloneView(*d_BV, curIndices);

        // Last block_size rows of VBV (minus final entry)
        Teuchos::SerialDenseMatrix<int,ScalarType> VBV_lastrow( Teuchos::View, *d_VBV, d_expansionSize, d_curDim, d_curDim, 0 );
        MVT::MvTransMv( ST::one(), *W_new, *BV_current, VBV_lastrow );

        // Last block_size columns of VBV
        Teuchos::SerialDenseMatrix<int,ScalarType> VBV_lastcol( Teuchos::View, *d_VBV, d_curDim+d_expansionSize, d_expansionSize, 0, d_curDim );
        MVT::MvTransMv( ST::one(), *W_all, *BV_new, VBV_lastcol );
    }

    // All bases are expanded, increase current subspace dimension
    d_curDim += d_expansionSize;

    d_ritzIndexValid   = false;
    d_ritzVectorsValid = false;
}

//---------------------------------------------------------------------------//
// Solve low dimensional eigenproblem using LAPACK
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
void GeneralizedDavidson<ScalarType,MV,OP>::solveProjectedEigenproblem()
{
    if( d_haveB )
    {
        // VAV and VBV need to stay unchanged, GGES will overwrite
        //  S and T with the triangular matrices from the generalized
        //  Schur form
        d_S->assign(*d_VAV);
        d_T->assign(*d_VBV);

        // Get QZ Decomposition of Projected Problem
        char leftVecs  = 'V'; // compute left vectors
        char rightVecs = 'V'; // compute right vectors
        char sortVals  = 'N'; // don't sort
        int sdim;
        // int work_size = 10*d_curDim+16;
        Teuchos::LAPACK<int,ScalarType> lapack;
        int info;
        // workspace query
        int work_size = -1;
        std::vector<ScalarType> work(1);
        lapack.GGES( leftVecs, rightVecs, sortVals, NULL, d_curDim, d_S->values(), d_S->stride(),
                     d_T->values(), d_T->stride(), &sdim, &d_alphar[0], &d_alphai[0], &d_betar[0],
                     d_Q->values(), d_Q->stride(), d_Z->values(), d_Z->stride(), &work[0], work_size, 0, &info );
        // actual call
        work_size = work[0];
        work.resize(work_size);
        lapack.GGES( leftVecs, rightVecs, sortVals, NULL, d_curDim, d_S->values(), d_S->stride(),
                     d_T->values(), d_T->stride(), &sdim, &d_alphar[0], &d_alphai[0], &d_betar[0],
                     d_Q->values(), d_Q->stride(), d_Z->values(), d_Z->stride(), &work[0], work_size, 0, &info );

        d_ritzIndexValid   = false;
        d_ritzVectorsValid = false;

        std::stringstream ss;
        ss << "Anasazi::GeneralizedDavidson: GGES returned error code " << info << std::endl;
        TEUCHOS_TEST_FOR_EXCEPTION( info!=0, std::runtime_error, ss.str() );
    }
    else
    {
        // VAV needs to stay unchanged, GGES will overwrite
        //  S with the triangular matrix from the Schur form
        d_S->assign(*d_VAV);

        // Get QR Decomposition of Projected Problem
        char vecs = 'V'; // compute Schur vectors
        int sdim;
        int work_size = 3*d_curDim;
        std::vector<ScalarType>  work(work_size);
        int info;

        Teuchos::LAPACK<int,ScalarType> lapack;
        lapack.GEES( vecs, d_curDim, d_S->values(), d_S->stride(), &sdim, &d_alphar[0], &d_alphai[0],
                     d_Z->values(), d_Z->stride(), &work[0], work_size, 0, 0, &info);

        std::fill( d_betar.begin(), d_betar.end(), ST::one() );

        d_ritzIndexValid   = false;
        d_ritzVectorsValid = false;

        std::stringstream ss;
        ss << "Anasazi::GeneralizedDavidson: GEES returned error code " << info << std::endl;
        TEUCHOS_TEST_FOR_EXCEPTION( info!=0, std::runtime_error, ss.str() );
    }
}

//---------------------------------------------------------------------------//
/*
 * Get index vector into current Ritz values/vectors
 *
 * The current ordering of d_alphar, d_alphai, d_betar will be used.
 * Reordering those vectors will invalidate the vector returned here.
 */
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
void GeneralizedDavidson<ScalarType,MV,OP>::computeRitzIndex()
{
    if( d_ritzIndexValid )
        return;

    d_ritzIndex.resize( d_curDim );
    int i=0;
    while( i < d_curDim )
    {
        if( d_alphai[i] == ST::zero() )
        {
            d_ritzIndex[i] = 0;
            i++;
        }
        else
        {
            d_ritzIndex[i]   =  1;
            d_ritzIndex[i+1] = -1;
            i+=2;
        }
    }
    d_ritzIndexValid = true;
}

//---------------------------------------------------------------------------//
/*
 * Compute current Ritz vectors
 *
 * The current ordering of d_alphar, d_alphai, d_betar will be used.
 * Reordering those vectors will invalidate the vector returned here.
 */
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
void GeneralizedDavidson<ScalarType,MV,OP>::computeRitzVectors()
{
    if( d_ritzVectorsValid )
        return;

    // Make Ritz indices current
    computeRitzIndex();

    // Get indices of converged vector
    std::vector<int> checkedIndices(d_residualSize);
    for( int ii=0; ii<d_residualSize; ++ii )
        checkedIndices[ii] = ii;

    // Get eigenvectors of projected system
    Teuchos::SerialDenseMatrix<int,ScalarType> X(Teuchos::Copy,*d_Z,d_curDim,d_curDim);
    computeProjectedEigenvectors( X );

    // Get view of wanted vectors
    Teuchos::SerialDenseMatrix<int,ScalarType> X_wanted(Teuchos::View,X,d_curDim,d_residualSize);

    // Get views of relevant portion of V, evecs
    d_ritzVecs = MVT::CloneViewNonConst( *d_ritzVecSpace, checkedIndices );

    std::vector<int> curIndices(d_curDim);
    for( int i=0; i<d_curDim; ++i )
        curIndices[i] = i;

    RCP<const MV> V_current = MVT::CloneView( *d_V, curIndices );

    // Now form Ritz vector
    MVT::MvTimesMatAddMv(ST::one(),*V_current,X_wanted,ST::zero(),*d_ritzVecs);

    // Normalize vectors, conjugate pairs get normalized together
    std::vector<MagnitudeType> scale(d_residualSize);
    MVT::MvNorm( *d_ritzVecs, scale );
    Teuchos::LAPACK<int,ScalarType> lapack;
    for( int i=0; i<d_residualSize; ++i )
    {
        if( d_ritzIndex[i] == 0 )
        {
            scale[i] = 1.0/scale[i];
        }
        else if( d_ritzIndex[i] == 1 )
        {
            MagnitudeType nrm = lapack.LAPY2(scale[i],scale[i+1]);
            scale[i]   = 1.0/nrm;
            scale[i+1] = 1.0/nrm;
        }
    }
    MVT::MvScale( *d_ritzVecs, scale );

    d_ritzVectorsValid = true;

}

//---------------------------------------------------------------------------//
// Use sort manager to sort generalized eigenvalues
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
void GeneralizedDavidson<ScalarType,MV,OP>::sortValues( std::vector<MagnitudeType> &realParts,
                                             std::vector<MagnitudeType> &imagParts,
                                             std::vector<int>    &permVec,
                                             int N)
{
    permVec.resize(N);

    TEUCHOS_TEST_FOR_EXCEPTION( (int) realParts.size()<N, std::runtime_error,
        "Anasazi::GeneralizedDavidson::SortValues: Number of requested sorted values greater than vector length." );
    TEUCHOS_TEST_FOR_EXCEPTION( (int) imagParts.size()<N, std::runtime_error,
        "Anasazi::GeneralizedDavidson::SortValues: Number of requested sorted values greater than vector length." );

    RCP< std::vector<int> > rcpPermVec = Teuchos::rcpFromRef(permVec);

    d_sortMan->sort( realParts, imagParts, rcpPermVec, N );

    d_ritzIndexValid = false;
    d_ritzVectorsValid = false;
}

//---------------------------------------------------------------------------//
/*
 * Compute (right) scaled eigenvectors of a pair of dense matrices
 *
 * This routine computes the eigenvectors for the generalized eigenvalue
 * problem \f$ \beta A x = \alpha B x \f$.  The input matrices are the upper
 * quasi-triangular matrices S and T from a real QZ decomposition,
 * the routine dtgevc will back-calculate the eigenvectors of the original
 * pencil (A,B) using the orthogonal matrices Q and Z.
 */
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
void GeneralizedDavidson<ScalarType,MV,OP>::computeProjectedEigenvectors(
        Teuchos::SerialDenseMatrix<int,ScalarType> &X )
{
    int N = X.numRows();
    if( d_haveB )
    {
        Teuchos::SerialDenseMatrix<int,ScalarType>  S(Teuchos::Copy, *d_S, N, N);
        Teuchos::SerialDenseMatrix<int,ScalarType>  T(Teuchos::Copy, *d_T, N, N);
        Teuchos::SerialDenseMatrix<int,ScalarType> VL(Teuchos::Copy, *d_Q, N, N);

        char whichVecs = 'R'; // only need right eigenvectors
        char howMany   = 'B'; // back-compute eigenvectors of original A,B (we have Schur decomposition here)
        int work_size = 6*d_maxSubspaceDim;
        std::vector<ScalarType> work(work_size,ST::zero());
        int info;
        int M;

        Teuchos::LAPACK<int,ScalarType> lapack;
        lapack.TGEVC( whichVecs, howMany, 0, N, S.values(), S.stride(), T.values(), T.stride(),
                      VL.values(), VL.stride(), X.values(), X.stride(), N, &M, &work[0], &info );

        std::stringstream ss;
        ss << "Anasazi::GeneralizedDavidson: TGEVC returned error code " << info << std::endl;
        TEUCHOS_TEST_FOR_EXCEPTION( info!=0, std::runtime_error, ss.str() );
    }
    else
    {
        Teuchos::SerialDenseMatrix<int,ScalarType>  S(Teuchos::Copy, *d_S, N, N);
        Teuchos::SerialDenseMatrix<int,ScalarType> VL(Teuchos::Copy, *d_Z, N, N);

        char whichVecs = 'R'; // only need right eigenvectors
        char howMany   = 'B'; // back-compute eigenvectors of original A (we have Schur decomposition here)
        int sel = 0;
        std::vector<ScalarType> work(3*N);
        int m;
        int info;

        Teuchos::LAPACK<int,ScalarType> lapack;

        lapack.TREVC( whichVecs, howMany, &sel, N, S.values(), S.stride(), VL.values(), VL.stride(),
                      X.values(), X.stride(), N, &m, &work[0], &info );

        std::stringstream ss;
        ss << "Anasazi::GeneralizedDavidson: TREVC returned error code " << info << std::endl;
        TEUCHOS_TEST_FOR_EXCEPTION( info!=0, std::runtime_error, ss.str() );
    }
}

//---------------------------------------------------------------------------//
// Scale eigenvectors by quasi-diagonal matrices alpha and beta
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
void GeneralizedDavidson<ScalarType,MV,OP>::scaleEigenvectors(
        const Teuchos::SerialDenseMatrix<int,ScalarType> &X,
              Teuchos::SerialDenseMatrix<int,ScalarType> &X_alpha,
              Teuchos::SerialDenseMatrix<int,ScalarType> &X_beta )
{
    int Nr = X.numRows();
    int Nc = X.numCols();

    TEUCHOS_TEST_FOR_EXCEPTION( Nr>d_curDim, std::logic_error,
        "Anasazi::GeneralizedDavidson::ScaleEigenvectors: Matrix size exceeds current dimension");
    TEUCHOS_TEST_FOR_EXCEPTION( Nc>d_curDim, std::logic_error,
        "Anasazi::GeneralizedDavidson::ScaleEigenvectors: Matrix size exceeds current dimension");
    TEUCHOS_TEST_FOR_EXCEPTION( X_alpha.numRows()!=Nr, std::logic_error,
        "Anasazi::GeneralizedDavidson::ScaleEigenvectors: number of rows in Xalpha does not match X");
    TEUCHOS_TEST_FOR_EXCEPTION( X_alpha.numCols()!=Nc, std::logic_error,
        "Anasazi::GeneralizedDavidson::ScaleEigenvectors: number of cols in Xalpha does not match X");
    TEUCHOS_TEST_FOR_EXCEPTION( X_beta.numRows()!=Nr, std::logic_error,
        "Anasazi::GeneralizedDavidson::ScaleEigenvectors: number of rows in Xbeta does not match X");
    TEUCHOS_TEST_FOR_EXCEPTION( X_beta.numCols()!=Nc, std::logic_error,
        "Anasazi::GeneralizedDavidson::ScaleEigenvectors: number of cols in Xbeta does not match X");

    // Now form quasi-diagonal matrices
    //  containing alpha and beta
    Teuchos::SerialDenseMatrix<int,ScalarType> Alpha(Nc,Nc,true);
    Teuchos::SerialDenseMatrix<int,ScalarType> Beta(Nc,Nc,true);

    computeRitzIndex();

    for( int i=0; i<Nc; ++i )
    {
        if( d_ritzIndex[i] == 0 )
        {
            Alpha(i,i) = d_alphar[i];
            Beta(i,i)  = d_betar[i];
        }
        else if( d_ritzIndex[i] == 1 )
        {
            Alpha(i,i)   = d_alphar[i];
            Alpha(i,i+1) = d_alphai[i];
            Beta(i,i)    = d_betar[i];
        }
        else
        {
            Alpha(i,i-1) = d_alphai[i];
            Alpha(i,i)   = d_alphar[i];
            Beta(i,i)    = d_betar[i];
        }
    }

    int err;

    // Multiply the eigenvectors by alpha
    err = X_alpha.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, ST::one(), X, Alpha, ST::zero() );
    std::stringstream astream;
    astream << "GeneralizedDavidson::ScaleEigenvectors: multiply returned error code " << err;
    TEUCHOS_TEST_FOR_EXCEPTION( err!=0, std::runtime_error, astream.str() );

    // Multiply the eigenvectors by beta
    err = X_beta.multiply( Teuchos::NO_TRANS, Teuchos::NO_TRANS, ST::one(), X, Beta, ST::zero() );
    std::stringstream bstream;
    bstream << "GeneralizedDavidson::ScaleEigenvectors: multiply returned error code " << err;
    TEUCHOS_TEST_FOR_EXCEPTION( err!=0, std::runtime_error, bstream.str() );
}

//---------------------------------------------------------------------------//
// Compute residual
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
void GeneralizedDavidson<ScalarType,MV,OP>::computeResidual()
{
    computeRitzIndex();

    // Determine how many residual vectors need to be computed
    d_residualSize = std::max( d_blockSize, d_NEV );
    if( d_curDim < d_residualSize )
    {
        d_residualSize = d_curDim;
    }
    else if( d_ritzIndex[d_residualSize-1] == 1 )
    {
        d_residualSize++;
    }

    // Get indices of all valid residual vectors
    std::vector<int> residualIndices(d_residualSize);
    for( int i=0; i<d_residualSize; ++i )
        residualIndices[i] = i;

    // X will store (right) eigenvectors of projected system
    Teuchos::SerialDenseMatrix<int,ScalarType> X(Teuchos::Copy,*d_Z,d_curDim,d_curDim);

    // Get eigenvectors of projected problem -- computed from previous Schur decomposition
    computeProjectedEigenvectors( X );

    // X_alpha and X_beta will be eigenvectors right-multiplied by alpha, beta (which are quasi-diagonal portions of S,T)
    Teuchos::SerialDenseMatrix<int,ScalarType> X_alpha(d_curDim,d_residualSize);
    Teuchos::SerialDenseMatrix<int,ScalarType>  X_beta(d_curDim,d_residualSize);

    // X_wanted is the wanted portion of X
    Teuchos::SerialDenseMatrix<int,ScalarType> X_wanted(Teuchos::View, X, d_curDim, d_residualSize);

    // Scale Eigenvectors by alpha or beta
    scaleEigenvectors( X_wanted, X_alpha, X_beta );

    // Get view of residual vector(s)
    RCP<MV> R_active = MVT::CloneViewNonConst( *d_R, residualIndices );

    // View of active portion of AV,BV
    std::vector<int> activeIndices(d_curDim);
    for( int i=0; i<d_curDim; ++i )
        activeIndices[i]=i;

    // Compute residual
    RCP<const MV> AV_active = MVT::CloneView( *d_AV, activeIndices );
    MVT::MvTimesMatAddMv(ST::one(),*AV_active, X_beta,  ST::zero(),*R_active);

    if( d_haveB )
    {
        RCP<const MV> BV_active = MVT::CloneView( *d_BV, activeIndices );
        MVT::MvTimesMatAddMv(ST::one(),*BV_active, X_alpha,-ST::one(), *R_active);
    }
    else
    {
        RCP<const MV> V_active = MVT::CloneView( *d_V, activeIndices );
        MVT::MvTimesMatAddMv(ST::one(),*V_active, X_alpha,-ST::one(), *R_active);
    }

    /* Apply a scaling to the residual
     * For generalized eigenvalue problems, LAPACK scales eigenvectors
     *  to have unit length in the infinity norm, we want them to have unit
     *  length in the 2-norm.  For conjugate pairs, the scaling is such that
     *  |xr|^2 + |xi|^2 = 1
     * Additionally, the residual is currently computed as r=beta*A*x-alpha*B*x
     *  but the "standard" residual is r=A*x-(alpha/beta)*B*x, or if we want
     *  to scale the residual by the Ritz value then it is r=(beta/alpha)*A*x-B*x
     *  Performing the scaling this way allows us to avoid the possibility of
     *  diving by infinity or zero if the StatusTest were allowed to handle the
     *  scaling.
     */
    Teuchos::LAPACK<int,ScalarType> lapack;
    Teuchos::BLAS<int,ScalarType> blas;
    std::vector<MagnitudeType> resScaling(d_residualSize);
    for( int icol=0; icol<d_residualSize; ++icol )
    {
        if( d_ritzIndex[icol] == 0 )
        {
            MagnitudeType Xnrm = blas.NRM2( d_curDim, X_wanted[icol], 1);
            MagnitudeType ABscaling = d_relativeConvergence ? d_alphar[icol] : d_betar[icol];
            resScaling[icol] = MT::one() / (Xnrm * ABscaling);
        }
        else if( d_ritzIndex[icol] == 1 )
        {
            MagnitudeType Xnrm1 = blas.NRM2( d_curDim, X_wanted[icol],   1 );
            MagnitudeType Xnrm2 = blas.NRM2( d_curDim, X_wanted[icol+1], 1 );
            MagnitudeType Xnrm  = lapack.LAPY2(Xnrm1,Xnrm2);
            MagnitudeType ABscaling = d_relativeConvergence ? lapack.LAPY2(d_alphar[icol],d_alphai[icol])
                                                            : d_betar[icol];
            resScaling[icol]   = MT::one() / (Xnrm * ABscaling);
            resScaling[icol+1] = MT::one() / (Xnrm * ABscaling);
        }
    }
    MVT::MvScale( *R_active, resScaling );

    // Compute residual norms
    d_resNorms.resize(d_residualSize);
    MVT::MvNorm(*R_active,d_resNorms);

    // If Ritz value i is real, then the corresponding residual vector
    //  is the true residual
    // If Ritz values i and i+1 form a conjugate pair, then the
    //  corresponding residual vectors are the real and imaginary components
    //  of the residual.  Adjust the residual norms appropriately...
    for( int i=0; i<d_residualSize; ++i )
    {
        if( d_ritzIndex[i] == 1 )
        {
            MagnitudeType nrm = lapack.LAPY2(d_resNorms[i],d_resNorms[i+1]);
            d_resNorms[i]   = nrm;
            d_resNorms[i+1] = nrm;
        }
    }

    // Evaluate with status test
    d_tester->checkStatus(this);

    // Determine which residual vectors should be used for subspace expansion
    if( d_useLeading || d_blockSize >= d_NEV )
    {
        d_expansionSize=d_blockSize;
        if( d_ritzIndex[d_blockSize-1]==1 )
            d_expansionSize++;

        d_expansionIndices.resize(d_expansionSize);
        for( int i=0; i<d_expansionSize; ++i )
            d_expansionIndices[i] = i;
    }
    else
    {
        std::vector<int> convergedVectors = d_tester->whichVecs();

        // Get index of first unconverged vector
        int startVec;
        for( startVec=0; startVec<d_residualSize; ++startVec )
        {
            if( std::find(convergedVectors.begin(),convergedVectors.end(),startVec)==convergedVectors.end() )
                break;
        }

        // Now get a contiguous block of indices starting at startVec
        // If this crosses the end of our residual vectors, take the final d_blockSize vectors
        int endVec = startVec + d_blockSize - 1;
        if( endVec > (d_residualSize-1) )
        {
            endVec   = d_residualSize-1;
            startVec = d_residualSize-d_blockSize;
        }

        // Don't split conjugate pairs on either end of the range
        if( d_ritzIndex[startVec]==-1 )
        {
            startVec--;
            endVec--;
        }

        if( d_ritzIndex[endVec]==1 )
            endVec++;

        d_expansionSize = 1+endVec-startVec;
        d_expansionIndices.resize(d_expansionSize);
        for( int i=0; i<d_expansionSize; ++i )
            d_expansionIndices[i] = startVec+i;
    }
}

//---------------------------------------------------------------------------//
// Print current status.
//---------------------------------------------------------------------------//
template <class ScalarType, class MV, class OP>
void GeneralizedDavidson<ScalarType,MV,OP>::currentStatus( std::ostream &myout )
{
    using std::endl;

    myout.setf(std::ios::scientific, std::ios::floatfield);
    myout.precision(6);
    myout <<endl;
    myout <<"================================================================================" << endl;
    myout << endl;
    myout <<"                    GeneralizedDavidson Solver Solver Status" << endl;
    myout << endl;
    myout <<"The solver is "<<(d_initialized ? "initialized." : "not initialized.") << endl;
    myout <<"The number of iterations performed is " << d_iteration << endl;
    myout <<"The number of operator applies performed is " << d_opApplies << endl;
    myout <<"The block size is         " << d_expansionSize << endl;
    myout <<"The current basis size is " << d_curDim << endl;
    myout <<"The number of requested eigenvalues is " << d_NEV << endl;
    myout <<"The number of converged values is " << d_tester->howMany() << endl;
    myout << endl;

    myout.setf(std::ios_base::right, std::ios_base::adjustfield);

    if( d_initialized )
    {
        myout << "CURRENT RITZ VALUES" << endl;

        myout << std::setw(24) << "Ritz Value"
              << std::setw(30) << "Residual Norm" << endl;
        myout << "--------------------------------------------------------------------------------" << endl;
        if( d_residualSize > 0 )
        {
            std::vector<MagnitudeType> realRitz(d_curDim), imagRitz(d_curDim);
            std::vector< Value<ScalarType> > ritzVals = getRitzValues();
            for( int i=0; i<d_curDim; ++i )
            {
                realRitz[i] = ritzVals[i].realpart;
                imagRitz[i] = ritzVals[i].imagpart;
            }
            std::vector<int> permvec;
            sortValues( realRitz, imagRitz, permvec, d_curDim );

            int numToPrint = std::max( d_numToPrint, d_NEV );
            numToPrint = std::min( d_curDim, numToPrint );

            // Because the sort manager does not use a stable sort, occasionally
            //  the portion of a conjugate pair with negative imaginary part will be placed
            //  first...in that case the following will not give the usual expected behavior
            //  and an extra value will be printed.  This is only an issue with the output
            //  format because the actually sorting of Schur forms is guaranteed to be stable.
            if( d_ritzIndex[permvec[numToPrint-1]] != 0 )
                numToPrint++;

            int i=0;
            while( i<numToPrint )
            {
                if( imagRitz[i] == ST::zero() )
                {
                    myout << std::setw(15) << realRitz[i];
                    myout << " + i" << std::setw(15) << ST::magnitude( imagRitz[i] );
                    if( i < d_residualSize )
                        myout << std::setw(20) << d_resNorms[permvec[i]] << endl;
                    else
                        myout << "        Not Computed" << endl;

                    i++;
                }
                else
                {
                    // Positive imaginary part
                    myout << std::setw(15) << realRitz[i];
                    myout << " + i" << std::setw(15) << ST::magnitude( imagRitz[i] );
                    if( i < d_residualSize )
                        myout << std::setw(20) << d_resNorms[permvec[i]] << endl;
                    else
                        myout << "        Not Computed" << endl;

                    // Negative imaginary part
                    myout << std::setw(15) << realRitz[i];
                    myout << " - i" << std::setw(15) << ST::magnitude( imagRitz[i] );
                    if( i < d_residualSize )
                        myout << std::setw(20) << d_resNorms[permvec[i]] << endl;
                    else
                        myout << "        Not Computed" << endl;

                    i+=2;
                }
            }
        }
        else
        {
            myout << std::setw(20) << "[ NONE COMPUTED ]" << endl;
        }
    }
    myout << endl;
    myout << "================================================================================" << endl;
    myout << endl;
}

} // namespace Anasazi

#endif // ANASAZI_GENERALIZED_DAVIDSON_HPP

