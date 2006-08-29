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


/*! \file AnasaziLOBPCG.hpp
  \brief Implementation of the locally-optimal block preconditioned conjugate gradient (LOBPCG) method
*/

/*
    LOBPCG contains local storage of up to 10*blockSize_ vectors, representing 10 entities
      X,H,P,R
      KX,KH,KP  (product of K and the above)
      MX,MH,MP  (product of M and the above, not allocated if we don't have an M matrix)
    If full orthogonalization is enabled, one extra multivector of blockSize_ vectors is required to 
    compute the local update of X and P.
    
    A solver is bound to an eigenproblem at declaration.
    Other solver parameters (e.g., block size, auxiliary vectors) can be changed dynamically.
    
    The orthogonalization manager is used to project away from the auxiliary vectors.
    If full orthogonalization is enabled, the orthogonalization manager is also used to construct an M orthonormal basis.
    The orthogonalization manager is subclass of MatOrthoManager, which LOBPCG assumes to be defined by the M inner product.
    LOBPCG will not work correctly if the orthomanager uses a different inner product.
 */


#ifndef ANASAZI_LOBPCG_HPP
#define ANASAZI_LOBPCG_HPP

#include "AnasaziTypes.hpp"

#include "AnasaziEigensolver.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "AnasaziMatOrthoManager.hpp"
#include "AnasaziModalSolverUtils.hpp"

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

/*!     \class Anasazi::LOBPCG
        
        \brief This class implements the Locally-Optimal Block Preconditioned
        Conjugate Gradient (%LOBPCG) method for solving Hermitian eigenvalue problems.
        
        This implementation is a modification of the one found in 
        A. Knyazev, "Toward the optimal preconditioned eigensolver:
        Locally optimal block preconditioner conjugate gradient method",
        SIAM J. Sci. Comput., vol 23, n 2, pp. 517-541.
        
        The modification consists of the orthogonalization steps recommended in
        U. Hetmaniuk and R. Lehoucq, "Basis Selection in LOBPCG", Journal of Computational Physics. 
        
        These modifcation are referred to as full orthogonalization, and consist of also conducting
        the local optimization using an orthonormal basis.
        

        \ingroup anasazi_solver_framework

        \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist
*/

namespace Anasazi {

  //! @name LOBPCG Structures
  //@{ 

  /** \brief Structure to contain pointers to Anasazi state variables.
   *
   * This struct is utilized by LOBPCG::initialize() and LOBPCG::getState().
   */
  template <class ScalarType, class MV>
  struct LOBPCGState {
    //! The current eigenvectors.
    Teuchos::RefCountPtr<const MV> X; 
    //! The image of the current eigenvectors under K.
    Teuchos::RefCountPtr<const MV> KX; 
    //! The image of the current eigenvectors under M, or Teuchos::null if M was not specified.
    Teuchos::RefCountPtr<const MV> MX;
    //! The current search direction.
    Teuchos::RefCountPtr<const MV> P; 
    //! The image of the current search direction under K.
    Teuchos::RefCountPtr<const MV> KP; 
    //! The image of the current search direction under M, or Teuchos::null if M was not specified.
    Teuchos::RefCountPtr<const MV> MP;
    /*! \brief The current preconditioned residual vectors.
     *
     *  H is only useful when LOBPCG::iterate() throw a LOBPCGRitzFailure exception.
     */
    Teuchos::RefCountPtr<const MV> H; 
    //! The image of the current preconditioned residual vectors under K.
    Teuchos::RefCountPtr<const MV> KH; 
    //! The image of the current preconditioned residual vectors under M, or Teuchos::null if M was not specified.
    Teuchos::RefCountPtr<const MV> MH;
    //! The current residual vectors.
    Teuchos::RefCountPtr<const MV> R;
    //! The current Ritz values.
    Teuchos::RefCountPtr<const std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> > T;
    LOBPCGState() : X(Teuchos::null),KX(Teuchos::null),MX(Teuchos::null),
                    P(Teuchos::null),KP(Teuchos::null),MP(Teuchos::null),
                    H(Teuchos::null),KH(Teuchos::null),MH(Teuchos::null),
                    R(Teuchos::null),T(Teuchos::null) {};
  };

  //@}

  //! @name LOBPCG Exceptions
  //@{ 

  /** \brief LOBPCGRitzFailure is thrown when the LOBPCG solver is unable to
   *  continue a call to LOBPCG::iterate() due to a failure of the algorithm.
   *
   *  This signals that the Rayleigh-Ritz analysis over the subspace \c
   *  colsp([X H P]) detected ill-conditioning of the projected mass matrix
   *  and the inability to generate a set of orthogonal eigenvectors for 
   *  the projected problem.
   *
   *  This exception is only thrown from the LOBPCG::iterate() routine. After
   *  catching this exception, the user can recover the subspace via
   *  LOBPCG::getState(). This information can be used to restart the solver.
   *
   */
  class LOBPCGRitzFailure : public AnasaziError {public:
    LOBPCGRitzFailure(const std::string& what_arg) : AnasaziError(what_arg)
    {}};

  /** \brief LOBPCGInitFailure is thrown when the LOBPCG solver is unable to
   * generate an initial iterate in the LOBPCG::initialize() routine. 
   *
   * This exception is thrown from the LOBPCG::initialize() method, which is
   * called by the user or from the LOBPCG::iterate() method when isInitialized()
   * == \c false.
   *
   * In the case that this exception is thrown, LOBPCG::hasP() and
   * LOBPCG::isInitialized() will be \c false and the user will need to provide
   * a new initial iterate to the solver.
   *
   */
  class LOBPCGInitFailure : public AnasaziError {public:
    LOBPCGInitFailure(const std::string& what_arg) : AnasaziError(what_arg)
    {}};

  /** \brief LOBPCGOrthoFailure is thrown when an orthogonalization attempt 
   * fails.
   *
   * This is thrown in one of two scenarios. After preconditioning the residual,
   * the orthogonalization manager is asked to orthogonalize the preconditioned
   * residual (H) against the auxiliary vectors. If full orthogonalization
   * is enabled, H is also orthogonalized against X and P and normalized.
   *
   * The second scenario involves the generation of new X and P from the
   * basis [X H P]. When full orthogonalization is enabled, an attempt is
   * made to select coefficients for X and P so that they will be
   * mutually orthogonal and orthonormal.
   *
   * If either of these attempts fail, the solver throws an LOBPCGOrthoFailure
   * exception.
   */
  class LOBPCGOrthoFailure : public AnasaziError {public:
    LOBPCGOrthoFailure(const std::string& what_arg) : AnasaziError(what_arg)
    {}};

  //@}


  template <class ScalarType, class MV, class OP>
  class LOBPCG : public Eigensolver<ScalarType,MV,OP> { 
  public:
    
    //! @name Constructor/Destructor
    //@{ 
    
    /*! \brief LOBPCG constructor with eigenproblem, solver utilities, and parameter list of solver options.
     *
     * This constructor takes pointers required by the eigensolver, in addition
     * to a parameter list of options for the eigensolver. These options include the following:
     *   - "Block Size" - an \c int specifying the block size used by the algorithm. This can also be specified using the setBlockSize() method.
     *   - "Full Ortho" - a \c bool specifying whether the solver should employ a full orthogonalization technique. This can also be specified using the setFullOrtho() method.
     */
    LOBPCG( const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem, 
            const Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > &sorter,
            const Teuchos::RefCountPtr<OutputManager<ScalarType> > &printer,
            const Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > &tester,
            const Teuchos::RefCountPtr<MatOrthoManager<ScalarType,MV,OP> > &ortho,
            Teuchos::ParameterList &params 
          );
    
    //! %LOBPCG destructor
    virtual ~LOBPCG() {};

    //@}

    //! @name Solver methods
    //@{
    
    /*! \brief This method performs %LOBPCG iterations until the status test
     * indicates the need to stop or an error occurs (in which case, an
     * exception is thrown).
     *
     * iterate() will first determine whether the solver is initialized; if
     * not, it will call initialize() using default arguments.  After
     * initialization, the solver performs %LOBPCG iterations until the status
     * test evaluates as Passed, at which point the method returns to the
     * caller.
     *
     * The %LOBPCG iteration proceeds as follows:
     * -# The current residual (R) is preconditioned to form H
     * -# H is orthogonalized against the auxiliary vectors and, if full orthogonalization\n
     *    is enabled, against X and P. 
     * -# The basis [X H P] is used to project the problem matrices.
     * -# The projected eigenproblem is solved, and the desired eigenvectors and eigenvalues are selected.
     * -# These are used to form the new eigenvector estimates (X) and the search directions (P).\n
     *    If full orthogonalization is enabled, these are generated to be mutually orthogonal and with orthonormal columns.
     * -# The new residual (R) is formed.
     *
     * The status test is queried at the beginning of the iteration.
     *
     * Possible exceptions thrown include std::logic_error, std::invalid_argument or
     * one of the LOBPCG-specific exceptions.
     *
    */
    void iterate();

    /*! \brief Initialize the solver to an iterate, optionally providing the
     * Ritz values, residual, and search direction.
     *
     * The %LOBPCG eigensolver contains a certain amount of state relating to
     * the current iterate, including the current residual, the current search
     * direction, and the images of these spaces under the operators.
     *
     * initialize() gives the user the opportunity to manually set these,
     * although this must be done with caution, abiding by the rules
     * given below. All notions of orthogonality and orthonormality are derived
     * from the inner product specified by the orthogonalization manager.
     *
     * \post 
     *   - isInitialized() == true (see post-conditions of isInitialize())
     *   - If newstate.P != Teuchos::null, hasP() == true.\n
     *     Otherwise, hasP() == false
     *
     * The user has the option of specifying any component of the state using
     * initialize(). However, these arguments are assumed to match the
     * post-conditions specified under isInitialized(). Any component of the
     * state (i.e., KX) not given to initialize() will be generated.
     *
     */
    void initialize(LOBPCGState<ScalarType,MV> newstate);

    /*! \brief Initialize the solver with the initial vectors from the eigenproblem
     *  or random data.
     */
    void initialize();

    /*! \brief Indicates whether the solver has been initialized or not.
     *
     * \return bool indicating the state of the solver.
     * \post
     * If isInitialized() == \c true:
     *   - X is orthogonal to auxiliary vectors and has orthonormal columns
     *   - KX == Op*X
     *   - MX == M*X if M != Teuchos::null\n
     *     Otherwise, MX == Teuchos::null
     *   - getRitzValues() returns the sorted Ritz values with respect to X
     *   - getResidualVecs() returns the residual vectors with respect to X
     *   - If hasP() == \c true,
     *      - P orthogonal to auxiliary vectors
     *      - If getFullOrtho() == \c true,
     *        - P is orthogonal to X and has orthonormal columns
     *      - KP == Op*P
     *      - MP == M*P if M != Teuchos::null\n
     *        Otherwise, MP == Teuchos::null
     */
    bool isInitialized() { return initialized_; }

    /*! \brief Get the current state of the eigensolver.
     * 
     * The data is only valid if isInitialized() == \c true. The
     * data for the search directions P is only meaningful if hasP() == \c
     * true. Finally, the data for the preconditioned residual (H) is only meaningful in the situation where
     * the solver throws an ::LOBPCGRitzFailure exception during iterate().
     *
     * \returns An LOBPCGState object containing const pointers to the current
     * solver state.
     */
    LOBPCGState<ScalarType,MV> getState() const {
      LOBPCGState<ScalarType,MV> state;
      state.X = X_;
      state.KX = KX_;
      state.P = P_;
      state.KP = KP_;
      state.H = H_;
      state.KH = KH_;
      state.R = R_;
      state.T = Teuchos::rcp(new std::vector<MagnitudeType>(theta_));
      if (hasM_) {
        state.MX = MX_;
        state.MP = MP_;
        state.MH = MH_;
      }
      else {
        state.MX = Teuchos::null;
        state.MP = Teuchos::null;
        state.MH = Teuchos::null;
      }
      return state;
    }

    //@}

    //! @name Status methods
    //@{

    //! \brief Get the current iteration count.
    int getNumIters() const { return(iter_); };

    //! \brief Reset the iteration count.
    void resetNumIters() { iter_=0; };

    /*! \brief Get the Ritz vectors from the previous iteration.
      
        \return A multivector with getBlockSize() vectors containing 
        the sorted Ritz vectors corresponding to the most significant Ritz values.
        The i-th vector of the return corresponds to the i-th Ritz vector; there is no need to use
        getRitzIndex().
     */
    Teuchos::RefCountPtr<const MV> getRitzVectors() {return X_;}

    /*! \brief Get the Ritz values from the previous iteration.
     *
     *  \return A vector of length getCurSubspaceDim() containing the Ritz values from the
     *  previous projected eigensolve.
     */
    std::vector<Value<ScalarType> > getRitzValues() { 
      std::vector<Value<ScalarType> > ret(nevLocal_);
      for (int i=0; i<nevLocal_; i++) {
        ret[i].realpart = theta_[i];
        ret[i].imagpart = ZERO;
      }
      return ret;
    }

    /*! \brief Get the index used for extracting Ritz vectors from getRitzVectors().
     *
     * Because BlockDavidson is a Hermitian solver, all Ritz values are real and all Ritz vectors can be represented in a 
     * single column of a multivector. Therefore, getRitzIndex() is not needed when using the output from getRitzVectors().
     *
     * \return An \c int vector of size getCurSubspaceDim() composed of zeros.
     */
    std::vector<int> getRitzIndex() {
      std::vector<int> ret(nevLocal_,0);
      return ret;
    }




    /*! \brief Get the current residual norms
     *
     *  \return A vector of length getCurSubspaceDim() containing the norms of the
     *  residuals, with respect to the orthogonalization manager norm() method.
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getResNorms();


    /*! \brief Get the current residual 2-norms
     *
     *  \return A vector of length getCurSubspaceDim() containing the 2-norms of the
     *  residuals. 
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getRes2Norms();


    /*! \brief Get the 2-norms of the Ritz residuals.
     *
     *  \return A vector of length getCurSubspaceDim() containing the 2-norms of the Ritz residuals.
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getRitzRes2Norms() {
      std::vector<MagnitudeType> ret = ritz2norms_;
      ret.resize(nevLocal_);
      return ret;
    }


    /*! \brief Get the dimension of the search subspace used to generate the current eigenvectors and eigenvalues.
     *
     *  %LOBPCG employs a sequential subspace iteration, maintaining a fixed-rank basis, as opposed to an expanding subspace
     *  mechanism employed by Krylov-subspace solvers like BlockKrylovSchur and BlockDavidson.
     *  
     *  \return An integer specifying the rank of the subspace generated by the eigensolver. If isInitialized() == \c false, 
     *  the return is 0. Otherwise, the return will be 2*getBlockSize() or 3*getBlockSize().
     */
    int getCurSubspaceDim() const {
      if (!initialized_) return 0;
      return nevLocal_;
    }

    /*! \brief Get the maximum dimension allocated for the search subspace. For %LOBPCG, this always returns 3*getBlockSize(), the dimension of the 
     *   subspace colspan([X H P]).
     */
    int getMaxSubspaceDim() const {return 3*blockSize_;}

    //@}

    //!  @name Accessor routines from Eigensolver
    //@{


    //! Get a constant reference to the eigenvalue problem.
    const Eigenproblem<ScalarType,MV,OP>& getProblem() const { return(*problem_); };


    /*! \brief Set the blocksize to be used by the iterative solver in solving
     * this eigenproblem.
     *  
     *  If the block size is reduced, then the new iterate (and residual and
     *  search direction) are chosen as the subset of the current iterate
     *  preferred by the sort manager.  Otherwise, the solver state is set to
     *  uninitialized.
     */
    void setBlockSize(int blockSize);


    //! Get the blocksize to be used by the iterative solver in solving this eigenproblem.
    int getBlockSize() const { return(blockSize_); }


    /*! \brief Set the auxiliary vectors for the solver.
     *
     *  Because the current iterate X and search direction P cannot be assumed
     *  orthogonal to the new auxiliary vectors, a call to setAuxVecs() with a
     *  non-empty argument will reset the solver to the uninitialized state.
     *
     *  In order to preserve the current state, the user will need to extract
     *  it from the solver using getState(), orthogonalize it against the new
     *  auxiliary vectors, and manually reinitialize the solver using
     *  initialize().
     */
    void setAuxVecs(const Teuchos::Array<Teuchos::RefCountPtr<const MV> > &auxvecs);

    //! Get the current auxiliary vectors.
    Teuchos::Array<Teuchos::RefCountPtr<const MV> > getAuxVecs() const {return auxVecs_;}

    //@}

    //!  @name %LOBPCG-specific accessor routines
    //@{

    /*! \brief Instruct the LOBPCG iteration to use full orthogonality.
     *
     *  If the getFullOrtho() == \c false and isInitialized() == \c true and hasP() == \c true, then
     *  P will be invalidated by setting full orthogonalization to \c true.
     */
    void setFullOrtho(bool fullOrtho);

    //! Determine if the LOBPCG iteration is using full orthogonality.
    bool getFullOrtho() const { return(fullOrtho_); }
    
    //! Indicates whether the search direction given by getState() is valid.
    bool hasP() {return hasP_;}

    //@}
    
    //!  @name Output methods
    //@{

    //! This method requests that the solver print out its current status to screen.
    void currentStatus(ostream &os);

    //@}

  private:
    //
    // Convenience typedefs
    //
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename SCT::magnitudeType MagnitudeType;
    const MagnitudeType ONE;  
    const MagnitudeType ZERO; 
    const MagnitudeType NANVAL;
    //
    // Internal structs
    //
    struct CheckList {
      bool checkX, checkMX, checkKX;
      bool checkH, checkMH;
      bool checkP, checkMP, checkKP;
      bool checkR, checkQ;
      CheckList() : checkX(false),checkMX(false),checkKX(false),
                    checkH(false),checkMH(false),
                    checkP(false),checkMP(false),checkKP(false),
                    checkR(false),checkQ(false) {};
    };
    //
    // Internal methods
    //
    string accuracyCheck(const CheckList &chk, const string &where) const;
    //
    // Classes inputed through constructor that define the eigenproblem to be solved.
    //
    const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> >     problem_;
    const Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> >      sm_;
    const Teuchos::RefCountPtr<OutputManager<ScalarType> >          om_;
    const Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> >       tester_;
    const Teuchos::RefCountPtr<MatOrthoManager<ScalarType,MV,OP> >  orthman_;
    //
    // Information obtained from the eigenproblem
    //
    Teuchos::RefCountPtr<OP> Op_;
    Teuchos::RefCountPtr<OP> MOp_;
    Teuchos::RefCountPtr<OP> Prec_;
    bool hasM_;
    //
    // Internal utilities class required by eigensolver.
    //
    ModalSolverUtils<ScalarType,MV,OP> MSUtils_;
    //
    // Internal timers
    //
    Teuchos::RefCountPtr<Teuchos::Time> timerOp_, timerMOp_, timerPrec_,
                                        timerSort_, 
                                        timerLocalProj_, timerDS_,
                                        timerLocalUpdate_, timerCompRes_,
                                        timerOrtho_, timerInit_;
    //
    // Counters
    //
    // Number of operator applications
    int count_ApplyOp_, count_ApplyM_, count_ApplyPrec_;

    //
    // Algorithmic parameters.
    //
    // blockSize_ is the solver block size
    int blockSize_;
    //
    // fullOrtho_ dictates whether the orthogonalization procedures specified by Hetmaniuk and Lehoucq should
    // be activated (see citations at the top of this file)
    bool fullOrtho_;

    //
    // Current solver state
    //
    // initialized_ specifies that the basis vectors have been initialized and the iterate() routine
    // is capable of running; _initialize is controlled  by the initialize() member method
    // For the implications of the state of initialized_, please see documentation for initialize()
    bool initialized_;
    //
    // nevLocal_ reflects how much of the current basis is valid (0 <= nevLocal_ <= 3*blockSize_)
    // this tells us how many of the values in theta_ are valid Ritz values
    int nevLocal_;
    //
    // hasP_ tells us whether there is valid data in P (and KP,MP)
    bool hasP_;
    //
    // State Multivecs
    Teuchos::RefCountPtr<MV> X_, KX_, MX_, R_,
                             H_, KH_, MH_,
                             P_, KP_, MP_;
    // tmpMV is needed only if fullOrtho_ == true
    // because is depends on fullOrtho_, which is easily toggled by the user, we will allocate it 
    // and deallocate it inside of iterate()
    Teuchos::RefCountPtr<MV> tmpMV_;        
    // 
    // auxiliary vectors
    Teuchos::Array<Teuchos::RefCountPtr<const MV> > auxVecs_;
    int numAuxVecs_;
    //
    // Number of iterations that have been performed.
    int iter_;
    // 
    // Current eigenvalues, residual norms
    std::vector<MagnitudeType> theta_, Rnorms_, R2norms_, ritz2norms_;
    // 
    // are the residual norms current with the residual?
    bool Rnorms_current_, R2norms_current_;

  };




  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor
  template <class ScalarType, class MV, class OP>
  LOBPCG<ScalarType,MV,OP>::LOBPCG(
        const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem, 
        const Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > &sorter,
        const Teuchos::RefCountPtr<OutputManager<ScalarType> > &printer,
        const Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > &tester,
        const Teuchos::RefCountPtr<MatOrthoManager<ScalarType,MV,OP> > &ortho,
        Teuchos::ParameterList &params
        ) :
    ONE(Teuchos::ScalarTraits<MagnitudeType>::one()),
    ZERO(Teuchos::ScalarTraits<MagnitudeType>::zero()),
    NANVAL(Teuchos::ScalarTraits<MagnitudeType>::nan()),
    // problem, tools
    problem_(problem), 
    sm_(sorter),
    om_(printer),
    tester_(tester),
    orthman_(ortho),
    MSUtils_(om_),
    // timers, counters
    timerOp_(Teuchos::TimeMonitor::getNewTimer("Operation Op*x")),
    timerMOp_(Teuchos::TimeMonitor::getNewTimer("Operation M*x")),
    timerPrec_(Teuchos::TimeMonitor::getNewTimer("Operation Prec*x")),
    timerSort_(Teuchos::TimeMonitor::getNewTimer("Sorting eigenvalues")),
    timerLocalProj_(Teuchos::TimeMonitor::getNewTimer("Local projection")),
    timerDS_(Teuchos::TimeMonitor::getNewTimer("Direct solve")),
    timerLocalUpdate_(Teuchos::TimeMonitor::getNewTimer("Local update")),
    timerCompRes_(Teuchos::TimeMonitor::getNewTimer("Computing residuals")),
    timerOrtho_(Teuchos::TimeMonitor::getNewTimer("Orthogonalization")),
    timerInit_(Teuchos::TimeMonitor::getNewTimer("Initialization")),
    count_ApplyOp_(0),
    count_ApplyM_(0),
    count_ApplyPrec_(0),
    // internal data
    blockSize_(0),
    fullOrtho_(params.get("Full Ortho", true)),
    initialized_(false),
    nevLocal_(0),
    hasP_(false),
    auxVecs_( Teuchos::Array<Teuchos::RefCountPtr<const MV> >(0) ), 
    numAuxVecs_(0),
    iter_(0),
    Rnorms_current_(false),
    R2norms_current_(false)
  {     
    TEST_FOR_EXCEPTION(problem_ == Teuchos::null,std::invalid_argument,
                       "Anasazi::LOBPCG::constructor: user passed null problem pointer.");
    TEST_FOR_EXCEPTION(sm_ == Teuchos::null,std::invalid_argument,
                       "Anasazi::LOBPCG::constructor: user passed null sort manager pointer.");
    TEST_FOR_EXCEPTION(om_ == Teuchos::null,std::invalid_argument,
                       "Anasazi::LOBPCG::constructor: user passed null output manager pointer.");
    TEST_FOR_EXCEPTION(tester_ == Teuchos::null,std::invalid_argument,
                       "Anasazi::LOBPCG::constructor: user passed null status test pointer.");
    TEST_FOR_EXCEPTION(orthman_ == Teuchos::null,std::invalid_argument,
                       "Anasazi::LOBPCG::constructor: user passed null orthogonalization manager pointer.");
    TEST_FOR_EXCEPTION(problem_->isProblemSet() == false, std::invalid_argument,
                       "Anasazi::LOBPCG::constructor: problem is not set.");
    TEST_FOR_EXCEPTION(problem_->isHermitian() == false, std::invalid_argument,
                       "Anasazi::LOBPCG::constructor: problem is not hermitian.");

    // get the problem operators
    Op_   = problem_->getOperator();
    TEST_FOR_EXCEPTION(Op_ == Teuchos::null, std::invalid_argument,
                       "Anasazi::LOBPCG::constructor: problem provides no operator.");
    MOp_  = problem_->getM();
    Prec_ = problem_->getPrec();
    hasM_ = (MOp_ != Teuchos::null);

    // set the block size and allocate data
    int bs = params.get("Block Size", problem_->getNEV());
    setBlockSize(bs);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the block size and make necessary adjustments.
  template <class ScalarType, class MV, class OP>
  void LOBPCG<ScalarType,MV,OP>::setBlockSize (int blockSize) 
  {
    // time spent here counts towards timerInit_
    Teuchos::TimeMonitor lcltimer( *timerInit_ );


    // This routine only allocates space; it doesn't not perform any computation
    // if size is decreased, take the first blockSize vectors of all and leave state as is
    // otherwise, grow/allocate space and set solver to unitialized

    TEST_FOR_EXCEPTION(blockSize <= 0, std::invalid_argument, "Anasazi::LOBPCG::setBlockSize was passed a non-positive block size");
    if (blockSize == blockSize_) {
      // do nothing
      return;
    }
    else if (blockSize < blockSize_) {
      // shrink vectors
      blockSize_ = blockSize;

      theta_.resize(3*blockSize_);
      ritz2norms_.resize(3*blockSize_);
      Rnorms_.resize(blockSize_);
      R2norms_.resize(blockSize_);

      if (initialized_) {
        // shrink multivectors with copy
        // create ind = {0, 1, ..., blockSize-1}
        std::vector<int> ind(blockSize_);
        for (int i=0; i<blockSize_; i++) ind[i] = i;
        
        X_  = MVT::CloneCopy(*X_,ind);
        KX_ = MVT::CloneCopy(*KX_,ind);
        if (hasM_) {
          MX_ = MVT::CloneCopy(*MX_,ind);
        }
        else {
          MX_ = X_;
        }
        R_  = MVT::CloneCopy(*R_,ind);
        P_  = MVT::CloneCopy(*P_,ind);
        KP_ = MVT::CloneCopy(*KP_,ind);
        if (hasM_) {
          MP_ = MVT::CloneCopy(*MP_,ind);
        }
        else {
          MP_ = P_;
        }
      }
      else {
        // shrink multivectors without copying
        X_ = MVT::Clone(*X_,blockSize_);
        KX_ = MVT::Clone(*KX_,blockSize_);
        if (hasM_) {
          MX_ = MVT::Clone(*MX_,blockSize_);
        }
        else {
          MX_ = X_;
        }
        R_ = MVT::Clone(*R_,blockSize_);
        P_ = MVT::Clone(*P_,blockSize_);
        KP_ = MVT::Clone(*KP_,blockSize_);
        if (hasM_) {
          MP_ = MVT::Clone(*MP_,blockSize_);
        }
        else {
          MP_ = P_;
        }
      }
      // shrink H
      H_ = MVT::Clone(*H_,blockSize_);
      KH_ = MVT::Clone(*KH_,blockSize_);
      if (hasM_) {
        MH_ = MVT::Clone(*MH_,blockSize_);
      }
      else {
        MH_ = H_;
      }
    } 
    else {  // blockSize > blockSize_
      // this is also the scenario for our initial call to setBlockSize(), in the constructor
      initialized_ = false;

      Teuchos::RefCountPtr<const MV> tmp;
      // grab some Multivector to Clone
      // in practice, getInitVec() should always provide this, but it is possible to use a 
      // Eigenproblem with nothing in getInitVec() by manually initializing with initialize(); 
      // in case of that strange scenario, we will try to Clone from X_
      if (blockSize_ > 0) {
        tmp = X_;
      }
      else {
        tmp = problem_->getInitVec();
        TEST_FOR_EXCEPTION(tmp == Teuchos::null,std::logic_error,
                           "Anasazi::LOBPCG::setBlockSize(): Eigenproblem did not specify initial vectors to clone from");
      }
      // grow/allocate vectors
      theta_.resize(3*blockSize,NANVAL);
      ritz2norms_.resize(3*blockSize_,NANVAL);
      Rnorms_.resize(blockSize,NANVAL);
      R2norms_.resize(blockSize,NANVAL);
      
      // clone multivectors off of tmp
      X_ = MVT::Clone(*tmp,blockSize);
      KX_ = MVT::Clone(*tmp,blockSize);
      if (hasM_) {
        MX_ = MVT::Clone(*tmp,blockSize);
      }
      else {
        MX_ = X_;
      }
      R_ = MVT::Clone(*tmp,blockSize);
      H_ = MVT::Clone(*tmp,blockSize);
      KH_ = MVT::Clone(*tmp,blockSize);
      if (hasM_) {
        MH_ = MVT::Clone(*tmp,blockSize);
      }
      else {
        MH_ = H_;
      }
      hasP_ = false;
      P_ = MVT::Clone(*tmp,blockSize);
      KP_ = MVT::Clone(*tmp,blockSize);
      if (hasM_) {
        MP_ = MVT::Clone(*tmp,blockSize);
      }
      else {
        MP_ = P_;
      }
      blockSize_ = blockSize;
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the auxiliary vectors
  template <class ScalarType, class MV, class OP>
  void LOBPCG<ScalarType,MV,OP>::setAuxVecs(const Teuchos::Array<Teuchos::RefCountPtr<const MV> > &auxvecs) {
    typedef typename Teuchos::Array<Teuchos::RefCountPtr<const MV> >::iterator tarcpmv;

    // set new auxiliary vectors
    auxVecs_ = auxvecs;
    
    numAuxVecs_ = 0;
    for (tarcpmv i=auxVecs_.begin(); i != auxVecs_.end(); i++) {
      numAuxVecs_ += MVT::GetNumberVecs(**i);
    }
    
    // If the solver has been initialized, X and P are not necessarily orthogonal to new auxiliary vectors
    if (numAuxVecs_ > 0 && initialized_) {
      initialized_ = false;
      hasP_ = false;
    }

    if (om_->isVerbosity( Debug ) ) {
      // Check almost everything here
      CheckList chk;
      chk.checkQ = true;
      om_->print( Debug, accuracyCheck(chk, ": in setAuxVecs()") );
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  /* Initialize the state of the solver
   * 
   * POST-CONDITIONS:
   *
   * initialized_ == true
   * X is orthonormal, orthogonal to auxVecs_
   * KX = Op*X
   * MX = M*X if hasM_
   * theta_ contains Ritz values of X
   * R = KX - MX*diag(theta_)
   * if hasP() == true,
   *   P orthogonal to auxVecs_
   *   if fullOrtho_ == true,
   *     P orthonormal and orthogonal to X
   *   KP = Op*P
   *   MP = M*P
   */
  template <class ScalarType, class MV, class OP>
  void LOBPCG<ScalarType,MV,OP>::initialize(LOBPCGState<ScalarType,MV> newstate)
  {
    // NOTE: memory has been allocated by setBlockSize(). Use SetBlock below; do not Clone
    // NOTE: Time spent in this routine is allotted to timerInit_, in addition to the respective sections.

    hasP_ = false;  // this will be set to true below if appropriate

    std::vector<int> bsind(blockSize_);
    for (int i=0; i<blockSize_; i++) bsind[i] = i;

    // set up X: if the user doesn't specify X, ignore the rest
    if (newstate.X != Teuchos::null && MVT::GetNumberVecs(*newstate.X) >= blockSize_ && MVT::GetVecLength(*newstate.X) == MVT::GetVecLength(*X_) ) {

      Teuchos::TimeMonitor lcltimer( *timerInit_ );

      // put data in X,MX,KX
      MVT::SetBlock(*newstate.X,bsind,*X_);
      if (hasM_) {
        if (newstate.MX != Teuchos::null && MVT::GetNumberVecs(*newstate.MX) >= blockSize_ && MVT::GetVecLength(*newstate.MX) == MVT::GetVecLength(*MX_) ) {
          MVT::SetBlock(*newstate.MX,bsind,*MX_);
        }
        else {
          Teuchos::TimeMonitor lcltimer( *timerMOp_ );
          OPT::Apply(*MOp_,*X_,*MX_);
          count_ApplyM_ += blockSize_;
        }
      }
      else {
        // an assignment would be redundant; take advantage of this opportunity to debug a little
        TEST_FOR_EXCEPTION(MX_ != X_, std::logic_error, "Anasazi::LOBPCG::initialize(): solver invariant not satisfied");
      }
      if (newstate.KX != Teuchos::null && MVT::GetNumberVecs(*newstate.KX) >= blockSize_ && MVT::GetVecLength(*newstate.KX) == MVT::GetVecLength(*KX_) ) {
        MVT::SetBlock(*newstate.KX,bsind,*KX_);
      }
      else {
        Teuchos::TimeMonitor lcltimer( *timerOp_ );
        OPT::Apply(*Op_,*X_,*KX_);
        count_ApplyOp_ += blockSize_;
      }

      // set up Ritz values
      theta_.resize(3*blockSize_,NANVAL);
      ritz2norms_.resize(3*blockSize_,NANVAL);
      if (newstate.T != Teuchos::null && (signed int)(newstate.T->size()) >= blockSize_) {
        for (int i=0; i<blockSize_; i++) {
          theta_[i] = (*newstate.T)[i];
        }
      }
      else {
        // get ritz vecs/vals
        Teuchos::SerialDenseMatrix<int,ScalarType> KK(blockSize_,blockSize_),
                                                   MM(blockSize_,blockSize_),
                                                    S(blockSize_,blockSize_);
        {
          Teuchos::TimeMonitor lcltimer( *timerLocalProj_ );
          // project K
          MVT::MvTransMv(ONE,*X_,*KX_,KK);
          // project M
          MVT::MvTransMv(ONE,*X_,*MX_,MM);
          nevLocal_ = blockSize_;
        }

        // solve the projected problem
        {
          Teuchos::TimeMonitor lcltimer( *timerDS_ );
          MSUtils_.directSolver(blockSize_, KK, &MM, &S, &theta_, &nevLocal_, 1);
          TEST_FOR_EXCEPTION(nevLocal_ != blockSize_,LOBPCGInitFailure,
                             "Anasazi::LOBPCG::initialize(): Not enough Ritz vectors to initialize algorithm.");
        }

        // We only have blockSize_ ritz pairs, but we still want them in the correct order
        {
          Teuchos::TimeMonitor lcltimer( *timerSort_ );

          std::vector<int> order(blockSize_);
          // 
          // sort the first blockSize_ values in theta_
          sm_->sort( this, blockSize_, theta_, &order );   // don't catch exception
          //
          // apply the same ordering to the primitive ritz vectors
          MSUtils_.permuteVectors(order,S);
        }

        // compute ritz residual norms
        {
          Teuchos::BLAS<int,ScalarType> blas;
          Teuchos::SerialDenseMatrix<int,ScalarType> R(blockSize_,blockSize_);
          // R = MM*S*diag(theta) - KK*S
          R.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,ONE,MM,S,ZERO);
          for (int i=0; i<blockSize_; i++) {
            blas.SCAL(blockSize_,theta_[i],R[i],1);
          }
          R.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,-ONE,KK,S,ONE);
          for (int i=0; i<blockSize_; i++) {
            ritz2norms_[i] = blas.NRM2(blockSize_,R[i],1);
          }
        }

        // update the solution
        {
          Teuchos::TimeMonitor lcltimer( *timerLocalUpdate_ );
          // X <- X*S
          MVT::MvAddMv( ONE, *X_, ZERO, *X_, *R_ );        
          MVT::MvTimesMatAddMv( ONE, *R_, S, ZERO, *X_ );
          // KX <- KX*S
          MVT::MvAddMv( ONE, *KX_, ZERO, *KX_, *R_ );        
          MVT::MvTimesMatAddMv( ONE, *R_, S, ZERO, *KX_ );
          if (hasM_) {
            // MX <- MX*S
            MVT::MvAddMv( ONE, *MX_, ZERO, *MX_, *R_ );        
            MVT::MvTimesMatAddMv( ONE, *R_, S, ZERO, *MX_ );
          }
        }
      }

  
      // set up R
      if (newstate.R != Teuchos::null && MVT::GetNumberVecs(*newstate.R) >= blockSize_ && MVT::GetVecLength(*newstate.R) == MVT::GetVecLength(*R_) ) {
        MVT::SetBlock(*newstate.R,bsind,*R_);
      }
      else {
        Teuchos::TimeMonitor lcltimer( *timerCompRes_ );
        // form R <- KX - MX*T
        MVT::MvAddMv(ZERO,*KX_,ONE,*KX_,*R_);
        Teuchos::SerialDenseMatrix<int,ScalarType> T(blockSize_,blockSize_);
        T.putScalar(ZERO);
        for (int i=0; i<blockSize_; i++) T(i,i) = theta_[i];
        MVT::MvTimesMatAddMv(-ONE,*MX_,T,ONE,*R_);
      }

      // R has been updated; mark the norms as out-of-date
      Rnorms_current_ = false;
      R2norms_current_ = false;
  
      // put data in P,KP,MP: P is not used to set theta
      if (newstate.P != Teuchos::null && MVT::GetNumberVecs(*newstate.P) >= blockSize_ && MVT::GetVecLength(*newstate.P) == MVT::GetVecLength(*P_) ) {
        hasP_ = true;

        MVT::SetBlock(*newstate.P,bsind,*P_);

        if (newstate.KP != Teuchos::null && MVT::GetNumberVecs(*newstate.KP) >= blockSize_ && MVT::GetVecLength(*newstate.KP) == MVT::GetVecLength(*KP_) ) {
          MVT::SetBlock(*newstate.KP,bsind,*KP_);
        }
        else {
          Teuchos::TimeMonitor lcltimer( *timerOp_ );
          OPT::Apply(*Op_,*P_,*KP_);
          count_ApplyOp_ += blockSize_;
        }

        if (hasM_) {
          if (newstate.MP != Teuchos::null && MVT::GetNumberVecs(*newstate.MP) >= blockSize_ && MVT::GetVecLength(*newstate.MP) == MVT::GetVecLength(*MP_) ) {
            MVT::SetBlock(*newstate.MP,bsind,*MP_);
          }
          else {
            Teuchos::TimeMonitor lcltimer( *timerMOp_ );
            OPT::Apply(*MOp_,*P_,*MP_);
            count_ApplyM_ += blockSize_;
          }
        }
      }

      // finally, we are initialized
      initialized_ = true;

      if (om_->isVerbosity( Debug ) ) {
        // Check almost everything here
        CheckList chk;
        chk.checkX = true;
        chk.checkKX = true;
        chk.checkMX = true;
        chk.checkP = true;
        chk.checkKP = true;
        chk.checkMP = true;
        chk.checkR = true;
        chk.checkQ = true;
        om_->print( Debug, accuracyCheck(chk, ": after initialize()") );
      }

    }
    else {

      LOBPCGState<ScalarType,MV> newstate;
      { // begin timer scope
        Teuchos::TimeMonitor lcltimer( *timerInit_ );

        // generate something, projectAndNormalize, call myself recursively
        Teuchos::RefCountPtr<const MV> ivec = problem_->getInitVec();
        TEST_FOR_EXCEPTION(ivec == Teuchos::null,std::logic_error,
                           "Anasazi::LOBPCG::initialize(): Eigenproblem did not specify initial vectors to clone from");

        int initSize = MVT::GetNumberVecs(*ivec);
        if (initSize > blockSize_) {
          // we need only the first blockSize_ vectors from ivec; get a view of them
          initSize = blockSize_;
          std::vector<int> ind(blockSize_);
          for (int i=0; i<blockSize_; i++) ind[i] = i;
          ivec = MVT::CloneView(*ivec,ind);
        }

        // alloc newX
        Teuchos::RefCountPtr<MV> newMX, newX = MVT::Clone(*ivec,blockSize_);
        // assign ivec to first part of newX
        std::vector<int> ind(initSize);
        if (initSize > 0) {
          for (int i=0; i<initSize; i++) ind[i] = i;
          MVT::SetBlock(*ivec,ind,*newX);
        }
        // fill the rest of newX with random
        if (blockSize_ > initSize) {
          ind.resize(blockSize_ - initSize);
          for (int i=0; i<blockSize_ - initSize; i++) ind[i] = initSize + i;
          Teuchos::RefCountPtr<MV> rX = MVT::CloneView(*newX,ind);
          MVT::MvRandom(*rX);
          rX = Teuchos::null;
        }

        // compute newMX if hasM_
        if (hasM_) {
          newMX = MVT::Clone(*ivec,blockSize_);
          {
            Teuchos::TimeMonitor lcltimer( *timerMOp_ );
            OPT::Apply(*MOp_,*newX,*newMX);
            count_ApplyM_ += blockSize_;
          }
        }
        else {
          newMX = Teuchos::null;
        }

        // remove auxVecs from newX and normalize newX
        if (auxVecs_.size() > 0) {
          Teuchos::TimeMonitor lcltimer( *timerOrtho_ );
          Teuchos::Array<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > > dummy;
          int rank = orthman_->projectAndNormalize(*newX,newMX,dummy,Teuchos::null,auxVecs_);
          TEST_FOR_EXCEPTION(rank != blockSize_,LOBPCGInitFailure,
                             "Anasazi::LOBPCG::initialize(): Couldn't generate initial basis of full rank.");
        }
        else {
          Teuchos::TimeMonitor lcltimer( *timerOrtho_ );
          int rank = orthman_->normalize(*newX,newMX,Teuchos::null);
          TEST_FOR_EXCEPTION(rank != blockSize_,LOBPCGInitFailure,
                             "Anasazi::LOBPCG::initialize(): Couldn't generate initial basis of full rank.");
        }

        // call myself recursively
        newstate.X = newX;
        newstate.MX = newMX;
      } // end of timer scope; we needed this because the following recursive call to initialize contains its own call to timerInit_
      initialize(newstate);
    }
  }

  template <class ScalarType, class MV, class OP>
  void LOBPCG<ScalarType,MV,OP>::initialize()
  {
    LOBPCGState<ScalarType,MV> empty;
    initialize(empty);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Instruct the solver to use full orthogonalization
  template <class ScalarType, class MV, class OP>
  void LOBPCG<ScalarType,MV,OP>::setFullOrtho (bool fullOrtho)
  {
    if ( fullOrtho_ == true || initialized_ == false || fullOrtho == fullOrtho_ ) {
      // state is already orthogonalized or solver is not initialized
      fullOrtho_ = fullOrtho;
      return;
    }

    // solver is initialized, state is not fully orthogonalized, and user has requested full orthogonalization
    fullOrtho_ = true;
    // throw away data in P
    hasP_ = false;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Perform LOBPCG iterations until the StatusTest tells us to stop.
  template <class ScalarType, class MV, class OP>
  void LOBPCG<ScalarType,MV,OP>::iterate () 
  {
    //
    // Allocate/initialize data structures
    //
    if (initialized_ == false) {
      initialize();
    }

    // if fullOrtho_ == true, then we must produce the following on every iteration:
    // [newX newP] = [X H P] [CX;CP]
    // the structure of [CX;CP] when using full orthogonalization does not allow us to 
    // do this in place, and R_ does not have enough storage for newX and newP. therefore, 
    // we must allocate additional storage for this.
    // otherwise, when not using full orthogonalization, the structure
    // [newX newP] = [X H P] [CX1  0 ]
    //                       [CX2 CP2]  allows us to work using only R as work space
    //                       [CX3 CP3] 
    if (fullOrtho_) {
      if (tmpMV_ == Teuchos::null || MVT::GetNumberVecs(*tmpMV_) != blockSize_) {
        tmpMV_ = MVT::Clone(*X_,blockSize_);
      }
    }
    else {
      tmpMV_ = Teuchos::null;
    }

    //
    // Miscellaneous definitions
    const int oneBlock    =   blockSize_;
    const int twoBlocks   = 2*blockSize_;
    const int threeBlocks = 3*blockSize_;
    
    //
    // Define dense projected/local matrices
    //   KK = Local stiffness matrix               (size: 3*blockSize_ x 3*blockSize_)
    //   MM = Local mass matrix                    (size: 3*blockSize_ x 3*blockSize_)
    //    S = Local eigenvectors                   (size: 3*blockSize_ x 3*blockSize_)
    Teuchos::SerialDenseMatrix<int,ScalarType> KK( threeBlocks, threeBlocks ), 
                                               MM( threeBlocks, threeBlocks ),
                                                S( threeBlocks, threeBlocks );

    while (tester_->checkStatus(this) != Passed) {

      // Print information on current status
      if (om_->isVerbosity(Debug)) {
        currentStatus( om_->stream(Debug) );
      }
      else if (om_->isVerbosity(IterationDetails)) {
        currentStatus( om_->stream(IterationDetails) );
      }

      iter_++;
      
      // Apply the preconditioner on the residuals: H <- Prec*R
      if (Prec_ != Teuchos::null) {
        Teuchos::TimeMonitor lcltimer( *timerPrec_ );
        OPT::Apply( *Prec_, *R_, *H_ );   // don't catch the exception
        count_ApplyPrec_ += blockSize_;
      }
      else {
        std::vector<int> ind(blockSize_);
        for (int i=0; i<blockSize_; i++) { ind[i] = i; }
        MVT::SetBlock(*R_,ind,*H_);
      }

      // Apply the mass matrix on H
      if (hasM_) {
        Teuchos::TimeMonitor lcltimer( *timerMOp_ );
        OPT::Apply( *MOp_, *H_, *MH_);    // don't catch the exception
        count_ApplyM_ += blockSize_;
      }

      // orthogonalize H against the auxiliary vectors
      // optionally: orthogonalize H against X and P ([X P] is already orthonormal)
      Teuchos::Array<Teuchos::RefCountPtr<const MV> > Q;
      Teuchos::Array<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > > C = 
        Teuchos::tuple<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > >(Teuchos::null);
      Q = auxVecs_;
      if (fullOrtho_) {
        Q.push_back(X_);
        if (hasP_) {
          Q.push_back(P_);
        }
      }
      {
        Teuchos::TimeMonitor lcltimer( *timerOrtho_ );
        int rank = orthman_->projectAndNormalize(*H_,MH_,C,Teuchos::null,Q);
        TEST_FOR_EXCEPTION(rank != blockSize_,LOBPCGOrthoFailure,
                           "Anasazi::LOBPCG::iterate(): unable to compute full basis for H");
      }

      if (om_->isVerbosity( Debug ) ) {
        CheckList chk;
        chk.checkH = true;
        chk.checkMH = true;
        om_->print( Debug, accuracyCheck(chk, ": after ortho H") );
      }
      else if (om_->isVerbosity( OrthoDetails ) ) {
        CheckList chk;
        chk.checkH = true;
        chk.checkMH = true;
        om_->print( OrthoDetails, accuracyCheck(chk,": after ortho H") );
      }

      // Apply the stiffness matrix to H
      {
        Teuchos::TimeMonitor lcltimer( *timerOp_ );
        OPT::Apply( *Op_, *H_, *KH_);   // don't catch the exception
        count_ApplyOp_ += blockSize_;
      }

      int localSize;
      if (hasP_) {
        localSize = threeBlocks;
      }
      else {
        localSize = twoBlocks;
      }

      // Form "local" mass and stiffness matrices
      {
        Teuchos::TimeMonitor lcltimer( *timerLocalProj_ );
        /* We will construct (block) upper triangular parts only.
                 (X^H)             (KK11 KK12 KK13)
            KK = (H^H) K [X H P] = ( --  KK22 KK23)
                 (P^H)             ( --   --  KK33)
                 (X^H)             (MM11 MM12 MM13)
            MM = (H^H) M [X H P] = ( --  MM22 MM23) 
                 (P^H)             ( --   --  MM33)
        */
        Teuchos::SerialDenseMatrix<int,ScalarType> KK11( Teuchos::View, KK, blockSize_, blockSize_ );
        MVT::MvTransMv( ONE, *X_, *KX_, KK11 );
        Teuchos::SerialDenseMatrix<int,ScalarType> KK12( Teuchos::View, KK, blockSize_, blockSize_, 0, blockSize_ );
        MVT::MvTransMv( ONE, *X_, *KH_, KK12 );
        Teuchos::SerialDenseMatrix<int,ScalarType> KK22( Teuchos::View, KK, blockSize_, blockSize_, blockSize_, blockSize_ );
        MVT::MvTransMv( ONE, *H_, *KH_, KK22 );
        
        Teuchos::SerialDenseMatrix<int,ScalarType> MM11( Teuchos::View, MM, blockSize_, blockSize_ );
        MVT::MvTransMv( ONE, *X_, *MX_, MM11 );
        Teuchos::SerialDenseMatrix<int,ScalarType> MM12( Teuchos::View, MM, blockSize_, blockSize_, 0, blockSize_ );
        MVT::MvTransMv( ONE, *X_, *MH_, MM12 );
        Teuchos::SerialDenseMatrix<int,ScalarType> MM22( Teuchos::View, MM, blockSize_, blockSize_, blockSize_, blockSize_ );
        MVT::MvTransMv( ONE, *H_, *MH_, MM22 );

        if (hasP_) {
          Teuchos::SerialDenseMatrix<int,ScalarType> KK13( Teuchos::View, KK, blockSize_, blockSize_, 0, twoBlocks );
          MVT::MvTransMv( ONE, *X_, *KP_, KK13 );
          Teuchos::SerialDenseMatrix<int,ScalarType> KK23( Teuchos::View, KK, blockSize_, blockSize_, blockSize_, twoBlocks );
          MVT::MvTransMv( ONE, *H_, *KP_, KK23 );
          Teuchos::SerialDenseMatrix<int,ScalarType> KK33( Teuchos::View, KK, blockSize_, blockSize_, twoBlocks, twoBlocks );
          MVT::MvTransMv( ONE, *P_, *KP_, KK33 );
          
          Teuchos::SerialDenseMatrix<int,ScalarType> MM13( Teuchos::View, MM, blockSize_, blockSize_, 0, twoBlocks );
          MVT::MvTransMv( ONE, *X_, *MP_, MM13 );
          Teuchos::SerialDenseMatrix<int,ScalarType> MM23( Teuchos::View, MM, blockSize_, blockSize_, blockSize_, twoBlocks );
          MVT::MvTransMv( ONE, *H_, *MP_, MM23 );
          Teuchos::SerialDenseMatrix<int,ScalarType> MM33( Teuchos::View, MM, blockSize_, blockSize_, twoBlocks, twoBlocks );
          MVT::MvTransMv( ONE, *P_, *MP_, MM33 );
        }

        // make MM and KK symmetric in memory: we need this so we can form ritz residual vectors and use MM for inner product below (with ease)
        for (int j=0; j<localSize; j++) {
          for (int i=j+1; i<localSize; i++) {
            KK(i,j) = KK(j,i);
            MM(i,j) = MM(j,i);
          }
        }
      } // end timing block

      // Perform a spectral decomposition
      {
        Teuchos::TimeMonitor lcltimer( *timerDS_ );
        nevLocal_ = localSize;
        MSUtils_.directSolver(localSize, KK, &MM, &S, &theta_, &nevLocal_, 0);    // don't catch the exception
        // localSize tells directSolver() how big KK,MM are
        // however, directSolver() may choose to use only the principle submatrices of KK,MM because of loss of 
        // MM-orthogonality in the projected eigenvectors
        // nevLocal_ tells us how much it used, in effect dictating back to us how big localSize is, 
        // and therefore telling us which of [X H P] to use in computing the new iterates below.
        // we will not tolerate this ill-conditioning, and will throw an exception.
      }
      om_->stream(Debug) << " After directSolve: localSize == " << localSize << " \tnevLocal == " << nevLocal_ << endl;
      TEST_FOR_EXCEPTION(nevLocal_ != localSize, LOBPCGRitzFailure, "Ill-conditioning detected in projecteded mass matrix." );

      Teuchos::LAPACK<int,ScalarType> lapack;
      Teuchos::BLAS<int,ScalarType> blas;
      //
      //---------------------------------------------------
      // Sort the ritz values using the sort manager
      //---------------------------------------------------
      {
        Teuchos::TimeMonitor lcltimer( *timerSort_ );

        std::vector<int> order(nevLocal_);
        // 
        // sort the first nevLocal_ values in theta_
        sm_->sort( this, nevLocal_, theta_, &order );   // don't catch exception
        //
        // Sort the primitive ritz vectors
        Teuchos::SerialDenseMatrix<int,ScalarType> curS(Teuchos::View,S,nevLocal_,nevLocal_);
        MSUtils_.permuteVectors(order,curS);
      }

      // compute ritz residual norms
      {
        Teuchos::SerialDenseMatrix<int,ScalarType> R(nevLocal_,nevLocal_);
        Teuchos::SerialDenseMatrix<int,ScalarType> lclKK(Teuchos::View,KK,nevLocal_,nevLocal_),
                                                   lclMM(Teuchos::View,MM,nevLocal_,nevLocal_);
        // R = MM*S*diag(theta) - KK*S
        R.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,ONE,lclMM,S,ZERO);
        for (int i=0; i<nevLocal_; i++) {
          blas.SCAL(nevLocal_,theta_[i],R[i],1);
        }
        R.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,-ONE,lclKK,S,ONE);
        for (int i=0; i<nevLocal_; i++) {
          ritz2norms_[i] = blas.NRM2(nevLocal_,R[i],1);
        }
      }


      // before computing X,P: perform second orthogonalization per Ulrich,Rich paper
      // CX will be the coefficients of [X,H,P] for new X, CP for new P
      // The paper suggests orthogonalizing CP against CX and orthonormalizing CP, w.r.t. MM
      // Here, we will also orthonormalize CX.
      // This is accomplished using the Cholesky factorization of [CX CP]^H MM [CX CP]
      Teuchos::SerialDenseMatrix<int,ScalarType> CX(threeBlocks,oneBlock), CP(threeBlocks,oneBlock);
      if (fullOrtho_ && localSize >= twoBlocks) {
        // build orthonormal basis for (  0  ) that is orthogonal to ( S11 )
        //                             ( S21 )                       ( S21 )
        //                             ( S31 )                       ( S31 )
        // Do this using Cholesky factorization of ( S11  0  )
        //                                         ( S21 S21 )
        //                                         ( S31 S31 )
        //           ( S11  0  )
        // Build C = ( S21 S21 )
        //           ( S31 S31 )
        Teuchos::SerialDenseMatrix<int,ScalarType> C(threeBlocks,twoBlocks),
                                                tmp1(threeBlocks,twoBlocks),
                                                tmp2(twoBlocks  ,twoBlocks);

        // first block of rows
        for (int i=0; i<oneBlock; i++) {
          // CX
          for (int j=0; j<oneBlock; j++) {
            C(i,j) = S(i,j);
          }
          // CP
          for (int j=oneBlock; j<twoBlocks; j++) {
            C(i,j) = ZERO;
          }
        }
        // second block of rows
        for (int j=0; j<oneBlock; j++) {
          for (int i=oneBlock; i<twoBlocks; i++) {
            // CX
            C(i,j)          = S(i,j);
            // CP
            C(i,j+oneBlock) = S(i,j);
          }
        }
        // third block of rows
        if (localSize == threeBlocks) {
          for (int j=0; j<oneBlock; j++) {
            for (int i=twoBlocks; i<threeBlocks; i++) {
              // CX
              C(i,j)          = S(i,j);
              // CP
              C(i,j+oneBlock) = S(i,j);
            }
          }
        }
        else {
          for (int j=0; j<twoBlocks; j++) {
            for (int i=twoBlocks; i<threeBlocks; i++) {
              C(i,j) = ZERO;
            }
          }
        }

        // compute tmp1 = MM*C
        int teuchosret;
        teuchosret = tmp1.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,ONE,MM,C,ZERO);
        TEST_FOR_EXCEPTION(teuchosret != 0,std::logic_error,"Logic error calling SerialDenseMatrix::multiply");

        // compute tmp2 = C^H*tmp1 == C^H*MM*C
        teuchosret = tmp2.multiply(Teuchos::CONJ_TRANS,Teuchos::NO_TRANS,ONE,C,tmp1,ZERO);
        TEST_FOR_EXCEPTION(teuchosret != 0,std::logic_error,"Logic error calling SerialDenseMatrix::multiply");

        // compute R (cholesky) of tmp2
        int info;
        lapack.POTRF('U',twoBlocks,tmp2.values(),tmp2.stride(),&info);
        TEST_FOR_EXCEPTION(info != 0, LOBPCGOrthoFailure, 
                           "Anasazi::LOBPCG::iterate(): Could not perform full orthogonalization.");
        // compute C = C inv(R)
        blas.TRSM(Teuchos::RIGHT_SIDE,Teuchos::UPPER_TRI,Teuchos::NO_TRANS,Teuchos::NON_UNIT_DIAG,
                  threeBlocks,twoBlocks,ONE,tmp2.values(),tmp2.stride(),C.values(),C.stride());
        // put C(:,0:oneBlock-1) into CX
        for (int j=0; j<oneBlock; j++) {
          for (int i=0; i<threeBlocks; i++) {
            CX(i,j) = C(i,j);
          }
        }
        // put C(:,oneBlock:twoBlocks-1) into CP
        for (int j=0; j<oneBlock; j++) {
          for (int i=0; i<threeBlocks; i++) {
            CP(i,j) = C(i,oneBlock+j);
          }
        }

        // check the results
        if (om_->isVerbosity( Debug ) ) {
          Teuchos::SerialDenseMatrix<int,ScalarType> tmp1(threeBlocks,oneBlock),
                                                     tmp2(oneBlock,oneBlock);
          MagnitudeType tmp;
          int teuchosret;
          stringstream os;
          os.precision(2);
          os.setf(ios::scientific, ios::floatfield);

          os << " Checking Full Ortho: iteration " << iter_ << endl;

          // check CX^T MM CX == I
          // compute tmp1 = MM*CX
          teuchosret = tmp1.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,ONE,MM,CX,ZERO);
          TEST_FOR_EXCEPTION(teuchosret != 0,std::logic_error,"Logic error calling SerialDenseMatrix::multiply");
          // compute tmp2 = CX^H*tmp1 == CX^H*MM*CX
          teuchosret = tmp2.multiply(Teuchos::CONJ_TRANS,Teuchos::NO_TRANS,ONE,CX,tmp1,ZERO);
          TEST_FOR_EXCEPTION(teuchosret != 0,std::logic_error,"Logic error calling SerialDenseMatrix::multiply");
          // subtrace tmp2 - I == CX^H * MM * CX - I
          for (int i=0; i<oneBlock; i++) tmp2(i,i) -= ONE;
          tmp = tmp2.normFrobenius();          
          os << " >> Error in CX^H MM CX == I : " << tmp << endl;

          // check CP^T MM CP == I
          // compute tmp1 = MM*CP
          teuchosret = tmp1.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,ONE,MM,CP,ZERO);
          TEST_FOR_EXCEPTION(teuchosret != 0,std::logic_error,"Logic error calling SerialDenseMatrix::multiply");
          // compute tmp2 = CP^H*tmp1 == CP^H*MM*CP
          teuchosret = tmp2.multiply(Teuchos::CONJ_TRANS,Teuchos::NO_TRANS,ONE,CP,tmp1,ZERO);
          TEST_FOR_EXCEPTION(teuchosret != 0,std::logic_error,"Logic error calling SerialDenseMatrix::multiply");
          // subtrace tmp2 - I == CP^H * MM * CP - I
          for (int i=0; i<oneBlock; i++) tmp2(i,i) -= ONE;
          tmp = tmp2.normFrobenius();          
          os << " >> Error in CP^H MM CP == I : " << tmp << endl;

          // check CX^T MM CP == 0
          // compute tmp1 = MM*CP
          teuchosret = tmp1.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,ONE,MM,CP,ZERO);
          TEST_FOR_EXCEPTION(teuchosret != 0,std::logic_error,"Logic error calling SerialDenseMatrix::multiply");
          // compute tmp2 = CX^H*tmp1 == CX^H*MM*CP
          teuchosret = tmp2.multiply(Teuchos::CONJ_TRANS,Teuchos::NO_TRANS,ONE,CX,tmp1,ZERO);
          TEST_FOR_EXCEPTION(teuchosret != 0,std::logic_error,"Logic error calling SerialDenseMatrix::multiply");
          // subtrace tmp2 == CX^H * MM * CP
          tmp = tmp2.normFrobenius();          
          os << " >> Error in CX^H MM CP == 0 : " << tmp << endl;

          os << endl;
          om_->print(Debug,os.str());
        }
      }
      else {
        //      [S11]
        // S1 = [S21]
        //      [S31]
        //
        // CX = S1
        //
        //      [ 0 ]
        // CP = [S21]
        //      [S31]
        //
        Teuchos::SerialDenseMatrix<int,ScalarType> S1(Teuchos::View,S,threeBlocks,oneBlock);
        CX.assign(S1);
        CP.assign(S1);
        // remove first block of rows from CP
        for (int j=0; j<oneBlock; j++) {
          for (int i=0; i<oneBlock; i++) {
            CP(i,j) = ZERO;
          }
        }
      }

      //
      // Update the spaces: compute new X and new P
      // P is only computed if we have localSize >= twoBlocks
      // Note: Use R as a temporary work space and (if full ortho) tmpMV
      hasP_ = false;
      {
        Teuchos::TimeMonitor lcltimer( *timerLocalUpdate_ );

        TEST_FOR_EXCEPTION(localSize==oneBlock,std::logic_error,"Anasazi::LOBPCG::iterate(): Logic error in local update.");
        if (localSize == twoBlocks) {
          Teuchos::SerialDenseMatrix<int,ScalarType> CX1( Teuchos::View, CX, blockSize_, blockSize_ ),
                                                     CX2( Teuchos::View, CX, blockSize_, blockSize_, oneBlock, 0 ),
                                                     CP1( Teuchos::View, CP, blockSize_, blockSize_ ),
                                                     CP2( Teuchos::View, CP, blockSize_, blockSize_, oneBlock, 0 );

          hasP_ = true;

          // X = [X H][CX1]
          //          [CX2]
          //
          // P = [X H][CP1]
          //          [CP2]
          //

          // R <- X
          MVT::MvAddMv( ONE, *X_, ZERO, *X_, *R_ );        
          // X <- R*CX1 + H*CX2 == X*CX1 + H*CX2
          MVT::MvTimesMatAddMv( ONE, *R_, CX1, ZERO, *X_ );
          MVT::MvTimesMatAddMv( ONE, *H_, CX2, ONE , *X_ );
          // P <- R*CP1 + H*CP2 == X*CP1 + H*CP2
          // however, if fullOrtho_ == false, CP1 == ZERO
          if (fullOrtho_) { 
            MVT::MvTimesMatAddMv( ONE, *R_, CP1, ZERO, *P_ );
            MVT::MvTimesMatAddMv( ONE, *H_, CP2, ONE, *P_ );
          }
          else {
            MVT::MvTimesMatAddMv( ONE, *H_, CP2, ZERO, *P_ );
          }

          // R  <- KX
          MVT::MvAddMv( ONE, *KX_, ZERO, *KX_, *R_ );        
          // KX <- R*CX1 + KH*CX2 == KX*CX1 + KH*CX2
          MVT::MvTimesMatAddMv( ONE, *R_, CX1, ZERO, *KX_ );
          MVT::MvTimesMatAddMv( ONE, *KH_, CX2, ONE , *KX_ );
          // KP <- R*CP1 + KH*CP2 == KX*CP1 + KH*CP2
          // however, if fullOrtho_ == false, CP1 == ZERO
          if (fullOrtho_) { 
            MVT::MvTimesMatAddMv( ONE, *R_, CP1, ZERO, *KP_ );
            MVT::MvTimesMatAddMv( ONE, *KH_, CP2, ONE, *KP_ );
          }
          else {
            MVT::MvTimesMatAddMv( ONE, *KH_, CP2, ZERO, *KP_ );
          }

          if (hasM_) {
            // R  <- MX
            MVT::MvAddMv( ONE, *MX_, ZERO, *MX_, *R_ );        
            // MX <- R*CX1 + MH*CX2 == MX*CX1 + MH*CX2
            MVT::MvTimesMatAddMv( ONE, *R_, CX1, ZERO, *MX_ );
            MVT::MvTimesMatAddMv( ONE, *MH_, CX2, ONE , *MX_ );
            // MP <- R*CP1 + MH*CP2 == MX*CP1 + MH*CP2
            // however, if fullOrtho_ == false, CP1 == ZERO
            if (fullOrtho_) { 
              MVT::MvTimesMatAddMv( ONE, *R_, CP1, ZERO, *MP_ );
              MVT::MvTimesMatAddMv( ONE, *MH_, CP2, ONE, *MP_ );
            }
            else {
              MVT::MvTimesMatAddMv( ONE, *MH_, CP2, ZERO, *MP_ );
            }
          }

        } // if (localSize == twoBlocks)
        else if (localSize == threeBlocks) {
          Teuchos::SerialDenseMatrix<int,ScalarType> CX1( Teuchos::View, CX, blockSize_, blockSize_ ),
                                                     CX2( Teuchos::View, CX, blockSize_, blockSize_, oneBlock, 0 ),
                                                     CX3( Teuchos::View, CX, blockSize_, blockSize_, twoBlocks, 0 ),
                                                     CP1( Teuchos::View, CP, blockSize_, blockSize_ ),
                                                     CP2( Teuchos::View, CP, blockSize_, blockSize_, oneBlock, 0 ),
                                                     CP3( Teuchos::View, CP, blockSize_, blockSize_, twoBlocks, 0 );

          hasP_ = true;

          // X <- X*CX1 + P*CX3
          // P <- X*CP1 + P*CP3 (note, CP1 == ZERO if fullOrtho_==false)
          if (fullOrtho_) {
            // copy X,P
            MVT::MvAddMv(ONE,*X_, ZERO,*X_, *R_);
            MVT::MvAddMv(ONE,*P_, ZERO,*P_, *tmpMV_);
            // perform [X P][CX1 CP1]
            //              [CX3 CP3]
            MVT::MvTimesMatAddMv( ONE, *R_,     CX1, ZERO, *X_ );
            MVT::MvTimesMatAddMv( ONE, *tmpMV_, CX3,  ONE, *X_ );
            MVT::MvTimesMatAddMv( ONE, *R_,     CP1, ZERO, *P_ );
            MVT::MvTimesMatAddMv( ONE, *tmpMV_, CP3,  ONE, *P_ );
          }
          else {
            // copy X
            MVT::MvAddMv(ONE,*X_, ZERO,*X_, *R_);
            // perform [X P][CX1  0 ]
            //              [CX3 CP3]
            MVT::MvTimesMatAddMv( ONE, *R_, CX1, ZERO, *X_ );
            MVT::MvTimesMatAddMv( ONE, *P_, CX3,  ONE, *X_ );
            // copy P
            MVT::MvAddMv(ONE,*P_, ZERO,*P_, *R_);
            MVT::MvTimesMatAddMv( ONE, *R_, CP3, ZERO, *P_ );
          }
          // X <- X + H*CX2
          // P <- P + H*CP2
          MVT::MvTimesMatAddMv( ONE, *H_, CX2,  ONE, *X_ );
          MVT::MvTimesMatAddMv( ONE, *H_, CP2,  ONE, *P_ );


          // KX <- KX*CX1 + KP*CX3
          // KP <- KX*CP1 + KP*CP3 (note, CP1 == ZERO if fullOrtho_==false)
          if (fullOrtho_) {
            // copy KX,KP
            MVT::MvAddMv(ONE,*KX_, ZERO,*KX_, *R_);
            MVT::MvAddMv(ONE,*KP_, ZERO,*KP_, *tmpMV_);
            // perform [KX KP][CX1 CP1]
            //                [CX3 CP3]
            MVT::MvTimesMatAddMv( ONE, *R_,     CX1, ZERO, *KX_ );
            MVT::MvTimesMatAddMv( ONE, *tmpMV_, CX3,  ONE, *KX_ );
            MVT::MvTimesMatAddMv( ONE, *R_,     CP1, ZERO, *KP_ );
            MVT::MvTimesMatAddMv( ONE, *tmpMV_, CP3,  ONE, *KP_ );
          }
          else {
            // copy KX
            MVT::MvAddMv(ONE,*KX_, ZERO,*KX_, *R_);
            // perform [KX KP][CX1  0 ]
            //                [CX3 CP3]
            MVT::MvTimesMatAddMv( ONE, *R_ , CX1, ZERO, *KX_ );
            MVT::MvTimesMatAddMv( ONE, *KP_, CX3,  ONE, *KX_ );
            // copy KP
            MVT::MvAddMv(ONE,*KP_, ZERO,*KP_, *R_);
            MVT::MvTimesMatAddMv( ONE, *R_, CP3, ZERO, *KP_ );
          }
          // KX <- KX + KH*CX2
          // KP <- KP + KH*CP2
          MVT::MvTimesMatAddMv( ONE, *KH_, CX2,  ONE, *KX_ );
          MVT::MvTimesMatAddMv( ONE, *KH_, CP2,  ONE, *KP_ );


          if (hasM_) {
            // MX <- MX*CX1 + MP*CX3
            // MP <- MX*CP1 + MP*CP3 (note, CP1 == ZERO if fullOrtho_==false)
            if (fullOrtho_) {
              // copy MX,MP
              MVT::MvAddMv(ONE,*MX_, ZERO,*MX_, *R_);
              MVT::MvAddMv(ONE,*MP_, ZERO,*MP_, *tmpMV_);
              // perform [MX MP][CX1 CP1]
              //                [CX3 CP3]
              MVT::MvTimesMatAddMv( ONE, *R_,     CX1, ZERO, *MX_ );
              MVT::MvTimesMatAddMv( ONE, *tmpMV_, CX3,  ONE, *MX_ );
              MVT::MvTimesMatAddMv( ONE, *R_,     CP1, ZERO, *MP_ );
              MVT::MvTimesMatAddMv( ONE, *tmpMV_, CP3,  ONE, *MP_ );
            }
            else {
              // copy MX
              MVT::MvAddMv(ONE,*MX_, ZERO,*MX_, *R_);
              // perform [MX MP][CX1  0 ]
              //                [CX3 CP3]
              MVT::MvTimesMatAddMv( ONE, *R_ , CX1, ZERO, *MX_ );
              MVT::MvTimesMatAddMv( ONE, *MP_, CX3,  ONE, *MX_ );
              // copy MP
              MVT::MvAddMv(ONE,*MP_, ZERO,*MP_, *R_);
              MVT::MvTimesMatAddMv( ONE, *R_, CP3, ZERO, *MP_ );
            }
            // MX <- MX + MH*CX2
            // MP <- MP + MH*CP2
            MVT::MvTimesMatAddMv( ONE, *MH_, CX2,  ONE, *MX_ );
            MVT::MvTimesMatAddMv( ONE, *MH_, CP2,  ONE, *MP_ );
          }

        } // if (localSize == threeBlocks)
      } // end timing block

      // Compute the new residuals, explicitly
      {
        Teuchos::TimeMonitor lcltimer( *timerCompRes_ );
        MVT::MvAddMv( ONE, *KX_, ZERO, *KX_, *R_ );
        Teuchos::SerialDenseMatrix<int,ScalarType> T( blockSize_, blockSize_ );
        for (int i = 0; i < blockSize_; i++) {
          T(i,i) = theta_[i];
        }
        MVT::MvTimesMatAddMv( -ONE, *MX_, T, ONE, *R_ );
      }

      // R has been updated; mark the norms as out-of-date
      Rnorms_current_ = false;
      R2norms_current_ = false;


      // When required, monitor some orthogonalities
      if (om_->isVerbosity( Debug ) ) {
        // Check almost everything here
        CheckList chk;
        chk.checkX = true;
        chk.checkKX = true;
        chk.checkMX = true;
        chk.checkP = true;
        chk.checkKP = true;
        chk.checkMP = true;
        chk.checkR = true;
        om_->print( Debug, accuracyCheck(chk, ": after local update") );
      }
      else if (om_->isVerbosity( OrthoDetails )) {
        CheckList chk;
        chk.checkX = true;
        chk.checkP = true;
        chk.checkR = true;
        om_->print( OrthoDetails, accuracyCheck(chk, ": after local update") );
      }
    } // end while (statusTest == false)
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // compute/return residual M-norms
  template <class ScalarType, class MV, class OP>
  std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> 
  LOBPCG<ScalarType,MV,OP>::getResNorms() {
    if (Rnorms_current_ == false) {
      // Update the residual norms
      orthman_->norm(*R_,&Rnorms_);
      Rnorms_current_ = true;
    }
    return Rnorms_;
  }

  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // compute/return residual 2-norms
  template <class ScalarType, class MV, class OP>
  std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> 
  LOBPCG<ScalarType,MV,OP>::getRes2Norms() {
    if (R2norms_current_ == false) {
      // Update the residual 2-norms 
      MVT::MvNorm(*R_,&R2norms_);
      R2norms_current_ = true;
    }
    return R2norms_;
  }




  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Check accuracy, orthogonality, and other debugging stuff
  // 
  // bools specify which tests we want to run (instead of running more than we actually care about)
  //
  // we don't bother checking the following because they are computed explicitly:
  //    H == Prec*R
  //   KH == K*H
  //
  // 
  // checkX : X orthonormal
  //          orthogonal to auxvecs
  // checkMX: check MX == M*X
  // checkKX: check KX == K*X
  // checkP : if fullortho P orthonormal and orthogonal to X
  //          orthogonal to auxvecs
  // checkMP: check MP == M*P
  // checkKP: check KP == K*P
  // checkH : if fullortho H orthonormal and orthogonal to X and P
  //          orthogonal to auxvecs
  // checkMH: check MH == M*H
  // checkR : check R orthogonal to X
  // checkQ : check that auxiliary vectors are actually orthonormal
  //
  // TODO: 
  //  add checkTheta 
  //
  template <class ScalarType, class MV, class OP>
  std::string LOBPCG<ScalarType,MV,OP>::accuracyCheck( const CheckList &chk, const string &where ) const 
  {
    stringstream os;
    os.precision(2);
    os.setf(ios::scientific, ios::floatfield);
    MagnitudeType tmp;

    os << " Debugging checks: iteration " << iter_ << where << endl;

    // X and friends
    if (chk.checkX && initialized_) {
      tmp = orthman_->orthonormError(*X_);
      os << " >> Error in X^H M X == I : " << tmp << endl;
      for (unsigned int i=0; i<auxVecs_.size(); i++) {
        tmp = orthman_->orthogError(*X_,*auxVecs_[i]);
        os << " >> Error in X^H M Q[" << i << "] == 0 : " << tmp << endl;
      }
    }
    if (chk.checkMX && hasM_ && initialized_) {
      tmp = MSUtils_.errorEquality(X_.get(), MX_.get(), MOp_.get());
      os << " >> Error in MX == M*X    : " << tmp << endl;
    }
    if (chk.checkKX && initialized_) {
      tmp = MSUtils_.errorEquality(X_.get(), KX_.get(), Op_.get());
      os << " >> Error in KX == K*X    : " << tmp << endl;
    }

    // P and friends
    if (chk.checkP && hasP_ && initialized_) {
      if (fullOrtho_) {
        tmp = orthman_->orthonormError(*P_);
        os << " >> Error in P^H M P == I : " << tmp << endl;
        tmp = orthman_->orthogError(*P_,*X_);
        os << " >> Error in P^H M X == 0 : " << tmp << endl;
      }
      for (unsigned int i=0; i<auxVecs_.size(); i++) {
        tmp = orthman_->orthogError(*P_,*auxVecs_[i]);
        os << " >> Error in P^H M Q[" << i << "] == 0 : " << tmp << endl;
      }
    }
    if (chk.checkMP && hasM_ && hasP_ && initialized_) {
      tmp = MSUtils_.errorEquality(P_.get(), MP_.get(), MOp_.get());
      os << " >> Error in MP == M*P    : " << tmp << endl;
    }
    if (chk.checkKP && hasP_ && initialized_) {
      tmp = MSUtils_.errorEquality(P_.get(), KP_.get(), Op_.get());
      os << " >> Error in KP == K*P    : " << tmp << endl;
    }

    // H and friends
    if (chk.checkH && initialized_) {
      if (fullOrtho_) {
        tmp = orthman_->orthonormError(*H_);
        os << " >> Error in H^H M H == I : " << tmp << endl;
        tmp = orthman_->orthogError(*H_,*X_);
        os << " >> Error in H^H M X == 0 : " << tmp << endl;
        if (hasP_) {
          tmp = orthman_->orthogError(*H_,*P_);
          os << " >> Error in H^H M P == 0 : " << tmp << endl;
        }
      }
      for (unsigned int i=0; i<auxVecs_.size(); i++) {
        tmp = orthman_->orthogError(*H_,*auxVecs_[i]);
        os << " >> Error in H^H M Q[" << i << "] == 0 : " << tmp << endl;
      }
    }
    if (chk.checkMH && hasM_ && initialized_) {
      tmp = MSUtils_.errorEquality(H_.get(), MH_.get(), MOp_.get());
      os << " >> Error in MH == M*H    : " << tmp << endl;
    }

    // R: this is not M-orthogonality, but standard euclidean orthogonality
    if (chk.checkR && initialized_) {
      Teuchos::SerialDenseMatrix<int,ScalarType> xTx(blockSize_,blockSize_);
      MVT::MvTransMv(ONE,*X_,*R_,xTx);
      tmp = xTx.normFrobenius();
      os << " >> Error in X^H R == 0   : " << tmp << endl;
    }

    // Q
    if (chk.checkQ) {
      for (unsigned int i=0; i<auxVecs_.size(); i++) {
        tmp = orthman_->orthonormError(*auxVecs_[i]);
        os << " >> Error in Q[" << i << "]^H M Q[" << i << "] == I : " << tmp << endl;
        for (unsigned int j=i+1; j<auxVecs_.size(); j++) {
          tmp = orthman_->orthogError(*auxVecs_[i],*auxVecs_[j]);
          os << " >> Error in Q[" << i << "]^H M Q[" << j << "] == 0 : " << tmp << endl;
        }
      }
    }

    os << endl;

    return os.str();
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Print the current status of the solver
  template <class ScalarType, class MV, class OP>
  void 
  LOBPCG<ScalarType,MV,OP>::currentStatus(ostream &os) 
  {
    os.setf(ios::scientific, ios::floatfield);  
    os.precision(6);
    os <<endl;
    os <<"================================================================================" << endl;
    os << endl;
    os <<"                              LOBPCG Solver Status" << endl;
    os << endl;
    os <<"The solver is "<<(initialized_ ? "initialized." : "not initialized.") << endl;
    os <<"The number of iterations performed is " << iter_       << endl;
    os <<"The current block size is             " << blockSize_  << endl;
    os <<"The number of auxiliary vectors is    " << numAuxVecs_ << endl;
    os <<"The number of operations Op*x   is " << count_ApplyOp_   << endl;
    os <<"The number of operations M*x    is " << count_ApplyM_    << endl;
    os <<"The number of operations Prec*x is " << count_ApplyPrec_ << endl;

    os.setf(ios_base::right, ios_base::adjustfield);

    if (initialized_) {
      os << endl;
      os <<"CURRENT EIGENVALUE ESTIMATES             "<<endl;
      os << std::setw(20) << "Eigenvalue" 
         << std::setw(20) << "Residual(M)"
         << std::setw(20) << "Residual(2)"
         << endl;
      os <<"--------------------------------------------------------------------------------"<<endl;
      for (int i=0; i<blockSize_; i++) {
        os << std::setw(20) << theta_[i];
        if (Rnorms_current_) os << std::setw(20) << Rnorms_[i];
        else os << std::setw(20) << "not current";
        if (R2norms_current_) os << std::setw(20) << R2norms_[i];
        else os << std::setw(20) << "not current";
        os << endl;
      }
    }
    os <<"================================================================================" << endl;
    os << endl;
  }

  
} // end Anasazi namespace

#endif // ANASAZI_LOBPCG_HPP
