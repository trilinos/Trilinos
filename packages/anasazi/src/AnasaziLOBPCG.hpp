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
    LOBPCG contains local storage of up to 10*_blockSize vectors, representing 10 entities
      X,H,P,R
      KX,KH,KP  (product of K and the above)
      MX,MH,MP  (product of M and the above, not allocated if we don't have an M matrix)
    If full orthogonalization is enabled, one extra multivector of _blockSize vectors is required to 
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
        

        \ingroup anasazi_solvers

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
     *   - getEigenvalues() returns the Ritz values with respect to X
     *   - getResidualVecs() returns the residual vectors with respect to X
     *   - If hasP() == \c true,
     *      - P orthogonal to auxiliary vectors
     *      - If getFullOrtho() == \c true,
     *        - P is orthogonal to X and has orthonormal columns
     *      - KP == Op*P
     *      - MP == M*P if M != Teuchos::null\n
     *        Otherwise, MP == Teuchos::null
     */
    bool isInitialized() { return _initialized; }

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
      state.X = _X;
      state.KX = _KX;
      state.P = _P;
      state.KP = _KP;
      state.H = _H;
      state.KH = _KH;
      state.R = _R;
      state.T = Teuchos::rcp(new std::vector<MagnitudeType>(_theta));
      if (_hasM) {
        state.MX = _MX;
        state.MP = _MP;
        state.MH = _MH;
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
    int getNumIters() const { return(_iter); };

    //! \brief Reset the iteration count.
    void resetNumIters() { _iter=0; };

    //! \brief Get the current approximate eigenvectors.
    Teuchos::RefCountPtr<const MV> getEvecs() {return _X;}

    //! \brief Get the residual vectors.
    Teuchos::RefCountPtr<const MV> getResidualVecs() {return _R;}

    /*! \brief Get the current eigenvalue estimates.
     *
     *  \return A vector of length getBlockSize() containing the eigenvalue
     *  estimates associated with the current iterate.
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getEigenvalues() { 
      std::vector<MagnitudeType> ret = _theta;
      ret.resize(_blockSize);
      return ret;
    }

    /*! \brief Get the Ritz values from the previous iteration.
     *
     *  \return A vector of length getCurSubspaceDim() containing the Ritz values from the
     *  previous projected eigensolve.
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getRitzValues() { 
      std::vector<MagnitudeType> ret = _theta;
      ret.resize(_nevLocal);
      return ret;
    }

    /*! \brief Get the current residual norms
     *
     *  \return A vector of length getCurSubspaceDim() containing the norms of the
     *  residuals, with respect to the orthogonalization manager norm() method.
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getResNorms()    {return _Rnorms;}


    /*! \brief Get the current residual 2-norms
     *
     *  \return A vector of length getCurSubspaceDim() containing the 2-norms of the
     *  residuals. 
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getRes2Norms()   {return _R2norms;}


    /*! \brief Get the 2-norms of the Ritz residuals.
     *
     *  \return A vector of length getCurSubspaceDim() containing the 2-norms of the Ritz residuals.
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getRitzRes2Norms() {
      std::vector<MagnitudeType> ret = _ritz2norms;
      ret.resize(_nevLocal);
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
    int getCurSubspaceDim() {
      if (!_initialized) return 0;
      return _nevLocal;
    }

    /*! \brief Get the maximum dimension allocated for the search subspace. For %LOBPCG, this always returns 3*getBlockSize(), the dimension of the 
     *   subspace colspan([X H P]).
     */
    int getMaxSubspaceDim() {return 3*_blockSize;}

    //@}

    //!  @name Accessor routines from Eigensolver
    //@{


    //! Get a constant reference to the eigenvalue problem.
    const Eigenproblem<ScalarType,MV,OP>& getProblem() const { return(*_problem); };


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
    int getBlockSize() const { return(_blockSize); }


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
    Teuchos::Array<Teuchos::RefCountPtr<const MV> > getAuxVecs() const {return _auxVecs;}

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
    bool getFullOrtho() const { return(_fullOrtho); }
    
    //! Indicates whether the search direction given by getState() is valid.
    bool hasP() {return _hasP;}

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
    const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> >     _problem;
    const Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> >      _sm;
    const Teuchos::RefCountPtr<OutputManager<ScalarType> >          _om;
    const Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> >       _tester;
    const Teuchos::RefCountPtr<MatOrthoManager<ScalarType,MV,OP> >  _orthman;
    //
    // Information obtained from the eigenproblem
    //
    Teuchos::RefCountPtr<OP> _Op;
    Teuchos::RefCountPtr<OP> _MOp;
    Teuchos::RefCountPtr<OP> _Prec;
    bool _hasM;
    //
    // Internal utilities class required by eigensolver.
    //
    ModalSolverUtils<ScalarType,MV,OP> _MSUtils;
    //
    // Internal timers
    //
    Teuchos::RefCountPtr<Teuchos::Time> _timerOp, _timerMOp, _timerPrec,
                                        _timerSort, 
                                        _timerLocalProj, _timerDS,
                                        _timerLocalUpdate, _timerCompRes,
                                        _timerOrtho, _timerInit;
    //
    // Counters
    //
    // Number of operator applications
    int _count_ApplyOp, _count_ApplyM, _count_ApplyPrec;

    //
    // Algorithmic parameters.
    //
    // _blockSize is the solver block size
    int _blockSize;
    //
    // _fullOrtho dictates whether the orthogonalization procedures specified by Hetmaniuk and Lehoucq should
    // be activated (see citations at the top of this file)
    bool _fullOrtho;

    //
    // Current solver state
    //
    // _initialized specifies that the basis vectors have been initialized and the iterate() routine
    // is capable of running; _initialize is controlled  by the initialize() member method
    // For the implications of the state of _initialized, please see documentation for initialize()
    bool _initialized;
    //
    // _nevLocal reflects how much of the current basis is valid (0 <= _nevLocal <= 3*_blockSize)
    // this tells us how many of the values in _theta are valid Ritz values
    int _nevLocal;
    //
    // _hasP tells us whether there is valid data in P (and KP,MP)
    bool _hasP;
    //
    // State Multivecs
    Teuchos::RefCountPtr<MV> _X, _KX, _MX, _R,
                             _H, _KH, _MH,
                             _P, _KP, _MP;
    // tmpMV is needed only if _fullOrtho == true
    // because is depends on _fullOrtho, which is easily toggled by the user, we will allocate it 
    // and deallocate it inside of iterate()
    Teuchos::RefCountPtr<MV> _tmpMV;        
    // 
    // auxiliary vectors
    Teuchos::Array<Teuchos::RefCountPtr<const MV> > _auxVecs;
    int _numAuxVecs;
    //
    // Number of iterations that have been performed.
    int _iter;
    // 
    // Current eigenvalues, residual norms
    std::vector<MagnitudeType> _theta, _Rnorms, _R2norms, _ritz2norms;

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
    _problem(problem), 
    _sm(sorter),
    _om(printer),
    _tester(tester),
    _orthman(ortho),
    _Op(_problem->getOperator()),
    _MOp(_problem->getM()),
    _Prec(_problem->getPrec()),
    _hasM(_MOp != Teuchos::null),
    _MSUtils(_om),
    // timers, counters
    _timerOp(Teuchos::TimeMonitor::getNewTimer("Operation Op*x")),
    _timerMOp(Teuchos::TimeMonitor::getNewTimer("Operation M*x")),
    _timerPrec(Teuchos::TimeMonitor::getNewTimer("Operation Prec*x")),
    _timerSort(Teuchos::TimeMonitor::getNewTimer("Sorting eigenvalues")),
    _timerLocalProj(Teuchos::TimeMonitor::getNewTimer("Local projection")),
    _timerDS(Teuchos::TimeMonitor::getNewTimer("Direct solve")),
    _timerLocalUpdate(Teuchos::TimeMonitor::getNewTimer("Local update")),
    _timerCompRes(Teuchos::TimeMonitor::getNewTimer("Computing residuals")),
    _timerOrtho(Teuchos::TimeMonitor::getNewTimer("Orthogonalization")),
    _timerInit(Teuchos::TimeMonitor::getNewTimer("Initialization")),
    _count_ApplyOp(0),
    _count_ApplyM(0),
    _count_ApplyPrec(0),
    // internal data
    _fullOrtho(params.get("Full Ortho", true)),
    _initialized(false),
    _nevLocal(0),
    _hasP(false),
    _auxVecs( Teuchos::Array<Teuchos::RefCountPtr<const MV> >(0) ), 
    _numAuxVecs(0),
    _iter(0)
  {     
    TEST_FOR_EXCEPTION(_problem == Teuchos::null,std::invalid_argument,
                       "Anasazi::LOBPCG::constructor: user specified null problem pointer.");
    TEST_FOR_EXCEPTION(_problem->isProblemSet() == false, std::invalid_argument,
                       "Anasazi::LOBPCG::constructor: user specified problem is not set.");
    TEST_FOR_EXCEPTION(_problem->isHermitian() == false, std::invalid_argument,
                       "Anasazi::LOBPCG::constructor: user specified problem is not hermitian.");

    // set the block size and allocate data
    _blockSize = 0;
    int bs = params.get("Block Size", _problem->getNEV());
    setBlockSize(bs);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the block size and make necessary adjustments.
  template <class ScalarType, class MV, class OP>
  void LOBPCG<ScalarType,MV,OP>::setBlockSize (int blockSize) 
  {
    // time spent here counts towards _timerInit
    Teuchos::TimeMonitor lcltimer( *_timerInit );


    // This routine only allocates space; it doesn't not perform any computation
    // if size is decreased, take the first blockSize vectors of all and leave state as is
    // otherwise, grow/allocate space and set solver to unitialized

    TEST_FOR_EXCEPTION(blockSize <= 0, std::invalid_argument, "Anasazi::LOBPCG::setBlockSize was passed a non-positive block size");
    if (blockSize == _blockSize) {
      // do nothing
      return;
    }
    else if (blockSize < _blockSize) {
      // shrink vectors
      _blockSize = blockSize;

      _theta.resize(3*_blockSize);
      _ritz2norms.resize(3*_blockSize);
      _Rnorms.resize(_blockSize);
      _R2norms.resize(_blockSize);

      if (_initialized) {
        // shrink multivectors with copy
        // create ind = {0, 1, ..., blockSize-1}
        std::vector<int> ind(_blockSize);
        for (int i=0; i<_blockSize; i++) ind[i] = i;
        
        _X  = MVT::CloneCopy(*_X,ind);
        _KX = MVT::CloneCopy(*_KX,ind);
        if (_hasM) {
          _MX = MVT::CloneCopy(*_MX,ind);
        }
        else {
          _MX = _X;
        }
        _R  = MVT::CloneCopy(*_R,ind);
        _P  = MVT::CloneCopy(*_P,ind);
        _KP = MVT::CloneCopy(*_KP,ind);
        if (_hasM) {
          _MP = MVT::CloneCopy(*_MP,ind);
        }
        else {
          _MP = _P;
        }
      }
      else {
        // shrink multivectors without copying
        _X = MVT::Clone(*_X,_blockSize);
        _KX = MVT::Clone(*_KX,_blockSize);
        if (_hasM) {
          _MX = MVT::Clone(*_MX,_blockSize);
        }
        else {
          _MX = _X;
        }
        _R = MVT::Clone(*_R,_blockSize);
        _P = MVT::Clone(*_P,_blockSize);
        _KP = MVT::Clone(*_KP,_blockSize);
        if (_hasM) {
          _MP = MVT::Clone(*_MP,_blockSize);
        }
        else {
          _MP = _P;
        }
      }
      // shrink H
      _H = MVT::Clone(*_H,_blockSize);
      _KH = MVT::Clone(*_KH,_blockSize);
      if (_hasM) {
        _MH = MVT::Clone(*_MH,_blockSize);
      }
      else {
        _MH = _H;
      }
    } 
    else {  // blockSize > _blockSize
      // this is also the scenario for our initial call to setBlockSize(), in the constructor
      _initialized = false;

      Teuchos::RefCountPtr<const MV> tmp;
      // grab some Multivector to Clone
      // in practice, getInitVec() should always provide this, but it is possible to use a 
      // Eigenproblem with nothing in getInitVec() by manually initializing with initialize(); 
      // in case of that strange scenario, we will try to Clone from _X
      if (_blockSize > 0) {
        tmp = _X;
      }
      else {
        tmp = _problem->getInitVec();
        TEST_FOR_EXCEPTION(tmp == Teuchos::null,std::logic_error,
                           "Anasazi::LOBPCG::setBlockSize(): Eigenproblem did not specify initial vectors to clone from");
      }
      // grow/allocate vectors
      _theta.resize(3*blockSize,NANVAL);
      _ritz2norms.resize(3*_blockSize,NANVAL);
      _Rnorms.resize(blockSize,NANVAL);
      _R2norms.resize(blockSize,NANVAL);
      
      // clone multivectors off of tmp
      _X = MVT::Clone(*tmp,blockSize);
      _KX = MVT::Clone(*tmp,blockSize);
      if (_hasM) {
        _MX = MVT::Clone(*tmp,blockSize);
      }
      else {
        _MX = _X;
      }
      _R = MVT::Clone(*tmp,blockSize);
      _H = MVT::Clone(*tmp,blockSize);
      _KH = MVT::Clone(*tmp,blockSize);
      if (_hasM) {
        _MH = MVT::Clone(*tmp,blockSize);
      }
      else {
        _MH = _H;
      }
      _hasP = false;
      _P = MVT::Clone(*tmp,blockSize);
      _KP = MVT::Clone(*tmp,blockSize);
      if (_hasM) {
        _MP = MVT::Clone(*tmp,blockSize);
      }
      else {
        _MP = _P;
      }
      _blockSize = blockSize;
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the auxiliary vectors
  template <class ScalarType, class MV, class OP>
  void LOBPCG<ScalarType,MV,OP>::setAuxVecs(const Teuchos::Array<Teuchos::RefCountPtr<const MV> > &auxvecs) {
    typedef typename Teuchos::Array<Teuchos::RefCountPtr<const MV> >::iterator tarcpmv;

    // set new auxiliary vectors
    _auxVecs = auxvecs;
    
    if (_om->isVerbosity( Debug ) ) {
      // Check almost everything here
      CheckList chk;
      chk.checkQ = true;
      _om->print( Debug, accuracyCheck(chk, ": in setAuxVecs()") );
    }

    _numAuxVecs = 0;
    for (tarcpmv i=_auxVecs.begin(); i != _auxVecs.end(); i++) {
      _numAuxVecs += MVT::GetNumberVecs(**i);
    }
    
    // If the solver has been initialized, X and P are not necessarily orthogonal to new auxiliary vectors
    if (_numAuxVecs > 0 && _initialized) {
      _initialized = false;
      _hasP = false;
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  /* Initialize the state of the solver
   * 
   * POST-CONDITIONS:
   *
   * _initialized == true
   * X is orthonormal, orthogonal to _auxVecs
   * KX = Op*X
   * MX = M*X if _hasM
   * _theta contains Ritz values of X
   * R = KX - MX*diag(_theta)
   * if hasP() == true,
   *   P orthogonal to _auxVecs
   *   if _fullOrtho == true,
   *     P orthonormal and orthogonal to X
   *   KP = Op*P
   *   MP = M*P
   */
  template <class ScalarType, class MV, class OP>
  void LOBPCG<ScalarType,MV,OP>::initialize(LOBPCGState<ScalarType,MV> newstate)
  {
    // NOTE: memory has been allocated by setBlockSize(). Use SetBlock below; do not Clone
    // NOTE: Time spent in this routine is allotted to _timerInit, in addition to the respective sections.

    _hasP = false;  // this will be set to true below if appropriate

    std::vector<int> bsind(_blockSize);
    for (int i=0; i<_blockSize; i++) bsind[i] = i;

    // set up X: if the user doesn't specify X, ignore the rest
    if (newstate.X != Teuchos::null && MVT::GetNumberVecs(*newstate.X) >= _blockSize && MVT::GetVecLength(*newstate.X) == MVT::GetVecLength(*_X) ) {

      Teuchos::TimeMonitor lcltimer( *_timerInit );

      // put data in X,MX,KX
      MVT::SetBlock(*newstate.X,bsind,*_X);
      if (_hasM) {
        if (newstate.MX != Teuchos::null && MVT::GetNumberVecs(*newstate.MX) >= _blockSize && MVT::GetVecLength(*newstate.MX) == MVT::GetVecLength(*_MX) ) {
          MVT::SetBlock(*newstate.MX,bsind,*_MX);
        }
        else {
          Teuchos::TimeMonitor lcltimer( *_timerMOp );
          OPT::Apply(*_MOp,*_X,*_MX);
          _count_ApplyM += _blockSize;
        }
      }
      else {
        // an assignment would be redundant; take advantage of this opportunity to debug a little
        TEST_FOR_EXCEPTION(_MX != _X, std::logic_error, "Anasazi::LOBPCG::initialize(): solver invariant not satisfied");
      }
      if (newstate.KX != Teuchos::null && MVT::GetNumberVecs(*newstate.KX) >= _blockSize && MVT::GetVecLength(*newstate.KX) == MVT::GetVecLength(*_KX) ) {
        MVT::SetBlock(*newstate.KX,bsind,*_KX);
      }
      else {
        Teuchos::TimeMonitor lcltimer( *_timerOp );
        OPT::Apply(*_Op,*_X,*_KX);
        _count_ApplyOp += _blockSize;
      }

      // set up Ritz values
      _theta.resize(3*_blockSize,NANVAL);
      _ritz2norms.resize(3*_blockSize,NANVAL);
      if (newstate.T != Teuchos::null && (signed int)(newstate.T->size()) >= _blockSize) {
        for (int i=0; i<_blockSize; i++) {
          _theta[i] = (*newstate.T)[i];
        }
      }
      else {
        // get ritz vecs/vals
        Teuchos::SerialDenseMatrix<int,ScalarType> KK(_blockSize,_blockSize),
                                                   MM(_blockSize,_blockSize),
                                                    S(_blockSize,_blockSize);
        // project the problem matrices
        {
          Teuchos::TimeMonitor lcltimer( *_timerLocalProj );
          // project K
          MVT::MvTransMv(ONE,*_X,*_KX,KK);
          // project M
          MVT::MvTransMv(ONE,*_X,*_MX,MM);
          _nevLocal = _blockSize;
        }

        // solve the projected problem
        {
          Teuchos::TimeMonitor lcltimer( *_timerDS );
          _MSUtils.directSolver(_blockSize, KK, &MM, &S, &_theta, &_nevLocal, 1);
          TEST_FOR_EXCEPTION(_nevLocal != _blockSize,LOBPCGInitFailure,
                             "Anasazi::LOBPCG::initialize(): Not enough Ritz vectors to initialize algorithm.");
        }

        // We only have _blockSize ritz pairs, but we still want them in the correct order
        {
          Teuchos::TimeMonitor lcltimer( *_timerSort );
          Teuchos::BLAS<int,ScalarType> blas;
          // The sort manager is templated on ScalarType
          // Make a ScalarType copy of _theta for sorting

          std::vector<int> _order(_blockSize);

          std::vector<ScalarType> _theta_st(_blockSize);
          std::copy(_theta.begin(),_theta.begin()+_blockSize,_theta_st.begin());

          _sm->sort( this, _blockSize, &(_theta_st[0]), &_order );   // don't catch exception
          
          // Put the sorted ritz values back into _theta
          for (int i=0; i<_blockSize; i++) {
            _theta[i] = SCT::real(_theta_st[i]);
          }

          // Sort the primitive ritz vectors
          Teuchos::SerialDenseMatrix<int,ScalarType> copyS( S );
          for (int i=0; i<_blockSize; i++) {
            blas.COPY(_blockSize, copyS[_order[i]], 1, S[i], 1);
          }
        }

        // compute ritz residual norms
        {
          Teuchos::BLAS<int,ScalarType> blas;
          Teuchos::SerialDenseMatrix<int,ScalarType> R(_blockSize,_blockSize);
          // R = MM*S*diag(theta) - KK*S
          R.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,ONE,MM,S,ZERO);
          for (int i=0; i<_blockSize; i++) {
            blas.SCAL(_blockSize,_theta[i],R[i],1);
          }
          R.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,-ONE,KK,S,ONE);
          for (int i=0; i<_blockSize; i++) {
            _ritz2norms[i] = blas.NRM2(_blockSize,R[i],1);
          }
        }

        // update the solution
        {
          Teuchos::TimeMonitor lcltimer( *_timerLocalUpdate );
          // X <- X*S
          MVT::MvAddMv( ONE, *_X, ZERO, *_X, *_R );        
          MVT::MvTimesMatAddMv( ONE, *_R, S, ZERO, *_X );
          // KX <- KX*S
          MVT::MvAddMv( ONE, *_KX, ZERO, *_KX, *_R );        
          MVT::MvTimesMatAddMv( ONE, *_R, S, ZERO, *_KX );
          if (_hasM) {
            // MX <- MX*S
            MVT::MvAddMv( ONE, *_MX, ZERO, *_MX, *_R );        
            MVT::MvTimesMatAddMv( ONE, *_R, S, ZERO, *_MX );
          }
        }
      }

  
      // set up R
      if (newstate.R != Teuchos::null && MVT::GetNumberVecs(*newstate.R) >= _blockSize && MVT::GetVecLength(*newstate.R) == MVT::GetVecLength(*_R) ) {
        MVT::SetBlock(*newstate.R,bsind,*_R);
      }
      else {
        Teuchos::TimeMonitor lcltimer( *_timerCompRes );
        // form R <- KX - MX*T
        MVT::MvAddMv(ZERO,*_KX,ONE,*_KX,*_R);
        Teuchos::SerialDenseMatrix<int,ScalarType> T(_blockSize,_blockSize);
        T.putScalar(ZERO);
        for (int i=0; i<_blockSize; i++) T(i,i) = _theta[i];
        MVT::MvTimesMatAddMv(-ONE,*_MX,T,ONE,*_R);
      }
      // Update the residual norms
      _orthman->norm(*_R,&_Rnorms);
      // Update the residual 2-norms 
      MVT::MvNorm(*_R,&_R2norms);

  
      // put data in P,KP,MP: P is not used to set theta
      if (newstate.P != Teuchos::null && MVT::GetNumberVecs(*newstate.P) >= _blockSize && MVT::GetVecLength(*newstate.P) == MVT::GetVecLength(*_P) ) {
        _hasP = true;

        MVT::SetBlock(*newstate.P,bsind,*_P);

        if (newstate.KP != Teuchos::null && MVT::GetNumberVecs(*newstate.KP) >= _blockSize && MVT::GetVecLength(*newstate.KP) == MVT::GetVecLength(*_KP) ) {
          MVT::SetBlock(*newstate.KP,bsind,*_KP);
        }
        else {
          Teuchos::TimeMonitor lcltimer( *_timerOp );
          OPT::Apply(*_Op,*_P,*_KP);
          _count_ApplyOp += _blockSize;
        }

        if (_hasM) {
          if (newstate.MP != Teuchos::null && MVT::GetNumberVecs(*newstate.MP) >= _blockSize && MVT::GetVecLength(*newstate.MP) == MVT::GetVecLength(*_MP) ) {
            MVT::SetBlock(*newstate.MP,bsind,*_MP);
          }
          else {
            Teuchos::TimeMonitor lcltimer( *_timerMOp );
            OPT::Apply(*_MOp,*_P,*_MP);
            _count_ApplyM += _blockSize;
          }
        }
      }

      _initialized = true;

      if (_om->isVerbosity( Debug ) ) {
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
        _om->print( Debug, accuracyCheck(chk, ": after initialize()") );
      }

    }
    else {

      LOBPCGState<ScalarType,MV> newstate;
      { // begin timer scope
        Teuchos::TimeMonitor lcltimer( *_timerInit );

        // generate something, projectAndNormalize, call myself recursively
        Teuchos::RefCountPtr<const MV> ivec = _problem->getInitVec();
        TEST_FOR_EXCEPTION(ivec == Teuchos::null,std::logic_error,
                           "Anasazi::LOBPCG::initialize(): Eigenproblem did not specify initial vectors to clone from");

        int initSize = MVT::GetNumberVecs(*ivec);
        if (initSize > _blockSize) {
          // we need only the first _blockSize vectors from ivec; get a view of them
          initSize = _blockSize;
          std::vector<int> ind(_blockSize);
          for (int i=0; i<_blockSize; i++) ind[i] = i;
          ivec = MVT::CloneView(*ivec,ind);
        }

        // alloc newX
        Teuchos::RefCountPtr<MV> newMX, newX = MVT::Clone(*ivec,_blockSize);
        // assign ivec to first part of newX
        std::vector<int> ind(initSize);
        if (initSize > 0) {
          for (int i=0; i<initSize; i++) ind[i] = i;
          MVT::SetBlock(*ivec,ind,*newX);
        }
        // fill the rest of newX with random
        if (_blockSize > initSize) {
          ind.resize(_blockSize - initSize);
          for (int i=0; i<_blockSize - initSize; i++) ind[i] = initSize + i;
          Teuchos::RefCountPtr<MV> rX = MVT::CloneView(*newX,ind);
          MVT::MvRandom(*rX);
          rX = Teuchos::null;
        }

        // compute newMX if _hasM
        if (_hasM) {
          newMX = MVT::Clone(*ivec,_blockSize);
          {
            Teuchos::TimeMonitor lcltimer( *_timerMOp );
            OPT::Apply(*_MOp,*newX,*newMX);
            _count_ApplyM += _blockSize;
          }
        }
        else {
          newMX = Teuchos::null;
        }

        // remove auxVecs from newX and normalize newX
        if (_auxVecs.size() > 0) {
          Teuchos::TimeMonitor lcltimer( *_timerOrtho );
          Teuchos::Array<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > > dummy;
          int rank = _orthman->projectAndNormalize(*newX,newMX,dummy,Teuchos::null,_auxVecs);
          TEST_FOR_EXCEPTION(rank != _blockSize,LOBPCGInitFailure,
                             "Anasazi::LOBPCG::initialize(): Couldn't generate initial basis of full rank.");
        }
        else {
          Teuchos::TimeMonitor lcltimer( *_timerOrtho );
          int rank = _orthman->normalize(*newX,newMX,Teuchos::null);
          TEST_FOR_EXCEPTION(rank != _blockSize,LOBPCGInitFailure,
                             "Anasazi::LOBPCG::initialize(): Couldn't generate initial basis of full rank.");
        }

        // call myself recursively
        newstate.X = newX;
        newstate.MX = newMX;
      } // end of timer scope; we needed this because the following recursive call to initialize contains its own call to _timerInit
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
    if ( _fullOrtho == true || _initialized == false || fullOrtho == _fullOrtho ) {
      // state is already orthogonalized or solver is not initialized
      _fullOrtho = fullOrtho;
      return;
    }

    // solver is initialized, state is not fully orthogonalized, and user has requested full orthogonalization
    _fullOrtho = true;
    // throw away data in P
    _hasP = false;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Perform LOBPCG iterations until the StatusTest tells us to stop.
  template <class ScalarType, class MV, class OP>
  void LOBPCG<ScalarType,MV,OP>::iterate () 
  {
    //
    // Allocate/initialize data structures
    //
    if (_initialized == false) {
      initialize();
    }

    // if _fullOrtho == true, then we must produce the following on every iteration:
    // [newX newP] = [X H P] [CX;CP]
    // the structure of [CX;CP] when using full orthogonalization does not allow us to 
    // do this in place, and _R does not have enough storage for newX and newP. therefore, 
    // we must allocate additional storage for this.
    // otherwise, when not using full orthogonalization, the structure
    // [newX newP] = [X H P] [CX1  0 ]
    //                       [CX2 CP2]  allows us to work using only R as work space
    //                       [CX3 CP3] 
    if (_fullOrtho) {
      if (_tmpMV == Teuchos::null || MVT::GetNumberVecs(*_tmpMV) != _blockSize) {
        _tmpMV = MVT::Clone(*_X,_blockSize);
      }
    }
    else {
      _tmpMV = Teuchos::null;
    }

    //
    // Miscellaneous definitions
    const int oneBlock    =   _blockSize;
    const int twoBlocks   = 2*_blockSize;
    const int threeBlocks = 3*_blockSize;
    
    //
    // Define dense projected/local matrices
    //   KK = Local stiffness matrix               (size: 3*_blockSize x 3*_blockSize)
    //   MM = Local mass matrix                    (size: 3*_blockSize x 3*_blockSize)
    //    S = Local eigenvectors                   (size: 3*_blockSize x 3*_blockSize)
    Teuchos::SerialDenseMatrix<int,ScalarType> KK( threeBlocks, threeBlocks ), 
                                               MM( threeBlocks, threeBlocks ),
                                                S( threeBlocks, threeBlocks );

    while (_tester->checkStatus(this) != Passed) {

      _iter++;
      
      // Apply the preconditioner on the residuals: H <- Prec*R
      if (_Prec != Teuchos::null) {
        Teuchos::TimeMonitor lcltimer( *_timerPrec );
        OPT::Apply( *_Prec, *_R, *_H );   // don't catch the exception
        _count_ApplyPrec += _blockSize;
      }
      else {
        std::vector<int> ind(_blockSize);
        for (int i=0; i<_blockSize; i++) { ind[i] = i; }
        MVT::SetBlock(*_R,ind,*_H);
      }

      // Apply the mass matrix on H
      if (_hasM) {
        Teuchos::TimeMonitor lcltimer( *_timerMOp );
        OPT::Apply( *_MOp, *_H, *_MH);    // don't catch the exception
        _count_ApplyM += _blockSize;
      }

      // orthogonalize H against the auxiliary vectors
      // optionally: orthogonalize H against X and P ([X P] is already orthonormal)
      Teuchos::Array<Teuchos::RefCountPtr<const MV> > Q;
      Teuchos::Array<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > > C = 
        Teuchos::tuple<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > >(Teuchos::null);
      Q = _auxVecs;
      if (_fullOrtho) {
        Q.push_back(_X);
        if (_hasP) {
          Q.push_back(_P);
        }
      }
      {
        Teuchos::TimeMonitor lcltimer( *_timerOrtho );
        int rank = _orthman->projectAndNormalize(*_H,_MH,C,Teuchos::null,Q);
        TEST_FOR_EXCEPTION(rank != _blockSize,LOBPCGOrthoFailure,
                           "Anasazi::LOBPCG::iterate(): unable to compute full basis for H");
      }

      if (_om->isVerbosity( Debug ) ) {
        CheckList chk;
        chk.checkH = true;
        chk.checkMH = true;
        _om->print( Debug, accuracyCheck(chk, ": after ortho H") );
      }
      else if (_om->isVerbosity( OrthoDetails ) ) {
        CheckList chk;
        chk.checkH = true;
        chk.checkMH = true;
        _om->print( OrthoDetails, accuracyCheck(chk,": after ortho H") );
      }

      // Apply the stiffness matrix to H
      {
        Teuchos::TimeMonitor lcltimer( *_timerOp );
        OPT::Apply( *_Op, *_H, *_KH);   // don't catch the exception
        _count_ApplyOp += _blockSize;
      }

      int localSize;
      if (_hasP) {
        localSize = threeBlocks;
      }
      else {
        localSize = twoBlocks;
      }

      // Form "local" mass and stiffness matrices
      {
        Teuchos::TimeMonitor lcltimer( *_timerLocalProj );
        /* We will construct (block) upper triangular parts only.
                 (X^H)             (KK11 KK12 KK13)
            KK = (H^H) K [X H P] = ( --  KK22 KK23)
                 (P^H)             ( --   --  KK33)
                 (X^H)             (MM11 MM12 MM13)
            MM = (H^H) M [X H P] = ( --  MM22 MM23) 
                 (P^H)             ( --   --  MM33)
        */
        Teuchos::SerialDenseMatrix<int,ScalarType> KK11( Teuchos::View, KK, _blockSize, _blockSize );
        MVT::MvTransMv( ONE, *_X, *_KX, KK11 );
        Teuchos::SerialDenseMatrix<int,ScalarType> KK12( Teuchos::View, KK, _blockSize, _blockSize, 0, _blockSize );
        MVT::MvTransMv( ONE, *_X, *_KH, KK12 );
        Teuchos::SerialDenseMatrix<int,ScalarType> KK22( Teuchos::View, KK, _blockSize, _blockSize, _blockSize, _blockSize );
        MVT::MvTransMv( ONE, *_H, *_KH, KK22 );
        
        Teuchos::SerialDenseMatrix<int,ScalarType> MM11( Teuchos::View, MM, _blockSize, _blockSize );
        MVT::MvTransMv( ONE, *_X, *_MX, MM11 );
        Teuchos::SerialDenseMatrix<int,ScalarType> MM12( Teuchos::View, MM, _blockSize, _blockSize, 0, _blockSize );
        MVT::MvTransMv( ONE, *_X, *_MH, MM12 );
        Teuchos::SerialDenseMatrix<int,ScalarType> MM22( Teuchos::View, MM, _blockSize, _blockSize, _blockSize, _blockSize );
        MVT::MvTransMv( ONE, *_H, *_MH, MM22 );

        if (_hasP) {
          Teuchos::SerialDenseMatrix<int,ScalarType> KK13( Teuchos::View, KK, _blockSize, _blockSize, 0, twoBlocks );
          MVT::MvTransMv( ONE, *_X, *_KP, KK13 );
          Teuchos::SerialDenseMatrix<int,ScalarType> KK23( Teuchos::View, KK, _blockSize, _blockSize, _blockSize, twoBlocks );
          MVT::MvTransMv( ONE, *_H, *_KP, KK23 );
          Teuchos::SerialDenseMatrix<int,ScalarType> KK33( Teuchos::View, KK, _blockSize, _blockSize, twoBlocks, twoBlocks );
          MVT::MvTransMv( ONE, *_P, *_KP, KK33 );
          
          Teuchos::SerialDenseMatrix<int,ScalarType> MM13( Teuchos::View, MM, _blockSize, _blockSize, 0, twoBlocks );
          MVT::MvTransMv( ONE, *_X, *_MP, MM13 );
          Teuchos::SerialDenseMatrix<int,ScalarType> MM23( Teuchos::View, MM, _blockSize, _blockSize, _blockSize, twoBlocks );
          MVT::MvTransMv( ONE, *_H, *_MP, MM23 );
          Teuchos::SerialDenseMatrix<int,ScalarType> MM33( Teuchos::View, MM, _blockSize, _blockSize, twoBlocks, twoBlocks );
          MVT::MvTransMv( ONE, *_P, *_MP, MM33 );
        }
      } // end timing block

      // Perform a spectral decomposition
      {
        Teuchos::TimeMonitor lcltimer( *_timerDS );
        _nevLocal = localSize;
        _MSUtils.directSolver(localSize, KK, &MM, &S, &_theta, &_nevLocal, 0);    // don't catch the exception
        // localSize tells directSolver() how big KK,MM are
        // however, directSolver() may choose to use only the principle submatrices of KK,MM because of loss of 
        // MM-orthogonality in the projected eigenvectors
        // _nevLocal tells us how much it used, in effect dictating back to us how big localSize is, 
        // and therefore telling us which of [X H P] to use in computing the new iterates below.
        // we will not tolerate this ill-conditioning, and will throw an exception.
      }
      _om->stream(Debug) << " After directSolve: localSize == " << localSize << " \tnevLocal == " << _nevLocal << endl;
      TEST_FOR_EXCEPTION(_nevLocal != localSize, LOBPCGRitzFailure, "Ill-conditioning detected in projecteded mass matrix." );

      Teuchos::LAPACK<int,ScalarType> lapack;
      Teuchos::BLAS<int,ScalarType> blas;
      //
      //---------------------------------------------------
      // Sort the ritz values using the sort manager
      //---------------------------------------------------
      // The sort manager is templated on ScalarType
      // Make a ScalarType copy of _theta for sorting
      {
        Teuchos::TimeMonitor lcltimer( *_timerSort );

        std::vector<int> _order(_nevLocal);

        std::vector<ScalarType> _theta_st(_nevLocal);
        std::copy(_theta.begin(),_theta.begin()+_nevLocal,_theta_st.begin());

        _sm->sort( this, _nevLocal, &(_theta_st[0]), &_order );   // don't catch exception
        
        // Put the sorted ritz values back into _theta
        for (int i=0; i<_nevLocal; i++) {
          _theta[i] = SCT::real(_theta_st[i]);
        }

        // Sort the primitive ritz vectors
        // We need the first _blockSize vectors ordered to generate the next
        // columns immediately below, as well as later, when/if we restart.
        Teuchos::SerialDenseMatrix<int,ScalarType> copyS( S );
        for (int i=0; i<_nevLocal; i++) {
          blas.COPY(_nevLocal, copyS[_order[i]], 1, S[i], 1);
        }
      }

      // compute ritz residual norms
      {
        Teuchos::SerialDenseMatrix<int,ScalarType> R(_nevLocal,_nevLocal);
        Teuchos::SerialDenseMatrix<int,ScalarType> lclKK(Teuchos::View,KK,_nevLocal,_nevLocal),
                                                   lclMM(Teuchos::View,MM,_nevLocal,_nevLocal);
        // R = MM*S*diag(theta) - KK*S
        R.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,ONE,lclMM,S,ZERO);
        for (int i=0; i<_nevLocal; i++) {
          blas.SCAL(_nevLocal,_theta[i],R[i],1);
        }
        R.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,-ONE,lclKK,S,ONE);
        for (int i=0; i<_nevLocal; i++) {
          _ritz2norms[i] = blas.NRM2(_nevLocal,R[i],1);
        }
      }


      // before computing X,P: perform second orthogonalization per Ulrich,Rich paper
      // CX will be the coefficients of [X,H,P] for new X, CP for new P
      // The paper suggests orthogonalizing CP against CX and orthonormalizing CP, w.r.t. MM
      // Here, we will also orthonormalize CX.
      // This is accomplished using the Cholesky factorization of [CX CP]^H MM [CX CP]
      Teuchos::SerialDenseMatrix<int,ScalarType> CX(threeBlocks,oneBlock), CP(threeBlocks,oneBlock);
      if (_fullOrtho && localSize >= twoBlocks) {
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
        if (_om->isVerbosity( Debug ) ) {
          Teuchos::SerialDenseMatrix<int,ScalarType> tmp1(threeBlocks,oneBlock),
                                                     tmp2(oneBlock,oneBlock);
          MagnitudeType tmp;
          int teuchosret;
          stringstream os;
          os.precision(2);
          os.setf(ios::scientific, ios::floatfield);

          os << " Checking Full Ortho: iteration " << _iter << endl;

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
          _om->print(Debug,os.str());
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
      _hasP = false;
      {
        Teuchos::TimeMonitor lcltimer( *_timerLocalUpdate );

        if (localSize == oneBlock) {

          Teuchos::SerialDenseMatrix<int,ScalarType> CX1( Teuchos::View, CX, _blockSize, _blockSize );

          // X <- R*CX1 == X*CX1 
          MVT::MvAddMv( ONE, *_X, ZERO, *_X, *_R );        
          MVT::MvTimesMatAddMv( ONE, *_R, CX1, ZERO, *_X );
          // KX <- R*CX1 == KX*CX1 
          MVT::MvAddMv( ONE, *_KX, ZERO, *_KX, *_R );
          MVT::MvTimesMatAddMv( ONE, *_R, CX1, ZERO, *_KX );
          if (_hasM) {
            // MX <- R*CX1 == MX*CX1 
            MVT::MvAddMv( ONE, *_MX, ZERO, *_MX, *_R );
            MVT::MvTimesMatAddMv( ONE, *_R, CX1, ZERO, *_MX );
          }
        } // if (localSize == oneBlocks)
        else if (localSize == twoBlocks) {
          Teuchos::SerialDenseMatrix<int,ScalarType> CX1( Teuchos::View, CX, _blockSize, _blockSize ),
                                                     CX2( Teuchos::View, CX, _blockSize, _blockSize, oneBlock, 0 ),
                                                     CP1( Teuchos::View, CP, _blockSize, _blockSize ),
                                                     CP2( Teuchos::View, CP, _blockSize, _blockSize, oneBlock, 0 );

          _hasP = true;

          // X = [X H][CX1]
          //          [CX2]
          //
          // P = [X H][CP1]
          //          [CP2]
          //

          // R <- X
          MVT::MvAddMv( ONE, *_X, ZERO, *_X, *_R );        
          // X <- R*CX1 + H*CX2 == X*CX1 + H*CX2
          MVT::MvTimesMatAddMv( ONE, *_R, CX1, ZERO, *_X );
          MVT::MvTimesMatAddMv( ONE, *_H, CX2, ONE , *_X );
          // P <- R*CP1 + H*CP2 == X*CP1 + H*CP2
          // however, if _fullOrtho == false, CP1 == ZERO
          if (_fullOrtho) { 
            MVT::MvTimesMatAddMv( ONE, *_R, CP1, ZERO, *_P );
            MVT::MvTimesMatAddMv( ONE, *_H, CP2, ONE, *_P );
          }
          else {
            MVT::MvTimesMatAddMv( ONE, *_H, CP2, ZERO, *_P );
          }

          // R  <- KX
          MVT::MvAddMv( ONE, *_KX, ZERO, *_KX, *_R );        
          // KX <- R*CX1 + KH*CX2 == KX*CX1 + KH*CX2
          MVT::MvTimesMatAddMv( ONE, *_R, CX1, ZERO, *_KX );
          MVT::MvTimesMatAddMv( ONE, *_KH, CX2, ONE , *_KX );
          // KP <- R*CP1 + KH*CP2 == KX*CP1 + KH*CP2
          // however, if _fullOrtho == false, CP1 == ZERO
          if (_fullOrtho) { 
            MVT::MvTimesMatAddMv( ONE, *_R, CP1, ZERO, *_KP );
            MVT::MvTimesMatAddMv( ONE, *_KH, CP2, ONE, *_KP );
          }
          else {
            MVT::MvTimesMatAddMv( ONE, *_KH, CP2, ZERO, *_KP );
          }

          if (_hasM) {
            // R  <- MX
            MVT::MvAddMv( ONE, *_MX, ZERO, *_MX, *_R );        
            // MX <- R*CX1 + MH*CX2 == MX*CX1 + MH*CX2
            MVT::MvTimesMatAddMv( ONE, *_R, CX1, ZERO, *_MX );
            MVT::MvTimesMatAddMv( ONE, *_MH, CX2, ONE , *_MX );
            // MP <- R*CP1 + MH*CP2 == MX*CP1 + MH*CP2
            // however, if _fullOrtho == false, CP1 == ZERO
            if (_fullOrtho) { 
              MVT::MvTimesMatAddMv( ONE, *_R, CP1, ZERO, *_MP );
              MVT::MvTimesMatAddMv( ONE, *_MH, CP2, ONE, *_MP );
            }
            else {
              MVT::MvTimesMatAddMv( ONE, *_MH, CP2, ZERO, *_MP );
            }
          }

        } // if (localSize == twoBlocks)
        else if (localSize == threeBlocks) {
          Teuchos::SerialDenseMatrix<int,ScalarType> CX1( Teuchos::View, CX, _blockSize, _blockSize ),
                                                     CX2( Teuchos::View, CX, _blockSize, _blockSize, oneBlock, 0 ),
                                                     CX3( Teuchos::View, CX, _blockSize, _blockSize, twoBlocks, 0 ),
                                                     CP1( Teuchos::View, CP, _blockSize, _blockSize ),
                                                     CP2( Teuchos::View, CP, _blockSize, _blockSize, oneBlock, 0 ),
                                                     CP3( Teuchos::View, CP, _blockSize, _blockSize, twoBlocks, 0 );

          _hasP = true;

          // X <- X*CX1 + P*CX3
          // P <- X*CP1 + P*CP3 (note, CP1 == ZERO if _fullOrtho==false)
          if (_fullOrtho) {
            // copy X,P
            MVT::MvAddMv(ONE,*_X, ZERO,*_X, *_R);
            MVT::MvAddMv(ONE,*_P, ZERO,*_P, *_tmpMV);
            // perform [X P][CX1 CP1]
            //              [CX3 CP3]
            MVT::MvTimesMatAddMv( ONE, *_R,     CX1, ZERO, *_X );
            MVT::MvTimesMatAddMv( ONE, *_tmpMV, CX3,  ONE, *_X );
            MVT::MvTimesMatAddMv( ONE, *_R,     CP1, ZERO, *_P );
            MVT::MvTimesMatAddMv( ONE, *_tmpMV, CP3,  ONE, *_P );
          }
          else {
            // copy X
            MVT::MvAddMv(ONE,*_X, ZERO,*_X, *_R);
            // perform [X P][CX1  0 ]
            //              [CX3 CP3]
            MVT::MvTimesMatAddMv( ONE, *_R, CX1, ZERO, *_X );
            MVT::MvTimesMatAddMv( ONE, *_P, CX3,  ONE, *_X );
            // copy P
            MVT::MvAddMv(ONE,*_P, ZERO,*_P, *_R);
            MVT::MvTimesMatAddMv( ONE, *_R, CP3, ZERO, *_P );
          }
          // X <- X + H*CX2
          // P <- P + H*CP2
          MVT::MvTimesMatAddMv( ONE, *_H, CX2,  ONE, *_X );
          MVT::MvTimesMatAddMv( ONE, *_H, CP2,  ONE, *_P );


          // KX <- KX*CX1 + KP*CX3
          // KP <- KX*CP1 + KP*CP3 (note, CP1 == ZERO if _fullOrtho==false)
          if (_fullOrtho) {
            // copy KX,KP
            MVT::MvAddMv(ONE,*_KX, ZERO,*_KX, *_R);
            MVT::MvAddMv(ONE,*_KP, ZERO,*_KP, *_tmpMV);
            // perform [KX KP][CX1 CP1]
            //                [CX3 CP3]
            MVT::MvTimesMatAddMv( ONE, *_R,     CX1, ZERO, *_KX );
            MVT::MvTimesMatAddMv( ONE, *_tmpMV, CX3,  ONE, *_KX );
            MVT::MvTimesMatAddMv( ONE, *_R,     CP1, ZERO, *_KP );
            MVT::MvTimesMatAddMv( ONE, *_tmpMV, CP3,  ONE, *_KP );
          }
          else {
            // copy KX
            MVT::MvAddMv(ONE,*_KX, ZERO,*_KX, *_R);
            // perform [KX KP][CX1  0 ]
            //                [CX3 CP3]
            MVT::MvTimesMatAddMv( ONE, *_R , CX1, ZERO, *_KX );
            MVT::MvTimesMatAddMv( ONE, *_KP, CX3,  ONE, *_KX );
            // copy KP
            MVT::MvAddMv(ONE,*_KP, ZERO,*_KP, *_R);
            MVT::MvTimesMatAddMv( ONE, *_R, CP3, ZERO, *_KP );
          }
          // KX <- KX + KH*CX2
          // KP <- KP + KH*CP2
          MVT::MvTimesMatAddMv( ONE, *_KH, CX2,  ONE, *_KX );
          MVT::MvTimesMatAddMv( ONE, *_KH, CP2,  ONE, *_KP );


          if (_hasM) {
            // MX <- MX*CX1 + MP*CX3
            // MP <- MX*CP1 + MP*CP3 (note, CP1 == ZERO if _fullOrtho==false)
            if (_fullOrtho) {
              // copy MX,MP
              MVT::MvAddMv(ONE,*_MX, ZERO,*_MX, *_R);
              MVT::MvAddMv(ONE,*_MP, ZERO,*_MP, *_tmpMV);
              // perform [MX MP][CX1 CP1]
              //                [CX3 CP3]
              MVT::MvTimesMatAddMv( ONE, *_R,     CX1, ZERO, *_MX );
              MVT::MvTimesMatAddMv( ONE, *_tmpMV, CX3,  ONE, *_MX );
              MVT::MvTimesMatAddMv( ONE, *_R,     CP1, ZERO, *_MP );
              MVT::MvTimesMatAddMv( ONE, *_tmpMV, CP3,  ONE, *_MP );
            }
            else {
              // copy MX
              MVT::MvAddMv(ONE,*_MX, ZERO,*_MX, *_R);
              // perform [MX MP][CX1  0 ]
              //                [CX3 CP3]
              MVT::MvTimesMatAddMv( ONE, *_R , CX1, ZERO, *_MX );
              MVT::MvTimesMatAddMv( ONE, *_MP, CX3,  ONE, *_MX );
              // copy MP
              MVT::MvAddMv(ONE,*_MP, ZERO,*_MP, *_R);
              MVT::MvTimesMatAddMv( ONE, *_R, CP3, ZERO, *_MP );
            }
            // MX <- MX + MH*CX2
            // MP <- MP + MH*CP2
            MVT::MvTimesMatAddMv( ONE, *_MH, CX2,  ONE, *_MX );
            MVT::MvTimesMatAddMv( ONE, *_MH, CP2,  ONE, *_MP );
          }

        } // if (localSize == threeBlocks)
      } // end timing block

      // Compute the new residuals, explicitly
      {
        Teuchos::TimeMonitor lcltimer( *_timerCompRes );
        MVT::MvAddMv( ONE, *_KX, ZERO, *_KX, *_R );
        Teuchos::SerialDenseMatrix<int,ScalarType> T( _blockSize, _blockSize );
        for (int i = 0; i < _blockSize; i++) {
          T(i,i) = _theta[i];
        }
        MVT::MvTimesMatAddMv( -ONE, *_MX, T, ONE, *_R );
      }
      
      // When required, monitor some orthogonalities
      if (_om->isVerbosity( Debug ) ) {
        // Check almost everything here
        CheckList chk;
        chk.checkX = true;
        chk.checkKX = true;
        chk.checkMX = true;
        chk.checkP = true;
        chk.checkKP = true;
        chk.checkMP = true;
        chk.checkR = true;
        _om->print( Debug, accuracyCheck(chk, ": after local update") );
      }
      else if (_om->isVerbosity( OrthoDetails )) {
        CheckList chk;
        chk.checkX = true;
        chk.checkP = true;
        chk.checkR = true;
        _om->print( OrthoDetails, accuracyCheck(chk, ": after local update") );
      }


      // Update the residual norms
      _orthman->norm(*_R,&_Rnorms);

      // Update the residual 2-norms 
      MVT::MvNorm(*_R,&_R2norms);

      // Print information on current status
      if (_om->isVerbosity(Debug)) {
        currentStatus( _om->stream(Debug) );
      }
      else if (_om->isVerbosity(IterationDetails)) {
        currentStatus( _om->stream(IterationDetails) );
      }
    } // end while (statusTest == false)
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

    os << " Debugging checks: iteration " << _iter << where << endl;

    // X and friends
    if (chk.checkX) {
      tmp = _orthman->orthonormError(*_X);
      os << " >> Error in X^H M X == I : " << tmp << endl;
      for (unsigned int i=0; i<_auxVecs.size(); i++) {
        tmp = _orthman->orthogError(*_X,*_auxVecs[i]);
        os << " >> Error in X^H M Q[" << i << "] == 0 : " << tmp << endl;
      }
    }
    if (chk.checkMX && _hasM) {
      tmp = _MSUtils.errorEquality(_X.get(), _MX.get(), _MOp.get());
      os << " >> Error in MX == M*X    : " << tmp << endl;
    }
    if (chk.checkKX) {
      tmp = _MSUtils.errorEquality(_X.get(), _KX.get(), _Op.get());
      os << " >> Error in KX == K*X    : " << tmp << endl;
    }

    // P and friends
    if (chk.checkP && _hasP) {
      if (_fullOrtho) {
        tmp = _orthman->orthonormError(*_P);
        os << " >> Error in P^H M P == I : " << tmp << endl;
        tmp = _orthman->orthogError(*_P,*_X);
        os << " >> Error in P^H M X == 0 : " << tmp << endl;
      }
      for (unsigned int i=0; i<_auxVecs.size(); i++) {
        tmp = _orthman->orthogError(*_P,*_auxVecs[i]);
        os << " >> Error in P^H M Q[" << i << "] == 0 : " << tmp << endl;
      }
    }
    if (chk.checkMP && _hasM && _hasP) {
      tmp = _MSUtils.errorEquality(_P.get(), _MP.get(), _MOp.get());
      os << " >> Error in MP == M*P    : " << tmp << endl;
    }
    if (chk.checkKP && _hasP) {
      tmp = _MSUtils.errorEquality(_P.get(), _KP.get(), _Op.get());
      os << " >> Error in KP == K*P    : " << tmp << endl;
    }

    // H and friends
    if (chk.checkH) {
      if (_fullOrtho) {
        tmp = _orthman->orthonormError(*_H);
        os << " >> Error in H^H M H == I : " << tmp << endl;
        tmp = _orthman->orthogError(*_H,*_X);
        os << " >> Error in H^H M X == 0 : " << tmp << endl;
        if (_hasP) {
          tmp = _orthman->orthogError(*_H,*_P);
          os << " >> Error in H^H M P == 0 : " << tmp << endl;
        }
      }
      for (unsigned int i=0; i<_auxVecs.size(); i++) {
        tmp = _orthman->orthogError(*_H,*_auxVecs[i]);
        os << " >> Error in H^H M Q[" << i << "] == 0 : " << tmp << endl;
      }
    }
    if (chk.checkMH && _hasM) {
      tmp = _MSUtils.errorEquality(_H.get(), _MH.get(), _MOp.get());
      os << " >> Error in MH == M*H    : " << tmp << endl;
    }

    // R: this is not M-orthogonality, but standard euclidean orthogonality
    if (chk.checkR) {
      Teuchos::SerialDenseMatrix<int,ScalarType> xTx(_blockSize,_blockSize);
      MVT::MvTransMv(ONE,*_X,*_R,xTx);
      tmp = xTx.normFrobenius();
      os << " >> Error in X^H R == 0   : " << tmp << endl;
    }

    // Q
    if (chk.checkQ) {
      for (unsigned int i=0; i<_auxVecs.size(); i++) {
        tmp = _orthman->orthonormError(*_auxVecs[i]);
        os << " >> Error in Q[" << i << "]^H M Q[" << i << "] == I : " << tmp << endl;
        for (unsigned int j=i+1; j<_auxVecs.size(); j++) {
          tmp = _orthman->orthogError(*_auxVecs[i],*_auxVecs[j]);
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
    os <<"The solver is "<<(_initialized ? "initialized." : "not initialized.") << endl;
    os <<"The number of iterations performed is " << _iter       << endl;
    os <<"The current block size is             " << _blockSize  << endl;
    os <<"The number of auxiliary vectors is    " << _numAuxVecs << endl;
    os <<"The number of operations Op*x   is " << _count_ApplyOp   << endl;
    os <<"The number of operations M*x    is " << _count_ApplyM    << endl;
    os <<"The number of operations Prec*x is " << _count_ApplyPrec << endl;

    os.setf(ios_base::right, ios_base::adjustfield);

    if (_initialized) {
      os << endl;
      os <<"CURRENT EIGENVALUE ESTIMATES             "<<endl;
      os << std::setw(20) << "Eigenvalue" 
         << std::setw(20) << "Residual(M)"
         << std::setw(20) << "Residual(2)"
         << endl;
      os <<"--------------------------------------------------------------------------------"<<endl;
      for (int i=0; i<_blockSize; i++) {
        os << std::setw(20) << _theta[i] 
           << std::setw(20) << _Rnorms[i] 
           << std::setw(20) << _R2norms[i] 
           << endl;
      }
    }
    os <<"================================================================================" << endl;
    os << endl;
  }

  
} // end Anasazi namespace

#endif // ANASAZI_LOBPCG_HPP
