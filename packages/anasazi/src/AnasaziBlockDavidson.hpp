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

/*! \file AnasaziBlockDavidson.hpp
  \brief Implementation of the block Davidson method
*/

// TODO: the documentation here needs to be made rigorous
// in particular, getState() and initialize() need to exactly describe their 
// input and output

#ifndef ANASAZI_BLOCK_DAVIDSON_HPP
#define ANASAZI_BLOCK_DAVIDSON_HPP

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

/*!     \class Anasazi::BlockDavidson
  
        \brief This class implements the block Davidson method, an iterative
        method for solving symmetric eigenvalue problems.

        \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist
*/

namespace Anasazi {

  //@{ \name BlockDavidson Structures 

  /** \brief Structure to contain pointers to BlockDavidson state variables.
   *
   * This struct is utilized by BlockDavidson::initialize() and BlockDavidson::getState().
   */
  template <class ScalarType, class MulVec>
  struct BlockDavidsonState {
    Teuchos::RefCountPtr<const MulVec> V, KV, MV;
    Teuchos::RefCountPtr<const MulVec> X, KX, MX;
    Teuchos::RefCountPtr<const MulVec> H, KH, MH;
    Teuchos::RefCountPtr<const MulVec> R;
    Teuchos::RefCountPtr<const Teuchos::SerialDenseMatrix<int,ScalarType> > KK, MM;
    Teuchos::RefCountPtr<const std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> > T;
    BlockDavidsonState() : V(Teuchos::null),KV(Teuchos::null),MV(Teuchos::null),
                           X(Teuchos::null),KX(Teuchos::null),MX(Teuchos::null),
                           KK(Teuchos::null),MM(Teuchos::null),
                           R(Teuchos::null),T(Teuchos::null) {};
  };

  //@}

  //@{ \name BlockDavidson Exceptions

  /** \brief BlockDavidsonInitFailure is thrown when the BlockDavidson solver is unable to
   * generate an initial iterate in the BlockDavidson::initialize() routine. 
   *
   * This exception is thrown from the BlockDavidson::initialize() method, which is
   * called by the user or from the BlockDavidson::iterate() method if isInitialized()
   * == \c false.
   *
   * In the case that this exception is thrown, BlockDavidson::hasP() and
   * BlockDavidson::isInitialized() will be \c false and the user will need to provide
   * a new initial iterate to the solver.
   *
   * \relates BlockDavidson, BlockDavidson::initialize()
   */
  class BlockDavidsonInitFailure : public AnasaziError {public:
    BlockDavidsonInitFailure(const std::string& what_arg) : AnasaziError(what_arg)
    {}};

  //@}
 
  
  template <class ScalarType, class MV, class OP>
  class BlockDavidson : public Eigensolver<ScalarType,MV,OP> { 
  public:
    //@{ \name Constructor/Destructor.
    
    //! %Anasazi::BlockDavidson constructor.
    BlockDavidson( const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem, 
                   const Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > &sorter,
                   const Teuchos::RefCountPtr<OutputManager<ScalarType> > &printer,
                   const Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > &tester,
                   const Teuchos::RefCountPtr<MatOrthoManager<ScalarType,MV,OP> > &ortho,
                   Teuchos::ParameterList &params 
                 );
    
    //! %Anasazi::BlockDavidson destructor.
    virtual ~BlockDavidson() {};
    //@}
    
    //@{ \name Solver methods.
    
    /*! \brief This method performs %BlockDavidson iterations until the status
     * test indicates the need to stop or an error occurs (in which case, an
     * exception is thrown).
     *
     * iterate() will first determine whether the solver is inintialized; if
     * not, it will call initialize() using default arguments. After
     * initialization, the solver performs %BlockDavidson iterations until the
     * status test evaluates as Passed, at which point the method returns to
     * the caller. 
     *
     * Possible exceptions thrown include: finish
     */
    void iterate();
    //@}

    /*! \brief Initialize the solver to an iterate, optionally providing the
     * Ritz values, residual, and search direction.
     *
     * The %BlockDavidson eigensolver contains a certain amount of state
     * relating to the current eigenvectors, including the current residual,
     * the current Krylov basis, and the images of these under the operators.
     *
     * initialize() gives the user the opportunity to manually set these,
     * although this must be done with caution, abiding by the rules given
     * below. All notions of orthogonality and orthonormality are derived from
     * the inner product specified by the orthogonalization manager.
     *
     * \post 
     * <li>isInitialized() == true (see post-conditions of isInitialize())
     *
     * The user has the option of specifying any component of the state using
     * initialize(). However, these arguments are assumed to match the
     * post-conditions specified under isInitialized(). Any component of the
     * state (i.e., KX) not given to initialize() will be generated.
     */
    void initialize(BlockDavidsonState<ScalarType,MV> state);
    void initialize();

    /*! \brief Indicates whether the solver has been initialized or not.
     *
     * \return bool indicating the state of the solver.
     * \post
     * <ul>
     * <li> finish
     * </ul>
     */
    bool isInitialized() { return _initialized; }

    /*! \brief Get the current state of the eigensolver.
     * 
     * The data is only valid if isInitialized() == \c true. 
     *
     * The data for the preconditioned residual is only meaningful in the
     * scenario that the solver throws an ::BlockDavidsonRitzFailure exception
     * during iterate().
     *
     * \returns A BlockDavidsonState object containing const pointers to the current
     * solver state.
     */
    BlockDavidsonState<ScalarType,MV> getState() const {
      BlockDavidsonState<ScalarType,MV> state;
      state.V = _V;
      state.KV = _KV;
      state.X = _X;
      state.KX = _KX;
      state.H = _H;
      state.KH = _KH;
      state.R = _R;
      state.KK = _KK;
      state.MM = _MM;
      state.T = Teuchos::rcp(new std::vector<MagnitudeType>(_theta));
      if (_hasM) {
        state.MV = _MV;
        state.MX = _MX;
        state.MH = _MH;
      }
      else {
        state.MV = Teuchos::null;
        state.MX = Teuchos::null;
        state.MH = Teuchos::null;
      }
      return state;
    }

    //@}

    //@{ \name Status methods.

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

    /*! \brief Get the Ritz values for the previous iteration.
     *
     *  \return A vector of length not exceeding 3*getBlockSize() containing the Ritz values from the
     *  previous projected eigensolve.
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getRitzValues() { 
      std::vector<MagnitudeType> ret = _theta;
      ret.resize(_curDim);
      return ret;
    }

    /*! \brief Get the current residual norms
     *
     *  \return A vector of length blockSize containing the norms of the
     *  residuals, with respect to the orthogonalization manager norm() method.
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getResNorms()    {return _Rnorms;}


    /*! \brief Get the current residual 2-norms
     *
     *  \return A vector of length blockSize containing the 2-norms of the
     *  residuals. 
     */
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> getRes2Norms()   {return _R2norms;}

    //@{ \name Accessor routines


    //! Get a constant reference to the eigenvalue problem.
    const Eigenproblem<ScalarType,MV,OP>& getProblem() const { return(*_problem); };


    /*! \brief Set the blocksize and number of blocks to be used by the
     * iterative solver in solving this eigenproblem.
     *  
     *  If the block size is reduced, then the new iterate (and residual and
     *  search direction) are chosen as the subset of the current iterate
     *  preferred by the sort manager.  Otherwise, the solver state is set to
     *  uninitialized.
     */
    void setSize(int blockSize, int numBlocks);

    //! \brief Set the blocksize. 
    void setBlockSize(int blockSize);

    //! Get the size of the basis: blockSize * numBlocks
    int getSize() const { return(_blockSize*_numBlocks); }


    //! Get the blocksize to be used by the iterative solver in solving this eigenproblem.
    int getBlockSize() const { return(_blockSize); }



    /*! \brief Set the auxilliary vectors for the solver.
     *
     *  Because the current iterate X and search direction P cannot be assumed
     *  orthogonal to the new auxilliary vectors, a call to setAuxVecs() may
     *  reset the solver to the uninitialized state. This happens only in the
     *  case where the user requests full orthogonalization when the solver is
     *  in an initialized state with full orthogonalization disabled.
     *
     *  In order to preserve the current iterate, the user will need to extract
     *  it from the solver using getEvecs(), orthogonalize it against the new
     *  auxilliary vectors, and manually reinitialize the solver using
     *  initialize().
     */
    void setAuxVecs(const Teuchos::Array<Teuchos::RefCountPtr<const MV> > &auxvecs);

    //! Get the auxilliary vectors for the solver.
    Teuchos::Array<Teuchos::RefCountPtr<const MV> > getAuxVecs() const {return _auxVecs;}

    //@}

    //@{ \name Output methods.
    
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
    typedef typename std::vector<ScalarType>::iterator STiter;
    typedef typename std::vector<MagnitudeType>::iterator MTiter;
    const MagnitudeType ONE;  
    const MagnitudeType ZERO; 
    const MagnitudeType NANVAL;
    //
    // Internal structs
    //
    struct CheckList {
      bool checkV, checkMV, checkKV;
      bool checkX, checkMX, checkKX;
      bool checkH, checkMH;
      bool checkR, checkQ;
      CheckList() : checkV(false),checkMV(false),checkKV(false),
                    checkX(false),checkMX(false),checkKX(false),
                    checkH(false),checkMH(false),
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
                                        _timerSortEval, _timerDS,
                                        _timerOrtho, _timerTotal;
    //
    // Counters
    //
    int _count_ApplyOp, _count_ApplyM, _count_ApplyPrec;

    //
    // Algorithmic parameters.
    //
    // _blockSize is the solver block size; it controls the number of eigenvectors that 
    // we compute, the number of residual vectors that we compute, and therefore the number
    // of vectors added to the basis on each iteration.
    int _blockSize;
    // _numBlocks is the size of the allocated space for the Krylov basis, in blocks.
    int _numBlocks; 
    
    // 
    // Current solver state
    //
    // _initialized specifies that the basis vectors have been initialized and the iterate() routine
    // is capable of running; _initialize is controlled  by the initialize() member method
    // For the implications of the state of _initialized, please see documentation for initialize()
    bool _initialized;
    //
    // _curDim reflects how much of the current basis is valid 
    // NOTE: 0 <= _curDim <= _blockSize*_numBlocks
    // this also tells us how many of the values in _theta are valid Ritz values
    int _curDim;
    //
    // State Multivecs
    Teuchos::RefCountPtr<MV> _X, _KX, _MX, _R,
                             _H, _KH, _MH,
                             _V, _KV, _MV;
    //
    // Projected matrices
    //
    Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > _KK, _MM;
    // 
    // Auxilliary vectors
    Teuchos::Array<Teuchos::RefCountPtr<const MV> > _auxVecs;
    int _numAuxVecs;
    //
    // Number of iterations that have been performed.
    int _iter;
    // 
    // Current eigenvalues, residual norms
    std::vector<MagnitudeType> _theta, _Rnorms, _R2norms;

  };


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructor
  template <class ScalarType, class MV, class OP>
  BlockDavidson<ScalarType,MV,OP>::BlockDavidson(
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
    _timerSortEval(Teuchos::TimeMonitor::getNewTimer("Sorting eigenvalues")),
    _timerDS(Teuchos::TimeMonitor::getNewTimer("Direct solve")),
    _timerOrtho(Teuchos::TimeMonitor::getNewTimer("Orthogonalization")),
    _timerTotal(Teuchos::TimeMonitor::getNewTimer("Total time")),
    _count_ApplyOp(0),
    _count_ApplyM(0),
    _count_ApplyPrec(0),
    // internal data
    _iter(0), 
    _initialized(false),
    _curDim(0),
    _auxVecs( Teuchos::Array<Teuchos::RefCountPtr<const MV> >(0) ), 
    _numAuxVecs(0)
  {     
    TEST_FOR_EXCEPTION(_problem == Teuchos::null,std::logic_error,
                       "Anasazi::BlockDavidson::constructor: user specified null problem pointer.");
    TEST_FOR_EXCEPTION(_problem->isProblemSet() == false, std::logic_error,
                       "Anasazi::BlockDavidson::constructor: user specified problem is not set.");
    TEST_FOR_EXCEPTION(_problem->isHermitian() == false, std::logic_error,
                       "Anasazi::BlockDavidson::constructor: user specified problem is not hermitian.");

    // set the block size and allocate data
    _blockSize = 0;
    _numBlocks = 0;
    int bs = params.get("Block Size", _problem->getNEV());
    int nb = params.get("Num Blocks", 1);
    setSize(bs,nb);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the block size
  // This simply calls setSize(), modifying the block size while retaining the number of blocks.
  template <class ScalarType, class MV, class OP>
  void BlockDavidson<ScalarType,MV,OP>::setBlockSize (int blockSize) 
  {
    setSize(blockSize,_numBlocks);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the block size and make necessary adjustments.
  template <class ScalarType, class MV, class OP>
  void BlockDavidson<ScalarType,MV,OP>::setSize (int blockSize, int numBlocks) 
  {
    // This routine only allocates space; it doesn't not perform any computation
    // if size is decreased, take the first blockSize vectors of all and leave state as is
    // otherwise, grow/allocate space and set solver to unitialized

    TEST_FOR_EXCEPTION(numBlocks <= 0 || blockSize <= 0, std::logic_error, "Anasazi::BlockDavidson::setSize was passed a non-positive argument.");
    if (blockSize == _blockSize && numBlocks == _numBlocks) {
      // do nothing
      return;
    }

    // We will perform no significant computation in this routine.

    // We are bound by the requirement that X contains the preferred Ritz
    // vectors with respect to the basis V. Therefore, if we modify V
    // (via truncation) we must invalidate the state of the solver.
    // Also, if we must enlarge X, we must invalidate the state of the 
    // solver.

    // handle those things dependant only on blockSize:
    // X,KX,MX, R, H,KH,MH, _Rnorms,_R2norms
    if (blockSize < _blockSize) {
      //
      // shrink vectors
      //
      // H,KH,MH have no state; just shrink them
      _H = MVT::Clone(*_H,blockSize);
      _KH = MVT::Clone(*_H,blockSize);
      if (_hasM) {
        _MH = MVT::Clone(*_H,blockSize);
      }
      else {
        _MH = _H;
      }

      // shrink the vectors for norms and values
      _Rnorms.resize(blockSize);
      _R2norms.resize(blockSize);

      // handle X,KX,MX,R
      if (_initialized) {
        // shrink multivectors with copy
        // create ind = {0, 1, ..., blockSize-1}
        std::vector<int> ind(blockSize);
        for (int i=0; i<blockSize; i++) ind[i] = i;
        
        _X  = MVT::CloneCopy(*_X,ind);
        _KX = MVT::CloneCopy(*_KX,ind);
        if (_hasM) {
          _MX = MVT::CloneCopy(*_MX,ind);
        }
        else {
          _MX = _X;
        }
        _R  = MVT::CloneCopy(*_R,ind);
      }
      else {
        // shrink multivectors without copying
        _X = MVT::Clone(*_X,blockSize);
        _KX = MVT::Clone(*_KX,blockSize);
        if (_hasM) {
          _MX = MVT::Clone(*_MX,blockSize);
        }
        else {
          _MX = _X;
        }
        _R = MVT::Clone(*_R,blockSize);
      }
    }
    else {  // blockSize > _blockSize
      // this is also the scenario for our initial call to setSize(), in the constructor
      // in this case, _blockSize == 0 and none of the multivecs have been allocated
      _initialized = false;

      Teuchos::RefCountPtr<const MV> tmp;
      // grab some Multivector to Clone
      // in practice, getInitVec() should always provide this, but it is possible to use a 
      // Eigenproblem with nothing in getInitVec() by manually initializing with initialize(); 
      // in case of that strange scenario, we will try to Clone from _X
      if (_X != Teuchos::null) { // this is equivalent to _blockSize > 0
        tmp = _X;
      }
      else {
        tmp = _problem->getInitVec();
        TEST_FOR_EXCEPTION(tmp == Teuchos::null,std::logic_error,
                           "Anasazi::BlockDavidson::setSize(): Eigenproblem did not specify initial vectors to clone from");
      }
      // grow/allocate vectors
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
    }

    // now, handle those things dependant on blockSize and numBlocks
    // V,KV,MV, KK,MM, theta
    if (blockSize*numBlocks > _blockSize*_numBlocks) {
      // grow the basis
      _theta.resize(blockSize*numBlocks);
      int newsd = blockSize*numBlocks;

      if (_initialized) {
        // copy the old data to the new basis
        std::vector<int> ind(_curDim);
        for (int i=0; i<_curDim; i++) ind[i] = i;
        Teuchos::RefCountPtr<MV> newV;
        // V
        newV = MVT::Clone(*_X,newsd);
        MVT::SetBlock(*_V,ind,*newV);
        _V = newV;
        // KV
        newV = MVT::Clone(*_X,newsd);
        MVT::SetBlock(*_KV,ind,*newV);
        _KV = newV;
        if (_hasM) {
          // MV
          newV = MVT::Clone(*_X,newsd);
          MVT::SetBlock(*_MV,ind,*newV);
          _MV = newV;
        }
        else {
          _MV = _V;
        }
        // KK,MM
        Teuchos::RefCountPtr< Teuchos::SerialDenseMatrix<int,ScalarType> > tmp1,tmp2,tmp3;
        // create new KK and a submatrix view; then assign from old KK
        tmp1 = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(newsd,newsd) );
        tmp2 = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( 
                  Teuchos::View, *tmp1, _curDim, _curDim
               ) );
        tmp3 = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(
                  Teuchos::View, *_KK, _curDim, _curDim
               ) );
        tmp2->assign(*tmp3);
        tmp3 = Teuchos::null;
        tmp2 = Teuchos::null;
        _KK = tmp1;
        tmp1 = Teuchos::null;
        // create new MM and a submatrix view; then assign from old MM
        tmp1 = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(newsd,newsd) );
        tmp2 = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( 
                  Teuchos::View, *tmp1, _curDim, _curDim
               ) );
        tmp3 = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(
                  Teuchos::View, *_MM, _curDim, _curDim
               ) );
        tmp2->assign(*tmp3);
        tmp3 = Teuchos::null;
        tmp2 = Teuchos::null;
        _MM = tmp1;
        tmp1 = Teuchos::null;
      }
      else {
        // just allocate space
        _V = MVT::Clone(*_X,newsd);
        _KV = MVT::Clone(*_X,newsd);
        if (_hasM) {
          _MV = MVT::Clone(*_X,newsd);
        }
        else {
          _MV = _V;
        }
        _KK = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(newsd,newsd) );
        _MM = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(newsd,newsd) );
      }
    }
    else if (blockSize*numBlocks < _blockSize*_numBlocks) {
      // the space has shrunk: if we have to truncate vectors, then reset to uninitialized
      int newsd = blockSize*numBlocks;
      if (newsd < _curDim) {
        _initialized == false;
      }

      // if we are still initialized, reallocate and copy
      // otherwise, just allocate
      if (_initialized) {
        _theta.resize(blockSize*numBlocks);

        Teuchos::RefCountPtr<MV> newV;
        std::vector<int> ind(_curDim);
        for (int i=0; i<_curDim; i++) ind[i] = i;
        // V
        newV = MVT::Clone(*_X,newsd);
        _V = MVT::CloneView(*_V,ind);
        MVT::SetBlock(*_V,ind,*newV);
        _V = newV;
        // KV
        newV = MVT::Clone(*_X,newsd);
        _KV = MVT::CloneView(*_KV,ind);
        MVT::SetBlock(*_KV,ind,*newV);
        _KV = newV;
        if (_hasM) {
          // MV
          newV = MVT::Clone(*_X,newsd);
          _MV = MVT::CloneView(*_MV,ind);
          MVT::SetBlock(*_MV,ind,*newV);
          _MV = newV;
        }
        else {
          _MV = _V;
        }
        // create new KK and a submatrix view; then assign from old KK
        Teuchos::RefCountPtr< Teuchos::SerialDenseMatrix<int,ScalarType> > tmp1,tmp2,tmp3;
        tmp1 = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(newsd,newsd) );
        tmp2 = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( 
                  Teuchos::View, *tmp1, _curDim, _curDim
               ) );
        tmp3 = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(
                  Teuchos::View, *_KK, _curDim, _curDim
               ) );
        tmp2->assign(*tmp3);
        tmp3 = Teuchos::null;
        tmp2 = Teuchos::null;
        _KK = tmp1;
        tmp1 = Teuchos::null;
        // create new MM and a submatrix view; then assign from old MM
        tmp1 = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(newsd,newsd) );
        tmp2 = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( 
                  Teuchos::View, *tmp1, _curDim, _curDim
               ) );
        tmp3 = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(
                  Teuchos::View, *_MM, _curDim, _curDim
               ) );
        tmp2->assign(*tmp3);
        tmp3 = Teuchos::null;
        tmp2 = Teuchos::null;
        _MM = tmp1;
        tmp1 = Teuchos::null;
      }
      else {
        _theta.resize(blockSize*numBlocks,NANVAL);
        // just allocate space
        _V = MVT::Clone(*_X,newsd);
        _KV = MVT::Clone(*_X,newsd);
        if (_hasM) {
          _MV = MVT::Clone(*_X,newsd);
        }
        else {
          _MV = _V;
        }
        _KK = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(newsd,newsd) );
        _MM = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(newsd,newsd) );
      }
    }

    // change the local numbers
    if (_initialized) {
      // check that everything is still okay
      CheckList chk;
      chk.checkX = true;
      chk.checkKX = true;
      chk.checkMX = true;
      chk.checkV = true;
      chk.checkKV = true;
      chk.checkMV = true;
      chk.checkR = true;
      _om->print(Debug, accuracyCheck(chk, ": after setSize()") );
    }
    else {
      _curDim = 0;
    }
    _blockSize = blockSize;
    _numBlocks = numBlocks;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the auxilliary vectors
  template <class ScalarType, class MV, class OP>
  void BlockDavidson<ScalarType,MV,OP>::setAuxVecs(const Teuchos::Array<Teuchos::RefCountPtr<const MV> > &auxvecs) {
    typedef typename Teuchos::Array<Teuchos::RefCountPtr<const MV> >::iterator tarcpmv;

    // set new auxilliary vectors
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
    
    // If the solver has been initialized, X and P are not necessarily orthogonal to new auxilliary vectors
    if (_auxVecs.size() > 0 && _initialized) {
      _initialized = false;
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  /* Initialize the state of the solver
   * 
   * POST-CONDITIONS:
   *
   * V is orthonormal, orthogonal to _auxVecs, for first _curDim vectors
   * KV = Op*V
   * MV = M*V if _hasM
   * _theta contains Ritz w.r.t. V
   * X is Ritz vectors w.r.t. V
   * KX = Op*X
   * MX = M*X if _hasM
   * R = KX - MX*diag(_theta)
   */
  template <class ScalarType, class MV, class OP>
  void BlockDavidson<ScalarType,MV,OP>::initialize(BlockDavidsonState<ScalarType,MV> state)
  {
    // NOTE: memory has been allocated by setBlockSize(). Use SetBlock below; do not Clone

    std::vector<int> bsind(_blockSize);
    for (int i=0; i<_blockSize; i++) bsind[i] = i;

    // in BlockDavidson, V is primary
    // the order of dependence follows like so.
    // 0> V
    //   1> KV,MV,KK,MM
    //     2> theta,X,KX,MX
    //        3> R
    // 
    // if the user specifies all data for a level, we will accept it.
    // otherwise, we will generate the whole level, and all subsequent levels.
    //
    // these levels are ordered based on data dependence and partitioned according
    // to the amount of work required to produce the items in a level.
    //
    // inconsitent multivectors widths and lengths will not be tolerated, and
    // will be treated with exceptions.
    std::string errstr("Anasazi::BlockDavidson::initialize(): multivectors must have a consistent length and width.");

    // set up V: if the user doesn't specify V, ignore the rest
    if (state.V != Teuchos::null) {
      TEST_FOR_EXCEPTION( MVT::GetVecLength(*state.V) != MVT::GetVecLength(*_V),
                          std::logic_error, errstr );
      TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*state.V) < _blockSize,
                          std::logic_error, errstr );
      TEST_FOR_EXCEPTION( MVT::GetNumberVecs(*state.V) > _blockSize*_numBlocks,
                          std::logic_error, errstr );

      _curDim = MVT::GetNumberVecs(*state.V);
      // pick an integral amount
      _curDim = (int)(_curDim / _blockSize)*_blockSize;
      std::vector<int> nevind(_curDim);
      for (int i=0; i<_curDim; i++) nevind[i] = i;

      // put data in V,MV,KV
      MVT::SetBlock(*state.V,nevind,*_V);

      // get local views of V,MV,KV: view of first _curDim vectors
      Teuchos::RefCountPtr<MV> lclV, lclKV, lclMV;
      // generate lclV in case we need it for KK,MM below
      lclV = MVT::CloneView(*_V,nevind);

      // M*V
      if (_hasM) {
        if (state.MV != Teuchos::null ) {
            
          TEST_FOR_EXCEPTION(MVT::GetVecLength(*state.MV) != MVT::GetVecLength(*_MV),
                             std::logic_error, errstr);
          TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.MV) != _curDim,
                             std::logic_error, errstr );
          // accept user's M*V
          MVT::SetBlock(*state.MV,nevind,*_MV);
          lclMV = MVT::CloneView(*_MV,nevind);
        }
        else {
          // generate our own M*V... 
          lclMV = MVT::CloneView(*_MV,nevind);
          OPT::Apply(*_MOp,*lclV,*lclMV);
          _count_ApplyM += _curDim;
          // ...and ignore other data from user
          state.MM = state.KK = Teuchos::null;
        }
      }
      else {
        // if _hasM == false, then (should) _MV == _V
        // an assignment would be redundant here
        // take advantage of this opportunity to debug a little
        TEST_FOR_EXCEPTION(_MV != _V, std::logic_error, "Anasazi::BlockDavidson::initialize(): invariant not satisfied");
        // generate lclMV in case we need it for KK,MM below
        lclMV = lclV;
      }

      // K*V
      if (state.KV != Teuchos::null) {
        TEST_FOR_EXCEPTION(MVT::GetVecLength(*state.KV) != MVT::GetVecLength(*_KV),
                           std::logic_error, errstr );
        TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.KV) != _curDim,
                            std::logic_error, errstr );
        // accept user's K*V... 
        MVT::SetBlock(*state.KV,nevind,*_KV);
        lclKV = MVT::CloneView(*_KV,nevind);
      }
      else {
        // generate our own K*V... 
        lclKV = MVT::CloneView(*_KV,nevind);
        OPT::Apply(*_Op,*_V,*_KV);
        _count_ApplyOp += _curDim;
        // ...and ignore other data from user
        state.MM = state.KK = Teuchos::null;
      }

      // KK and MM
      if (state.MM != Teuchos::null && state.KK != Teuchos::null) {
        TEST_FOR_EXCEPTION(state.MM->numRows() != _curDim 
                           || state.MM->numCols() != _curDim 
                           || state.KK->numRows() != _curDim 
                           || state.KK->numCols() != _curDim,
                           std::logic_error, errstr);
        // copy the part we want into _KK and _MM
        Teuchos::SerialDenseMatrix<int,ScalarType>
            lclKK(Teuchos::View,*_KK,_curDim,_curDim),
            lclMM(Teuchos::View,*_MM,_curDim,_curDim);
        lclKK.assign(*state.KK);
        lclMM.assign(*state.MM);
      }
      else {
        // generate MM and KK...
        Teuchos::SerialDenseMatrix<int,ScalarType> 
            lclKK(Teuchos::View,*_KK,_curDim,_curDim),
            lclMM(Teuchos::View,*_MM,_curDim,_curDim);
        MVT::MvTransMv(ONE,*lclV,*lclKV,lclKK);
        MVT::MvTransMv(ONE,*lclV,*lclMV,lclMM);
        // ...and don't accept X,theta
        state.X = Teuchos::null;
        state.T = Teuchos::null;
      }

      // X,MX,KX,theta require Ritz analisys; no point in accepting one without the rest
      if (state.X != Teuchos::null && state.KX != Teuchos::null 
          && (!_hasM || state.MX != Teuchos::null) && state.T != Teuchos::null) {
        TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.X) != _blockSize,
                            std::logic_error, errstr );
        TEST_FOR_EXCEPTION(MVT::GetVecLength(*state.X) != MVT::GetVecLength(*_X),
                            std::logic_error, errstr );
        TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.KX) != _blockSize,
                            std::logic_error, errstr );
        TEST_FOR_EXCEPTION(MVT::GetVecLength(*state.KX) != MVT::GetVecLength(*_X),
                            std::logic_error, errstr );
        if (_hasM) {
          TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.MX) != _blockSize,
                              std::logic_error, errstr );
          TEST_FOR_EXCEPTION(MVT::GetVecLength(*state.KX) != MVT::GetVecLength(*_X),
                              std::logic_error, errstr );
        }
        TEST_FOR_EXCEPTION((signed int)(state.T->size()) != _curDim,
                            std::logic_error, errstr );
        MVT::SetBlock(*state.X,bsind,*_X);
        MVT::SetBlock(*state.KX,bsind,*_KX);
        if (_hasM) {
          MVT::SetBlock(*state.MX,bsind,*_MX);
        }
        std::copy(state.T->begin(),state.T->end(),_theta.begin());
      }
      else {
        // compute ritz vecs/vals
        Teuchos::SerialDenseMatrix<int,ScalarType> S(_curDim,_curDim);
        //
        int rank = _curDim;
        _MSUtils.directSolver(_curDim, *_KK, _MM.get(), &S, &_theta, &rank, 1);
        // we want all ritz values back
        TEST_FOR_EXCEPTION(rank != _curDim,BlockDavidsonInitFailure,
                           "Anasazi::BlockDavidson::initialize(): Not enough Ritz vectors to initialize algorithm.");
        {
          std::vector<int> _order(_curDim);
          // make a ScalarType copy of theta
          std::vector<ScalarType> _theta_st(_curDim);
          std::copy(_theta.begin(),_theta.begin()+_curDim,_theta_st.begin());
          // sort it
          _sm->sort( NULL, _curDim, &(_theta_st[0]), &_order );   // don't catch exception
          // copy back to the MagnitudeType 
          std::copy(_theta_st.begin(),_theta_st.begin()+_curDim,_theta.begin());
          // Sort the primitive ritz vectors
          Teuchos::SerialDenseMatrix<int,ScalarType> copyS( S );
          Teuchos::BLAS<int,ScalarType> blas;
          for (int i=0; i<_curDim; i++) {
            blas.COPY(_curDim, copyS[_order[i]], 1, S[i], 1);
          }
        }
        Teuchos::SerialDenseMatrix<int,ScalarType> S1(Teuchos::View,S,_curDim,_blockSize);
        // X <- lclV*S
        MVT::MvTimesMatAddMv( ONE, *lclV, S1, ZERO, *_X );
        // KX <- lclKV*S
        MVT::MvTimesMatAddMv( ONE, *lclKV, S1, ZERO, *_KX );
        if (_hasM) {
          // MX <- lclMV*S
          MVT::MvTimesMatAddMv( ONE, *lclMV, S1, ZERO, *_MX );
        }
        // we generated theta,X,KX,MX so we don't want to use the user's R
        state.R = Teuchos::null;
      }

      // done with local pointers
      lclV = Teuchos::null;
      lclKV = Teuchos::null;
      lclMV = Teuchos::null;

      // set up R
      if (state.R != Teuchos::null) {
        TEST_FOR_EXCEPTION(MVT::GetNumberVecs(*state.R) != _blockSize ,
                           std::logic_error, errstr );
        TEST_FOR_EXCEPTION(MVT::GetVecLength(*state.R) != MVT::GetVecLength(*_X),
                           std::logic_error, errstr );
        MVT::SetBlock(*state.R,bsind,*_R);
      }
      else {
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

      _initialized = true;

      if (_om->isVerbosity( Debug ) ) {
        // Check almost everything here
        CheckList chk;
        chk.checkV = true;
        chk.checkKV = true;
        chk.checkMV = true;
        chk.checkX = true;
        chk.checkKX = true;
        chk.checkMX = true;
        chk.checkR = true;
        chk.checkQ = true;
        _om->print( Debug, accuracyCheck(chk, ": after initialize()") );
      }
    }
    else {
      // user did not specify a basis V
      // generate one
      // 
      // generate something, projectAndNormalize, call myself recursively
      Teuchos::RefCountPtr<const MV> ivec = _problem->getInitVec();
      TEST_FOR_EXCEPTION(ivec == Teuchos::null,std::logic_error,
                         "Anasazi::BlockDavdison::initialize(): Eigenproblem did not specify initial vectors to clone from");

      _curDim = MVT::GetNumberVecs(*ivec);
      // pick the largest multiple of _blockSize
      _curDim = (int)(_curDim / _blockSize)*_blockSize;
      bool userand = false;
      if (_curDim < _blockSize) {
        // we need at least _blockSize vectors
        // use a random multivec
        userand = true;
        _curDim = _blockSize;
      }

      // make an index
      std::vector<int> nevind(_curDim);
      for (int i=0; i<_curDim; i++) nevind[i] = i;

      // alloc newV
      Teuchos::RefCountPtr<MV> newMV, newV = MVT::Clone(*ivec,_curDim);
      if (userand) {
        MVT::MvRandom(*newV);
      }
      else {
        // assign ivec to first part of newV
        MVT::SetBlock(*ivec,nevind,*newV);
      }

      // compute newMV if _hasM
      if (_hasM) {
        newMV = MVT::Clone(*newV,_curDim);
        OPT::Apply(*_MOp,*newV,*newMV);
        _count_ApplyM += _blockSize;
      }
      else {
        newMV = Teuchos::null;
      }

      // remove auxVecs from newV and normalize newV
      if (_auxVecs.size() > 0) {
        Teuchos::Array<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > > dummy;
        int rank = _orthman->projectAndNormalize(*newV,newMV,dummy,Teuchos::null,_auxVecs);
        TEST_FOR_EXCEPTION(rank != _curDim,BlockDavidsonInitFailure,
                           "Anasazi::BlockDavidson::initialize(): Couldn't generate initial basis of full rank.");
      }
      else {
        int rank = _orthman->normalize(*newV,newMV,Teuchos::null);
        TEST_FOR_EXCEPTION(rank != _curDim,BlockDavidsonInitFailure,
                           "Anasazi::BlockDavidson::initialize(): Couldn't generate initial basis of full rank.");
      }

      // call myself recursively
      BlockDavidsonState<ScalarType,MV> newstate;
      newstate.V = newV;
      newstate.MV = newMV;
      initialize(newstate);
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // initialize the solver with default state
  template <class ScalarType, class MV, class OP>
  void BlockDavidson<ScalarType,MV,OP>::initialize()
  {
    BlockDavidsonState<ScalarType,MV> empty;
    initialize(empty);
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Perform BlockDavidson iterations until the StatusTest tells us to stop.
  template <class ScalarType, class MV, class OP>
  void BlockDavidson<ScalarType,MV,OP>::iterate () 
  {
    //
    // Allocate/initialize data structures
    //
    if (_initialized == false) {
      initialize();
    }

    // as a data member, this would be redundant and require synchronization with
    // _blockSize and _numBlocks; we'll use a constant here.
    const int searchDim = _blockSize*_numBlocks;
    // we use this often enough...
    std::vector<int> bsind(_blockSize);
    for (int i=0; i<_blockSize; i++) { ind[i] = i; }

    //
    // The projected matrices are part of the state, but the eigenvectors can be defined
    // locally.
    //    S = Local eigenvectors         (size: _searchDim * _searchDim
    Teuchos::SerialDenseMatrix<int,ScalarType> S( _searchDim, _searchDim );


    ////////////////////////////////////////////////////////////////
    // iterate until the status test tells us to stop.
    while (_tester->checkStatus(this) != Passed) {

      _iter++;
      
      // Apply the preconditioner on the residuals: H <- Prec*R
      if (_Prec != Teuchos::null) {
        Teuchos::TimeMonitor PrecTimer( *_timerPrec );
        OPT::Apply( *_Prec, *_R, *_H );   // don't catch the exception
        _count_ApplyPrec += _blockSize;
      }
      else {
        MVT::SetBlock(*_R,ind,*_H);
      }

      // Apply the mass matrix on H
      if (_hasM) {
        Teuchos::TimeMonitor MOpTimer( *_timerMOp );
        OPT::Apply( *_MOp, *_H, *_MH);    // don't catch the exception
        _count_ApplyM += _blockSize;
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
  // checkV : V orthonormal
  //          orthogonal to auxvecs
  // checkMV: check MV == M*V
  // checkKV: check KV == K*V
  // checkX : X orthonormal
  //          orthogonal to auxvecs
  // checkMX: check MX == M*X
  // checkKX: check KX == K*X
  // checkH : H orthonormal 
  //          orthogonal to V and H and auxvecs
  // checkMH: check MH == M*H
  // checkR : check R orthogonal to X
  // checkQ : check that auxilliary vectors are actually orthonormal
  //
  // TODO: 
  //  add checkTheta 
  //
  template <class ScalarType, class MV, class OP>
  std::string BlockDavidson<ScalarType,MV,OP>::accuracyCheck( const CheckList &chk, const string &where ) const 
  {
    stringstream os;
    os.precision(2);
    os.setf(ios::scientific, ios::floatfield);
    MagnitudeType tmp;

    os << " Debugging checks: iteration " << _iter << where << endl;

    // V and friends
    std::vector<int> lclind(_curDim);
    for (int i=0; i<_curDim; i++) lclind[i] = i;
    Teuchos::RefCountPtr<MV> lclV,lclMV,lclKV;
    lclV = MVT::CloneView(*_V,lclind);
    if (chk.checkV) {
      tmp = _orthman->orthonormError(*lclV);
      os << " >> Error in V^H M V == I : " << tmp << endl;
      for (unsigned int i=0; i<_auxVecs.size(); i++) {
        tmp = _orthman->orthogError(*lclV,*_auxVecs[i]);
        os << " >> Error in V^H M Q[" << i << "] == 0 : " << tmp << endl;
      }
    }
    if (chk.checkMV && _hasM) {
      lclMV = MVT::CloneView(*_MV,lclind);
      tmp = _MSUtils.errorEquality(lclV.get(), lclMV.get(), _MOp.get());
      os << " >> Error in MV == M*V    : " << tmp << endl;
    }
    if (chk.checkKV) {
      lclKV = MVT::CloneView(*_KV,lclind);
      tmp = _MSUtils.errorEquality(lclV.get(), lclKV.get(), _Op.get());
      os << " >> Error in KV == K*V    : " << tmp << endl;
    }

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

    // H and friends
    if (chk.checkH) {
      tmp = _orthman->orthonormError(*_H);
      os << " >> Error in H^H M H == I : " << tmp << endl;
      tmp = _orthman->orthogError(*_H,*lclV);
      os << " >> Error in H^H M V == 0 : " << tmp << endl;
      tmp = _orthman->orthogError(*_H,*_X);
      os << " >> Error in H^H M X == 0 : " << tmp << endl;
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
  BlockDavidson<ScalarType,MV,OP>::currentStatus(ostream &os) 
  {
    os.setf(ios::scientific, ios::floatfield);
    os.precision(6);
    os <<endl;
    os <<"******************* CURRENT STATUS *******************"<<endl;
    os <<"The number of iterations performed is " <<_iter<<endl;
    os <<"The block size is         " << _blockSize<<endl;
    os <<"The number of blocks is   " << _numBlocks<<endl;
    os <<"The current basis size is " << _curDim<<endl;
    os <<"The number of operations Op*x   is "<<_count_ApplyOp<<endl;
    os <<"The number of operations M*x    is "<<_count_ApplyM<<endl;
    os <<"The number of operations Prec*x is "<<_count_ApplyPrec<<endl;
    os <<endl;

    os.setf(ios_base::right, ios_base::adjustfield);

    os <<"CURRENT EIGENVALUE ESTIMATES             "<<endl;
    os << std::setw(16) << "Ritz value" 
        << std::setw(16) << "Residual"
        << endl;
    os <<"------------------------------------------------------"<<endl;
    if ( _iter > 0 || _curDim > 0 ) {
      for (int i=0; i<_blockSize; i++) {
        os << std::setw(16) << _theta[i] 
           << std::setw(16) << _Rnorms[i]
           << endl;
      }
    } 
    else {
      os <<"[none computed]" << endl;
    }
    os << "******************************************************" << endl;  
    os << endl; 
  }




  // finish





  /*

  template <class ScalarType, class MV, class OP>
  ReturnType 
  BlockDavidson<ScalarType,MV,OP>::solve () 
  {
    int info, nb, lwork;
    int bStart = 0, offSet = 0;
    int nFound = _blockSize;
    bool reStart = false;
    bool criticalExit = false;
    ReturnType ret;
    Teuchos::RefCountPtr<MV> Xcurrent, Xprev, Xtotal, Xnext;
    ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    Teuchos::BLAS<int,ScalarType> blas;
    Teuchos::LAPACK<int,ScalarType> lapack;
    std::vector<int> _order;

    // Work vector for GEQRF and ORGQR
    std::vector<ScalarType> tau( _dimSearch );
    
    //
    // Compute the required workspace for this factorization.
    //
    lwork = lapack.ILAENV( 1, "geqrf", "", _numBlocks*_blockSize, _numBlocks*_blockSize );
    lwork *= _blockSize;
    std::vector<ScalarType> work( lwork );

    for (nb = bStart; nb < maxBlock ; nb++) 
    {
      // 
      // Increment iteration counter
      //
      _iter++;
      //
      // Get a view of the current vectors being worked on.
      //
      int localSize = nb*_blockSize;
      std::vector<int> index( _blockSize );
      for (i=0; i < _blockSize; i++)
        index[i] = localSize + i;
      //
      Xcurrent = MVT::CloneView( *X, index );
      //
      if (localSize > 0) {
        index.resize( localSize );
        for (i=0; i < localSize; i++)
          index[i] = i;
        
        Xprev = MVT::CloneView( *X, index );
      }
      //
      // Apply the mass matrix.
      //        
      if (_MOp.get()) {
        {
          Teuchos::TimeMonitor MassTimer(*_timerMOp);
          OPT::Apply( *_MOp, *Xcurrent, *MX );
        }
        _count_ApplyM += MVT::GetNumberVecs( *Xcurrent );
      }
      //
      // Orthonormalize Xcurrent against the known eigenvectors and previous vectors.
      //
      if (nb == bStart) {
        if (nFound > 0) {
          if (_knownEV == 0) {
            {
              Teuchos::TimeMonitor OrtoTimer(*_timerOrtho);
              ret = _orthman->normalize(*Xcurrent,MX,Teuchos::null,info);
            }
          }
          else {
            {
              Teuchos::TimeMonitor OrthoTimer(*_timerOrtho);
              ret = _orthman->projectAndNormalize(*Xcurrent,MX,
                                                  Teuchos::tuple<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > >(Teuchos::null),
                                                  Teuchos::null,
                                                  Teuchos::tuple<Teuchos::RefCountPtr<const MV> >(Xprev),
                                                  info);
            }
          }
        }
        nFound = 0;
      } 
      else {
        {
          Teuchos::TimeMonitor OrthoTimer(*_timerOrtho);
          ret = _orthman->projectAndNormalize(*Xcurrent,MX,
                                              Teuchos::tuple<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > >(Teuchos::null),
                                              Teuchos::null,
                                              Teuchos::tuple<Teuchos::RefCountPtr<const MV> >(Xprev),
                                              info);
        }
      }
      //
      // Exit the code if there has been a problem.
      //
      if (ret != Ok) {
        return Failed;
      }
      //
      // Check orthogonality of X ( if required )
      //
      if (_om->isVerbosity( Anasazi::OrthoDetails ) ) {
        if (localSize > 0) {
          accuracyCheck( Xcurrent, MX, Xprev);
        }
        else {
          accuracyCheck( Xcurrent, MX, Teuchos::null );
        }
      }
      //
      // Apply the stiffness matrix.
      //
      {
        Teuchos::TimeMonitor OpTimer(*_timerOp);
        OPT::Apply( *_Op, *Xcurrent, *KX );
      }
      _count_ApplyOp += MVT::GetNumberVecs( *Xcurrent );
      //
      // Update the local stiffness matrix ( Xtotal^T * K * Xcurrent where Xtotal = [Xprev Xcurrent] )
      // Note:  Only the upper half of the matrix is stored in KK
      //
      index.resize( localSize + _blockSize );
      for (i=0; i < localSize + _blockSize; i++) {
        index[i] = _knownEV + i;
      }
      Xtotal = MVT::CloneView( *X, index );
      Teuchos::SerialDenseMatrix<int,ScalarType> subKK( Teuchos::View, KK, localSize+_blockSize, _blockSize, 0, localSize );
      MVT::MvTransMv( one, *Xtotal, *KX, subKK );
      //
      // Perform spectral decomposition
      //
      int nevLocal = localSize+_blockSize;
      {
        Teuchos::TimeMonitor DSTimer(*_timerDS);
        info = _MSUtils.directSolver(localSize+_blockSize, KK, 0, &S, &_theta, &nevLocal, 10);
      }
      //
      // Exit the code if there has been a problem.
      //
      if (info != 0) {
        //
        // Detect critical failure
        // 
        if (info<0) {
          criticalExit = true;
          break;
        }
        //
        // Restart:  Spectral decomposition failed.
        // Reinitialize all counters and randomize the starting block
        //
        reStart = true;
        _numRestarts++;
        index.resize( _blockSize );
        for (i=0; i<_blockSize; i++) {
          index[i] = _knownEV + i;
        }
        Teuchos::RefCountPtr<MV> Xinit = MVT::CloneView( *X, index );
        MVT::MvRandom( *Xinit );
        nFound = _blockSize;
        bStart = 0;
        break;
      }
      //
      //---------------------------------------------------
      // Sort the ritz values using the sort manager
      //---------------------------------------------------
      // The sort manager is templated on ScalarType
      // Make a ScalarType copy of _theta for sorting
      std::vector<ScalarType> _theta_st(_theta.size());
      {
        Teuchos::TimeMonitor SortTimer(*_timerSortEval);
        std::copy<MTiter,STiter>(_theta.begin(),_theta.begin()+localSize+_blockSize,_theta_st.begin());
        _order.resize(localSize+_blockSize);
        ret = _sm->sort( this, localSize+_blockSize, &(_theta_st[0]), &_order );

        // Reorder _theta according to sorting results from _theta_st
        std::vector<MagnitudeType> _theta_copy(_theta);
        for (i=0; i<localSize+_blockSize; i++) {
          _theta[i] = _theta_copy[_order[i]];
        }
      }
      if (ret != Ok) {
        _om->print(Anasazi::Error, "ERROR : Sorting in solve()\n" );
        criticalExit = true;
        break; // break out of for(nb) loop
      }
      // Sort the primitive ritz vectors
      // We need the first _blockSize vectors ordered to generate the next
      // columns immediately below, as well as when/if we restart.
      Teuchos::SerialDenseMatrix<int,ScalarType> copyS( S );
      for (i=0; i<localSize+_blockSize; i++) {
        blas.COPY(localSize+_blockSize, copyS[_order[i]], 1, S[i], 1);
      }
      // Create a view matrix of the first _blockSize vectors
      Teuchos::SerialDenseMatrix<int,ScalarType> subS( Teuchos::View, S, localSize+_blockSize, _blockSize );
      //
      //---------------------------------------------------
      // Update the search space 
      //---------------------------------------------------
      // KX = Xtotal * S where S is the eigenvectors of the projected problem.
      //
      MVT::MvTimesMatAddMv( one, *Xtotal, subS, zero, *KX );
      //
      // Apply the mass matrix for the next block
      // 
      if (_MOp.get()) {
        {
          Teuchos::TimeMonitor MOpTimer(*_timerMOp);
          OPT::Apply( *_MOp, *KX, *MX );
        }
        _count_ApplyM += MVT::GetNumberVecs( *KX );
      }
      //
      // Apply the stiffness matrix for the next block
      //
      {
        Teuchos::TimeMonitor OpTimer(*_timerOp);
        OPT::Apply( *_Op, *KX, *R );
      }
      _count_ApplyOp += MVT::GetNumberVecs( *KX );
      //
      // Compute the residual :
      // R = KX - MX*diag(theta)
      // 
      Teuchos::SerialDenseMatrix<int,ScalarType> D(_blockSize, _blockSize);
      for (i=0; i<_blockSize; i++ ) {
        D(i,i) = -_theta[i];
      }
      //
      if (_MOp.get()) {
        MVT::MvTimesMatAddMv( one, *MX, D, one, *R );
      }
      else {
        MVT::MvTimesMatAddMv( one, *KX, D, one, *R );
      }
      MVT::MvNorm( *R, &_normR );
      //
      // Scale the norms of residuals with the eigenvalues and check for converged eigenvectors.
      //
      nFound = 0;
      for (j=0; j<_blockSize; j++) {
        // Scale eigenvalues if _theta is non-zero.
        if ( _theta[j] != zero ) {
          _normR[j] /= SCT::magnitude(_theta[j]);
        }
        // Check for convergence
        if (_normR[j] < _residual_tolerance) {
          nFound ++;          
        }
      }
      // Print information on current iteration
      _om->stream( Anasazi::IterationDetails)
            << " Iteration " << _iter << " - Number of converged eigenvectors "
            << _knownEV + nFound << endl;
      
      if (_om->isVerbosity( Anasazi::IterationDetails )) {
        stringstream os;
        os << endl;
        os.precision(2);
        os.setf(ios::scientific, ios::floatfield);
        for (i=0; i<_blockSize; ++i) {
          os << " Iteration " << _iter << " - Scaled Norm of Residual " << i;
          os << " = " << _normR[i] << endl;
        }
        os << endl;
        os.precision(2);
        for (i=0; i<localSize + _blockSize; ++i) {
          os << " Iteration "<< _iter << " - Ritz eigenvalue " << i;
          os.setf((fabs(_theta[i]) < 0.01) ? ios::scientific : ios::fixed, ios::floatfield);
          os << " = " << _theta[i] << endl;
        }
        os << endl;
        
        _om->print(Anasazi::IterationDetails,os.str());
      }
      //
      // Exit the loop to treat the converged eigenvectors
      //
      if (nFound > 0) {
        nb += 1;
        offSet = 0;
        break;
      }
      //
      // Apply the preconditioner on the residuals
      //
      if (maxBlock == 1) {
        if (_Prec.get()) {
          {
            Teuchos::TimeMonitor PrecTimer(*_timerPrec);
            OPT::Apply( *_Prec, *R, *Xcurrent );
          }
          _count_ApplyPrec += MVT::GetNumberVecs( *R );
        }
        else
          MVT::MvAddMv( one, *R, zero, *R, *Xcurrent );
        //
        // Update the preconditioned residual 
        //
        MVT::MvAddMv( one, *KX, -one, *Xcurrent, *Xcurrent );
        break;
      } // if (maxBlock == 1)

      if (nb == maxBlock - 1) {
        nb++;
        break;
      }
      //
      // Prepare next block in factorization
      //
      index.resize( _blockSize );
      for( i=0; i<_blockSize; i++) {
        index[i] = _knownEV + localSize + _blockSize + i;
      }
      Xnext = MVT::CloneView( *X, index );
      if (_Prec.get()) {
        {
          Teuchos::TimeMonitor PrecTimer(*_timerPrec);
          OPT::Apply( *_Prec, *R, *Xnext );
        }
        _count_ApplyPrec += MVT::GetNumberVecs( *R );
      }
      else {
        MVT::MvAddMv( one, *R, zero, *R, *Xnext );
      }
    } // for (nb = bStart; nb < maxBlock; nb++)

  } // end solve()

  */

  
  } // End of namespace Anasazi

#endif

// End of file AnasaziBlockDavidson.hpp
