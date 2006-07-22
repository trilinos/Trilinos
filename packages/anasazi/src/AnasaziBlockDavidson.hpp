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
  template <class ScalarType, class MV>
  struct BlockDavidsonState {
    Teuchos::RefCountPtr<const MV> X, KX, MX;
    Teuchos::RefCountPtr<const MV> V, KV, MV;
    Teuchos::RefCountPtr<const MV> H, KH, MH;
    Teuchos::RefCountPtr<const MV> R;
    Teuchos::RefCountPtr<const std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType> > T;
    BlockDavidson() : X(Teuchos::null),KX(Teuchos::null),MX(Teuchos::null),
                      V(Teuchos::null),KV(Teuchos::null),MV(Teuchos::null),
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
    void solve();
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
    void initialize(LOBPCGState<ScalarType,MV> state);
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
     * scenario that the solver throws an ::LOBPCGRitzFailure exception
     * during iterate().
     *
     * \returns A BlockDavidsonState object containing const pointers to the current
     * solver state.
     */
    BlockDavidsonState<ScalarType,MV> getState() const {
      BlockDavidsonState<ScalarType,MV> state;
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
      ret.resize(_nevLocal);
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
    //
    // Internal structs
    //
    struct CheckList {
      bool checkX, checkMX, checkKX;
      bool checkH, checkMH;
      bool checkV, checkMV, checkKV;
      bool checkR, checkQ;
      CheckList() : checkX(false),checkMX(false),checkKX(false),
                    checkH(false),checkMH(false),
                    checkV(false),checkMV(false),checkKV(false),
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
    // _nevLocal reflects how much of the current basis is valid 
    // NOTE: 0 <= _nevLocal <= _blockSize*_numBlocks
    // this also tells us how many of the values in _theta are valid Ritz values
    int _nevLocal;
    //
    // State Multivecs
    Teuchos::RefCountPtr<MV> _X, _KX, _MX, _R,
                             _H, _KH, _MH,
                             _V, _KV, _MV;
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
    _fullOrtho(params.get("Full Ortho", true)),
    _initialized(false),
    _nevLocal(0),
    _hasP(false),
    _auxVecs( Teuchos::Array<Teuchos::RefCountPtr<const MV> >(0) ), 
    _numAuxVecs(0)
  {     
    // set the block size and allocate data
    _blockSize = 0;
    _numBlocks = 0;
    int bs = params.get("Block Size", _problem->getNEV());
    int nb = params.get("Num Blocks", 1);
    setSize(bs,nb);
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
      // this is also the scenario for our initial call to setBlockSize(), in the constructor
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
    // V,KV,MV, _theta
    if (blockSize*numBlocks > _blockSize*_numBlocks) {
      // grow the basis
      _theta.resize(blockSize*numBlocks);
      int newsd = blockSize*numBlocks;

      if (_initialized) {
        // copy the old data to the new basis
        std::vector<int> ind(_nevLocal);
        for (int i=0; i<_nevLocal; i++) ind[i] = i;
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
      }
    }
    else if (blockSize*numBlocks < _blockSize*_numBlocks) {
      // the space has shrunk: if we have to truncate vectors, then reset to uninitialized
      int newsd = blockSize*numBlocks;
      if (newsd < _nevLocal) {
        _initialized == false;
      }

      // if we are still initialized, reallocate and copy
      // otherwise, just allocate
      if (_initialized) {
        _theta.resize(blockSize*numBlocks);

        Techos::RefCountPtr<MV> newV;
        std::vector<int> ind(_nevLocal);
        for (int i=0; i<_nevLocal; i++) ind[i] = i;
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
      }
    }

    // change the local numbers
    if (_initialized) {
      // check that everything is still okay
      CheckList chk;
      chk.X = true;
      chk.KX = true;
      chk.MX = true;
      chk.V = true;
      chk.KV = true;
      chk.MV = true;
      chk.R = true;
      _om->print(Debug, accuracyCheck(chk, ": after setSize()") );
    }
    else {
      _nevLocal = 0;
    }
    _blockSize = blockSize;
    _numBlocks = numBlocks;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Set the auxilliary vectors
  template <class ScalarType, class MV, class OP>
  void LOBPCG<ScalarType,MV,OP>::setAuxVecs(const Teuchos::Array<Teuchos::RefCountPtr<const MV> > &auxvecs) {
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
      _hasP = false;
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  /* Initialize the state of the solver
   * 
   * POST-CONDITIONS:
   *
   * V is orthonormal, orthogonal to _auxVecs, for first _nevLocal vectors
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
  // finish: finish writing this method
  {
    // NOTE: memory has been allocated by setBlockSize(). Use SetBlock below; do not Clone

    std::vector<int> bsind(_blockSize);
    for (int i=0; i<_blockSize; i++) bsind[i] = i;

    // set up X: if the user doesn't specify X, ignore the rest
    if (state.X != Teuchos::null && MVT::GetNumberVecs(*state.X) >= _blockSize && MVT::GetVecLength(*state.X) == MVT::GetVecLength(*_X) ) {

      // put data in X,MX,KX
      MVT::SetBlock(*state.X,bsind,*_X);
      if (_hasM) {
        if (state.MX != Teuchos::null && MVT::GetNumberVecs(*state.MX) >= _blockSize && MVT::GetVecLength(*state.MX) == MVT::GetVecLength(*_MX) ) {
          MVT::SetBlock(*state.MX,bsind,*_MX);
        }
        else {
          OPT::Apply(*_MOp,*_X,*_MX);
          _count_ApplyM += _blockSize;
        }
      }
      else {
        // an assignment would be redundant; take advantage of this opportunity to debug a little
        TEST_FOR_EXCEPTION(_MX != _X, std::logic_error, "Anasazi::LOBPCG::initialize(): invariant not satisfied");
      }
      if (state.KX != Teuchos::null && MVT::GetNumberVecs(*state.KX) >= _blockSize && MVT::GetVecLength(*state.KX) == MVT::GetVecLength(*_KX) ) {
        MVT::SetBlock(*state.KX,bsind,*_KX);
      }
      else {
        OPT::Apply(*_Op,*_X,*_KX);
        _count_ApplyOp += _blockSize;
      }

      // set up Ritz values
      _theta.resize(3*_blockSize,NANVAL);
      if (state.T != Teuchos::null && (signed int)(state.T->size()) >= _blockSize) {
        for (int i=0; i<_blockSize; i++) {
          _theta[i] = (*state.T)[i];
        }
      }
      else {
        // get ritz vecs/vals
        Teuchos::SerialDenseMatrix<int,ScalarType> KK(_blockSize,_blockSize),
                                                   MM(_blockSize,_blockSize),
                                                    S(_blockSize,_blockSize);
        // project K
        MVT::MvTransMv(ONE,*_X,*_KX,KK);
        // project M
        MVT::MvTransMv(ONE,*_X,*_MX,MM);
        _nevLocal = _blockSize;
        _MSUtils.directSolver(_blockSize, KK, &MM, &S, &_theta, &_nevLocal, 1);
        TEST_FOR_EXCEPTION(_nevLocal < _blockSize,LOBPCGInitFailure,
                           "Anasazi::LOBPCG::initialize(): Not enough Ritz vectors to initialize algorithm.");
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
  
      // set up R
      if (state.R != Teuchos::null && MVT::GetNumberVecs(*state.R) >= _blockSize && MVT::GetVecLength(*state.R) == MVT::GetVecLength(*_R) ) {
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

  
      // put data in P,KP,MP: P is not used to set theta
      if (state.P != Teuchos::null && MVT::GetNumberVecs(*state.P) >= _blockSize && MVT::GetVecLength(*state.P) == MVT::GetVecLength(*_P) ) {
        _hasP = true;

        MVT::SetBlock(*state.P,bsind,*_P);

        if (state.KP != Teuchos::null && MVT::GetNumberVecs(*state.KP) >= _blockSize && MVT::GetVecLength(*state.KP) == MVT::GetVecLength(*_KP) ) {
          MVT::SetBlock(*state.KP,bsind,*_KP);
        }
        else {
          OPT::Apply(*_Op,*_P,*_KP);
          _count_ApplyOp += _blockSize;
        }

        if (_hasM) {
          if (state.MP != Teuchos::null && MVT::GetNumberVecs(*state.MP) >= _blockSize && MVT::GetVecLength(*state.MP) == MVT::GetVecLength(*_MP) ) {
            MVT::SetBlock(*state.MP,bsind,*_MP);
          }
          else {
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
        OPT::Apply(*_MOp,*newX,*newMX);
        _count_ApplyM += _blockSize;
      }
      else {
        newMX = Teuchos::null;
      }

      // remove auxVecs from newX and normalize newX
      if (_auxVecs.size() > 0) {
        Teuchos::Array<Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > > dummy;
        int rank = _orthman->projectAndNormalize(*newX,newMX,dummy,Teuchos::null,_auxVecs);
        TEST_FOR_EXCEPTION(rank != _blockSize,LOBPCGInitFailure,
                           "Anasazi::LOBPCG::initialize(): Couldn't generate initial basis of full rank.");
      }
      else {
        int rank = _orthman->normalize(*newX,newMX,Teuchos::null);
        TEST_FOR_EXCEPTION(rank != _blockSize,LOBPCGInitFailure,
                           "Anasazi::LOBPCG::initialize(): Couldn't generate initial basis of full rank.");
      }

      // call myself recursively
      LOBPCGState<ScalarType,MV> newstate;
      newstate.X = newX;
      newstate.MX = newMX;
      initialize(newstate);
    }
  }

  template <class ScalarType, class MV, class OP>
  void LOBPCG<ScalarType,MV,OP>::initialize()
  {
    LOBPCGState<ScalarType,MV> empty;
    initialize(empty);
  }












  template <class ScalarType, class MV, class OP>
  ReturnType 
  BlockDavidson<ScalarType,MV,OP>::solve () 
  {
    Teuchos::TimeMonitor LocalTimer(*_timerTotal,_restartTimers);
    
    if ( _restartTimers ) {
      _timerOp->reset();
      _timerMOp->reset();
      _timerPrec->reset();
      _timerSortEval->reset();
      _timerDS->reset();
      _timerOrtho->reset();
      _count_ApplyOp = 0;
      _count_ApplyM = 0;
      _count_ApplyPrec = 0;
    }

    //
    // Check the Anasazi::Eigenproblem was set by user, if not, return failed.
    //
    if ( !_problem->IsProblemSet() ) {
      _om->print(Anasazi::Error,"ERROR : Anasazi::Eigenproblem was not set, call Anasazi::Eigenproblem::SetProblem() before calling solve\n");
      return Failed;
    }
    //
    // Check the Anasazi::Eigenproblem is symmetric, if not, return failed.
    //
    if ( !_problem->IsSymmetric() ) {
      _om->print(Anasazi::Error, "ERROR : Anasazi::Eigenproblem is not symmetric\n" );
      return Failed;
    }
    //
    // Retrieve the initial vector and operator information from the Anasazi::Eigenproblem.
    //
    Teuchos::RefCountPtr<MV> iVec = _problem->GetInitVec();
    //
    if ( iVec.get() == 0 ) {
      _om->print( Anasazi::Error, "ERROR : Initial vector is not specified, set initial vector in eigenproblem\n");
      return Failed;
    }
    
    int dim = MVT::GetVecLength( *iVec );
    //
    // Check that the maximum number of blocks for the eigensolver is a positive number
    //    
    if ( _numBlocks<=0 ) {
      _om->stream(Anasazi::Error) << "ERROR : numBlocks = "<< _numBlocks <<" [ should be positive number ] " << endl;
      return Failed;
    } 
    //
    // Check that the maximum number of iterations is a positive number
    //    
    if ( _maxIter<=0 ) {
      _om->stream(Anasazi::Error) << "ERROR : maxIter = "<< _maxIter <<" [ should be positive number ] " << endl;
      return Failed;
    } 
    //
    // Check that the search subspace is larger than the number of eigenvalues requested
    //
    if ( _numBlocks*_blockSize < _nev ) {
      _om->stream( Anasazi::Error ) 
            << "ERROR : Search space dimension (numBlocks*blockSize) = "<< _numBlocks*_blockSize 
            << " [ should be greater than "<< _nev << " ] " << endl;
      return Failed;
    } 
    //
    // If the search subspace dimension is the same size as the number of requested eigenvalues,
    // then we must be computing all of them.
    //
    if ( (_numBlocks*_blockSize == _nev) && (_nev != dim) ) {
      _om->stream( Anasazi::Error )
            << "ERROR : Search space dimension (numBlocks*blockSize) = "<< _numBlocks*_blockSize 
            << " [ should be greater than "<< _nev << " ] " << endl;
      return Failed;
    }
    //
    // If the search subspace dimension is larger than the dimension of the operator, reset
    // the maximum number of blocks accordingly.
    //    
    if (_numBlocks*_blockSize > dim ) {
      _om->stream( Anasazi::Warning )
            << "WARNING : Search space dimension (numBlocks*blockSize) = "<< _numBlocks*_blockSize 
            <<" [ should not be greater than " << dim << " ] " << endl;
      
      // Set the maximum number of blocks in the factorization below the dimension of the space.
      _numBlocks = dim / _blockSize;
      
      _om->stream( Anasazi::Warning )
            << "WARNING : numBlocks reset to "<< _numBlocks << endl;
    }
    //
    // Reinitialize internal data and pointers, prepare for solve
    //
    _numRestarts = 0; 
    _iter = 0; 
    _knownEV = 0;
    //
    // Necessary variables
    //
    int i, j;
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
    //
    // Define local block vectors
    //
    // MX = Working vectors (storing M*X if M is specified, else pointing to X)
    // KX = Working vectors (storing K*X)
    //
    // R = Residuals
    //
    Teuchos::RefCountPtr<MV> X, MX, KX, R;
    X = MVT::Clone( *iVec, _dimSearch + _blockSize );
    KX = MVT::Clone( *iVec, _blockSize );
    R = MVT::Clone( *iVec, _blockSize );
    //
    // Check to see if there is a mass matrix, so we know how much space is required.
    // [ If there isn't a mass matrix we can use a view of X.]
    //
    if (_MOp.get()) {
      MX = MVT::Clone( *iVec, _blockSize );
    }
    else {
      std::vector<int> index( _blockSize );
      for (int i=0; i<_blockSize; i++) {
        index[i] = i;
      }
      MX = MVT::CloneView( *X, index );
    }
    //
    // KK = Local stiffness matrix               (size: dimSearch x dimSearch)
    //
    // S = Local eigenvectors                    (size: dimSearch x dimSearch)
    //
    Teuchos::SerialDenseMatrix<int,ScalarType> KK( _dimSearch, _dimSearch ), S( _dimSearch, _dimSearch );

    // Work vector for GEQRF and ORGQR
    std::vector<ScalarType> tau( _dimSearch );
    
    //
    // Initialize the workspace.
    //
    {
      // View vector into the basis
      Teuchos::RefCountPtr<MV> tmpXinit;
      // index vector
      std::vector<int> index;

      //
      // Determine how many init vectors were specified by the user
      //
      int numIVecs = MVT::GetNumberVecs( *iVec );
      if (numIVecs > _blockSize) {
        numIVecs = _blockSize;
      }

      //
      // Get a view into the basis
      //
      index.resize( numIVecs );
      for (i=0; i<numIVecs; i++) {
        index[i] = i;
      }
      tmpXinit = MVT::CloneView( *X, index );

      //
      // Copy the first numIVecs of the initial vectors into the first
      // numIVecs vectors of the basis (any additional vectors in iVec are ignored)
      //
      MVT::SetBlock( *iVec, index, *tmpXinit );

      //
      // Augment the initial vectors with random vectors if necessary
      //
      int leftOver = _blockSize - numIVecs;
      if (leftOver > 0) {
        index.resize(leftOver);
        for (i=0; i<leftOver; i++) {
          index[i] = numIVecs + i;
        }
        tmpXinit = MVT::CloneView( *X, index );
        MVT::MvRandom( *tmpXinit );
      }
    }

    //
    // Determine the maximum number of blocks for this factorization.
    // ( NOTE:  This will be the _numBlocks since we don't know about already converged vectors here )
    int maxBlock = _numBlocks;
    //
    // Compute the required workspace for this factorization.
    //
    lwork = lapack.ILAENV( 1, "geqrf", "", maxBlock*_blockSize, maxBlock*_blockSize );
    lwork *= _blockSize;
    std::vector<ScalarType> work( lwork );
    //
    while ( _iter <= _maxIter ) {
       
      for (nb = bStart; nb < maxBlock ; nb++) 
      {
        // 
        // Increment iteration counter
        //
        _iter++;
        if ( _iter > _maxIter )
          break;
        //
        // Get a view of the current vectors being worked on.
        // The offset is dependent on the number of iterations been taken so far (nb)
        // and the number of known eigenvalues (_knownEV).
        //
        int localSize = nb*_blockSize;
        std::vector<int> index( _blockSize );
        for (i=0; i < _blockSize; i++)
          index[i] = localSize + _knownEV + i;
        //
        Xcurrent = MVT::CloneView( *X, index );
        //
        if (_knownEV + localSize > 0) {
          index.resize( _knownEV + localSize );
          for (i=0; i < _knownEV + localSize; i++)
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

      //
      // Check if there is any reason we need to skip the rest of the while loop.
      //
      if (_iter > _maxIter) {
        break;
      }

      // Experienced unrecoverable error
      if ( criticalExit ) {
        break;
      }

      if (reStart == true) {
        reStart = false;
        continue;
      }

      //
      // Store the final converged eigenvectors
      //
      if (_knownEV + nFound >= _nev) {
        std::vector<int> index(1);
        for (j=0; j<_blockSize; j++) {
            if (_normR[j] < _residual_tolerance) {
              index[0] = j;
              Teuchos::RefCountPtr<MV> tmpKX = MVT::CloneView( *KX, index );
              index[0] = _knownEV;
              MVT::SetBlock( *tmpKX, index, *_evecs );
              (*_evals)[_knownEV] = _theta[j];
              _resids[_knownEV] = _normR[j];
              _knownEV++;
            }
            if (_knownEV == _nev) {
              break;
            }
        }
        break;
      } // if (_knownEV + nFound >= _nev)      
      //
      // Store the converged eigenvectors and define the restarting block
      // in the particular case of 1 block ( maxBlock == 1 ).
      //
      if (maxBlock == 1) {
        if (nFound > 0) {
          std::vector<int> index(1);
          int tmp_ptr = _knownEV + nFound;
          nFound = 0;
          for (j=0; j<_blockSize; j++) {
            //
            // Get a view of the current prospective eigenvector.
            //
            index[0] = j;
            Teuchos::RefCountPtr<MV> tmpKX = MVT::CloneView( *KX, index );
            if (_normR[j] < _residual_tolerance) {
              index[0] = _knownEV;
              MVT::SetBlock( *tmpKX, index, *X );              
              MVT::SetBlock( *tmpKX, index, *_evecs );              
              (*_evals)[_knownEV] = _theta[j];
              _resids[_knownEV] = _normR[j];
              _knownEV++;
              nFound++;              
            }
            else {
              index[0] = tmp_ptr + (j-nFound);
              MVT::SetBlock( *tmpKX, index, *X );
            }
          } // for (j=0; j<_blockSize; j++)
          //
          index.resize( nFound );
          for (i=0; i<nFound; i++) {
            index[i] = _knownEV + _blockSize - nFound + i;
          }
          Xnext = MVT::CloneView( *X, index );
          MVT::MvRandom( *Xnext );
        }
        else {
          nFound = _blockSize;
        }
        continue;
      } // if (maxBlock==1)
      //
      // Define the restarting block when maxBlock > 1
      //
      if (nFound > 0) {
        //
        // Copy the converged eigenvalues and residuals.
        //
        int count=0;
        std::vector<int> index(nFound);
        for (j=0; j<_blockSize; j++) {
          if (_normR[j] < _residual_tolerance) {
            index[count] = j;
            (*_evals)[_knownEV + count] = _theta[j];
            _resids[_knownEV + count] = _normR[j];
            count++;
          }
        }
        Teuchos::RefCountPtr<MV> tmpKX = MVT::CloneView( *KX, index );
        for (j=0; j<nFound; j++) {
          index[j] = _knownEV + j;
        }
        MVT::SetBlock( *tmpKX, index, *_evecs );

        int firstIndex = _blockSize;
        //
        // Find the first index where the unconverged eigenvectors start.
        //
        for (j=0; j<_blockSize; j++) {
          if (_normR[j] >= _residual_tolerance) {
            firstIndex = j;
            break;
          }
        }
        //
        // If the firstIndex comes before the number found, then we need
        // to move the converged eigenvectors to the front of the spectral
        // transformation.
        //
        std::vector<ScalarType> tmp_swap_vec(nb*_blockSize);
        while (firstIndex < nFound) {
          for (j=firstIndex; j<_blockSize; j++) {
            if (_normR[j] < _residual_tolerance) {
              //
              // Swap and j-th and firstIndex-th position
              //
              blas.COPY(nb*_blockSize, S[ j ], 1, &tmp_swap_vec[0], 1);
              blas.COPY(nb*_blockSize, S[ firstIndex ], 1, S[ j ], 1 );
              blas.COPY(nb*_blockSize, &tmp_swap_vec[0], 1, S[ firstIndex ], 1 );
              // Swap _theta
              std::swap<MagnitudeType>(_theta[j],_theta[firstIndex]);
              // Swap _normR
              std::swap<MagnitudeType>(_normR[j],_normR[firstIndex]);
              break;
            }
          } 
          for (j=0; j<_blockSize; j++) {
            if (_normR[j] >= _residual_tolerance) {
              firstIndex = j;
              break;
            }
          }
        } // while (firstIndex < nFound)
      } // if (nFound > 0)
      //
      // Define the restarting size
      //
      bStart = ((nb - offSet) > 2 ) ? (nb - offSet)/2 : 0;
      //
      // Define the restarting space and local stiffness matrix
      //
      KK.putScalar( zero );
      for (j=0; j<bStart*_blockSize; j++) {
        KK(j,j) = _theta[j + nFound];
      }
      //
      // Form the restarting space
      //
      int oldCol = nb*_blockSize;
      int newCol = nFound + (bStart+1)*_blockSize;
      newCol = (newCol > oldCol) ? oldCol : newCol;
      std::vector<int> index( oldCol );
      lapack.GEQRF(oldCol, newCol, S.values(), S.stride(), &tau[0], &work[0], lwork, &info);
      lapack.UNGQR(oldCol, newCol, newCol, S.values(), S.stride(), &tau[0], &work[0], lwork, &info);      
      for (i=0; i<oldCol; i++) {
        index[i] = _knownEV + i;
      }
      Teuchos::RefCountPtr<MV> oldX = MVT::CloneView( *X, index );
      index.resize( newCol );
      for (i=0; i<newCol; i++) {
        index[i] = _knownEV + i; 
      }
      Teuchos::RefCountPtr<MV> newX = MVT::CloneView( *X, index );
      Teuchos::RefCountPtr<MV> temp_newX = MVT::Clone( *X, newCol );
      Teuchos::SerialDenseMatrix<int,ScalarType> _Sview( Teuchos::View, S, oldCol, newCol );
      MVT::MvTimesMatAddMv( one, *oldX, _Sview, zero, *temp_newX );
      MVT::MvAddMv( one, *temp_newX, zero, *temp_newX, *newX ); 
      
      if (nFound == 0) {
        offSet++;
      }
      
      _knownEV += nFound;
      maxBlock = (_dimSearch/_blockSize) - (_knownEV/_blockSize);
      //
      // Put random vectors if the Rayleigh Ritz vectors are not enough
      // 
      newCol = nFound + (bStart+1)*_blockSize;
      if (newCol > oldCol) {
        index.resize( nFound );
        for (i=0; i<nFound; i++) {
          index[i] = _knownEV + _blockSize - nFound + i;
        }
        Xnext = MVT::CloneView( *X, index );
        MVT::MvRandom( *Xnext );
        continue;
      }
      //
      // Reset nFound counter.
      // 
      nFound = 0;
      //
    } // while (_iter < _maxIter)
    //
    // Sort the computed eigenvalues, eigenvectors, and residuals
    //
    if ((info==0) && (_knownEV > 0)) {
      // Sort the eigenvalues
      _order.resize(_knownEV);
      {
        Teuchos::TimeMonitor SortTimer(*_timerSortEval);
        ret = _sm->sort( this, _knownEV, &(*_evals)[0], &_order );
      }
      if (ret != Ok) {
        _om->print(Error,"ERROR : Sorting in solve()\n");
        criticalExit = true;
      }
      else {
        // use _order to permute _evecs and _resids
        _MSUtils.permuteVectors(_knownEV,_order,*_evecs,&_resids);
      }
    }

    //
    // Print out a final summary if necessary
    //
    if (_om->isVerbosity( FinalSummary )) {
      currentStatus();
    }

    // Print timing details 

    // Stop timer.
    _timerTotal->stop();

    // Reset format that will be used to print the summary
    Teuchos::TimeMonitor::format().setPageWidth(54);

    if (_om->isVerbosity( Anasazi::TimingDetails )) {
      _om->print(Anasazi::TimingDetails,"********************TIMING DETAILS********************\n");
      Teuchos::TimeMonitor::summarize( _om->stream(Anasazi::TimingDetails) );
      _om->print(Anasazi::TimingDetails,"******************************************************\n");
    }

    if (_knownEV == _nev) {
      return Ok;
    }
    else {
      return Unconverged;
    }

  } // end solve()

  template <class ScalarType, class MV, class OP>
  void 
  BlockDavidson<ScalarType,MV,OP>::accuracyCheck(Teuchos::RefCountPtr<const MV> X,  
                                                 Teuchos::RefCountPtr<const MV> MX, 
                                                 Teuchos::RefCountPtr<const MV> Q) const 
    {
      stringstream os;
      os.precision(2);
      os.setf(ios::scientific, ios::floatfield);
      ScalarType tmp;
      
      if (X != Teuchos::null) {
        if (MX != Teuchos::null) {
          tmp = _MSUtils.errorEquality(X.get(), MX.get(), _MOp.get());
          os << " >> Difference between MX and M*X = " << tmp << endl;
        }
        tmp = _orthman->orthonormError(*X,MX);
        os << " >> Error in X^T M X - I = " << tmp << endl;
      }
      
      if (Q == Teuchos::null) {
        os << endl;     
        return;
      }      

      tmp = _orthman->orthonormError(*Q);
      os << " >> Error in Q^T M Q - I = " << tmp << endl;
      if (X != Teuchos::null) {
        tmp = _orthman->orthogError(*X,MX,*Q);
        os << " >> Orthogonality Q^T M X up to " << tmp << endl;
      }
      os << endl;     
      
      _om->print(Error,os.str());
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
      os <<"The current basis size is " << _nevLocal<<end;
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
      if ( _iter > 0 || _nevLocal > 0 ) {
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

  
  } // End of namespace Anasazi

#endif

// End of file AnasaziBlockDavidson.hpp
