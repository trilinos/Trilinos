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
 *  \brief Implementation of the block Jacobi-Davidson method
 *
 *  \fixme: STILL TO BE COMPLETED!
 */

#ifndef ANASAZI_BLOCK_JACOBI_DAVIDSON_HPP
#define ANASAZI_BLOCK_JACOBI_DAVIDSON_HPP

#include "AnasaziConfigDefs.hpp"
#include "AnasaziEigensolver.hpp"
#include "AnasaziEigenproblem.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziEigensolver.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziProjectedEigenproblem.hpp"
#include "AnasaziOutputManager.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace Anasazi {

  template<class ScalarType, class MagnitudeType,  class MV, class OP>
  class BlockJacobiDavidson : public Eigensolver<ScalarType,MV,OP> 
  {
  public:
    //@{ \name Constructor/Destructor.

    //! %Anasazi::Block Jacobi-Davidson constructor.
    BlockJacobiDavidson (const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem,
                   const Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > &sm,
                   const Teuchos::RefCountPtr<OutputManager<ScalarType> > &om,
                   Teuchos::ParameterList &pl
                  );

    //! %Anasazi::BlockDavidson destructor.
    virtual ~BlockJacobiDavidson() {};
    //@}

    //@{ \name Solver application methods.
    /*! 
     * \brief This method uses iterate to compute approximate solutions to the
     * original problem.  It may return without converging if it has
     * taken the maximum number of iterations or numerical breakdown is
     * observed.
     */
    ReturnType solve();
    //@}

    //@{ \name Solver status methods.

    //! Get the current restart count of the iteration method.
    int GetNumRestarts() const 
    {
      return(_numRestarts); 
    }

    int GetNumIters() const 
    {
      return(_iters); 
    }

    /*! \brief Get a constant reference to the current linear problem, 
      which may include a current solution.
    */
    Eigenproblem<ScalarType,MV,OP>& GetEigenproblem() const { return(*_problem); };
    
    //! Get the blocksize to be used by the iterative solver in solving this eigenproblem.
    int GetBlockSize() const 
    { 
      return(_blockSize); 
    }

    //@}

    //@{ \name Output methods.

    //! This method requests that the solver print out its current status to screen.
    void currentStatus();

    //@}
  private:

    /*! \brief These methods will not be defined.     */
    BlockJacobiDavidson(const BlockJacobiDavidson<ScalarType,MagnitudeType,MV,OP> &method);
    BlockJacobiDavidson<ScalarType,MagnitudeType,MV,OP>& operator=
      (const BlockJacobiDavidson<ScalarType,MagnitudeType,MV,OP> &method);

    // 
    // Classes inputed through constructor that define the eigenproblem to be solved.
    //
    Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > _problem;
    Teuchos::RefCountPtr<OutputManager<ScalarType> > _om;
    Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > _sm;
    Teuchos::ParameterList _pl;
    //
    // Output stream from the output manager
    //
    std::ostream& _os;
    //
    // Information obtained from the eigenproblem
    //
    Teuchos::RefCountPtr<OP> _A;
    Teuchos::RefCountPtr<OP> _B;
    Teuchos::RefCountPtr<OP> _Prec;
    Teuchos::RefCountPtr<MV> _evecs;
    Teuchos::RefCountPtr<std::vector<ScalarType> > _evals;
    //
    // Internal data.
    // 
    ScalarType _residual_tolerance;
    ScalarType _TARGET;
    int _maxIters;
    int _blockSize;
    int _SMIN;
    int _SMAX;
    int _nev;
    int _numRestarts;
    int _iters;
    int _knownEV;
    std::vector<MagnitudeType> _normR;
    //
    // Convenience typedefs
    //
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
  };
  //
  // Implementation
  //
  template <class ScalarType, class MagnitudeType, class MV, class OP>
  BlockJacobiDavidson<ScalarType,MagnitudeType,MV,OP>::
  BlockJacobiDavidson(const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem,
                        const Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > &sm,
                        const Teuchos::RefCountPtr<OutputManager<ScalarType> > &om,
                        Teuchos::ParameterList &pl
                       ):
    _problem(problem), 
    _om(om),
    _sm(sm),
    _pl(pl),
    _os(_om->GetOStream()),
    _A(_problem->GetOperator()),
    _B(_problem->GetM()),
    _Prec(_problem->GetPrec()),
    _evecs(_problem->GetEvecs()),
    _evals(problem->GetEvals()), 
    _nev(problem->GetNEV()),
    _maxIters(_pl.get("Max Iters", 300)),
    _blockSize(_pl.get("Block Size", 1)),
    _residual_tolerance(_pl.get("Tol", 1.0e-6)),
    _numRestarts(0),
    _iters(0),
    _TARGET(_pl.get("Target", 0.0)),
    _SMIN(_pl.get("SMIN", 1)),
    _SMAX(_pl.get("SMAX", 16)),
    _knownEV(0)
  {
  }

  template <class ScalarType, class MagnitudeType, class MV, class OP>
  ReturnType BlockJacobiDavidson<ScalarType,MagnitudeType,MV,OP>::solve()
  {
    // ============= //
    // Sanity checks //
    // ============= //
    
    // Check the Anasazi::Eigenproblem was set by user, if not, return failed.
    if (!_problem->IsProblemSet()) 
    {
      if (_om->isVerbosityAndPrint(Anasazi::Error))
        _os << "ERROR : Anasazi::Eigenproblem was not set, call Anasazi::Eigenproblem::SetProblem() before calling solve"<< endl;
      return Failed;
    }
    
    // FIXME: Is this needed?
    // Check the Anasazi::Eigenproblem is symmetric, if not, return failed.
    if (!_problem->IsSymmetric()) 
    {
      if (_om->isVerbosityAndPrint( Anasazi::Error ))
        _os << "ERROR : Anasazi::Eigenproblem is not symmetric" << endl;
      return Failed;
    }
    
    // Retrieve the initial vector and operator information from the Anasazi::Eigenproblem.
    Teuchos::RefCountPtr<MV> iVec = _problem->GetInitVec();
    
    if (iVec.get() == 0) 
    {
      if (_om->isVerbosityAndPrint(Anasazi::Error)) 
        _os << "ERROR : Initial vector is not specified, set initial vector in eigenproblem "<<endl;
      return Failed;
    }
  
    if (iVec->GetNumberVecs() < 1) 
    {
      if (_om->isVerbosityAndPrint(Anasazi::Error)) 
        _os << "ERROR : Search space has incorrect dimension" << endl;
      return Failed;
    }
    
    // Check that the maximum number of blocks for the eigensolver is a positive number
    if (_blockSize <= 0) 
    {
      if (_om->isVerbosityAndPrint(Anasazi::Error)) 
        _os << "ERROR : blockSize = "<< _blockSize <<" [ should be positive number ] " << endl;
      return Failed;
    } 
    
    // Check that the maximum number of iterations is a positive number
    if (_maxIters <= 0) 
    {
      if (_om->isVerbosityAndPrint(Anasazi::Error)) 
        _os << "ERROR : maxIter = "<< _maxIters <<" [ should be positive number ] " << endl;
      return Failed;
    } 

    // Check that the maximum number of iterations is a positive number
    if (_evals->size() < _nev) 
    {
      if (_om->isVerbosityAndPrint(Anasazi::Error)) 
        _os << "ERROR : evals->size() is less than nev" << endl;
      return Failed;
    } 

    // ===================== //
    // Start with allocation //
    // ===================== //
    
    // Allocate the subspace U (and AU, BU), the Ritz-Galerkin matrix and the Ritz values ...
    int s;
    std::vector<int> sublist;

    // FIXME: is BU really necessary?
    Teuchos::RefCountPtr<MV> U, AU, BU;
    U  = MVT::Clone(*iVec, _SMAX);
    AU = MVT::Clone(*iVec, _SMAX);
    BU = MVT::Clone(*iVec, _SMAX);

    // Convenience multivectors, used to extract views of U, AU, BU
    MV *Utmp;
    MV *AUtmp;
    MV *BUtmp;
    MV *Rtmp;

    _normR.resize(_SMAX); // container for norms
    
    // working array
    std::vector<ScalarType> theta(_SMAX);
    std::vector<int> perm(_SMAX), iperm(_SMAX);
    Teuchos::SerialDenseMatrix<int, ScalarType> Z(_SMAX, _SMAX);
    Teuchos::SerialDenseMatrix<int, ScalarType> Ztmp;

    // Allocate space for the residuals ...
    Teuchos::RefCountPtr<MV> R;
    R = MVT::Clone(*iVec, _blockSize);
    R->MvRandom();

    // FIXME: are we using iVec ???
    //
    // Allocate space for the result ...
    int k;
    int conv;   // number of converged eigenpairs
    int purged; // number of deflated eigenpairs
    MV *Qtmp;
    MV *BQtmp;
    // FIXME: what is BQ??
    Teuchos::RefCountPtr<MV> BQ;
    BQ = MVT::Clone(*iVec, _nev);

    // Instantiate the ProjectedEigenproblem component and set some of the values ...
    ProjectedEigenproblem<int, ScalarType, MV> PE("Symmetric", _SMAX);
    PE.SetTheta(&theta);
    PE.SetZ(&Z);

    _knownEV = 0;
    s = 0;
    _iters = 0;

    if (_om->doPrint()) 
    {
      _os << "[Starting Solver]" << endl;
    }

    // ========================================= //
    // Start the (Block-) Davidson iteration ... //
    // ========================================= //
    
    while(_iters < _maxIters)
    {
      // ================================================== //
      // Add the (preconditioned) residuals to the search   //
      // space. To this end, orthogonalise w.r.t. foregoing //
      // subspace vectors and w.r.t. found eigenvectors Q   //
      // ================================================== //
      
      sublist.resize(1);
      std::vector<ScalarType> eta(1);
      std::vector<MagnitudeType> nrm(1);
      for(int i=0; i < _blockSize; i++)
      {
        sublist[0] = i;
        Rtmp = R->CloneView(sublist);

        // Gram-Schmidt reduce the vector ...
        for(int j=0; j < s; j++){
          sublist[0] = j;
          Utmp  = U->CloneView(sublist);
          BUtmp = BU->CloneView(sublist);
          Rtmp->MvDot((*BUtmp), &eta);
          Rtmp->MvAddMv(1.0, (*Rtmp), -eta[0], (*Utmp));
          delete Utmp, BUtmp;
        }

        for(int j=0; j < _knownEV; j++){
          sublist[0] = j;
          Qtmp = _evecs->CloneView(sublist);
          BQtmp = BQ->CloneView(sublist);
          Rtmp->MvDot((*BQtmp), &eta);
          Rtmp->MvAddMv(1.0, (*Rtmp), -eta[0], (*Qtmp));	      
          delete Qtmp, BQtmp;
        }

        // Now B-normalise the vector ...
        sublist[0] = _knownEV;
        Qtmp = _evecs->CloneView(sublist);
        _B->Apply((*Rtmp), (*Qtmp));
        Rtmp->MvDot((*Qtmp), &nrm);
        nrm[0] = sqrt(nrm[0]);
        Rtmp->MvAddMv(1.0/nrm[0], (*Rtmp), 0.0, (*Rtmp));

        sublist[0] = s;
        U->SetBlock((*Rtmp),sublist);

        delete Rtmp, Qtmp;
        s++;
      }

      // Update AU and BU ...
      sublist.resize(_blockSize);
      for(int i=0; i<_blockSize; ++i) sublist[i]=(s-_blockSize) + i;

      Utmp = U->CloneView(sublist);

      AUtmp = AU->CloneView(sublist);
      _A->Apply((*Utmp), (*AUtmp));

      BUtmp = BU->CloneView(sublist);
      _B->Apply((*Utmp), (*BUtmp));

      // Update the ProjectedEigenproblem component by telling it
      // about the space increase
      PE.Add((*Utmp), (*AUtmp), (*BUtmp));
      delete Utmp;
      delete AUtmp;
      delete BUtmp;

      // ====================================================== //
      // Extract eigenpair approximations from the space, then  //
      // Sort eigenpairs with respect to TARGET and compute the //
      // inverse permutation.                                   //
      // CAUTION: At return of the sorting manager, the thetas  //
      // are sorted!!                                           //
      // ====================================================== //
      
      PE.Extract();

      for(int i=0; i<s; ++i) {theta[i] -= _TARGET;}
      _sm->sort((Eigensolver<ScalarType,MV,OP>*)NULL, s, &theta[0], &perm);
      for(int i=0; i<s; ++i) {theta[i] += _TARGET;}
      for(int i=0; i<s; ++i) { iperm[perm[i]]=i;}

      // Permute the Z entries according to perm
      Ztmp.shape(s,1);

      for(int j=0; j<s; ++j){
        if (perm[j] != j){
          for(int i=0; i<s; ++i) {Ztmp[0][i] = Z[j][i];}
          for(int i=0; i<s; ++i) {Z[j][i] = Z[perm[j]][i];}
          for(int i=0; i<s; ++i) {Z[perm[j]][i] = Ztmp[0][i];}

          perm[iperm[j]] = perm[j];	  
        }  
      }

      // Check for convergence of the approximative eigenpairs. If
      // some of them are accurate enough, add them to the space _evecs
      Ztmp.shape(s,1);

      sublist.resize(s);
      for(int i=0; i<s; ++i) sublist[i]=i;	  

      Utmp = U->CloneView(sublist);
      AUtmp = AU->CloneView(sublist);
      BUtmp = BU->CloneView(sublist);

      sublist.resize(1);
      sublist[0] = 0;
      Rtmp = R->CloneView(sublist);

      for(int i=0; i<s; ++i) {Ztmp[0][i] = Z[0][i];}
      // Rtmp = AU*Z(:,1)
      Rtmp->MvTimesMatAddMv (1.0, (*AUtmp), Ztmp, 0.0); 
      // Rtmp = Rtmp - theta[0]*AU*Z(:,1)
      Rtmp->MvTimesMatAddMv (-theta[0], (*BUtmp), Ztmp, 1.0);
      // nrm  = ||Rtmp|| 
      Rtmp->MvNorm(&nrm);

#if 0
#define MIN(a,b) (a)<(b)?(a):(b)

      // FIXME: output manager
      cout << "It: " << _iters << " _knownEV: " << _knownEV;
      cout << " [" << nrm[0] << "]: ";
      cout << "(" << theta[0] << ")";
      int j = MIN(s,5);
      for(int i=1; i<j; ++i){
        cout << ", (" << theta[i] << ")";
      }
      cout << endl;
#endif

      conv = 0;
      while(nrm[0] < _residual_tolerance){

        sublist.resize(1);
        sublist[0] = _knownEV;
        Qtmp = _evecs->CloneView(sublist);
        // FIXME: the following was Q...
        BQtmp = BQ->CloneView(sublist);

        // Qtmp = U*Z(:,1)
        Qtmp->MvTimesMatAddMv (1.0, (*Utmp), Ztmp, 0.0);
        // BQtmp = B*Qtmp
        _B->Apply((*Qtmp), (*BQtmp));

        (*_evals)[_knownEV] = theta[conv];
        _normR[_knownEV] = nrm[0];
        _knownEV++;

        if (_knownEV == _nev) break;
        conv++;

        // Ztmp = Z(:,conv)
        for(int i=0; i<s; ++i) {Ztmp[0][i] = Z[conv][i];}
        // Rtmp = AU*Z(:,conv)
        Rtmp->MvTimesMatAddMv (1.0, (*AUtmp), Ztmp, 0.0);
        // Rtmp = Rtmp - theta[conv]*AU*Z(:,conv)
        Rtmp->MvTimesMatAddMv (-theta[conv], (*BUtmp), Ztmp, 1.0);
        // nrm  = ||Rtmp||
        Rtmp->MvNorm(&nrm);

        delete Qtmp;
      }

      delete Utmp, AUtmp, BUtmp;
      if (_knownEV == _nev) break;

      // ========================================================= //
      // conv of the approximations are saved since they are       //
      // accurate enough Perform a deflation (so far only on the Z //
      // matrix. The U, AU and BU matrices are adapted ONLY in the //
      // end)                                                      //
      // ========================================================= //
      
      // FIXME: Output manager everywhere!!!
      if (conv > 0)
      {
        if (_om->doPrint()) 
        {
          _os << "[Converged (" << conv << ")]" << endl;
        }
      }

      // ====================================================== //
      // Compute the residuals of the best blockSize eigenpair  //
      // approximations and precondition them, i.e. compute the //
      // new search space directions                            //
      // ====================================================== //
      
      Ztmp.shape(s,1);

      sublist.resize(s);
      for(int i=0; i<s; ++i) sublist[i]=i;	  

      AUtmp = AU->CloneView(sublist);
      BUtmp = BU->CloneView(sublist);

      sublist.resize(1);
      sublist[0]=_knownEV;	  
      Qtmp = _evecs->CloneView(sublist);

      for(int j=0; j<_blockSize; ++j){
        sublist[0] = j;
        Rtmp = R->CloneView(sublist);

        for(int i=0; i<s; ++i) {Ztmp[0][i] = Z[conv+j][i];}
        // Qtmp = AU*Z(:,1)
        Qtmp->MvTimesMatAddMv (1.0, (*AUtmp), Ztmp, 0.0);
        // Qtmp = Qtmp - theta[0]*AU*Z(:,1)
        Qtmp->MvTimesMatAddMv (-theta[conv+j], (*BUtmp), Ztmp, 1.0);
        // FIXME: XXXXXXXXXXXXXXXXXXX
        // if _Prec or _B are 0, the code crashes..
        _Prec->Apply((*Qtmp), (*Rtmp));

        delete Rtmp;
      }

      delete AUtmp, BUtmp, Qtmp;

      // ========================================================= //
      // Check the actual search space size and perform a restart, //
      // if adding the new directions would lead to an oversized   //
      // subspace                                                  //
      // ========================================================= //
      
      purged = 0;
      if ((s-conv+_blockSize)>_SMAX){
        // increment restart counter
        ++_numRestarts;

        // 
        purged = (s-conv) - (_SMIN-_blockSize);

        if (_om->doPrint()) 
        {
          _os << "[Restart (" << s - conv  << " -> "
              << s - conv - purged << ")]" << endl;
        }
      } 

      // Reshape the Z matrix according to the convergence
      // and/restart behaviour
      if ((conv + purged)>0)
      {
        sublist.resize(s);
        for(int i=0; i<s; ++i) sublist[i] = i;

        Utmp = U->CloneView(sublist);
        AUtmp = AU->CloneView(sublist);
        BUtmp = BU->CloneView(sublist);

        sublist.resize(s-conv-purged);
        for(int i=0; i<(s-conv-purged); ++i) sublist[i] = conv+i;
        Ztmp.shapeUninitialized(s,s-conv-purged);
        for(int i=0; i<s; ++i)
          for(int j=0; j<(s-conv-purged); ++j)
            Ztmp(i,j) = Z(i,j+conv);

        Utmp->MvTimesMatAddMv(1.0, (*Utmp), Ztmp, 0.0);
        AUtmp->MvTimesMatAddMv(1.0, (*AUtmp), Ztmp, 0.0);
        BUtmp->MvTimesMatAddMv(1.0, (*BUtmp), Ztmp, 0.0);

        PE.Rotate(Ztmp);	    

        for(int i=0; i<s; ++i){
          for(int j=0; j<(s-conv-purged); ++j){
            Z(i,j) = Z(i,j+conv);
          }
        }

        delete Utmp, AUtmp, BUtmp, Ztmp;

        // Finally, compute the new search space size
        s = s - conv - purged;
      }	  


      // Update the iteration step counter
      _iters++;
    }

    return(Ok);
  }

  // ==========================================================================
  template <class ScalarType, class MagnitudeType, class MV, class OP>
  void BlockJacobiDavidson<ScalarType,MagnitudeType,MV,OP>::
  currentStatus()
  {
    if (_om->doPrint()) 
    {
      _os.setf(ios::scientific, ios::floatfield);
      _os.precision(6);
      _os <<endl;
      _os << "******************* CURRENT STATUS *******************" << endl;
      _os << "Anasazi Block Jacobi-Davidson solver" << endl;
      _os << "The number of iterations performed is " <<_iters << endl;
      _os << "The number of restarts performed is " << _numRestarts << endl;
      _os << "The block size is " << _blockSize << endl;
      _os << "The number of eigenvalues requested is " << _nev << endl;
      _os << "The number of eigenvalues computed is "<< _knownEV << endl;
      _os << "The target is " << _TARGET << endl;
      _os << "The requested residual tolerance is "<< _residual_tolerance << endl;
      _os << endl;
      _os << "COMPUTED EIGENVALUES                 "<<endl;
      _os << std::setw(16) << std::right << "Eigenvalue" 
          << std::setw(16) << std::right << "Ritz Residual"
          << endl;
      _os << "------------------------------------------------------"<<endl;

      if (_knownEV > 0) 
      {
        for (int i = 0; i < _knownEV; i++) 
        {
          _os << std::setw(16) << std::right << (*_evals)[i]
              << std::setw(16) << std::right << _normR[i]
              << endl;
        }
      } 
      else 
      {
        _os << "[none computed]"<<endl;
      }
      _os << "******************************************************" << endl;  
      _os << endl; 
    }
  }
} // namespace Anasazi

#endif
