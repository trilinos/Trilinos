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

#include <sys/times.h>
#include <sys/timeb.h>

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziEigensolver.hpp"
#include "AnasaziEigenproblem.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperator.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziProjectedEigenproblem.hpp"
#include "AnasaziOutputManager.hpp"
#include "AnasaziUtils.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

//_SAB_Special: use only, if there exists a special fast MvTimesMatAddMv function, which can take a view of "this" as argument
//#define _SAB_Special

//_SAB_lateDefl: usefull, if _SAB_reorthog is also defined: late deflation is done and the SearchSpace will be reorthogonalized not only after restart, but also after conversion
//#define _SAB_lateDefl

//_SAB_modTargetSwitch: switch target separately for each search vector
//#define _SAB_modTargetSwitch

//_SAB_reorthog: to reorthogonalise search space after restart
//#define _SAB_reorthog

//_SAB_compare: some output and changes to compare with Oskar's Algo
#define _SAB_compare

//_SAB_DEBUG: for debugging purpose, write certain data to files
//#define _SAB_DEBUG
#define _SAB_DEBUG_ITER 30

//_SAB_TIMING: prints out needed time for the phases of the algorithm
//#define _SAB_TIMING

//_SAB_STATISTICS: prints out needed iterations and gmres steps and total time needed
//#define _SAB_STATISTICS

//_SAB_without_os: without output
//#define _SAB_without_os

/*!
\class Anasazi::BlockJacobiDavidson

\brief This class implements the block Jacobi-Davidson method, an iterative
method for solving eigenvalue problems.

The following discussion of the Jacobi-Davidson algorithm is derived from <I>
O. Chinellato, The Complex-Symmetric Jacobi-Davidson Algorithm and its
Application to the Computation of some Resonance Frequencies of Anisotropic
Lossy Axisymmetric Cavities, PhD Thesis, ETHZ, 2005</I>, Chapter 6.4.

The <I>Jacobi-Davidson Algorithm</I> was introduced by Sleijpen and van der
Vorst in 1996, and it is a typical subspace solver, designed to compute a
modest number of eigenpair approximations close to a given target \f$\tau\f$.
The eigenproblem to be solved is: Find <TT>nev</TT> eigenvectors of
\f[
A U = \lambda B U
\f]
that are close to \f$\tau\f$.

\fixme: something on A, B???

In an abstract fashion, the algorithm can be written as follows:
\f[
\begin{tabular}{p{1cm}{l}}
{\tt  1} &  Procedure JacobiDavidson$(A, B, \tau, U, \epsilon_{conv})$: \\
{\tt  2} &  Do: \\
{\tt  3} &  $\quad\quad$ $[\Theta, Z]$ := Extract$(A, B, U, \tau, \sigma)$ \\
{\tt  4} &  $\quad\quad$ $[\Lambda, Q, U, \Theta, Z, r]$ := Convergence/Deflation$(\Theta, Z, \epsilon_{conv})$ \\
{\tt  5} &  $\quad\quad$ $\tau$ := TargetSwitch$(r, \tau, \Theta)$ \\
{\tt  6} &  $\quad\quad$ Solve $(A - \sigma B) c \approx - r$ \\
{\tt  7} &  $\quad\quad$ $[U, \Theta, Z]$ := Restart($U, \Theta, Z, s_{min}, s_{max})$ \\
{\tt  8} &  $\quad\quad$ $U$ := Orthogonalize$(B, Q, U, C)$ \\
{\tt  9} &  EndProcedure 
\end{tabular}
\f]
This template is composed of several building blocks. We will now give some details on the implementation of each of them.

<P><B>Step 3: Subspace Extraction</B>. This step
finds eigenpairs residing in a given subspace \f$\mathcal{U}\f$,
which satify certain optimility conditions. Given the subspaces
\f$\mathcal{U}\f$, \f$A \mathcal{U}\f$ and \f$B \mathcal{U}\f$, several
strategies exist to extractthe (approximated) eigenpairs. Perhaps the best
known is the <I>Ritz-Galerkin approach</I>, which first defines the reduced 
operators:
\f[
\begin{tabular}{lcr}
$A_{\mathcal{U}}$ & = & $\mathcal{U}^H A \mathcal{U}$ \\
$B_{\mathcal{U}}$ & = & $\mathcal{U}^H B \mathcal{U}$, \\
\end{tabular}
\f]
then it computes the eigenpairs \f$(theta_i, z_i)\f$ of the reduced eigenproblem
\f[
A_{\mathcal{U}} z_i = \theta_i B_{\mathcal{U}} z_i.
\f]
Since this problem is dense, we can take advantage of LAPACK routines. This is done by the Anasazi::ProjectedEigenproblem class.

Given \f$z_i\f$, one can recover the eigenvalue corresponding to \f$\theta_i\f$
as
\f[
U_i = \mathcal{U} z_i.
\f]
The associated residuals can be expressed as \f$r(\theta_i, \mathcal{U} z_i) = A \mathcal{U} z_i - \theta_i B \mathcal{U} z_i\f$.


<P><B>Step 4: Convergence, Deflation</B>. Each eigenpair approximation obtained from the extract component is associated witha residual \f$r(\theta_i, \mathcal{U} z_i)\f$. If the corresponding <I>residual norm</I> is smaller than a required tolerence \f$\epsilon_{conv}\f$, i.e. if
\f[
\frac{
\| (\theta_i, \mathcal{U} z_i) \| 
}{
\| \mathcal{U} z_i  \|
} = 
\frac{
\| A \mathcal{U} z_i - \theta_i B \mathcal{U} z_i  \|
}{
\| \mathcal{U} z_i \|
} \leq \epsilon_{conv}
\f]
then the eigenpair is considered to be accurate enough and hence added to the
space of already converged eigenpairs \f$(\Lambda, Q)\f$. The Convergence/Deflation component iterates over the extracted eigenpair approximations one at a time and checks for thir convergence. Each accepted eigenpari is appended to the space of the already copnvernged ones and the subspace \f$\mathcal{U}\f$ is adapted appropriately. If the overlall number of converged eigenpairs reached the desired number of eigenpairs <TT>nev</TT>, the process terminates. 
 When an eigenpair is accepted, the search spaces are purged of all the components pointing in the directions of this pair. This guarantees that the same eigenvectors will not be found again, and it is called Deflation.
In class
Anasazi::BlockJacobiDavidson, \f$\Lambda\f$ and \f$Q\f$ are defined within the
input Anasazi::Eigenproblem class. 


<P><B>Step 5: Target Switch</B>. The Jacobi-Davidson algorithm is a typical shift-and-invert algorithm; however, in contrast to Krylov space methods, the shift is allowed to vary from iteration to iteration. Choosing the shift close to an actual eigenvalue of the matrix pencil \f$(A, B)\f$ typically yeilds good subspace extensions in the direction of the eigenvector associated with this eigenvalue, which in turn will speed up the convergence towards this eigenpair. A possible strategy to define the target switch is as follows. As soon as the relative norm falls below the requested threshold, the shift is chosen to equal the approximation eigenvalue \f$\theta_i\f$.


<P><B>Step 6: Correction Equation</B>. The Jacobi-Davidson algorithm requires,
  at each step, the expansion of the search space \f$\mathcal{U}\f$. Clearly,
  the quality of the approximations depends on the procedure adopted to expand
  this search space. A very elegant way to obtain good additional candidates
  consists in solving the so called <I>correction equation</I>.  This equation
  can be derived by considering the application of an eigenpair correction
  \f$(\delta, c)\f$ to the eigenpair approximation 
  \f$(\theta_1, \mathcal{U} z_1)\f$ in such a way
  that the new residual \f$r(\theta_1 + \delta, \mathcal{U} z_1 + c)\f$ will vanish. By expanding the "corrected" residual we have
\f[
A(\mathcal{U} z_1 + c) - (\theta_1 + \delta) B (\mathcal{U} z_1 + c) = 0.
\f]
We can replace \f$(\theta_1 + \delta)\f$ with either \f$theta_1\f$ 
(if the algorithm is close to converge to the desired eigenvalue, and
 therefore \f$\delta\f$ will be small) or \f$\tau\f$ 
(if the current approximation is too far from an actual eigenvalue). Therefore,
  the above equation can be rewritten as
\f[
(A - \chi B)c = - r(\theta_1, \mathcal{U} z_1).
\f]
with an appropriate \f$\chi\f$. Obviously, \f$c = \mathcal{U} z_1\f$ satisfies
this equation, yet without yeilding an additional direction. Therefore,
  the space where the additional correction vector is searched must be restricted. According to Sleijpen and van der Vorst (1996), the space being B-orthogonal to the matrix \f$\hat{Q} = (Q, q)\f$, where \f$q\f$ is the actual eigenvector approximation, represents a reasonable choice. Hence, the correction equation to be solved reads
\f[
(A - \chi B) c =  r(\theta_1, \mathcal{U} z_1) \mbox{ s.t. }
\hat{Q}^T B c = 0.
\f]
Note that this constraint offers the further advantage of inhibiting a potentially rapid condition number growth induced by the \f$\sigma's\f$ convergence towards a true eigenvalue. 


Note that the correction equation is the major advantage of the
Jacovi-Davidson algorthm over other iterative eigenvalue solvers. In contract
to other algorithms where the subspace is augmented following strict rules,
  the correction equation, which reflects the actual approximation status, is
  solved in each iteration. In this way the algorithm is guarantee to adapt to
  the actual situation which in turn promises a good performance.
            

<P><B>Step 7: Restart</B> 
In order to keep the storage requirements low, the dimension of the
subspace is reduced as soon as its dimension reaches an upper limit,
\f$s_{max}\f$. When doing so, the keep the eigenvectors approximations
associated with the \f$s_{min} -1\f$ best eigenpair approximations contained
in the subspace, where \f$s_{min}\f$ denotes the lower dimension limit of
the search space \f$\mathcal{U}\f$. This way, the algorithm maintains some
of the approximations gained during the previous iterations. In contrast to Krylov-based solvers, for which a restart constritutes an involved task, this step can be carried out easily and efficiently, and it is implemented as follows:
\f[
\begin{tabular}{p{1cm}{l}}
{\tt  1} &  Procedure Restart$(\mathcal{U}, \Theta, Z, s_{min}, s_{max})$ \\
{\tt  2} &  Do: \\
{\tt  3} &  $\quad\quad$If ($s = s_{max}$) Then \\
{\tt  4} &  $\quad\quad\quad\quad$\mathcal{U} := \mathcal{U} Z(1:s_{max}, 1:s_{min -1)$ \\
{\tt  5} &  $\quad\quad$End \\
{\tt  6} &  $\quad\quad$Return (\mathcal{U}, \Theta, Z)$ \\
{\tt  7} &  EndProcedure  \\
\end{tabular}
\f]

\author Oscar Chinellato (ETHZ/ICOS), Sabine Emch (ETHZ), Marzio Sala (ETHZ/COLAB)

\date Last updated on 01-Nov-05
*/
namespace Anasazi {

  template<class ScalarType, class MV, class OP>
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

    //! Get the number of iterations.
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
     typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
     
    /*! \brief These methods will not be defined.     */
    BlockJacobiDavidson(const BlockJacobiDavidson<ScalarType,MV,OP> &method);
    BlockJacobiDavidson<ScalarType,MV,OP>& operator=
      (const BlockJacobiDavidson<ScalarType,MV,OP> &method);

    inline Teuchos::RefCountPtr<MV> CV(MV& rhs, int which)
    {
      std::vector<int> list(1);
      list[0] = which;

      Teuchos::RefCountPtr<MV> res = Teuchos::rcp(rhs.CloneView(list));
      return(res);
    }
    
    inline Teuchos::RefCountPtr<MV> CC(MV& rhs, int which)
    {
      std::vector<int> list(1);
      list[0] = which;

      Teuchos::RefCountPtr<MV> res = Teuchos::rcp(rhs.CloneCopy(list));
      return(res);
    }

    inline Teuchos::RefCountPtr<MV> CV(MV& rhs, int first, int last)
    {
      std::vector<int> list(last - first);
      for (int i = first ; i < last ; ++i)
        list[i - first] = i;

      Teuchos::RefCountPtr<MV> res = Teuchos::rcp(rhs.CloneView(list));
      return(res);
    }

    inline Teuchos::RefCountPtr<MV> CC(MV& rhs, int first, int last)
    {
      std::vector<int> list(last - first);
      for (int i = first ; i < last ; ++i)
        list[i - first] = i;

      Teuchos::RefCountPtr<MV> res = Teuchos::rcp(rhs.CloneCopy(list));
      return(res);
    }
    
    inline Teuchos::RefCountPtr<MV> CV(MV& rhs, std::vector<int>& list)
    {
      Teuchos::RefCountPtr<MV> res = Teuchos::rcp(rhs.CloneView(list));
      return(res);
    }

     inline Teuchos::RefCountPtr<MV> CC(MV& rhs, std::vector<int>& list)
    {
      Teuchos::RefCountPtr<MV> res = Teuchos::rcp(rhs.CloneCopy(list));
      return(res);
    }
    
    inline Teuchos::RefCountPtr<MV> CV(Teuchos::RefCountPtr<MV>& rhs, int which)
    {
      std::vector<int> list(1);
      list[0] = which;

      Teuchos::RefCountPtr<MV> res = Teuchos::rcp(rhs->CloneView(list));
      return(res);
    }

    inline Teuchos::RefCountPtr<MV> CC(Teuchos::RefCountPtr<MV>& rhs, int which)
    {
      std::vector<int> list(1);
      list[0] = which;

      Teuchos::RefCountPtr<MV> res = Teuchos::rcp(rhs->CloneCopy(list));
      return(res);
    }
    
    inline Teuchos::RefCountPtr<MV> CV(Teuchos::RefCountPtr<MV>& rhs, 
                                   int first, int last)
    {
      std::vector<int> list(last - first);
      for (int i = first ; i < last ; ++i)
        list[i - first] = i;

      Teuchos::RefCountPtr<MV> res = Teuchos::rcp(rhs->CloneView(list));
      return(res);
    }

    inline Teuchos::RefCountPtr<MV> CC(Teuchos::RefCountPtr<MV>& rhs, 
                                   int first, int last)
    {
      std::vector<int> list(last - first);
      for (int i = first ; i < last ; ++i)
        list[i - first] = i;

      Teuchos::RefCountPtr<MV> res = Teuchos::rcp(rhs->CloneCopy(list));
      return(res);
    }
    
    inline Teuchos::RefCountPtr<MV> CV(Teuchos::RefCountPtr<MV>& rhs, 
                                       std::vector<int>& list)
    {
      Teuchos::RefCountPtr<MV> res = Teuchos::rcp(rhs->CloneView(list));
      return(res);
    }

    // 
    // Classes inputed through constructor that define the eigenproblem to be solved.
    //
    //! Problem to be solved.
    Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > _problem;
    //! Output manager to be used for any output operation.
    Teuchos::RefCountPtr<OutputManager<ScalarType> > _om;
    //! Sorting manager to be used to sort the eigenvalues.
    Teuchos::RefCountPtr<SortManager<ScalarType,MV,OP> > _sm;
    //! Parameter list containing options for this solver.
    Teuchos::ParameterList _pl;
    //! Output stream, reference from the output manager
    std::ostream& _os;
    
    // @}
    // @{ \name Information obtained from the eigenproblem
    
    //! Operator A
    Teuchos::RefCountPtr<OP> _A;
    //! Operator B
    Teuchos::RefCountPtr<OP> _B;
    //! Preconditioner
    Teuchos::RefCountPtr<OP> _Prec;
    //! Multi-vector that will contain the computed eigenvectors.
    Teuchos::RefCountPtr<MV> _evecs;
    //! Will contain the computed eigenvalues.
    Teuchos::RefCountPtr<std::vector<ScalarType> > _evals;
    
    // @}
    // @{ \name Internal data
    
    //! Tolerance for convergence
    MagnitudeType _residual_tolerance;
    //! Target, the solver will compute the eigenvalues close to this value.
    ScalarType _TARGET;
    //! Maximum number of allowed iterations.
    int _maxIters;
    //! Block size.
    int _blockSize;
    int _SMIN;
    int _SMAX;
    //! Number of requested eigenvalues.
    int _nev;
    //! Number of performed restarts.
    int _numRestarts;
    //! Number of performed iterations.
    int _iters;
    //! Number of computed and accepted eigenvalues.
    int _knownEV;
    //! Residual norm of computed and accepted eigenvalues.
    std::vector<MagnitudeType> _normR;
    bool _haveKinv;
    int _LSIterMax, _LSRestart, _gamma;
    double _elin, _eswitch;
    
    //
    // Convenience typedefs
    //
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;
  };
  
  //
  // Implementation
  //
  template <class ScalarType, class MV, class OP>
  BlockJacobiDavidson<ScalarType,MV,OP>::
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
    _maxIters(_pl.get("Max Iters", 100)),
    _blockSize(_pl.get("Block Size", 1)),
    _residual_tolerance(_pl.get("Tol", 1e-6)),
    _numRestarts(0),
    _iters(0),
   // _TARGET(_pl.get<ScalarType>("Target", Teuchos::ScalarTraits<ScalarType>::zero())),
    _TARGET((ScalarType)1),
    _SMIN(_pl.get("SMIN", 20)), 
    _SMAX(_pl.get("SMAX", 30)),
    _knownEV(0),
    _LSIterMax(_pl.get("KrylovSolver: MaxIters", 200)),
    _LSRestart(_pl.get("KrylovSolver: Restart",20)),
    _elin(_pl.get("KrylovSolver: Tolerance", 10e-9)),
    _gamma(_pl.get("gamma", 2)),
    _eswitch(_pl.get("Switch", 0.1))
  {
  }

  template <class ScalarType, class MV, class OP>
  ReturnType BlockJacobiDavidson<ScalarType,MV,OP>::solve()
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
    
    // Symmetric or Hermitian, but not general
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
 
    if (iVec->GetNumberVecs() != _blockSize) //Alternative: ds=iVec->GetNumberVecs setzen (ds>0!)
    {
      if (_om->isVerbosityAndPrint(Anasazi::Error)) 
        _os << "ERROR : Search space has incorrect dimension" << endl;
	_os << "iVec.numVecs: " << iVec->GetNumberVecs() << "blocksize: " << _blockSize << endl;
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

    // Check that the eigenvalue vector is allocated with correct size
    if (_evals->size() < _nev) 
    {
      if (_om->isVerbosityAndPrint(Anasazi::Error)) 
        _os << "ERROR : evals->size() is less than nev" << endl;
      return Failed;
    } 

    if (iVec->GetVecLength()-_nev<_SMAX) {
        if (_om->isVerbosityAndPrint(Anasazi::Error))
		_os << "ERROR : itialization vector length bigger than SMAX-nev" << endl;
	return Failed;
    }
    
    if (_SMIN>=_SMAX) {
        if (_om->isVerbosityAndPrint(Anasazi::Error))
		_os << "ERROR : SMIN >= SMAX" << endl;
	return Failed;
    }
    
    if (_SMIN<=_blockSize) {
    	if (_om->isVerbosityAndPrint(Anasazi::Error))
		_os << "ERROR : SMIN <= blocksize" << endl;
	return Failed;
    }
    
    // =============================== //
    // templated values of 1's and 0's //
    // =============================== //

    ScalarType ScalarOne  = Teuchos::ScalarTraits<ScalarType>::one();
    ScalarType ScalarZero = Teuchos::ScalarTraits<ScalarType>::zero();

    // ===================== //
    // Start with allocation //
    // ===================== //
    
    // Allocate the subspace U (and AU, BU), the Ritz-Galerkin matrix and the Ritz values ...
    int SearchSpaceSize = 0, outer=0, ds=_blockSize, restartsOld=0;
    std::vector<int> sublist;

    Teuchos::RefCountPtr<MV> U, AU, BU;
    U  = MVT::Clone(*iVec, _SMAX);
    AU = MVT::Clone(*iVec, _SMAX);
    BU = MVT::Clone(*iVec, _SMAX);

    #ifdef _SAB_TIMING
    struct timeb tb0, tbOrth, tbExt, tbConv, tbCorr, tbGMRES, tbRest;
    struct tms tm0, tmOrth, tmExt, tmConv, tmCorr, tmGMRES, tmRest;
    int tOrth1=0, tOrth2=0, tExt1=0, tExt2=0, tConv1=0, tConv2=0, tCorr1=0, tCorr2=0, tGMRES1=0, tGMRES2=0, tRest1=0, tRest2=0;
    #endif
    
    #ifdef _SAB_STATISTICS
    int numtotGMRES=0;
    struct timeb tbTot0, tbTot;
    struct tms tmTot0, tmTot;
    
    ftime(&tbTot0);
    times(&tmTot0);
    #endif
    
    _normR.resize(_nev); // container for norms
    
    // working array
    std::vector<ScalarType> theta(_SMAX), vc(1);
    std::vector<int> perm(_SMAX), iperm(_SMAX), cperm(_SMAX), vQ2;
    Teuchos::SerialDenseMatrix<int, ScalarType> Z(_SMAX, _SMAX), LkTrans(_nev,_nev), Rk(_nev,_nev);
    Teuchos::SerialDenseMatrix<int, ScalarType> Ztmp;

    // Allocate space for the residuals ...
    std::vector<Teuchos::RefCountPtr<MV> > R(_blockSize);
    //std::vector<MV*> R(_blockSize);
    sublist.resize(1);
    assert(iVec->GetNumberVecs()==_blockSize);
    for (int i=0; i<_blockSize; i++) {
    	sublist[0]=i;
    	R[i]=Teuchos::rcp(iVec->CloneCopy(sublist));
	//R[i]=iVec->CloneCopy(sublist);
    }
    Teuchos::RefCountPtr<MV> mvtmp;
    mvtmp = MVT::Clone(*iVec,1);

    // Allocate space for the result ...
    int conv = 0;   // number of converged eigenpairs per iteration
    int purged = 0; // number of deflated eigenpairs per iteration
    std::vector<ScalarType> sigma(_blockSize);
    MV *Qtmp;
    MV *BQtmp;
    MV *KBQtmp, *h;
    
    #ifdef _SAB_compare
    FILE *file;
    file=fopen("Res2.txt","w");
    #endif
    
    #ifdef _SAB_DEBUG
    FILE *file2, *fileC, *fileAU, *fileBU, *fileTheta, *fileZ, *fileKsys, *fileAsys, *fileR, *fileMisc;
    FILE *fileLk, *fileRk;
    
    fprintf(file,"Trilinos Algo\n");
    
    file2=fopen("EntriesU.txt","w");
    fileAU=fopen("EntriesAU.txt","w");
    fileBU=fopen("EntriesBU.txt","w");
    fileTheta=fopen("EntriesTheta.txt","w");
    fileZ=fopen("EntriesZ.txt","w");
    fileKsys=fopen("EntriesKsys.txt","w");
    fileAsys=fopen("EntriesAsys.txt","w");
    fileR=fopen("EntriesR.txt","w");
    fileMisc=fopen("EntriesMisc.txt","w");
    fileC=fopen("EntriesC.txt","w");
    fileLk=fopen("EntriesLk.txt","w");
    fileRk=fopen("EntriesRk.txt","w");
    #endif
    
     
    Teuchos::RefCountPtr<MV> Qtmp2, BQtmp2, KBQtmp2;

    Teuchos::RefCountPtr<MV> BQ, KBQ;
    Teuchos::BLAS<int,ScalarType> blas;
    Anasazi::ReturnType ok;
    
    BQ = MVT::Clone(*iVec, _nev);
    KBQ = MVT::Clone(*iVec, _nev);
    
    // Instantiate the ProjectedEigenproblem component and set some of the values ...
    // FIXME: for Hermitian problems...
    ProjectedEigenproblem<int, ScalarType, MV> PE("Symmetric", _SMAX);
    PE.SetTheta(&theta);
    PE.SetZ(&Z);

    _knownEV = 0;
    SearchSpaceSize = 0;
    _iters = 0;

    #ifndef _SAB_without_os
    if (_om->doPrint()) 
    {
      _os << "[Starting Solver]" << endl;
    }
    #endif

    // ================================================ //
    // Start the (Block-) Jacobi/Davidson iteration ... //
    // ================================================ //
    
    while(_iters < _maxIters)
    {
      // ================================================== //
      // Add the (preconditioned) residuals to the search   //
      // space. To this end, orthogonalise w.r.t. foregoing //
      // subspace vectors and w.r.t. found eigenvectors Q   //
      // ================================================== //
      
      #ifdef _SAB_TIMING
      ftime(&tb0);
      times(&tm0);
      #endif
      
      std::vector<ScalarType> eta(1);
      std::vector<MagnitudeType> nrm(1);
      std::vector<MagnitudeType> nrm_q(1);
      
      sublist.resize(1);
      
      for(int i=0; i < ds; i++)
      {    
	// B-normalise the vector
	if (_B.get() != 0) {
		ok=_B->Apply(*(R[i]), *mvtmp);
		assert(ok==Anasazi::Ok);
		R[i]->MvDot(*mvtmp, &eta, Anasazi::NO_CONJ);
	} else R[i]->MvDot(*R[i], &eta, Anasazi::NO_CONJ);
	
	eta[0]=sqrt(eta[0]); //FIXME: replace with corrected Teuchos-function        
        R[i]->MvAddMv(ScalarOne / eta[0], *(R[i]), ScalarZero, *mvtmp);
	
	
	MagnitudeType l=(MagnitudeType)0.001;
	int onemore = 0;
  	while(1){

    		if (onemore == 1) {break;}
    		if ((l < 1.1) && (l > 0.9)) {onemore = 1;} //Achtung: ist SMAX>DIM-nev, dann gibt's hier eine Endlosschleife
		#ifdef _SAB_DEBUG
		else printf("need more\n");
		#endif
				
		// Gram-Schmidt reduce the vector ...
		for(int j=0; j < SearchSpaceSize ; j++)
		{
			R[i]->MvDot(*CV(BU, j), &eta, Anasazi::NO_CONJ);
			R[i]->MvAddMv(ScalarOne, *(R[i]), -eta[0], *CV(U, j));
		}
		
		for(int j=0; j < _knownEV; j++)
		{
			R[i]->MvDot(*CV(BQ, j), &eta, Anasazi::NO_CONJ);
			R[i]->MvAddMv(ScalarOne, *(R[i]), -eta[0], *CV(_evecs, j));
		}
		
		// Now B-normalise the vector ...
		if (_B.get() != 0) {
			ok=_B->Apply(*(R[i]), *mvtmp);
			assert(ok==Anasazi::Ok);
			R[i]->MvDot(*mvtmp, &eta, Anasazi::NO_CONJ);
		} else R[i]->MvDot(*(R[i]), &eta, Anasazi::NO_CONJ);
		
		eta[0]=sqrt(eta[0]); //FIXME: replace with corrected Teuchos-function
		l=Teuchos::ScalarTraits<ScalarType>::magnitude(eta[0]);
		R[i]->MvAddMv(ScalarOne / eta[0], *(R[i]), ScalarZero, *mvtmp);
	}
	
	 #ifdef _SAB_DEBUG
      	//Teste orhogonalitaet
	/*
	if (_B.get() != 0) _B->Apply(*CV(R, i), *CV(_evecs, _knownEV));
	else {
		sublist[0] = i;
		CV(_evecs, _knownEV)->SetBlock(*CV(R, i), sublist);
	}
		
	CV(R, i)->MvDot(*CV(_evecs, _knownEV), &eta, Anasazi::NO_CONJ);
	cout << "eta: " << eta[0] << endl;
	*/
      	#endif
	
	//update U, BU
	sublist[0]=SearchSpaceSize;
        U->SetBlock(*(R[i]),sublist);
        if (_B.get()) {
		ok=_B->Apply(*(R[i]), *CV(BU, SearchSpaceSize));
		assert(ok==Anasazi::Ok);
	} else CV(BU, SearchSpaceSize)->SetBlock(*(R[i]), sublist);
	  
        SearchSpaceSize++;
      }
      
      // Update AU
      sublist.resize(_blockSize);
      for(int i=0; i<_blockSize; ++i) {
      	sublist[i]=(SearchSpaceSize - _blockSize) + i;
      }
      
      ok=_A->Apply(*CV(U, sublist), *CV(AU, sublist));
      assert(ok==Anasazi::Ok);
      
      #ifdef _SAB_DEBUG
      if ((_iters<_SAB_DEBUG_ITER) || (restartsOld!=_numRestarts)) {
	fprintf(file2,"%d: U:\n",_iters);
	CV(U,0,SearchSpaceSize)->MvPrintf_abs(file2);
	
	fprintf(fileAU,"%d: AU:\n",_iters);
	CV(AU,0,SearchSpaceSize)->MvPrintf_abs(fileAU);
	
	fprintf(fileBU,"%d: BU:\n",_iters);
	CV(BU,0,SearchSpaceSize)->MvPrintf_abs(fileBU);
      }
      #endif
      
      // Update the ProjectedEigenproblem component by telling it
      // about the space increase
      PE.Add(*CV(U, sublist), *CV(AU, sublist), *CV(BU, sublist));

      #ifdef _SAB_TIMING
      ftime(&tbOrth);
      times(&tmOrth);
      
      tOrth1+=(tbOrth.time-tb0.time)*1000+tbOrth.millitm-tb0.millitm;
      tOrth2+=tmOrth.tms_utime-tm0.tms_utime+tmOrth.tms_stime-tm0.tms_stime;
      #endif
      
      // ====================================================== //
      // Extract eigenpair approximations from the space, then  //
      // Sort eigenpairs with respect to TARGET and compute the //
      // inverse permutation.                                   //
      // CAUTION: At return of the sorting manager, the thetas  //
      // are sorted!!                                           //
      // ====================================================== //
      
      #ifdef _SAB_TIMING
      ftime(&tb0);
      times(&tm0);
      #endif
      
      PE.Extract();
      
      #ifdef _SAB_DEBUG
      if ((_iters<_SAB_DEBUG_ITER) || (restartsOld!=_numRestarts)) {
      	fprintf(fileTheta,"%d: theta:\n",_iters);
      	for (int i=0; i<SearchSpaceSize; i++) {
		fprintf(fileTheta,"%1.16f %1.16f\n",theta[i].real(), theta[i].imag());
	}
	
	fprintf(fileZ,"%d: Z:\n",_iters);
	for (int i=0; i<SearchSpaceSize; i++) {
		fprintf(fileZ,"[ ");
		for (int j=0; j<SearchSpaceSize; j++) {
			fprintf(fileZ,"(%1.16f, %1.16f) ",Z(i,j).real(), Z(i,j).imag());
		}
		fprintf(fileZ,"\n");
	}
      }
      #endif
      
      for(int i=0; i<SearchSpaceSize ; ++i) {theta[i] -= _TARGET;}
      _sm->sort((Eigensolver<ScalarType,MV,OP>*)NULL, SearchSpaceSize , &theta[0], &perm);
      
      for(int i=0; i<SearchSpaceSize ; ++i) {theta[i] += _TARGET;}
      for(int i=0; i<SearchSpaceSize ; ++i) { iperm[i]=i;}
      for(int i=0; i<SearchSpaceSize; i++) { cperm[i]=i;}

      // Permute the Z entries according to perm
      Ztmp.shape(SearchSpaceSize ,1);
      
      for(int j=0; j<SearchSpaceSize ; ++j){
        if (iperm[perm[j]] != j){
          for(int i=0; i<SearchSpaceSize ; ++i) {Ztmp[0][i] = Z[j][i];}
          for(int i=0; i<SearchSpaceSize ; ++i) {Z[j][i] = Z[iperm[perm[j]]][i];}
          for(int i=0; i<SearchSpaceSize ; ++i) {Z[iperm[perm[j]]][i] = Ztmp[0][i];}

          int h=iperm[perm[j]];
	  iperm[perm[j]]=j;
	  iperm[cperm[j]]=h;
	  
	  cperm[h]=cperm[j];
	  cperm[j]=perm[j];
        }
      }
      
      
      #ifdef _SAB_DEBUG
      if ((_iters<_SAB_DEBUG_ITER) || (restartsOld!=_numRestarts)) {
      	fprintf(fileTheta,"%dss: theta:\n",_iters);
      	for (int i=0; i<SearchSpaceSize; i++) {
		fprintf(fileTheta,"%1.16f %1.16f\n",theta[i].real(), theta[i].imag());
	}
	
	fprintf(fileZ,"%dss: Z:\n",_iters);
	for (int i=0; i<SearchSpaceSize; i++) {
		fprintf(fileZ,"[ ");
		for (int j=0; j<SearchSpaceSize; j++) {
			fprintf(fileZ,"(%1.16f, %1.16f) ",Z(i,j).real(), Z(i,j).imag());
		}
		fprintf(fileZ,"\n");
	}
      }
      #endif
      
      //Skalieren
      for(int j=0; j<SearchSpaceSize; j++) {
       Teuchos::SerialDenseMatrix<int,ScalarType> H(1,1), Zview(Teuchos::View, Z[j], Z.stride(), SearchSpaceSize, 1);
        H.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, ScalarOne, Zview, Zview, ScalarZero);
       H(0,0) = 1.0/sqrt(H(0,0)); //FIXME: replace later with corrected Teuchos-function
     
       blas.SCAL(SearchSpaceSize, H(0,0), Z[j], 1);
     }
     
     #ifdef _SAB_DEBUG
      if ((_iters<_SAB_DEBUG_ITER) || (restartsOld!=_numRestarts)) {
      	fprintf(fileTheta,"%dss: theta:\n",_iters);
      	for (int i=0; i<SearchSpaceSize; i++) {
		fprintf(fileTheta,"%1.16f %1.16f\n",theta[i].real(), theta[i].imag());
	}
	
	fprintf(fileZ,"%dss: Z:\n",_iters);
	for (int i=0; i<SearchSpaceSize; i++) {
		fprintf(fileZ,"[ ");
		for (int j=0; j<SearchSpaceSize; j++) {
			fprintf(fileZ,"(%1.16f, %1.16f) ",Z(i,j).real(), Z(i,j).imag());
		}
		fprintf(fileZ,"\n");
	}
      }
      #endif
     
     #ifdef _SAB_TIMING
      ftime(&tbExt);
      times(&tmExt);
      
      tExt1+=(tbExt.time-tb0.time)*1000+tbExt.millitm-tb0.millitm;
      tExt2+=tmExt.tms_utime-tm0.tms_utime+tmExt.tms_stime-tm0.tms_stime;
      #endif
      
      #ifdef _SAB_TIMING
      ftime(&tb0);
      times(&tm0);
      #endif
      
      // Check for convergence of the approximative eigenpairs. If
      // some of them are accurate enough, add them to the space _evecs
      
      // Ztmp = Z(:,0)
      Ztmp.shape(SearchSpaceSize, 1);
      for(int i=0; i<SearchSpaceSize ; ++i) {Ztmp[0][i] = Z[0][i];}

      // Qtmp = U*Z(:,0), BQtmp = BU*Z(:,0), nrm_q = ||Qtmp||
      CV(_evecs, _knownEV)->MvTimesMatAddMv (ScalarOne, *CV(U, 0, SearchSpaceSize), Ztmp, ScalarZero);
      CV(_evecs, _knownEV)->MvNorm(&nrm_q);
      CV(BQ, _knownEV)->MvTimesMatAddMv(ScalarOne, *CV(BU, 0, SearchSpaceSize), Ztmp, ScalarZero);
      
      // Rtmp = AU*Z(:,0)-theta[0]*BU*Z(:,0), nrm = ||Rtmp||
      R[0]->MvTimesMatAddMv (ScalarOne, *CV(AU, 0, SearchSpaceSize), Ztmp, ScalarZero);
      R[0]->MvAddMv (ScalarOne, *(R[0]), -theta[0], *CV(BQ, _knownEV));
      R[0]->MvNorm(&nrm);
       
      //relative norm     
      nrm[0] /= nrm_q[0];

      #ifndef _SAB_without_os
      if (_om->doPrint()) 
      {
        _os << "It: " << _iters << ", knownEV: " << _knownEV;
       // _os << ", s: " << s << ", rel_normr = " << nrm[0];
        _os << ", SearchSpaceSize : " << SearchSpaceSize  << ", rel_normr = " << nrm[0];
        _os << ", theta_i = ";
        _os << "(" << theta[0] << ")";
        int j = ANASAZI_MIN(SearchSpaceSize ,5);
        for(int i=1; i<j; ++i){
          _os << ", (" << theta[i] << ")";
        }
        _os << "(Tol=" << _residual_tolerance << ")" << endl;
      }
      #endif
      
      #ifdef _SAB_compare
      fprintf(file,"it=%d nrm/nrm_q=%1.16f\n",_iters, nrm[0]);
      #endif
      
      conv = 0;

      #ifndef _SAB_compare      
      while ((conv<SearchSpaceSize-1)&&(nrm[0] <= _residual_tolerance)) 
      #else
      if (nrm[0] <= _residual_tolerance) //nur fuer Vergleich mit Oskars Code
      #endif
      {
 	#ifdef _SAB_compare
        fprintf(file,"converged EP\n");
	#endif
	
        //new eigenvector already saved in _evecs
	//save new eigenvalue in _evals
        (*_evals)[_knownEV] = theta[conv];
	
        _normR[_knownEV] = nrm[0]; // relative norm
        _knownEV++;

        if (_knownEV == _nev) break;

        conv++;
	
	outer=0;
      
	//BQ already updated (BQ := [BQ BQtmp])
	
	if (_Prec.get()) 
        {
		//KBQ := [KBQ _Prec*BQtmp]
		ok=_Prec->Apply(*CV(BQ, _knownEV-1), *CV(KBQ, _knownEV-1));
		assert(ok==Anasazi::Ok);
	
		//update Lk an Rk	
		if (_knownEV==1) {
			LkTrans(0,0)=ScalarOne;
			CV(BQ, _knownEV-1) ->MvDot(*CV(KBQ, _knownEV-1), &vc, Anasazi::NO_CONJ);
			Rk(0,0)=vc[0];
		} else {
			Teuchos::SerialDenseMatrix<int,ScalarType> lk(Teuchos::View, LkTrans[_knownEV-1], LkTrans.stride(), _knownEV-1, 1), rk(Teuchos::View, Rk[_knownEV-1], Rk.stride(), _knownEV-1, 1), H(1,1);
			
			CV(BQ, _knownEV-1)->MvTransMv(ScalarOne, *CV(KBQ, 0, _knownEV - 1), lk, Anasazi::NO_CONJ);			
			blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, Teuchos::NON_UNIT_DIAG, lk.numRows(), lk.numCols(), ScalarOne, Rk.values(), Rk.stride(), lk.values(), lk.stride());
			
			CV(KBQ,_knownEV-1)->MvTransMv(ScalarOne, *CV(BQ,0,_knownEV - 1), rk, Anasazi::NO_CONJ);
			blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, Teuchos::UNIT_DIAG, rk.numRows(), rk.numCols(), ScalarOne, LkTrans.values(), LkTrans.stride(), rk.values(), rk.stride());
			
                        CV(BQ, _knownEV-1)->MvDot(*CV(KBQ, _knownEV-1), &vc, Anasazi::NO_CONJ);
			ScalarType rho=vc[0];
			H.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, ScalarOne, lk, rk, ScalarZero); //rho+=lk*rk
			rho-=H(0,0);
			
			LkTrans(_knownEV-1,_knownEV-1)=ScalarOne; //_Lk:=[Lk 0; lk^T 1]			
			Rk(_knownEV-1,_knownEV-1)=rho; //Rk:=[Rk rk; 0 rho]
		}
		
		#ifdef _SAB_DEBUG
		fprintf(fileLk,"%dconv: Lk:\n",_iters);
		for (int i=0; i<_knownEV; i++) {
			fprintf(fileLk,"[ ");
			for (int j=0; j<_knownEV; j++) {
				fprintf(fileLk,"(%1.16f, %1.16f) ",LkTrans(j,i).real(), LkTrans(j,i).imag());
			}
			fprintf(fileLk,"]\n");
		}
		
		fprintf(fileRk,"%dconv: Rk:\n",_iters);
		for (int i=0; i<_knownEV; i++) {
			fprintf(fileRk,"[ ");
			for (int j=0; j<_knownEV; j++) {
				fprintf(fileRk,"(%1.16f, %1.16f) ",Rk(i,j).real(), Rk(i,j).imag());
			}
			fprintf(fileRk,"]\n");
		}
		#endif
	}
	
	//Ztmp = Z(:,conv)
        for(int i=0; i<SearchSpaceSize ; ++i) {Ztmp[0][i] = Z[conv][i];}

	// Qtmp = U*Z(:,conv), BQtmp = BU*Z(:,conv), nrm_q = ||Qtmp||
        CV(_evecs, _knownEV)->MvTimesMatAddMv (ScalarOne, *CV(U, 0, SearchSpaceSize), Ztmp, ScalarZero);
        CV(_evecs, _knownEV)->MvNorm(&nrm_q);
        CV(BQ, _knownEV)->MvTimesMatAddMv(ScalarOne, *CV(BU, 0, SearchSpaceSize), Ztmp, ScalarZero);

        // Rtmp = AU*Z(:,conv) - theta[conv]*BU*Z(:,conv), nrm = ||Rtmp||
        R[0]->MvTimesMatAddMv (ScalarOne, *CV(AU, 0, SearchSpaceSize), Ztmp, ScalarZero);
        R[0]->MvAddMv (ScalarOne, *R[0], -theta[conv], *CV(BQ, _knownEV));
        R[0]->MvNorm(&nrm);
        
	//relative norm
        nrm[0] /= nrm_q[0];
	#ifdef _SAB_compare
	fprintf(file,"it=%d nrm/nrm_q=%1.16f\n",_iters, nrm[0]);
	#endif
      }
      
      if (_knownEV == _nev) {
      	break;
      }

      // ========================================================= //
      // conv of the approximations are saved since they are       //
      // accurate enough Perform a deflation (so far only on the Z //
      // matrix. The U, AU and BU matrices are adapted ONLY in the //
      // end)                                                      //
      // ========================================================= //
      if (conv > 0)
      {
      	#ifndef _SAB_without_os
        if (_om->doPrint()) 
        {
          _os << "[Converged (" << conv << ")]" << endl;
	}
	#endif
	
	ds=ANASAZI_MIN(SearchSpaceSize-conv,_blockSize);
	
	#ifndef _SAB_lateDefl
	assert(SearchSpaceSize-conv>0);
        Teuchos::SerialDenseMatrix<int,ScalarType> ZZ(Teuchos::View, Z[conv], Z.stride(), SearchSpaceSize, SearchSpaceSize - conv), ZZ2(Teuchos::View, Z[0], Z.stride(), SearchSpaceSize-conv, SearchSpaceSize-conv);
	
	#ifdef _SAB_Special
	
	//update U, AU, BU	
	CV(U,0,SearchSpaceSize-conv)->MvTimesMatAddMv(ScalarOne, *CV(U,0,SearchSpaceSize), ZZ, ScalarZero);
	CV(AU,0,SearchSpaceSize-conv)->MvTimesMatAddMv(ScalarOne, *CV(AU,0,SearchSpaceSize), ZZ, ScalarZero);
	CV(BU,0,SearchSpaceSize-conv)->MvTimesMatAddMv(ScalarOne, *CV(BU,0,SearchSpaceSize), ZZ, ScalarZero);
	
	#else
	
	MV *mvtmp2=U->Clone(SearchSpaceSize-conv);
	
	sublist.resize(SearchSpaceSize-conv);	
	for (int i=0; i<SearchSpaceSize-conv; i++) sublist[i]=i;
	
	//update U, AU, BU
	mvtmp2->MvTimesMatAddMv(ScalarOne, *CV(U,0,SearchSpaceSize), ZZ, ScalarZero);
	U->SetBlock(*mvtmp2, sublist);
	mvtmp2->MvTimesMatAddMv(ScalarOne, *CV(AU,0,SearchSpaceSize), ZZ, ScalarZero);
	AU->SetBlock(*mvtmp2, sublist);
	mvtmp2->MvTimesMatAddMv(ScalarOne, *CV(BU,0,SearchSpaceSize), ZZ, ScalarZero);
	BU->SetBlock(*mvtmp2, sublist);
	
	delete mvtmp2;
	
	#endif

	//update G (projected eigenvalue problem)
        //PE.Rotate(ZZ);
	PE.UpdateG(SearchSpaceSize-conv,conv);

	//update Z and theta
	ZZ2.putScalar(ScalarZero);
	for(int i=0; i<(SearchSpaceSize-conv); ++i){
	  theta[i] = theta[i+conv];
	  Z(i,i) = ScalarOne;
        }

        // Finally, compute the new search space size
        SearchSpaceSize = SearchSpaceSize - conv;
	
	#ifdef _SAB_DEBUG
	if ((_iters<_SAB_DEBUG_ITER) || (restartsOld!=_numRestarts)) {
		fprintf(fileTheta,"%dconv: theta:\n",_iters);
		for (int i=0; i<SearchSpaceSize; i++) {
			fprintf(fileTheta,"%1.16f %1.16f\n",theta[i].real(), theta[i].imag());
		}
		
		fprintf(fileZ,"%dconv: Z:\n",_iters);
		for (int i=0; i<SearchSpaceSize; i++) {
			fprintf(fileZ,"[ ");
			for (int j=0; j<SearchSpaceSize; j++) {
				fprintf(fileZ,"(%1.16f, %1.16f) ",Z(i,j).real(), Z(i,j).imag());
			}
			fprintf(fileZ,"\n");
		}
	}
	#endif
	
	#endif
      }	  
      
      #ifdef _SAB_TIMING
      ftime(&tbConv);
      times(&tmConv);
      
      tConv1+=(tbConv.time-tb0.time)*1000+tbConv.millitm-tb0.millitm;
      tConv2+=tmConv.tms_utime-tm0.tms_utime+tmConv.tms_stime-tm0.tms_stime;
      #endif
      
      #ifndef _SAB_modTargetSwitch
      // Switch target
      for(int i=0; i<ds; ++i)
      {
 	if (nrm[0]<_eswitch) sigma[i] = theta[0];
 	else sigma[i] = _TARGET;
      }
      #endif

      // ====================================================== //
      // Compute the residuals of the best blockSize eigenpair  //
      // approximations and precondition them, i.e. compute the //
      // new search space directions                            //
      // ====================================================== //
      
      #ifdef _SAB_TIMING
      ftime(&tb0);
      times(&tm0);
      #endif

      //views needed for correction
      Qtmp2 = CV(_evecs,0,_knownEV+1);
      BQtmp2 = CV(BQ,0,_knownEV+1);
      KBQtmp2 = CV(KBQ,0,_knownEV+1);
      
      sublist.resize(1);
      for (int j=0; j<ds ; j++){
      	#ifdef _SAB_lateDefl
	Teuchos::SerialDenseMatrix<int,ScalarType> ZZ(Teuchos::View, Z[j+conv], Z.stride(), SearchSpaceSize, 1);
	#else
	Teuchos::SerialDenseMatrix<int,ScalarType> ZZ(Teuchos::View, Z[j], Z.stride(), SearchSpaceSize, 1);
	#endif
	
	sublist[0]=_knownEV;
	
	// Qtmp = U*ZZ, BQtmp = BU*ZZ
	CV(_evecs, _knownEV)->MvTimesMatAddMv (ScalarOne, *CV(U, 0, SearchSpaceSize), ZZ, ScalarZero);
	if (_B.get()) {
		ok=_B->Apply(*CV(_evecs, _knownEV), *CV(BQ, _knownEV));
		assert(ok==Anasazi::Ok);
	} else CV(BQ, _knownEV)->SetBlock(*CV(_evecs, _knownEV), sublist);
	
	//Rtmp = AU*ZZ - theta[j(+conv)]*BU*ZZ
 	R[j]->MvTimesMatAddMv (ScalarOne, *CV(AU, 0, SearchSpaceSize), ZZ, ScalarZero);
	#ifdef _SAB_lateDefl
	R[j]->MvAddMv (ScalarOne, *(R[j]), -theta[j+conv], *CV(BQ, _knownEV));
	#else
	R[j]->MvAddMv (ScalarOne, *(R[j]), -theta[j], *CV(BQ, _knownEV));
	#endif
	
	#ifdef _SAB_modTargetSwitch
	// nrm = ||Rtmp||
        R[j]->MvNorm(&nrm);

	// nrm_q = ||Qtmp||
        CV(_evecs,_knownEV)->MvNorm(&nrm_q);
	
	if (nrm[0]<nrm_q[0]*_eswitch) {
		#ifndef _SAB_lateDefl
		sigma[j]=theta[j];
		#else
		sigma[j]=theta[j+conv];
		#endif
	} else sigma[j]=_TARGET;
	#endif
	
	Teuchos::RefCountPtr<Anasazi::OpBase<ScalarType, MV> > Ksys, Asys;
	
	//operator Ksys
	if (!_Prec.get()) {
		//Ksys:=I-Q*Qb^(T)
		Ksys = Teuchos::rcp(new OpKsysSimple<ScalarType,MV>(Qtmp2, BQtmp2));
	} else {
		sublist[0]=_knownEV;
		
		//KBQtmp = _Prec*BQtmp
		ok=_Prec->Apply(*CV(BQ,_knownEV), *CV(KBQ,_knownEV));
		assert(ok==Anasazi::Ok);
		
		if (_knownEV==0) {
			LkTrans(0,0)=ScalarOne; //Lk(0,0):=1
			CV(BQ,_knownEV)->MvDot(*CV(KBQ,_knownEV), &vc, Anasazi::NO_CONJ); //vc[0]:=BQtmp*KBQtmp
			Rk(0,0)=vc[0];//Rk(0,0)=vc[0]
		} else {
			Teuchos::SerialDenseMatrix<int,ScalarType> lk(Teuchos::View, LkTrans[_knownEV], LkTrans.stride(), _knownEV, 1);
			
			CV(BQ,_knownEV)->MvTransMv(ScalarOne, *CV(KBQ,0,_knownEV), lk, Anasazi::NO_CONJ);
			blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, Teuchos::NON_UNIT_DIAG, lk.numRows(), lk.numCols(), ScalarOne, Rk.values(), Rk.stride(), lk.values(), lk.stride());
			
			Teuchos::SerialDenseMatrix<int,ScalarType> rk(Teuchos::View, Rk[_knownEV], Rk.stride(), _knownEV, 1);
			
			CV(KBQ,_knownEV)->MvTransMv(ScalarOne, *CV(BQ,0,_knownEV), rk, Anasazi::NO_CONJ);
			blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::TRANS, Teuchos::UNIT_DIAG, rk.numRows(), rk.numCols(), ScalarOne, LkTrans.values(), LkTrans.stride(), rk.values(), rk.stride());
			
			CV(BQ,_knownEV)->MvDot(*CV(KBQ,_knownEV), &vc, Anasazi::NO_CONJ);
			ScalarType rho=vc[0];
			Teuchos::SerialDenseMatrix<int,ScalarType> H(1,1);
			H.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, ScalarOne, lk, rk, ScalarZero); //rho+=lk*rk
			rho-=H(0,0);

			LkTrans(_knownEV,_knownEV)=ScalarOne; //_Lk:=[Lk 0; lk^T 1]
			Rk(_knownEV,_knownEV)=rho; //Rk:=[Rk rk; 0 rho]
		}
		Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > LkTransView, RkView;
		LkTransView = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int,ScalarType> (Teuchos::View, LkTrans.values(), LkTrans.stride(), _knownEV+1, _knownEV+1));
		
		RkView = Teuchos::rcp(new Teuchos::SerialDenseMatrix<int,ScalarType> (Teuchos::View, Rk.values(), Rk.stride(), _knownEV+1, _knownEV+1));
	
		Ksys=Teuchos::rcp(new OpKsysSoph<ScalarType,MV,OP>(KBQtmp2, LkTransView, RkView, BQtmp2, _Prec));
	}
	
	//operator Asys:=(I-Qb*Q^(T))*(_A-_sigma*_B)
	Asys=Teuchos::rcp(new OpAsys<ScalarType, MV, OP>(BQtmp2, Qtmp2, _A, sigma[j], _B)); 
	
	Teuchos::SerialDenseMatrix<int,ScalarType> H; 
	H.shapeUninitialized(_knownEV+1,1);

	R[j]->MvTransMv(ScalarOne, *Qtmp2, H, Anasazi::NO_CONJ);
	Teuchos::RefCountPtr<MV> r=Teuchos::rcp(R[j]->CloneCopy());
	r->MvTimesMatAddMv(ScalarOne, *BQtmp2, H, -ScalarOne); //r:=-(I-Qb*Q^(T))*_r

	GMRES<ScalarType,MV,Anasazi::OpBase<ScalarType,MV> > gmres;
	
	gmres.setMaxIter(_LSIterMax);
	gmres.setTolerance(ANASAZI_MAX(pow(_gamma,-(double)outer),_elin));
	gmres.setRestart(_LSRestart);
	
	R[j]->MvInit(ScalarZero);
	
	#ifdef _SAB_DEBUG
	fprintf(fileLk,"%d: Lk:\n",_iters);
	for (int i=0; i<=_knownEV; i++) {
		fprintf(fileLk,"[ ");
		for (int j=0; j<=_knownEV; j++) {
			fprintf(fileLk,"(%1.16f, %1.16f) ",LkTrans(j,i).real(), LkTrans(j,i).imag());
		}
		fprintf(fileLk,"]\n");
	}
	
	fprintf(fileRk,"%d: Rk:\n",_iters);
	for (int i=0; i<=_knownEV; i++) {
		fprintf(fileRk,"[ ");
		for (int j=0; j<=_knownEV; j++) {
			fprintf(fileRk,"(%1.16f, %1.16f) ",Rk(i,j).real(), Rk(i,j).imag());
		}
		fprintf(fileRk,"]\n");
	}
		
	if ((_iters<_SAB_DEBUG_ITER) || (restartsOld!=_numRestarts)) {
		MV *bla, *bla2;
		
		bla=R[j]->CloneCopy();
		bla2=R[j]->CloneCopy();
		
		bla->MvInit(std::complex<double>(0.0,0.0));
		for (int ibla=0; ibla<=_knownEV; ibla++) {
			(*bla)[0][ibla]=std::complex<double>(1.0,0.0);
			Ksys->Apply(*bla,*bla2);
			fprintf(fileKsys, "%d: Ksys[%d]:\n",_iters,ibla);
			bla2->MvPrintf(fileKsys);
			(*bla)[0][ibla]=std::complex<double>(0.0,0.0);
		}
		
		for (int ibla=0; ibla<=_knownEV; ibla++) {
			(*bla)[0][ibla]=std::complex<double>(1.0,0.0);
			Asys->Apply(*bla,*bla2);
			fprintf(fileAsys, "%d: Asys:\n",_iters);
			bla2->MvPrintf(fileAsys);
			(*bla)[0][ibla]=std::complex<double>(0.0,0.0);
		}
		
		fprintf(fileR, "%d: r\n",_iters);
		r->MvPrintf(fileR);
		
		fprintf(fileMisc, "%d: Q:\n",_iters);
		Qtmp2->MvPrintf(fileMisc);
		
		fprintf(fileMisc, "%d: Qb:\n",_iters);
		BQtmp2->MvPrintf(fileMisc);
		
		if (_Prec.get()) {
			fprintf(fileMisc, "%d: Qk:\n",_iters);
			KBQtmp2->MvPrintf(fileMisc);
		}
		
		delete bla;
		delete bla2;
	}
	#endif
	
	#ifdef _SAB_TIMING
        ftime(&tbOrth);
        times(&tmOrth);
        #endif
	
	
	ok=gmres.solve(*Asys, *Ksys, *r, *(R[j])); 

	#ifdef _SAB_TIMING
	ftime(&tbGMRES);
	times(&tmGMRES);
	
	tGMRES1+=(tbGMRES.time-tbOrth.time)*1000+tbGMRES.millitm-tbOrth.millitm;
	tGMRES2+=tmGMRES.tms_utime-tmOrth.tms_utime+tmGMRES.tms_stime-tmOrth.tms_stime;
	#endif

	//printf("GMRES: %d %d\n",gmres.getIterations(), gmres.getRestarts());
	
	#ifdef _SAB_STATISTICS
	numtotGMRES+=gmres.getIterations();
	#endif
	
      #ifdef _SAB_DEBUG      
      if ((_iters<_SAB_DEBUG_ITER) || (restartsOld!=_numRestarts)) {
        //printf("nach GMRES\n");
      	fprintf(fileC,"%d: c:\n",_iters);
      	R[j]->MvPrintf(fileC);
      }
      #endif
     
	//printf("GMRES (it=%d, restarts=%d)\n",gmres.getIterations(), gmres.getRestarts());
      }
      
      #ifdef _SAB_TIMING
      ftime(&tbCorr);
      times(&tmCorr);
      
      tCorr1+=(tbCorr.time-tb0.time)*1000+tbCorr.millitm-tb0.millitm;
      tCorr2+=tmCorr.tms_utime-tm0.tms_utime+tmCorr.tms_stime-tm0.tms_stime;
      #endif

      // ========================================================= //
      // Check the actual search space size and perform a restart, //
      // if adding the new directions would lead to an oversized   //
      // subspace                                                  //
      // ========================================================= //
      
      #ifdef _SAB_TIMING
      ftime(&tb0);
      times(&tm0);
      #endif
    
      restartsOld=_numRestarts;
      
      purged = 0;
      
      #ifdef _SAB_lateDefl
      if ((SearchSpaceSize -conv + _blockSize)>_SMAX){
      //if ((SearchSpaceSize - conv + _blockSize)>ANASAZI_MIN(iVec->GetVecLength()-_knownEV,_SMAX) ) {
      	purged = SearchSpaceSize - conv - (_SMIN - _blockSize);
      #else
       if ((SearchSpaceSize + _blockSize)>_SMAX){
       //if ((SearchSpaceSize + _blockSize)> ANASAZI_MIN(iVec->GetVecLength()-_knownEV,_SMAX) ) {
      	purged = SearchSpaceSize - (_SMIN - _blockSize);
      #endif

      #ifndef _SAB_without_os
        if (_om->doPrint()) 
        {
          _os << "[Restart (" << SearchSpaceSize  << " -> "
              << SearchSpaceSize - purged << ")]" << endl;
        }
	#endif
	
        // increment restart counter
        _numRestarts++;
	
	assert(SearchSpaceSize-purged>0);
	
	#ifdef _SAB_lateDefl
	}
	if ((conv + purged)>0) {
	#endif
	
	#ifdef _SAB_lateDefl
	Teuchos::SerialDenseMatrix<int, ScalarType> ZZ(Teuchos::View, Z[conv],
                                                       Z.stride(), SearchSpaceSize,
                                                       SearchSpaceSize - conv - purged);
	
	MV *mvtmp2=U->Clone(SearchSpaceSize-conv-purged);
	
	sublist.resize(SearchSpaceSize-conv-purged);	
	for (int i=0; i<SearchSpaceSize-conv-purged; i++) sublist[i]=i;
	
	mvtmp2->MvTimesMatAddMv(ScalarOne, *CV(U,0,SearchSpaceSize), ZZ, ScalarZero);
	U->SetBlock(*mvtmp2, sublist);
	
	#ifndef _SAB_reorthog
	// wird sonst bei Reorthogonalisierung berechnet
	mvtmp2->MvTimesMatAddMv(ScalarOne, *CV(AU,0,SearchSpaceSize), ZZ, ScalarZero);
	AU->SetBlock(*mvtmp2, sublist);
	mvtmp2->MvTimesMatAddMv(ScalarOne, *CV(BU,0,SearchSpaceSize), ZZ, ScalarZero);
	BU->SetBlock(*mvtmp2, sublist);
	#endif
	
	delete mvtmp2;
	
	#else
	
	// Reshape the Z matrix according to the restart behaviour	
        Teuchos::SerialDenseMatrix<int, ScalarType> ZZ(Teuchos::View, Z[0],
                                                       Z.stride(), SearchSpaceSize,
                                                       SearchSpaceSize - purged);
		
	#ifdef _SAB_Special
	
	CV(U,0,SearchSpaceSize-purged)->MvTimesMatAddMv(ScalarOne, *CV(U,0,SearchSpaceSize), ZZ, ScalarZero);
	
	#ifndef _SAB_reorthog
	// wird sonst bei Reorthogonalisierung berechnet
	CV(AU,0,SearchSpaceSize-purged)->MvTimesMatAddMv(ScalarOne, *CV(AU,0,SearchSpaceSize), ZZ, ScalarZero);
	CV(BU,0,SearchSpaceSize-purged)->MvTimesMatAddMv(ScalarOne, *CV(BU,0,SearchSpaceSize), ZZ, ScalarZero);
	#endif
	
	#else				       
	
	MV *mvtmp2=U->Clone(SearchSpaceSize-purged);
	
	sublist.resize(SearchSpaceSize-purged);	
	for (int i=0; i<SearchSpaceSize-purged; i++) sublist[i]=i;
	
	mvtmp2->MvTimesMatAddMv(ScalarOne, *CV(U,0,SearchSpaceSize), ZZ, ScalarZero);
	U->SetBlock(*mvtmp2, sublist);
	
	#ifndef _SAB_reorthog
	// wird sonst bei Reorthogonalisierung berechnet
	mvtmp2->MvTimesMatAddMv(ScalarOne, *CV(AU,0,SearchSpaceSize), ZZ, ScalarZero);
	AU->SetBlock(*mvtmp2, sublist);
	mvtmp2->MvTimesMatAddMv(ScalarOne, *CV(BU,0,SearchSpaceSize), ZZ, ScalarZero);
	BU->SetBlock(*mvtmp2, sublist);
	#endif
	
	delete mvtmp2;
	
	#endif
	#endif
	
	// Finally, compute the new search space size
	#ifdef _SAB_lateDefl
	SearchSpaceSize = SearchSpaceSize - conv - purged;
	#else
	SearchSpaceSize = SearchSpaceSize - purged;
	#endif
	
	#ifndef _SAB_reorthog
	
	#ifndef _SAB_lateDefl
	PE.UpdateG(SearchSpaceSize,0);
	#else
	PE.UpdateG(SearchSpaceSize,conv);
	#endif
        //PE.Rotate(ZZ);
	#else
	
	//Reorthogonalisation
	sublist.resize(1);
		
	for(int i=0; i < SearchSpaceSize; i++) {
		sublist[0]=i;
		MV *utmp=U->CloneView(sublist);
		
		if (_B.get() != 0) {
			ok=_B->Apply(*utmp, *mvtmp);
			assert(ok==Anasazi::Ok);
			utmp->MvDot(*mvtmp, &eta, Anasazi::NO_CONJ);
		} else utmp->MvDot(*utmp, &eta, Anasazi::NO_CONJ);
		
		eta[0]=sqrt(eta[0]); //FIXME: replace with corrected Teuchos-function        
		utmp->MvAddMv(ScalarOne / eta[0], *utmp, ScalarZero, *mvtmp);
		
		MagnitudeType l=(MagnitudeType)0.001;
		int onemore = 0;
		while(1){
	
			if (onemore == 1) {break;}
			if ((l < 1.1) && (l > 0.9)) {onemore = 1;} //Achtung: ist SMAX>DIM-nev, dann gibt's hier eine Endlosschleife
			#ifdef _SAB_DEBUG
			else printf("need more (l=%f)\n",l);
			#endif
				
			// Gram-Schmidt reduce the vector ...
			for(int j=0; j < i ; j++)
			{
				utmp->MvDot(*CV(BU, j), &eta, Anasazi::NO_CONJ);
				utmp->MvAddMv(ScalarOne, *utmp, -eta[0], *CV(U,j));
			}
			
			for(int j=0; j < _knownEV; j++)
			{
				utmp->MvDot(*CV(BQ, j), &eta, Anasazi::NO_CONJ);
				utmp->MvAddMv(ScalarOne, *utmp, -eta[0], *CV(_evecs, j));
			}
			
			// Now B-normalise the vector ...
			if (_B.get() != 0) {
				ok=_B->Apply(*utmp, *mvtmp);
				assert(ok==Anasazi::Ok);
				utmp->MvDot(*mvtmp, &eta, Anasazi::NO_CONJ);
			} else utmp->MvDot(*utmp, &eta, Anasazi::NO_CONJ);
			
			eta[0]=sqrt(eta[0]); //FIXME: replace with corrected Teuchos-function
			
			l=Teuchos::ScalarTraits<ScalarType>::magnitude(eta[0]);
			utmp->MvAddMv(ScalarOne / eta[0], *utmp, ScalarZero, *mvtmp);
		}
		
		//update BU
		sublist[0]=i;
		if (_B.get()) {
			ok=_B->Apply(*utmp, *CV(BU, i));
			assert(ok==Anasazi::Ok);
        	else CV(BU, i)->SetBlock(*utmp, sublist);
		
		delete utmp;
	}
	
	//update AU
	sublist.resize(SearchSpaceSize);
	for(int i=0; i<SearchSpaceSize; ++i) sublist[i]=i;
      
	ok=_A->Apply(*CV(U, sublist), *CV(AU, sublist));
	assert(ok==Anasazi::Ok);
	
	//update G (projected eigenvalue problem)
	PE.Restart();
	PE.Add(*CV(U, sublist), *CV(AU, sublist), *CV(BU, sublist));
	#endif
	
      }
      
      #ifdef _SAB_TIMING
      ftime(&tbRest);
      times(&tmRest);
      
      tRest1+=(tbRest.time-tb0.time)*1000+tbRest.millitm-tb0.millitm;
      tRest2+=tmRest.tms_utime-tm0.tms_utime+tmRest.tms_stime-tm0.tms_stime;
      #endif

      // Update the iteration step counter
      _iters++;
      
      outer++;
    }

    #ifdef _SAB_compare
    fclose(file);
    #endif
    
    #ifdef _SAB_DEBUG
    fclose(file2);
    fclose(fileAU);
    fclose(fileBU);
    fclose(fileTheta);
    fclose(fileZ);
    fclose(fileKsys);
    fclose(fileAsys);
    fclose(fileR);
    fclose(fileMisc);
    fclose(fileC);
    fclose(fileLk);
    fclose(fileRk);
    #endif
    
    #ifdef _SAB_STATISTICS
    ftime(&tbTot);
    times(&tmTot);
    
    if (_om->doPrint()) {
    	_os << "Time: " << (tbTot.time-tbTot0.time)*1000+tbTot.millitm-tbTot0.millitm << " (" << tmTot.tms_utime-tmTot0.tms_utime+tmTot.tms_stime-tmTot0.tms_stime << ")" << endl;
    
    	_os << "Iterations: " << _iters << " (GMRES: " << numtotGMRES << ")" << endl;
    }
    #endif
      
    #ifdef _SAB_TIMING
    if (_om->doPrint()) 
    {      
    	_os << "Orthog: " << tOrth1 << " ms  (" << tOrth2 << ")" << endl;
    	_os << "Extract: " << tExt1 << " ms  (" << tExt2 << ")" << endl;
    	_os << "Conv: " << tConv1 << " ms  (" << tConv2 << ")" << endl;
    	_os << "Corr: " << tCorr1 << " ms  (" << tCorr2 << ")" << endl;
    	_os << "GMRES: " << tGMRES1 << " ms  (" << tGMRES2 << ")" << endl;
    	_os << "Rest: " << tRest1 << " ms (" << tRest2 << ")" << endl;
    }
    #endif
    
    //return(Ok);
    if (_knownEV==_nev) return Ok;
    else return Unconverged;
    
  } // solve()

  // ==========================================================================
  template <class ScalarType, class MV, class OP>
  void BlockJacobiDavidson<ScalarType,MV,OP>::
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
          << std::setw(24) << std::right << "Relative Residual"
          << endl;
      _os << "------------------------------------------------------"<<endl;

      if (_knownEV > 0) 
      {
        for (int i = 0; i < _knownEV; i++) 
        {
          _os << std::setw(16) << std::right << (*_evals)[i]
              << std::setw(24) << std::right << _normR[i]
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
