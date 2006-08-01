
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

#ifndef ANASAZI_SIMPLE_LOBPCG_SOLMGR_HPP
#define ANASAZI_SIMPLE_LOBPCG_SOLMGR_HPP

/*! \file AnasaziSimpleLOBPCGSolMgr.hpp
  \brief The Anasazi::SimpleLOBPCGSolMgr provides a simple solver manager over the LOBPCG 
  eigensolver.
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"

#include "AnasaziEigenproblem.hpp"
#include "AnasaziSolverManager.hpp"

#include "AnasaziLOBPCG.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziSVQBOrthoManager.hpp"
#include "AnasaziStatusTestMaxIters.hpp"
#include "AnasaziStatusTestResNorm.hpp"
#include "AnasaziStatusTestCombo.hpp"
#include "AnasaziStatusTestOutput.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "AnasaziModalSolverUtils.hpp"

#include "Teuchos_TimeMonitor.hpp"

/*! \class Anasazi::SimpleLOBPCGSolMgr
  \brief The Anasazi::SimpleLOBPCGSolMgr provides a simple solver
  manager over the LOBPCG eigensolver. 

  Anasazi::SimpleLOBPCGSolMgr allows the user to specify convergence
  tolerance, verbosity level and block size. When block size is less than the
  number of requested eigenvalues specified in the eigenproblem, checkpointing
  is activated.

  The purpose of this solver manager was to provide an example of a simple
  solver manager, useful for demonstration as well as a jumping-off point for
  solvermanager development. Also, the solver manager is useful for testing
  some of the features of the Anasazi::LOBPCG eigensolver, principally the use
  of auxiliary vectors.

  This solver manager does not verify before quitting that the nev eigenvectors
  that have converged are also the smallest nev eigenvectors that are known.

  \ingroup anasazi_solvermanagers

  \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, Heidi Thornquist
*/

namespace Anasazi {

template<class ScalarType, class MV, class OP>
class SimpleLOBPCGSolMgr : public SolverManager<ScalarType,MV,OP> {

  private:
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    
  public:

  //!@name Constructors/Destructor 
  //@{ 

  /*! \brief Basic constructor for SimpleLOBPCGSolMgr.
   *
   * This constructor accepts the Eigenproblem to be solved in addition
   * to a parameter list of options for the solver manager. These options include the following:
   *   - "Which" - a \c string specifying the desired eigenvalues: SM, LM, SR or LR. Default: SR
   *   - "Block Size" - a \c int specifying the block size to be used by the underlying LOBPCG solver. Default: problem->getNEV()
   *   - "Maximum Iterations" - a \c int specifying the maximum number of iterations the underlying solver is allowed to perform. Default: 100
   *   - "Verbosity" - a sum of MsgType specifying the verbosity. Default: Anasazi::Errors
   *   - "Convergence Tolerance" - a \c MagnitudeType specifying the level that residual norms must reach to decide convergence. Default: machine precision
   */
  SimpleLOBPCGSolMgr( const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem,
                             Teuchos::ParameterList &pl );

  //! Destructor.
  virtual ~SimpleLOBPCGSolMgr() {};
  //@}
  
  //! @name Accessor methods
  //@{ 

  Eigenproblem<ScalarType,MV,OP>& getProblem() const {
    return *_problem;
  }

  //@}

  //! @name Solver application methods
  //@{ 
    
  /*! \brief This method performs possibly repeated calls to the underlying eigensolver's iterate() routine
   * until the problem has been solved (as decided by the solver manager) or the solver manager decides to 
   * quit.
   *
   * \returns ::ReturnType specifying:
   *    - ::Converged: the eigenproblem was solved to the specification required by the solver manager.
   *    - ::Unconverged: the eigenproblem was not solved to the specification desired by the solver manager
  */
  ReturnType solve();
  //@}

  private:
  Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > _problem;
  string _whch; 
  MagnitudeType _tol;
  int _verb;
  int _blockSize;
  int _maxIters;
};


////////////////////////////////////////////////////////////////////////////////////////
template<class ScalarType, class MV, class OP>
SimpleLOBPCGSolMgr<ScalarType,MV,OP>::SimpleLOBPCGSolMgr( 
        const Teuchos::RefCountPtr<Eigenproblem<ScalarType,MV,OP> > &problem,
        Teuchos::ParameterList &pl ) : 
  _problem(problem),
  _whch("LM"),
  _tol(1e-6),
  _verb(Anasazi::Errors),
  _blockSize(0),
  _maxIters(100)
{
  TEST_FOR_EXCEPTION(_problem == Teuchos::null,              std::invalid_argument, "Problem not given to solver manager.");
  TEST_FOR_EXCEPTION(!_problem->isProblemSet(),              std::invalid_argument, "Problem not set.");
  TEST_FOR_EXCEPTION(!_problem->isHermitian(),               std::invalid_argument, "Problem not symmetric.");
  TEST_FOR_EXCEPTION(_problem->getInitVec() == Teuchos::null,std::invalid_argument, "Problem does not contain initial vectors to clone from.");

  _whch = pl.get("Which","SR");
  TEST_FOR_EXCEPTION(_whch != "SM" && _whch != "LM" && _whch != "SR" && _whch != "LR",
                     AnasaziError,
                     "SimpleLOBPCGSolMgr: \"Which\" parameter must be SM, LM, SR or LR.");

  _tol = pl.get("Convergence Tolerance",_tol);
  TEST_FOR_EXCEPTION(_tol <= 0,
                     AnasaziError,
                     "SimpleLOBPCGSolMgr: \"Tolerance\" parameter must be strictly postiive.");

  _verb = pl.get("Verbosity",_verb);

  _blockSize= pl.get("Block Size",_problem->getNEV());
  TEST_FOR_EXCEPTION(_blockSize <= 0,
                     AnasaziError,
                     "SimpleLOBPCGSolMgr: \"Block Size\" parameter must be strictly positive.");

  _maxIters = pl.get("Maximum Iterations",_maxIters);
}



////////////////////////////////////////////////////////////////////////////////////////
template<class ScalarType, class MV, class OP>
ReturnType 
SimpleLOBPCGSolMgr<ScalarType,MV,OP>::solve() {

  // sort manager
  Teuchos::RefCountPtr<BasicSort<ScalarType,MV,OP> > sorter = Teuchos::rcp( new BasicSort<ScalarType,MV,OP>(_whch) );
  // output manager
  Teuchos::RefCountPtr<BasicOutputManager<ScalarType> > printer = Teuchos::rcp( new BasicOutputManager<ScalarType>(_verb) );
  // status tests
  Teuchos::RefCountPtr<StatusTestMaxIters<ScalarType,MV,OP> > max;
  if (_maxIters > 0) {
    max = Teuchos::rcp( new StatusTestMaxIters<ScalarType,MV,OP>(_maxIters) );
  }
  else {
    max = Teuchos::null;
  }
  Teuchos::RefCountPtr<StatusTestResNorm<ScalarType,MV,OP> > norm 
      = Teuchos::rcp( new StatusTestResNorm<ScalarType,MV,OP>(_tol) );
  Teuchos::Array< Teuchos::RefCountPtr<StatusTest<ScalarType,MV,OP> > > alltests;
  alltests.push_back(norm);
  if (max != Teuchos::null) alltests.push_back(max);
  Teuchos::RefCountPtr<StatusTestCombo<ScalarType,MV,OP> > combo 
      = Teuchos::rcp( new StatusTestCombo<ScalarType,MV,OP>(
              StatusTestCombo<ScalarType,MV,OP>::OR, alltests
        ));
  // printing StatusTest
  Teuchos::RefCountPtr<StatusTestOutput<ScalarType,MV,OP> > outputtest
    = Teuchos::rcp( new StatusTestOutput<ScalarType,MV,OP>( printer,combo,1,Passed ) );
  // orthomanager
  Teuchos::RefCountPtr<SVQBOrthoManager<ScalarType,MV,OP> > ortho 
    = Teuchos::rcp( new SVQBOrthoManager<ScalarType,MV,OP>(_problem->getM()) );
  // parameter list
  Teuchos::ParameterList plist;
  plist.set("Block Size",_blockSize);
  plist.set("Full Ortho",true);

  // create an LOBPCG solver
  Teuchos::RefCountPtr<LOBPCG<ScalarType,MV,OP> > lobpcg_solver 
    = Teuchos::rcp( new LOBPCG<ScalarType,MV,OP>(_problem,sorter,printer,outputtest,ortho,plist) );
  // add the auxillary vecs from the eigenproblem to the solver
  if (_problem->getAuxVecs() != Teuchos::null) {
    lobpcg_solver->setAuxVecs( Teuchos::tuple<Teuchos::RefCountPtr<const MV> >(_problem->getAuxVecs()) );
  }

  int numfound = 0;
  int nev = _problem->getNEV();
  Teuchos::Array< Teuchos::RefCountPtr<MV> > foundvecs;
  Teuchos::Array< Teuchos::RefCountPtr< std::vector<MagnitudeType> > > foundvals;
  while (numfound < nev) {
    // reduce the strain on norm test, if we are almost done
    if (nev - numfound < _blockSize) {
      norm->setQuorum(nev-numfound);
    }

    // tell the solver to iterate
    try {
      lobpcg_solver->iterate();
    }
    catch (std::exception e) {
      // we are a simple solver manager. we don't catch exceptions. set solution empty, then rethrow.
      printer->stream(Anasazi::Errors) << "Exception: " << e.what() << endl;
      Eigensolution<ScalarType,MV> sol;
      sol.numVecs = 0;
      _problem->setSolution(sol);
      throw;
    }

    // check the status tests
    if (norm->getStatus() == Passed) {

      int num = norm->howMany();
      // if num < _blockSize, it is because we are on the last iteration: num+numfound>=nev
      TEST_FOR_EXCEPTION(num < _blockSize && num+numfound < nev,
                         std::logic_error,
                         "Anasazi::SimpleLOBPCGSolMgr::solve(): logic error.");
      std::vector<int> ind = norm->whichVecs();
      // just grab the ones that we need
      if (num + numfound > nev) {
        num = nev - numfound;
        ind.resize(num);
      }

      // copy the converged eigenvectors
      Teuchos::RefCountPtr<MV> newvecs = MVT::CloneCopy(*lobpcg_solver->getEvecs(),ind);
      // store them
      foundvecs.push_back(newvecs);
      // add them as auxiliary vectors
      Teuchos::Array<Teuchos::RefCountPtr<const MV> > auxvecs = lobpcg_solver->getAuxVecs();
      auxvecs.push_back(newvecs);
      // setAuxVecs() will reset the solver to uninitialized, without messing with numIters()
      lobpcg_solver->setAuxVecs(auxvecs);

      // copy the converged eigenvalues
      Teuchos::RefCountPtr<std::vector<MagnitudeType> > newvals = Teuchos::rcp( new std::vector<MagnitudeType>(num) );
      std::vector<MagnitudeType> all = lobpcg_solver->getEigenvalues();
      for (int i=0; i<num; i++) {
        (*newvals)[i] = all[ind[i]];
      }
      foundvals.push_back(newvals);

      numfound += num;
    }
    else if (max != Teuchos::null && max->getStatus() == Passed) {

      int num = norm->howMany();
      std::vector<int> ind = norm->whichVecs();
      
      if (num > 0) {
        // copy the converged eigenvectors
        Teuchos::RefCountPtr<MV> newvecs = MVT::CloneCopy(*lobpcg_solver->getEvecs(),ind);
        // orthornormalize to be safe
        ortho->normalize(*newvecs,Teuchos::null,Teuchos::null);
        // store them
        foundvecs.push_back(newvecs);
        // don't bother adding them as auxiliary vectors; we have reached maxiters and are going to quit
        
        // copy the converged eigenvalues
        Teuchos::RefCountPtr<std::vector<MagnitudeType> > newvals = Teuchos::rcp( new std::vector<MagnitudeType>(num) );
        std::vector<MagnitudeType> all = lobpcg_solver->getEigenvalues();
        for (int i=0; i<num; i++) {
          (*newvals)[i] = all[ind[i]];
        }
        foundvals.push_back(newvals);
  
        numfound += num;
      }
      break;  // while(numfound < nev)
    }
    else {
      TEST_FOR_EXCEPTION(true,std::logic_error,"Anasazi::SimpleLOBPCGSolMgr::solve(): solver returned without satisfy status test.");
    }
  } // end of while(numfound < nev)

  TEST_FOR_EXCEPTION(foundvecs.size() != foundvals.size(),std::logic_error,"Anasazi::SimpleLOBPCGSolMgr::solve(): inconsistent array sizes");

  // create contiguous storage for all eigenvectors, eigenvalues
  Eigensolution<ScalarType,MV> sol;
  sol.numVecs = numfound;
  if (numfound > 0) {
    // allocate space for eigenvectors
    sol.Evecs = MVT::Clone(*_problem->getInitVec(),numfound);
  }
  else {
    sol.Evecs = Teuchos::null;
  }
  sol.Espace = sol.Evecs;
  // allocate space for eigenvalues
  sol.Evals.resize(numfound);
  // all real eigenvalues: set index vectors [0,...,numfound-1]
  sol.index.resize(numfound,0);
  // store eigenvectors, eigenvalues
  int curttl = 0;
  for (unsigned int i=0; i<foundvals.size(); i++) {
    TEST_FOR_EXCEPTION((signed int)(foundvals[i]->size()) != MVT::GetNumberVecs(*foundvecs[i]), std::logic_error, "Anasazi::SimpleLOBPCGSolMgr::solve(): inconsistent sizes");
    unsigned int lclnum = foundvals[i]->size();
    std::vector<int> lclind(lclnum);
    for (unsigned int j=0; j<lclnum; j++) lclind[j] = curttl+j;
    // put the eigenvectors
    MVT::SetBlock(*foundvecs[i],lclind,*sol.Evecs);
    // put the eigenvalues
    copy( foundvals[i]->begin(), foundvals[i]->end(), sol.Evals.begin()+curttl );

    curttl += lclnum;
  }
  TEST_FOR_EXCEPTION( curttl != sol.numVecs, std::logic_error, "Anasazi::SimpleLOBPCGSolMgr::solve(): inconsistent sizes");

  // sort the eigenvalues and permute the eigenvectors appropriately
  // this requires making a ScalarType version of our MagnitudeType eigenvalues
  if (numfound > 0) {
    std::vector<int> order(sol.numVecs);
    std::vector<ScalarType> vals_st(sol.numVecs);
    std::copy(sol.Evals.begin(),sol.Evals.end(),vals_st.begin());
    sorter->sort( lobpcg_solver.get(), sol.numVecs, &vals_st[0], &order );
    for (int i=0; i<sol.numVecs; i++) {
      sol.Evals[i] = SCT::real( vals_st[i] );
    }
    // now permute the eigenvectors according to order
    ModalSolverUtils<ScalarType,MV,OP> msutils(printer);
    msutils.permuteVectors(sol.numVecs,order,*sol.Evecs);
  }

  // print final summary
  lobpcg_solver->currentStatus(printer->stream(FinalSummary));

  // print timing information
  Teuchos::TimeMonitor::summarize(printer->stream(TimingDetails));

  // send the solution to the eigenproblem
  _problem->setSolution(sol);
  printer->stream(Debug) << "Returning " << sol.numVecs << " eigenpairs to eigenproblem." << endl;

  // return from SolMgr::solve()
  if (sol.numVecs < nev) return Unconverged;
  return Converged;
}




} // end Anasazi namespace

#endif /* ANASAZI_SIMPLE_LOBPCG_SOLMGR_HPP */
