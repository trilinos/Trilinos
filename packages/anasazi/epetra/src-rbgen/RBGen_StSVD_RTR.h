// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RBGEN_STSVD_RTR_H
#define RBGEN_STSVD_RTR_H

#include "RBGen_PODMethod.hpp"
#include "RBGen_Method.hpp"
#include "AnasaziOrthoManager.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_LAPACK.hpp"
#include "AnasaziSolverUtils.hpp"


//
// This class encapsulates a dominant SVD computation via Riemannian Trust-Region manifold
// optimization.
// 
// We are given a snapshot matrix A of dimension m by n. We wish to find the rank-k 
// dominant SVD.
//
// The manifold is M = St(k,m) x St(k,n), where St(k,p) is the manifold of rank-k orthonormal
// bases of R^p (the orthogonal Stiefel manifold).
//
// The objective function is
// f : M --> R 
//   : (U,V) --> trace(U^T A V N), 
// where N is a diagonal matrix with distinct, ascending (or descending) strictly positive elements.
// For example, N = diag(5,4,3,2,1)
// N serves to order the singular vectors in U and V
// In the case where the dominant singular values are distinct, this function has a unique global maximizer,
// which is a strict global maximizer. Otherwise, there is a unique region where the global maximum is reached, 
// which corresponds to rotations of the singular vectors corresponding to the non-distinct singular values.
//
// This solver applies the Riemannian Trust-Region method (Absil, Baker and Gallivan) to maximize f 
// over M.
//

namespace RBGen {

  //! Class for producing a POD basis using a trust-region optimization on the Stiefel manifold.
  class StSVDRTR : public virtual Method<Epetra_MultiVector,Epetra_Operator>, public virtual PODMethod<double> {

  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    StSVDRTR();

    //! Destructor.
    virtual ~StSVDRTR() {};
    //@}

    //! @name Computation Methods
    //@{

    //! Computes bases for the left and (optionally) right singular subspaces, along with singular vaues.
    void computeBasis();

    //! Update the current basis using a new set of snapshots.
    void updateBasis( const Teuchos::RCP< Epetra_MultiVector >& update_ss );

    //@}

    //! @name Get Methods
    //@{

    //! Return a basis for the left singular subspace.
    Teuchos::RCP<const Epetra_MultiVector> getBasis() const;

    //! Return a basis for the right singular subspace.
    Teuchos::RCP<const Epetra_MultiVector> getRightBasis() const;

    //! Return the singular values.
    std::vector<double> getSingularValues() const;

    //! Return the cummulative wall-clock time.
    double getCompTime() const { return timerComp_.totalElapsedTime(); }

    //! Return the scaled residual norms.
    std::vector<double> getResNorms();

    //@}

    //! @name Set Methods
    //@{

    //! Initialize the method with the given parameter list and snapshot set.
    void Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params,
                     const Teuchos::RCP< const Epetra_MultiVector >& init,
                     const Teuchos::RCP< RBGen::FileIOHandler< Epetra_Operator > >& fileio = Teuchos::null );

    void Reset( const Teuchos::RCP<Epetra_MultiVector>& new_ss );

    //@}

    //! @name Status Methods
    //@{

    bool isInitialized() { return isInitialized_; }

    //@}

  protected:

    typedef Anasazi::SolverUtils<double,Epetra_MultiVector,Epetra_Operator> Utils;
    // debugging checklist
    struct CheckList {
      // outer checks
      // U,V are ortho
      bool checkUV;
      // R is residual
      bool checkRes;
      // UAVNsym, VAUNsym
      bool checkSyms;
      // check sigma
      bool checkSigma;
      // check f(U,V)
      bool checkF;

      // inner checks
      // tangency
      bool checkE, checkHE, checkD, checkHD, checkR;
      // length
      bool checkElen, checkRlen;
      // conjugacy
      bool checkEHR, checkDHR;

      CheckList() : checkUV(false), checkRes(false), checkSyms(false), 
                    checkSigma(false), checkF(false),
                    checkE(false), checkHE(false), checkD(false), checkHD(false), checkR(false),
                    checkElen(false), checkRlen(false),
                    checkEHR(false), checkDHR(false) {};
    };
    // subproblem return codes
    enum trRetType {
      UNINITIALIZED = 0,
      MAXIMUM_ITERATIONS,
      NEGATIVE_CURVATURE,
      EXCEEDED_TR,
      KAPPA_CONVERGENCE,
      THETA_CONVERGENCE,
      NOTHING
    };
    // these correspond to above
    std::vector<std::string> stopReasons_;
    typedef Teuchos::ScalarTraits<double> SCT;

    // private members 
    // Riemannian metric
    // g(eta,eta)
    double innerProduct( const Epetra_MultiVector &etaU, 
                         const Epetra_MultiVector &etaV ) const;
    // g(eta,zeta)
    double innerProduct( const Epetra_MultiVector &etaU, 
                         const Epetra_MultiVector &etaV, 
                         const Epetra_MultiVector &zetaU, 
                         const Epetra_MultiVector &zetaV ) const;
    // Retraction, in situ
    void retract( const Epetra_MultiVector &xU, 
                  const Epetra_MultiVector &xV, 
                  Epetra_MultiVector &etaU, 
                  Epetra_MultiVector &etaV ) const;
    // Apply Hessian
    void Hess( const Epetra_MultiVector &xU, 
               const Epetra_MultiVector &xV, 
               const Epetra_MultiVector &etaU, 
               const Epetra_MultiVector &etaV,
                     Epetra_MultiVector &HetaU,
                     Epetra_MultiVector &HetaV ) const;
    // Project back to tangent plane
    void Proj( const Epetra_MultiVector &xU, 
               const Epetra_MultiVector &xV, 
               Epetra_MultiVector &etaU, 
               Epetra_MultiVector &etaV ) const;
    // solve the trust-region subproblem
    void solveTRSubproblem();
    // private initialize routine
    void initialize();
    // compute sym(S) = 0.5*(S+S')
    void Sym(Epetra_MultiVector &S) const;
    // compute f(x) and other stuff
    void updateF();
    // compute residuals and their norms
    void updateResiduals();
    // debugging checks
    void Debug(const CheckList &chk, std::string where) const;
    // print status
    void printStatus() const;


    // Is this object initialized?
    bool isInitialized_;

    // Vector holding N
    std::vector<double> N_;

    // Pointer to the snapshots
    Teuchos::RCP<const Epetra_MultiVector> A_;
    // state multivecs
    Teuchos::RCP<Epetra_MultiVector> U_, V_,           // current iterate
                                     AU_, AV_,         // A*V and A'*U
                                     RU_, RV_,         // residual, gradient and residual of model minimization
                                     etaU_, etaV_,     // subproblem solution
                                     HeU_, HeV_,       // eta under Hessian
                                     deltaU_, deltaV_, // search directions
                                     HdU_, HdV_;       // search directions under Hessian

    // Vector holding singular values.
    std::vector<double> sigma_;

    // Convergence tolerance
    double tol_;

    // ortho manager: need this for qf()
    Teuchos::RCP< Anasazi::OrthoManager<double,Epetra_MultiVector> > ortho_;

    // cummulative timer
    Teuchos::Time timerComp_;

    // debug flag
    bool debug_;

    // verb level
    int verbLevel_;

    // residual norms: individual and combined
    std::vector<double> resNorms_;
    std::vector<double> resUNorms_, resVNorms_;
    double maxScaledNorm_;

    // trust-region state
    // initial and current trust-region radius
    double Delta0_, Delta_, Delta_bar_;
    // |eta|
    double etaLen_;
    // acceptance parameter
    double rhoPrime_;
    // norm of the initial gradient
    double normGrad0_;
    // number of outer iterations
    int iter_;
    // maximum number of outer iterations
    int maxOuterIters_;
    // maximum number of inner iterations
    int maxInnerIters_;
    // convergence parameters
    double conv_kappa_, conv_theta_;
    // most recent rho
    double rho_, rhonum_, rhoden_;
    bool tiny_rhonum_, zero_rhoden_, neg_rho_;
    // current objective function
    double fx_;
    // last inner stop
    trRetType innerStop_;
    // previous update was accepted or not
    bool accepted_;
    // num inner iterations
    int numInner_;
    // trust region adjustment
    std::string tradjust_;
    // dimensions of problem
    int m_, n_, rank_;
    // is V local or distributed
    bool localV_;
    // init with elementary vectors or random
    bool initElem_;

    // UAVNsym, VAUNsym
    // we need sym(U'*A*V*N) and sym(V'*A'*U*N)
    // for the Hessian
    Teuchos::RCP<Epetra_MultiVector> UAVNsym_, VAUNsym_;

    // DGESVD workspace
    Teuchos::LAPACK<int,double> lapack;
    Teuchos::RCP<Epetra_MultiVector> dgesvd_A_;
    std::vector<double> dgesvd_work_;
  };

} // end of RBGen namespace

#endif // RBGEN_STSVD_RTR_H
