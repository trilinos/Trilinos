// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RBGEN_INCSVD_POD_H
#define RBGEN_INCSVD_POD_H

#include "RBGen_PODMethod.hpp"
#include "RBGen_Method.hpp"
#include "RBGen_Filter.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziOrthoManager.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Time.hpp"


//
// computeBasis()
//       |
//   makePass()     ___
//   while () {    /    expand()
//      incStep()  ---- SVD()
//   }             \___ shrink()
//
// makePass(), expand() and shrink() are pure virtual
//
// makePass() is implemented in a base class that decides 
// if a method is Multipass, and if so, in what manner
// makePass() puts the data for the next pass into the proper columns
// of U_, then calls incstep (telling it how many new columns)
//
// incstep calls expand to construct expanded U_ and V_, and B_
// it computes the SVD of B_
// it passes this data to shrink(), which shrink according to the 
// base class.
//
// updateBasis() allows the current factorization to be updated
// by appending new snapshots. This allows a truly incremental 
// factorization. It is also pure virtual, and is implemented
// along side makePass(). It may not be implemented in all 
// IncSVDPOD methods; in fact, it is currently implemented 
// only by derived classes of ISVDSingle
//
// expand(),shrink() are implemented in a base class that 
// decides the representation: UDV, QRW, QBW
// 
//                                   IncSVDPOD
//                                       |
//      -------------------------------------------------------------------
//      |       |        |           |            |            |          |
//  ISVDUDV  ISVDQRW  ISVDQBW   ISVDMultiCD ISVDMultiSDA ISVDMultiSDB ISVDSingle
//      |       |        |           |            |            |          |
//      ------------------           --------------------------------------
//              \                                       /
//               \                                     /
//                \                                   /
//                 \                                 /
//                  \                               /
//                   \---  Concrete Base Class ----/
//
// Then a concrete base class (one of 3x4==12 varieties) is formed through
// inheritence. This is the Template Pattern type of Design Pattern.
//

namespace RBGen {

  //! Class for producing a basis using the Incremental SVD
  class IncSVDPOD : public virtual Method<Epetra_MultiVector,Epetra_Operator>, public virtual PODMethod<double> {

  public:
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    IncSVDPOD();

    //! Destructor.
    virtual ~IncSVDPOD() {};
    //@}

    //! @name Computation Methods
    //@{

    //! Computes bases for the left and (optionally) right singular subspaces, along with singular vaues.
    void computeBasis();

    //! Update the current basis by appending new snapshots.
    virtual void updateBasis( const Teuchos::RCP< Epetra_MultiVector >& update_ss ) = 0;

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
    const std::vector<double> & getResNorms();

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

    // private member for performing inc steps
    void incStep(int lup);
    virtual void makePass() = 0;
    virtual void expand(const int lup) = 0;
    virtual void shrink(const int down, std::vector<double> &S, Epetra_SerialDenseMatrix &U, Epetra_SerialDenseMatrix &V) = 0;

    // Is this object initialized?
    bool isInitialized_;

    // Singular value filter
    Teuchos::RCP< Filter<double> > filter_;

    // Max allocation size.
    // The maximum rank DURING any step:
    //    lup <= maxBasisSize_ - curRank_
    int maxBasisSize_;

    // Current rank of the factorization
    int curRank_;

    // Pointers to the snapshots and reduced basis.
    Teuchos::RCP<const Epetra_MultiVector> A_;
    Teuchos::RCP<Epetra_MultiVector> U_, V_;

    // SerialDenseMatrix holding current core matrix B
    Teuchos::RCP<Epetra_SerialDenseMatrix> B_;

    // Vector holding singular values.
    std::vector<double> sigma_;

    // Number of snapshots processed thus far
    int numProc_;

    // Maximum allowable number of passes through A
    int maxNumPasses_;

    // Current number of passes through A
    int curNumPasses_;

    // Convergence tolerance
    double tol_;

    // min,max number of update vectors
    int lmin_;
    int lmax_;
    int startRank_;

    // ortho manager
    Teuchos::RCP< Anasazi::OrthoManager<double,Epetra_MultiVector> > ortho_;

    // cummulative timer
    Teuchos::Time timerComp_;

    // debug flag
    bool debug_;

    // verb level
    int verbLevel_;

    // residual norms
    std::vector<double> resNorms_;

    // maximum number of data set columns
    // i.e., maximum number of rows in V
    int maxNumCols_;

    // is V locally replicated or globally distributed?
    bool Vlocal_;
  };

} // end of RBGen namespace

#endif // RBGEN_INCSVD_POD_H
