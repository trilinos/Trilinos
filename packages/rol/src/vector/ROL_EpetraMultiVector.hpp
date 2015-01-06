// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_EPETRAMULTIVECTOR_H
#define ROL_EPETRAMULTIVECTOR_H

#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "ROL_Vector.hpp"

/** \class ROL::EpetraMultiVector
    \brief Implements the ROL::Vector interface for an Epetra_MultiVector.
*/

namespace ROL {

template <class Real>
class EpetraMultiVector : public Vector<Real> {
private:

  Teuchos::RCP<Epetra_MultiVector>  epetra_vec_;

public:
  virtual ~EpetraMultiVector() {}

  EpetraMultiVector(const Teuchos::RCP<Epetra_MultiVector> & epetra_vec) : epetra_vec_(epetra_vec) {}

  /** \brief Compute \f$y \leftarrow x + y\f$ where \f$y = \mbox{*this}\f$.
  */
  void plus( const Vector<Real> &x ) {
    const EpetraMultiVector &ex = Teuchos::dyn_cast<const EpetraMultiVector>(x);
    epetra_vec_->Update( 1.0, *ex.getVector(), 1.0 );
  }

  /** \brief Compute \f$y \leftarrow \alpha y\f$ where \f$y = \mbox{*this}\f$.
  */
  void scale( const Real alpha ) { 
    epetra_vec_->Scale( (double)alpha );
  }

  /** \brief Returns \f$ \langle y,x \rangle \f$ where \f$y = \mbox{*this}\f$.
  */
  Real dot( const Vector<Real> &x ) const {
    double val;
    const EpetraMultiVector &ex = Teuchos::dyn_cast<const EpetraMultiVector>(x);
    epetra_vec_->Dot( *ex.getVector(), &val );
    return (Real)val;
  }

  /** \brief Returns \f$ \| y \| \f$ where \f$y = \mbox{*this}\f$.
  */
  Real norm() const {
    double val;
    epetra_vec_->Norm2(&val);
    return (Real) val;
  } 

  /** \brief Clone to make a new (uninitialized) vector.
  */
  Teuchos::RCP<Vector<Real> > clone() const{
    return Teuchos::rcp(new EpetraMultiVector( 
  	     Teuchos::rcp(new Epetra_MultiVector(epetra_vec_->Map(),epetra_vec_->NumVectors(),false)) ));
  }

  /** \brief Compute \f$y \leftarrow \alpha x + y\f$ where \f$y = \mbox{*this}\f$.
  */
  virtual void axpy( const Real alpha, const Vector<Real> &x ) {
    const EpetraMultiVector &ex = Teuchos::dyn_cast<const EpetraMultiVector>(x);
    epetra_vec_->Update( alpha, *ex.getVector(), 1.0 );
  }

  /**  \brief Set to zero vector.
  */
  virtual void zero() {
    epetra_vec_->PutScalar(0.0);
  }

  virtual void PutScalar(Real alpha) {
    epetra_vec_->PutScalar((double) alpha);
  }

  /**  \brief Set \f$y \leftarrow x\f$ where \f$y = \mbox{*this}\f$.
  */
  virtual void set( const Vector<Real> &x ) {
    const EpetraMultiVector &ex = Teuchos::dyn_cast<const EpetraMultiVector>(x);
    epetra_vec_->Scale(1.0,*ex.getVector());
  }

  Teuchos::RCP<const Epetra_MultiVector> getVector() const {
    return this->epetra_vec_;
  }

  Teuchos::RCP<Epetra_MultiVector> getVector() {
    return this->epetra_vec_;
  }

  Teuchos::RCP<Vector<Real> > basis( const int i ) const {
    Teuchos::RCP<EpetraMultiVector> e = Teuchos::rcp( new EpetraMultiVector( Teuchos::rcp(new Epetra_MultiVector(epetra_vec_->Map(),epetra_vec_->NumVectors(),true)) ));
    const Epetra_BlockMap & domainMap = e->getVector()->Map();

    Epetra_Map linearMap(domainMap.NumGlobalElements(), domainMap.NumMyElements(), 0, domainMap.Comm());
    int lid = linearMap.LID(i);
    if(lid >=0)
      (*e->getVector())[0][lid]= 1.0;

    return e;

    /*


    // Build IntVector of GIDs on all processors.
    const Epetra_Comm & comm = domainMap.Comm();
    int numMyElements = domainMap.NumMyElements();
    Epetra_BlockMap allGidsMap(-1, numMyElements, 1, 0, comm);
    Epetra_IntVector allGids(allGidsMap);
    for (int j=0; j<numMyElements; j++) {allGids[j] = domainMap.GID(j);}

    // Import my GIDs into an all-inclusive map. 
    int numGlobalElements = domainMap.NumGlobalElements();
    Epetra_LocalMap allGidsOnRootMap(numGlobalElements, 0, comm);
    Epetra_Import importer(allGidsOnRootMap, allGidsMap);
    Epetra_IntVector allGidsOnRoot(allGidsOnRootMap);
    allGidsOnRoot.Import(allGids, importer, Insert);
    Epetra_Map rootDomainMap(-1, allGidsOnRoot.MyLength(), allGidsOnRoot.Values(), domainMap.IndexBase(), comm);

    for (int j = 0; j < this->dimension(); j++) {
      // Put 1's in slots
      int curGlobalCol = rootDomainMap.GID(i); // Should return same value on all processors
      if (domainMap.MyGID(curGlobalCol)){
        int curCol = domainMap.LID(curGlobalCol);
        (*e->getVector())[0][curCol]= 1.0;
      }
    }

    return e;

    */
  }

  int dimension() const {return epetra_vec_->GlobalLength();}


}; // class EpetraMultiVector

} // namespace ROL

#endif

