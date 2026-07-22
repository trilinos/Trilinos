// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

  ROL::Ptr<Epetra_MultiVector>  epetra_vec_;

public:
  virtual ~EpetraMultiVector() {}

  EpetraMultiVector(const ROL::Ptr<Epetra_MultiVector> & epetra_vec) : epetra_vec_(epetra_vec) {}

  /** \brief Compute \f$y \leftarrow x + y\f$ where \f$y = \mbox{*this}\f$.
  */
  void plus( const Vector<Real> &x ) {
    const EpetraMultiVector &ex = dynamic_cast<const EpetraMultiVector&>(x);
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
    const EpetraMultiVector &ex = dynamic_cast<const EpetraMultiVector&>(x);
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
  ROL::Ptr<Vector<Real> > clone() const{
    return ROL::makePtr<EpetraMultiVector>( 
  	     ROL::makePtr<Epetra_MultiVector>(epetra_vec_->Map(),epetra_vec_->NumVectors(),false) );
  }

  /** \brief Compute \f$y \leftarrow \alpha x + y\f$ where \f$y = \mbox{*this}\f$.
  */
  void axpy( const Real alpha, const Vector<Real> &x ) {
    const EpetraMultiVector &ex = dynamic_cast<const EpetraMultiVector&>(x);
    epetra_vec_->Update( alpha, *ex.getVector(), 1.0 );
  }

  /**  \brief Set to zero vector.
  */
  void zero() {
    epetra_vec_->PutScalar(0.0);
  }

  void PutScalar(const Real C) {
    epetra_vec_->PutScalar(static_cast<double>(C));
  }

  void setScalar(const Real C) {
    epetra_vec_->PutScalar(static_cast<double>(C));
  }

  void randomize(const Real l=0.0, const Real u=1.0) {
    epetra_vec_->Random();                         // Random numbers between -1 and 1
    this->scale(static_cast<Real>(0.5)*(u-l));     // Random numbers between (l-u)/2 and (u-l)/2
    Ptr<Vector<Real>> copy = this->clone();
    copy->setScalar(static_cast<Real>(0.5)*(u+l)); // Constant vector with values (u+l)/2
    this->plus(*copy);                             // Random numbers between l and u
  }

  /**  \brief Set \f$y \leftarrow x\f$ where \f$y = \mbox{*this}\f$.
  */
  void set( const Vector<Real> &x ) {
    const EpetraMultiVector &ex = dynamic_cast<const EpetraMultiVector&>(x);
    epetra_vec_->Scale(1.0,*ex.getVector());
  }

  ROL::Ptr<const Epetra_MultiVector> getVector() const {
    return this->epetra_vec_;
  }

  ROL::Ptr<Epetra_MultiVector> getVector() {
    return this->epetra_vec_;
  }

  ROL::Ptr<Vector<Real> > basis( const int i ) const {
    ROL::Ptr<EpetraMultiVector> e = 
    ROL::makePtr<EpetraMultiVector>( ROL::makePtr<Epetra_MultiVector>(epetra_vec_->Map(),epetra_vec_->NumVectors(),true));
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

