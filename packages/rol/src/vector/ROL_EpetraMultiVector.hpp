//@HEADER
// ***********************************************************************
//
//                     Rapid Optimization Library
//
// Questions? Contact:    Drew Kouri (dpkouri@sandia.gov)
//                      Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef ROL_EPETRAMULTIVECTOR_H
#define ROL_EPETRAMULTIVECTOR_H

#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "ROL_Vector.hpp"

/** \class ROL::EpetraMultiVector
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
    EpetraMultiVector &ex = Teuchos::dyn_cast<EpetraMultiVector>(const_cast <Vector<Real>&>(x));
    Teuchos::RCP<const Epetra_MultiVector> xvalptr = ex.getVector();
    epetra_vec_->Update( 1.0, *xvalptr, 1.0 );
  }

  /** \brief Compute \f$y \leftarrow \alpha y\f$ where \f$y = \mbox{*this}\f$.
  */
  void scale( const Real alpha ) { 
    epetra_vec_->Scale( (double)alpha );
  }

  /** \brief Returns \f$ \langle y,x \rangle \f$ where \f$y = \mbox{*this}\f$.
  */
  Real dot( const Vector<Real> &x ) const {
    double val[1];
    EpetraMultiVector &ex = Teuchos::dyn_cast<EpetraMultiVector>(const_cast <Vector<Real>&>(x));
    Teuchos::RCP<const Epetra_MultiVector> xvalptr = ex.getVector();
    epetra_vec_->Dot( *xvalptr, val );
    return (Real)val[0];
  }

  /** \brief Returns \f$ \| y \| \f$ where \f$y = \mbox{*this}\f$.
  */
  Real norm() const {
    double val[1];
    epetra_vec_->Dot( *epetra_vec_, val );
    return (Real)sqrt(val[0]);
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
    EpetraMultiVector &ex = Teuchos::dyn_cast<EpetraMultiVector>(const_cast <Vector<Real>&>(x));
    Teuchos::RCP<const Epetra_MultiVector> xvalptr = ex.getVector();
    epetra_vec_->Update( alpha, *xvalptr, 1.0 );
  }

  /**  \brief Set to zero vector.
  */
  virtual void zero() {
    epetra_vec_->PutScalar(0.0);
  }

  /**  \brief Set \f$y \leftarrow x\f$ where \f$y = \mbox{*this}\f$.
  */
  virtual void set( const Vector<Real> &x ) {
    EpetraMultiVector &ex = Teuchos::dyn_cast<EpetraMultiVector>(const_cast <Vector<Real>&>(x));
    Teuchos::RCP<const Epetra_MultiVector> xvalptr = ex.getVector();
    epetra_vec_->Scale(1.0,*xvalptr);
  }

  Teuchos::RCP<const Epetra_MultiVector> getVector() const {
    return this->epetra_vec_;
  }

  Teuchos::RCP<Vector<Real> > basis( const int i ) const {
    Teuchos::RCP<EpetraMultiVector> e = Teuchos::rcp( new EpetraMultiVector( Teuchos::rcp(new Epetra_MultiVector(epetra_vec_->Map(),epetra_vec_->NumVectors(),true)) ));
    const Epetra_BlockMap & domainMap = const_cast <Epetra_BlockMap &> (const_cast <Epetra_MultiVector &> ((*e->getVector())).Map());

    // Build IntVector of GIDs on all processors.
    const Epetra_Comm & comm = domainMap.Comm();
    int numMyElements = domainMap.NumMyElements();
    Epetra_BlockMap allGidsMap(-1, numMyElements, 1, 0, comm);
    Epetra_IntVector allGids(allGidsMap);
    for (int j=0; j<numMyElements; j++) {allGids[j] = domainMap.GID64(j);}

    // Import my GIDs into an all-inclusive map. 
    int numGlobalElements = domainMap.NumGlobalElements64();
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
        (const_cast <Epetra_MultiVector &> (*e->getVector()))[0][curCol]= 1.0;
      }
    }

    return e;
  }

  int dimension() const {return epetra_vec_->GlobalLength();}


}; // class EpetraMultiVector

} // namespace ROL

#endif

