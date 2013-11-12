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

  Teuchos::RCP<const Epetra_MultiVector > getVector() const {
    return this->epetra_vec_;
  }

}; // class EpetraMultiVector

} // namespace ROL

#endif

