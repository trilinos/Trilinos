//@HEADER
// ***********************************************************************
//
//                     Rapid Optimization Library
//
// Questions? Contact   Drew Kouri (dpkouri@sandia.gov)
//                    Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef ROL_STDVECTOR_H
#define ROL_STDVECTOR_H

#include "ROL_Vector.hpp"

/** \class ROL::StdVector
    \brief Provides the std::vector implementation of the vector space interface.
*/


namespace ROL {

template <class Real, class Element=Real>
class StdVector : public Vector<Real> {

private:

  Teuchos::RCP<std::vector<Element> >  std_vec_;

public:

  StdVector(const Teuchos::RCP<std::vector<Element> > & std_vec) : std_vec_(std_vec) {}

  void plus( const Vector<Real> &x ) {
    StdVector &ex = Teuchos::dyn_cast<StdVector>(const_cast <Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
    int dimension  = (int)(this->std_vec_->size());
    for (int i=0; i<dimension; i++) {
      (*this->std_vec_)[i] += (*xvalptr)[i];
    }
  }

  void scale( const Real alpha ) {
    int dimension = (int)(this->std_vec_->size());
    for (int i=0; i<dimension; i++) {
      (*this->std_vec_)[i] *= alpha;
    }
  }

  Real dot( const Vector<Real> &x ) const {
    StdVector & ex = Teuchos::dyn_cast<StdVector>(const_cast <Vector<Real> &>(x));
    Teuchos::RCP<const std::vector<Element> > xvalptr = ex.getVector();
    int dimension  = (int)(this->std_vec_->size());
    Real val = 0;
    for (int i=0; i<dimension; i++) {
      val += (*this->std_vec_)[i]*(*xvalptr)[i];
    }
    return val;
  }

  Real norm () const {
    Real val = 0;
    val = sqrt( this->dot(*this) );
    return val;
  }

  Teuchos::RCP<Vector<Real> > clone() const {
    return Teuchos::rcp( new StdVector( Teuchos::rcp(new std::vector<Element>(this->std_vec_->size())) ));
  }

  Teuchos::RCP<const std::vector<Element> > getVector() const {
    return this->std_vec_;
  }

  Teuchos::RCP<Vector<Real> > basis( const int i ) const {
    Teuchos::RCP<StdVector> e = Teuchos::rcp( new StdVector( Teuchos::rcp(new std::vector<Element>(this->std_vec_->size(), 0.0)) ));
    (const_cast <std::vector<Element> &> (*e->getVector()))[i]= 1.0;
    return e;
  }

  int dimension() {return this->std_vec_->size();}

}; // class StdVector

} // namespace ROL

#endif
