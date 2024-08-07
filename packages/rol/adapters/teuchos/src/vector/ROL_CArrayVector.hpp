// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CARRAY_VECTOR_H
#define ROL_CARRAY_VECTOR_H

#include "ROL_Vector.hpp"
#include "Teuchos_ArrayRCP.hpp"

/** \class ROL::CArrayVector
    \brief Provides the C array implementation of the ROL::Vector interface
           for use with NumPy->C Array passing by pointer. This class in intended
           to be used with the Python ROL interface. 
 
    \details See https://github.com/cython/cython/wiki/tutorials-NumpyPointerToC
    \author Greg von Winckel (gvonwin@sandia.gov)
*/

namespace ROL {

template <class Real, class Element=Real>
class CArrayVector : public Vector<Real> {

    private:
        unsigned int dim_; 
        Teuchos::ArrayRCP<Element> array_;
    public:
        // Create from C array with raw ptr
        CArrayVector(Element* array, unsigned int dim) : 
            dim_(dim),array_(array,0,dim,false) {}
 
        // Create from Teuchos ArrayRCP
        CArrayVector(const Teuchos::ArrayRCP<Element> array) : 
            dim_(array.size()),array_(array) {}
 
        // Create an array of all zeros
        CArrayVector(unsigned int dim) : 
            dim_(dim),array_(dim) {}

        Teuchos::ArrayRCP<Element> getVector() const {
            return array_;
        }

        Teuchos::ArrayRCP<Element> getVector() {
            return array_;
        }
       
        void plus( const Vector<Real> &x ) {
            // Need to make sure object has a getVector method
            const CArrayVector &ex = dynamic_cast<const CArrayVector&>(x);            

            Teuchos::ArrayRCP<Element> xp(ex.getVector());
            for(unsigned int i=0; i<dim_; ++i) {
                (array_)[i] += xp[i];
            }
        }

        Real dot( const Vector<Real> &x ) const {
            const CArrayVector &ex = dynamic_cast<const CArrayVector&>(x);
            
            Teuchos::ArrayRCP<Element> xp(ex.getVector());
            Real val = 0;
            for(unsigned int i=0; i<dim_; ++i){
                val += (array_)[i]*(xp)[i];
            }            
             return val;
        }
      
        Real norm() const {
            Real val = 0;
            val = std::sqrt( dot(*this) );
            return val; 
        }

        void scale( const Real alpha ) {
           for(unsigned int i=0; i<dim_; ++i) {
               (array_)[i] *= alpha;
           }
        } 
 
        int dimension() const {
            return dim_;
        }

        ROL::Ptr<Vector<Real> > clone() const {
            return ROL::makePtr<CArrayVector>( Teuchos::ArrayRCP<Element>(dim_) );   
        }
        
        ROL::Ptr<Vector<Real> > basis (const int i ) const {
            ROL::Ptr<CArrayVector> e = 
                ROL::makePtr<CArrayVector>(dim_);
            (e->getVector())[i] = 1.0;
            return e; 
        }         

        void setScalar( const Real C ) {
           for(unsigned int i=0; i<dim_; ++i) {
               (array_)[i] = C;
           }
        } 

        void randomize( const Real l=0.0, const Real u=1.0 ) {
           Real a = (u-l);
           Real b = l;
           Real x(0);
           for(unsigned int i=0; i<dim_; ++i) {
               x = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
               (array_)[i] = a*x + b;
           }
        }

};    


}

#endif

