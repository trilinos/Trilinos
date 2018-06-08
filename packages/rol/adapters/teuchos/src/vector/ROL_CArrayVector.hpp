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

};    


}

#endif

