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

/** @ingroup la_group
 * \class ROL::MultiVectorDefault
 * \brief Default implementation of the ROL::MultiVector container class
 */


#ifndef ROL_MULTIVECTOR_DEFAULT_HPP
#define ROL_MULTIVECTOR_DEFAULT_HPP

#include "ROL_MultiVector.hpp"

namespace ROL {

template<class Real>
class MultiVectorDefault : public MultiVector<Real> {

    typedef Vector<Real>              V;       // Single vector 
    typedef Teuchos::RCP<V>           PV;      // Pointer to a vector
    typedef Teuchos::ArrayRCP<PV>     APV;     // Array of pointers to vectors
    typedef MultiVector<Real>         MV;      // Instance of this class
    typedef Teuchos::RCP<MV>          PMV;     // Pointer to an instance of this class 

    private:
        APV mvec_;             // Array of pointers to vectors 
        int numVectors_;       // number of vectors (elements in the array)
        int length_;           // number of elements in a vector
  
        virtual bool dimensionMismatch(const MV &A) const {
            
            bool equalWidth = ( A.getNumberOfVectors() != numVectors_ );
            bool equalLength = ( A.getLength() != length_ );
            return ( equalWidth | equalLength );
        } 

    public:

        // Create a MultiVector from an array of pointers to vectors
        MultiVectorDefault(APV mvec) : mvec_(mvec), 
                                       numVectors_(mvec.size()),
                                       length_(mvec[0]->dimension()) {}

        // Create a MultiVector from a pointer to a single vector
        MultiVectorDefault(PV vec) : mvec_(APV(1,vec)), 
                                     numVectors_(1),
                                     length_(vec->dimension()) {} 

        ~MultiVectorDefault() {}

        // Make a new MultiVector of the same dimensions 
        PMV clone() const {
            APV x(numVectors_);
            for(int i=0;i<numVectors_;++i) {
                x[i] = mvec_[i]->clone();
            }    
            return Teuchos::rcp(new MultiVectorDefault<Real>(x));
        } 
 
       // Make a new MultiVector of specified dimension
        PMV clone( const int numvecs ) const {
            APV x(numvecs);
            for(int i=0;i<numvecs;++i) {
                x[i] = mvec_[i]->clone();
            }    
            return Teuchos::rcp(new MultiVectorDefault<Real>(x));
        } 




        // Make a deep copy of this MultiVector
        PMV deepCopy() const {
            APV x(numVectors_);
            for(int i=0;i<numVectors_;++i) {
                x[i] = mvec_[i]->clone();
                x[i]->set(*mvec_[i]);
            }
            return Teuchos::rcp(new MultiVectorDefault<Real>(x));
        }

        // Make a deep copy specified vectors in the MultiVector
        PMV deepCopy(const std::vector<int> &index) const {
            int n = index.size();
            APV x(n);
            for(int i=0;i<n;++i) {
                int j = index[i]; 
                x[i] = mvec_[j]->clone();
                x[i]->set(*mvec_[j]);
            }
            return Teuchos::rcp(new MultiVectorDefault<Real>(x));
        } 

        // Make a shallow copy specified vectors in the MultiVector
        PMV shallowCopy(const std::vector<int> &index) {
            int n = index.size();
            APV x(n);
            for(int i=0;i<n;++i) {
                int j = index[i]; 
                x[i] = mvec_[j]; 
            }
            return Teuchos::rcp(new MultiVectorDefault<Real>(x));
        } 

        // Make a const shallow copy specified vectors in the MultiVector
        const PMV shallowCopyConst(const std::vector<int> &index) const {
            int n = index.size();
            APV x(n);
            for(int i=0;i<n;++i) {
                int j = index[i]; 
                x[i] = mvec_[j]; 
            }
            return Teuchos::rcp(new MultiVectorDefault<Real>(x));
        } 

        // Get the number of elements of a vector in the MultiVector
        ptrdiff_t getLength() const {
            return length_;
        }

        // Get the number of vectors in the MultiVector
        int getNumberOfVectors() const {
            return numVectors_;
        }

        void axpy(const Real alpha, const MV& x) {
            for(int i=0;i<numVectors_;++i) {
                mvec_[i]->axpy(alpha,*(x.getVector(i)));
            }  
        } 

        // Generic BLAS level 3 matrix multiplication
        // \f$\text{this}\leftarrow \alpha A B+\beta\text{this}\f$   
        void gemm(const Real alpha,
                  const MV& A,
                  const Teuchos::SerialDenseMatrix<int,Real> &B,
                  const Real beta) {

            TEUCHOS_TEST_FOR_EXCEPTION( this->dimensionMismatch(A),
                std::invalid_argument,
                "Error: MultiVectors must have the same dimensions.");

            // Scale this by beta
            this->scale(beta);

            for(int i=0;i<numVectors_;++i) {
                for(int j=0;j<numVectors_;++j) {
                    mvec_[i]->axpy(alpha*B(i,j),*A.getVector(j));  
                }
            }
        } 

        // Scale the MultiVector by a single scalar alpha 
        // \f$\text{this}\leftarrow\alpha\text{this}\f$
        void scale(const Real alpha) {
            for(int i=0;i<numVectors_;++i) {
                mvec_[i]->scale(alpha);  
            }
        }

        // Scale each vector in the MultiVector by a different alpha
        // \f$\text{this}[i]\leftarrow\alpha[i]\text{this}[i]\f$
        void scale(const std::vector<Real> &alpha) {

            TEUCHOS_TEST_FOR_EXCEPTION( alpha.size() != numVectors_,
                std::invalid_argument,
                "Error: alpha must have the same length as the number of vectors.");  
 
            for(int i=0;i<numVectors_;++i) {
                mvec_[i]->scale(alpha[i]);  
            }
        } 

        // Set the MultiVector equal to another MultiVector
        void set(const MV &A) {

            TEUCHOS_TEST_FOR_EXCEPTION( this->dimensionMismatch(A),
                std::invalid_argument,
                "Error: MultiVectors must have the same dimensions.");

            for(int i=0;i<numVectors_;++i) {
                mvec_[i]->set(*A.getVector(i));
            }
        }

        
        // Set some of the vectors in this MultiVector equal to corresponding 
        // vectors in another MultiVector
        void set(const MV &A, const std::vector<int> &index) {

            TEUCHOS_TEST_FOR_EXCEPTION( this->dimensionMismatch(A),
                std::invalid_argument,
                "Error: MultiVectors must have the same dimensions.");

            int n = index.size();
             
            for(int i=0;i<n;++i) {
                int k = index[i];
                mvec_[k]->set(*A.getVector(k));
            }
        }

        // Compute \f$\alpha A^\top \text{this}\f$ 
        void innerProducts(const Real alpha,
                           const MV &A,
                           Teuchos::SerialDenseMatrix<int,Real> &B) const {

            TEUCHOS_TEST_FOR_EXCEPTION( this->dimensionMismatch(A),
                std::invalid_argument,
                "Error: MultiVectors must have the same dimensions.");

            for(int i=0;i<numVectors_;++i) {
                for(int j=0;j<numVectors_;++j) {
                    B(i,j) = alpha*mvec_[j]->dot(*A.getVector(i));
                }  
            }
        }                  

        // Compute dot products of pairs of vectors
        void dots(const MV &A,
                  std::vector<Real> &b) const {

            TEUCHOS_TEST_FOR_EXCEPTION( this->dimensionMismatch(A),
                std::invalid_argument,
                "Error: MultiVectors must have the same dimensions.");

            for(int i=0;i<numVectors_;++i) {
                b[i] = mvec_[i]->dot(*A.getVector(i));
            }    
        } 

        // Compute the norm of each vector in the MultiVector
        void norms(std::vector<Real> &normvec) const {

            TEUCHOS_TEST_FOR_EXCEPTION( normvec.size()!=numVectors_,
                std::invalid_argument,
                "Error: normvec must have the same length as number of vectors.");

            for(int i=0;i<numVectors_;++i) {
                normvec[i] = mvec_[i]->norm(); 
            }    
        }

        // Zero each of the vectors in the MultiVector
        void zero() {
            for(int i=0;i<numVectors_;++i) {
                mvec_[i]->zero();
            }    
        }
         
        // Return a pointer to the ith vector
        PV getVector(int i) const {

            TEUCHOS_TEST_FOR_EXCEPTION( i>=numVectors_, 
                std::invalid_argument,
                "Error: index out of bounds");

            return mvec_[i]; 
        } 


};
}


#endif
