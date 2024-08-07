// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    typedef ROL::Ptr<V>               PV;      // Pointer to a vector
    typedef Teuchos::ArrayRCP<PV>     APV;     // Array of pointers to vectors
    typedef MultiVector<Real>         MV;      // Instance of this class
    typedef ROL::Ptr<MV>              PMV;     // Pointer to an instance of this class 

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

        // Create a MultiVector from a pointer to a constant vector

        ~MultiVectorDefault() {}

        // Make a new MultiVector of the same dimensions 
        PMV clone() const {
            APV x(numVectors_);
            for(int i=0;i<numVectors_;++i) {
                x[i] = mvec_[i]->clone();
            }    
            return ROL::makePtr<MultiVectorDefault<Real>>(x);
        } 
 
       // Make a new MultiVector of specified dimension
        PMV clone( const int numvecs ) const {
            APV x(numvecs);

            for(int i=0;i<numvecs;++i) {
                x[i] = mvec_[0]->clone();
            }    
            return ROL::makePtr<MultiVectorDefault<Real>>(x);
        } 




        // Make a deep copy of this MultiVector
        PMV deepCopy() const {
            APV x(numVectors_);
            for(int i=0;i<numVectors_;++i) {
                x[i] = mvec_[i]->clone();
                x[i]->set(*mvec_[i]);
            }
            return ROL::makePtr<MultiVectorDefault<Real>>(x);
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
            return ROL::makePtr<MultiVectorDefault<Real>>(x);
        } 

        // Make a shallow copy specified vectors in the MultiVector
        PMV shallowCopy(const std::vector<int> &index) {
            int n = index.size();
            APV x(n);
            for(int i=0;i<n;++i) {
                int j = index[i]; 
                x[i] = mvec_[j]; 
            }
            return ROL::makePtr<MultiVectorDefault<Real>>(x);
        } 

        // Make a const shallow copy specified vectors in the MultiVector
        const PMV shallowCopyConst(const std::vector<int> &index) const {
            int n = index.size();
            APV x(n);
            for(int i=0;i<n;++i) {
                int j = index[i]; 
                x[i] = mvec_[j]; 
            }
            return ROL::makePtr<MultiVectorDefault<Real>>(x);
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

           // Scale this by beta
            this->scale(beta);

            for(int i=0;i<B.numRows();++i) {
                for(int j=0;j<B.numCols();++j) {
                    mvec_[j]->axpy(alpha*B(i,j),*A.getVector(i));  
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

            ROL_TEST_FOR_EXCEPTION( static_cast<int>(alpha.size()) != numVectors_,
                std::invalid_argument,
                "Error: alpha must have the same length as the number of vectors.");  
 
            for(int i=0;i<numVectors_;++i) {
                mvec_[i]->scale(alpha[i]);  
            }
        } 

        // Set the MultiVector equal to another MultiVector
        void set(const MV &A) {

//            ROL_TEST_FOR_EXCEPTION( this->dimensionMismatch(A),
//                std::invalid_argument,
//                "Error: MultiVectors must have the same dimensions.");

            for(int i=0;i<numVectors_;++i) {
                mvec_[i]->set(*(A.getVector(i)));
            }
        }

        
        // Set some of the vectors in this MultiVector equal to corresponding 
        // vectors in another MultiVector
        void set(const MV &A, const std::vector<int> &index) {

//            ROL_TEST_FOR_EXCEPTION( this->dimensionMismatch(A),
//                std::invalid_argument,
//                "Error: MultiVectors must have the same dimensions.");

            int n = index.size();
            
            for(int i=0;i<n;++i) {
                int k = index[i];
                if(k<numVectors_ && i<A.getNumberOfVectors()) { 
                    mvec_[k]->set(*A.getVector(i));
                }
                
            }
        }

        // Compute \f$\alpha A^\top \text{this}\f$ 
        void innerProducts(const Real alpha,
                           const MV &A,
                           Teuchos::SerialDenseMatrix<int,Real> &B) const {

//            ROL_TEST_FOR_EXCEPTION( this->dimensionMismatch(A),
//                std::invalid_argument,
//                "Error: MultiVectors must have the same dimensions.");

            for(int i=0;i<A.getNumberOfVectors();++i) {
                for(int j=0;j<numVectors_;++j) {
                    B(i,j) = alpha*mvec_[j]->dot(*A.getVector(i));
                }  
            }
        }                  

        // Compute dot products of pairs of vectors
        void dots(const MV &A,
                  std::vector<Real> &b) const {

            ROL_TEST_FOR_EXCEPTION( this->dimensionMismatch(A),
                std::invalid_argument,
                "Error: MultiVectors must have the same dimensions.");

            for(int i=0;i<numVectors_;++i) {
                b[i] = mvec_[i]->dot(*A.getVector(i));
            }    
        } 

        // Compute the norm of each vector in the MultiVector
        void norms(std::vector<Real> &normvec) const {

            int min = numVectors_ < static_cast<int>(normvec.size()) ? numVectors_ : normvec.size();

            for(int i=0;i<min;++i) {
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

            ROL_TEST_FOR_EXCEPTION( i>=numVectors_, 
                std::invalid_argument,
                "Error: index out of bounds");

            return mvec_[i]; 
        } 

};
}


#endif
