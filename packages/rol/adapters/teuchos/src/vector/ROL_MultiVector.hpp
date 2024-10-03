// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** @ingroup la_group
 * \class ROL::MultiVector
 * \brief Provides a container and operations on multiple ROL vectors for
 *        use with other Trilinos packages which require multivectors
 */


#ifndef ROL_MULTIVECTOR_HPP
#define ROL_MULTIVECTOR_HPP

#include <vector>

#include "ROL_Vector.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace ROL {

template<class Real>
class MultiVector {

    typedef Vector<Real>          V;       // Single vector 
    typedef ROL::Ptr<V>           PV;      // Pointer to a vector
    typedef MultiVector<Real>     MV;      // Instance of the base class
    typedef ROL::Ptr<MV>          PMV;     // Pointer to an instance of the class 

    public:

        virtual ~MultiVector() {}

        /** \brief Make a new MultiVector of the same dimensions 

             @return A reference-counted pointer to the cloned MultiVector
        */
        virtual PMV clone() const = 0;

        /** \brief Make a new MultiVector of specified "width" 

             @return A reference-counted pointer to the cloned MultiVector
        */
        virtual PMV clone( const int numvecs ) const = 0;

        /** \brief Make a deep copy of this MultiVector

             @return A reference-counted pointer to a new MultiVector containing
                     deep copied values             
        */
        virtual PMV deepCopy() const  = 0;


        /** \brief Make a deep copy of this MultiVector
             @param[in] Array of indices of the vectors to copy
             @return A reference-counted pointer to a new MultiVector containing
                     deep copied values             
        */
        virtual PMV deepCopy(const std::vector<int> &index) const = 0; 


        /** \brief Make a shallow copy of this MultiVector
             @param[in] Array of indices of the vectors to copy
             @return A reference-counted pointer to a new MultiVector where the 
                     elements point to the data in *this
        */
        virtual PMV shallowCopy(const std::vector<int> &index) = 0;


        /** \brief Make a shallow copy of this MultiVector
             @param[in] Array of indices of the vectors to copy
             @return A reference-counted pointer to a new MultiVector where the 
                     elements point to the data in *this
        */
        virtual const PMV shallowCopyConst(const std::vector<int> &index) const = 0;


        /** \brief Get the number of elements of a vector in the MultiVector

            @return Number of elements in a vector
        */
        virtual ptrdiff_t getLength() const = 0;


        /** \brief Get the number of vectors in the MultiVector
 
             @return Number of vectors in the MultiVector
        */
        virtual int getNumberOfVectors() const = 0;

        /** \brief Generic BLAS level 3 matrix multiplication
            \f$\text{this}\leftarrow \alpha A B+\beta\text{*this}\f$   
            @param[in] alpha is a multiplicative factor of @b A
            @param[in] @b A is a MultiVector
            @param[in] @b B is a SerialDenseMatrix applied to A from the right
            @param[in] beta is a multiplicative factor of *this
        */
        virtual void gemm(const Real alpha,
                          const MV& A,
                          const Teuchos::SerialDenseMatrix<int,Real> &B,
                          const Real beta) = 0;

        /** \brief Perform the axpy operation columnwise on the MultiVector
            \f$ y_i\leftarrow y_i+\alpha x_i\f$ where \f$y\f$ is this MultiVector
            @param[in] alpha is the scaling factor
            @param[in] mv is the 
        */
        virtual void axpy(const Real alpha, const MV& x) = 0;


        /** \brief Scale the MultiVector by a single scalar alpha 
            \f$\text{this}\leftarrow\alpha\text{this}\f$
            @param[in] alpha is a scalar multiplicative factor
        */
        virtual void scale(const Real alpha) = 0;


        /** \brief Scale each vector in the MultiVector by a different alpha
            \f$\text{this}[i]\leftarrow\alpha[i]\text{*this}[i]\f$
            @param[in] alpha is a vector of multiplicative factors
        */
        virtual void scale(const std::vector<Real> &alpha) = 0;


        /** \brief Set the MultiVector equal to another MultiVector 
            @param[in] @b A is a MultiVector
        */
        virtual void set(const MV &A) = 0;

        
        /** \brief Set some of the vectors in this MultiVector equal to corresponding 
                   vectors in another MultiVector
            @param[in] @b A is a MultiVector
        */
        virtual void set(const MV &A, const std::vector<int> &index) = 0;


        /** \brief Compute \f$\alpha A^\top \text{*this}\f$ 
            @param[in] alpha is a multiplicative factor
            @param[in] @b A is a MultiVector
            @param[out] @b B is a SerialDenseMartrix in which to store the alpha-scaled
                        dot products of this MultiVector's vectors with those of @b A 
        */
        virtual void innerProducts(const Real alpha,
                                   const MV &A,
                                   Teuchos::SerialDenseMatrix<int,Real> &B) const = 0;


        /** \brief Compute dot products of pairs of vectors
            @param[in] @b A is a MultiVector
            @param[out] &b b is a vector containing the dot products between 
                           vectors contained in @b A and this MultiVector
        */
        virtual void dots(const MV &A,
                          std::vector<Real> &b) const = 0;


        /** \brief Compute the norm of each vector in the MultiVector
             @param[out] &b b is a vector containing the norms of the vectors 
              contained in this MultiVector
        */
        virtual void norms(std::vector<Real> &normvec) const = 0;


        /** \brief Zero each of the vectors in the MultiVector 
        */
        virtual void zero() = 0; 
        

        /** \brief Return a pointer to the ith vector 
            @param[in] i is the index of the desired vector
            @return A reference-counted pointer to the desired vector
        */
        virtual PV getVector(int i) const = 0;

};
}

#endif
