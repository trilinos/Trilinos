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


/** 
    \class Belos::MultiVecTraits<Scalar,ROL::Vector<Scalar>>
    \brief Provides interface for using ROL::Vector with Belos solvers
           (excluding block solvers).
                
    \author Created by Greg von Winckel
*/


#ifndef ROL_BELOS_VECTOR_HPP
#define ROL_BELOS_VECTOR_HPP

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosMultiVecTraits.hpp"
#include "BelosTypes.hpp"

#include "ROL_Vector.hpp"

namespace Belos {

    using ROL::Vector;
    using Teuchos::RCP;
    using Teuchos::rcp; 

    template<class Scalar>
    class MultiVecTraits<Scalar, ROL::Vector<Scalar>> {
    public:

        /// Not a multivector, but use the same symbol as the Belos examples
        typedef ROL::Vector<Scalar> MV; 

        /// \brief Create a new ROL Vector of the same dimension. The
        /// default constructor will be called on the elements of the 
        /// encapsulated data type 
        static Teuchos::RCP<MV> 
        Clone(const MV& mv, const int numVecs) {
            TEUCHOS_TEST_FOR_EXCEPTION(numVecs != 1, std::invalid_argument,
                "Belos::MultiVecTraits<Scalar,ROL::Vector<Scalar>>::Clone()\n" 
                "ROL Vectors have only one column")
            
            return mv.clone();
        }

        /// \brief Make a deep copy
        static Teuchos::RCP<MV> 
        CloneCopy(const MV& mv) {

            auto mv_clone = mv.clone();
            mv_clone->set(mv);

            return mv_clone;
        }

        /// \brief Make a deep copy
        static Teuchos::RCP<MV>
        CloneCopy(const MV& mv, const std::vector<int> &index) {
            TEUCHOS_TEST_FOR_EXCEPTION( (index.size() != 1) | (index[0] !=0), 
                std::invalid_argument,
                "Belos::MultiVecTraits<Scalar,ROL::Vector<Scalar>>::CloneCopy()\n"
                "ROL Vectors have only one column");

            auto mv_clone = mv.clone();
            mv_clone->set(mv);

            return mv_clone;
        }

        /// \brief Make a deep copy
        static Teuchos::RCP<MV>
        CloneCopy(const MV& mv, const Teuchos::Range1D& index) {
            TEUCHOS_TEST_FOR_EXCEPTION( (index.size() != 1) | (index.lbound() !=0), std::invalid_argument,
                "Belos::MultiVecTraits<Scalar,ROL::StdVector<Scalar>>::CloneCopy()\n"
                "ROL::StdVector has only one column");      

            auto mv_clone = mv.clone();
            mv_clone->set(mv);

            return mv_clone();
        }
        
        /// \brief Make a shallow copy
        static Teuchos::RCP<MV>
        CloneViewNonConst(MV& mv, const std::vector<int>& index) { 
            TEUCHOS_TEST_FOR_EXCEPTION( (index.size() != 1) | (index[0] !=0), 
                std::invalid_argument,
                "Belos::MultiVecTraits<Scalar,ROL::Vector<Scalar>>::CloneCopy()\n"
                "ROL Vectors have only one column");

              // Deep copy until ROL interface is updated
              auto mv_clone = mv.clone();
              mv_clone->set(mv); 
              return mv_clone;

//            return mv.cloneViewNonConst();
        }
 
        /// \brief Make a shallow copy
        static Teuchos::RCP<MV>
        CloneViewNonConst(MV& mv, const Teuchos::Range1D& index) {
            TEUCHOS_TEST_FOR_EXCEPTION( (index.size() != 1) | (index.lbound() !=0), std::invalid_argument,
                "Belos::MultiVecTraits<Scalar,ROL::StdVector<Scalar>>::CloneViewNonConst()\n"
                "ROL::StdVector has only one column");      
              // Deep copy until ROL interface is updated
              auto mv_clone = mv.clone();
              mv_clone->set(mv); 
              return mv_clone;

//            return mv.cloneViewNonConst();  
        }

        /// \brief Make a shallow copy
        static Teuchos::RCP<const MV>
        CloneView(const MV& mv, const std::vector<int>& index) {
            TEUCHOS_TEST_FOR_EXCEPTION( (index.size() != 1) | (index[0] !=0), 
                std::invalid_argument,
                "Belos::MultiVecTraits<Scalar,ROL::Vector<Scalar>>::CloneCopy()\n"
                "ROL Vectors have only one column");
 
            // Deep copy until ROL interface is updated
            auto mv_clone = mv.clone();
            mv_clone->set(mv); 
            return mv_clone;

 //           return mv.cloneView();
         }

        /// \brief Make a shallow copy
        static Teuchos::RCP<const MV>
        CloneView(const MV& mv, const Teuchos::Range1D& index) {
            TEUCHOS_TEST_FOR_EXCEPTION( (index.size() != 1) | (index.lbound() !=0), std::invalid_argument,
                "Belos::MultiVecTraits<Scalar,ROL::StdVector<Scalar>>::CloneCopy()\n"
                "ROL::StdVector has only one column");      
            // Deep copy until ROL interface is updated
            auto mv_clone = mv.clone();
            mv_clone->set(mv); 
            return mv_clone;

   //         return mv.cloneView();
        }

        /// \brief Return number of elements in vector
        static ptrdiff_t 
        GetGlobalLength(const MV& mv) {
            return static_cast<ptrdiff_t>(mv.dimension());
        }   

        /// \brief ROL Vectors have a single column
        static int 
        GetNumberVecs(const MV& mv) {
            return 1;
        }
                  
        /// \brief Stride is undefined for ROL::Vector
        static bool 
        HasConstantStride(const MV& mv) {
            return false;
        }

        /// \brief Replace mv with \f$\alpha AB+\beta mv\f$
        static void
        MvTimesMatAddMv(Scalar alpha, const MV& A, 
                        const Teuchos::SerialDenseMatrix<int, Scalar>& B,
                        Scalar beta, MV& mv) {
            TEUCHOS_TEST_FOR_EXCEPTION((B.numRows()!=1)|(B.numCols()!=1),
                std::invalid_argument,
                "Belos::MultiVecTraits<Scalar,ROL::Vector<Scalar>>::MvTimesMatAddMv()\n"
                "ROL::Vector has only column");  
                mv.scale(beta);
                mv.axpy(alpha*B[0][0],A);  
        }

        /// \brief \f$ mv\leftarrow mv + \alpha A + \beta B \f$
        static void
        MvAddMv(Scalar alpha, const MV& A, Scalar beta, const MV& B, MV& mv) {
            mv.axpy(alpha, A);
            mv.axpy(beta, B);
        }
         
        /// \brief Scale every element in the vector by alpha
        static void MvScale(MV& mv, Scalar alpha) {
            mv.scale(alpha); 
        }

        /// \brief Scale every element in the vector by alpha
        static void MvScale(MV& mv, const std::vector<Scalar>& alphas) {
            TEUCHOS_TEST_FOR_EXCEPTION( alphas.size()>1, std::invalid_argument,
                "Belos::MultiVecTraits<Scalar,ROL::Vector<Scalar>>::MvScale()\n"
                "ROL::Vector has only one column");
                mv.scale(alphas[0]); 
        }


        /// \brief Action of a matrix is undefined for ROL::Vector
        static void
        MvTransMv(const Scalar alpha, const MV& A,
                   const MV& B, Teuchos::SerialDenseMatrix<int,Scalar>& C) {
            TEUCHOS_TEST_FOR_EXCEPTION((C.numRows()!=1)|(C.numCols()!=1),
                std::invalid_argument,
                "Belos::MultiVecTraits<Scalar,ROL::Vector<Scalar>>::MvTransMv()\n"
                "Inner product of two ROL is a scalar");
            C[0][0] = alpha*A.dot(B);   
        }

        /// \brief Vector dot product
        static void
        MvDot(const MV& A, const MV& B, std::vector<Scalar> &dots) {
            TEUCHOS_TEST_FOR_EXCEPTION( dots.size()>1, std::invalid_argument,
                "Belos::MultiVecTraits<Scalar,ROL::Vector<Scalar>>::MvDot()\n"
                "ROL::Vector has only one column");
            
            dots[0] = A.dot(B); 
        }

        /// \brief Vector norm 
        static void
        MvNorm(const MV& mv, 
                std::vector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& normvec,
                NormType type=TwoNorm) {
             TEUCHOS_TEST_FOR_EXCEPTION( normvec.size()>1, std::invalid_argument,
                "Belos::MultiVecTraits<Scalar,ROL::Vector<Scalar>>::MvNorm()\n"
                "ROL::Vector has only one column");
             TEUCHOS_TEST_FOR_EXCEPTION( (type !=TwoNorm) , std::invalid_argument,
                "Belos::MultiVecTraits<Scalar,ROL::Vector<Scalar>>::MvNorm()\n"
                "ROL::Vector supports only Euclidean norm");
             normvec[0] = mv.norm(); 
        }

        /// \brief Deep copy into mv
        static void
        SetBlock(const MV& A, const std::vector<int>& index, MV& mv) {
            TEUCHOS_TEST_FOR_EXCEPTION( (index.size() != 1) | (index[0] !=0), 
                std::invalid_argument,
                "Belos::MultiVecTraits<Scalar,ROL::Vector<Scalar>>::CloneCopy()\n"
                "ROL Vectors have only one column");
            mv.set(A);
        }

        /// \brief Deep copy into mv
        static void
        SetBlock(const MV& A, const Teuchos::Range1D& index, MV& mv) {
            TEUCHOS_TEST_FOR_EXCEPTION( (index.size() != 1) | (index.lbound() !=0), 
                std::invalid_argument,
                "Belos::MultiVecTraits<Scalar,ROL::Vector<Scalar>>::CloneCopy()\n"
                "ROL Vectors have only one column");
            mv.set(A);
        }

        /// \brief Deep copy into mv
        static void
        Assign(const MV& A, MV& mv) {
            mv.set(A);
        }

        /// \brief Not implemented
        static void 
        MvRandom(MV& mv) {
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                "Belos::MultiVecTraits<Scalar,ROL::Vector<Scalar>>::Random()\n"
                "Random initialization not implemented for ROL::Vector");
        }

        /// \brief Initialize vector to all alpha           
        static void
        MvInit(MV& mv, const Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero ()) {
            TEUCHOS_TEST_FOR_EXCEPTION(alpha != 0,std::invalid_argument,
                 "Belos::MultiVecTraits<Scalar,ROL::Vector<Scalar>>::MvInit()\n");      
            mv.zero();
        }

        /// \brief Not implemented
        static void 
        MvPrint(const MV& mv, std::ostream& os) {
            TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                "Belos::MultiVecTraits<Scalar,ROL::Vector<Scalar>>::MvPrint()\n"
                "Print not implemented for ROL::Vector");
        }

    }; // END MultiVecTraits

}

#endif

