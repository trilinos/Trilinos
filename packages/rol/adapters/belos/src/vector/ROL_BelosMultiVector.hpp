// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
 * \class ROL::Belos::MultiVector
 * \brief Specializes the Belos::MultiVecTraits to use ROL::MultiVector
 */


#ifndef ROL_BELOS_MULTIVECTOR_HPP
#define ROL_BELOS_MULTIVECTOR_HPP

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosMultiVecTraits.hpp"
#include "BelosTypes.hpp"

#include "ROL_MultiVector.hpp"

#ifdef HAVE_BELOS_TSQR
namespace { // (anonymous)

  /// \brief Stub TSQR adapter for the ROL::MultiVector specialization of
  ///   Belos::MultiVecTraits.
  ///
  /// This is a stub.  All the methods throw.  The point is to make
  /// everything compile if TSQR is enabled.
  template<class Scalar>
  class RolStubTsqrAdapter {
  public:
    RolStubTsqrAdapter (const Teuchos::RCP<Teuchos::ParameterList>&)
    {
      ROL_TEST_FOR_EXCEPTION(false, std::logic_error, "Not implemented");
    }

    RolStubTsqrAdapter ()
    {
      ROL_TEST_FOR_EXCEPTION(false, std::logic_error, "Not implemented");
    }

    Teuchos::RCP<const Teuchos::ParameterList>
    getValidParameters () const
    {
      return Teuchos::rcp(new Teuchos::ParameterList>());
    }

    void
    setParameterList (const Teuchos::RCP<Teuchos::ParameterList>&)
    {}

    void
    factorExplicit (ROL::MultiVector<Scalar>& /* A */,
                    ROL::MultiVector<Scalar>& /* Q */,
                    Teuchos::SerialDenseMatrix<int, Scalar>& /* R */,
                    const bool forceNonnegativeDiagonal = false)
    {
      (void) forceNonnegativeDiagonal; // forestall "unused" warning
      ROL_TEST_FOR_EXCEPTION(false, std::logic_error, "Not implemented");
    }

    int
    revealRank (ROL::MultiVector<Scalar>& /* Q */,
                Teuchos::SerialDenseMatrix<int, Scalar>& /* R */,
                const typename ROL::ScalarTraits<Scalar>::magnitudeType& /* tol */)
    {
      return 0;
    }
  };

} // namespace (anonymous)
#endif // HAVE_BELOS_TSQR

namespace Belos {

    template<class Scalar>
    class MultiVecTraits<Scalar, ROL::MultiVector<Scalar> > {

        typedef ROL::MultiVector<Scalar>  MV;      // ROL::MultiVector object
        typedef Teuchos::RCP<MV>          PMV;     // Pointer to ROL::MultiVector object

        typedef Teuchos::SerialDenseMatrix<int,Scalar> Matrix;

        typedef typename ROL::ScalarTraits<Scalar>::magnitudeType magnitudeType;

        public:

            /** \brief Make a new MultiVector containing a specified number of vectors

                @return A reference-counted pointer to the cloned MultiVector
            */
            static PMV Clone(const MV& mv, const int numVecs) {
                return mv.clone(numVecs);
            }


            /** \brief Make a new MultiVector of the same size and containing a deep copy
                of the original MultiVector's data

                @return A reference-counted pointer to the cloned MultiVector
            */
            static PMV CloneCopy(const MV& mv) {
                return mv.deepCopy();
            }

            /** \brief Make a new MultiVector containing a deep copy of specified columns
                of the original MultiVector's data

                @return A reference-counted pointer to the cloned MultiVector
            */
            static PMV CloneCopy(const MV& mv, const std::vector<int>& index) {
                return mv.deepCopy(index);
            }

            /** \brief Make a new MultiVector containing a deep copy of specified
                contiguous block of columns of the original MultiVector's data

                @return A reference-counted pointer to the cloned MultiVector
            */
            static PMV CloneCopy(const MV& mv, const Teuchos::Range1D& index) {
                int n = index.size();
                int l = index.lbound();

                std::vector<int> stdex(n);
                for( int i=0; i<n; ++i ) {
                    stdex[i] = l+i;
                }

                return mv.deepCopy(stdex);
            }

            /** \brief Make a new MultiVector containing pointers to specified
                columns of the original MultiVector's data

                @return A reference-counted pointer to the cloned MultiVector
            */
            static PMV CloneViewNonConst(MV& mv, const std::vector<int>& index) {
                return mv.shallowCopy(index);
            }


            /** \brief Make a new MultiVector containing pointers to a contiguous
                block of columns of the original MultiVector's data

                @return A reference-counted pointer to the cloned MultiVector
            */
            static PMV CloneViewNonConst(MV& mv, const Teuchos::Range1D& index) {
                int n = index.size();
                int l = index.lbound();

                std::vector<int> stdex(n);
                for( int i=0; i<n; ++i ) {
                    stdex[i] = l+i;
                }

                return mv.shallowCopy(stdex);
            }

            /** \brief Make a new const MultiVector containing pointers to specified
                columns of the original MultiVector's data

                @return A reference-counted pointer to the cloned MultiVector
            */
            static PMV CloneView(const MV& mv, const std::vector<int>& index) {
                return mv.shallowCopyConst(index);
            }

            /** \brief Make a new const MultiVector containing pointers to a contiguous
                block of columns of the original MultiVector's data

                @return A reference-counted pointer to the cloned MultiVector
            */
            static PMV CloneView(const MV& mv, const Teuchos::Range1D& index) {
                int n = index.size();
                int l = index.lbound();

                std::vector<int> stdex(n);
                for( int i=0; i<n; ++i ) {
                    stdex[i] = l+i;
                }

                return mv.shallowCopyConst(stdex);
            }


            /** \brief Get length of the columns of the MultiVector

                @return Vector length
            */
            static ptrdiff_t GetGlobalLength(const MV& mv) {
                return mv.getLength();
            }

            /** \brief Return the number of vectors in the MultiVector
            */
            static int GetNumberVecs(const MV& mv) {
                return mv.getNumberOfVectors();
            }

            /// \brief Stride is undefined for ROL::Vector
            static bool HasConstantStride(const MV& mv) {
                return false;
            }

            /** \brief Replace mv with \f$\alpha AB+\beta mv\f$
            */
            static void MvTimesMatAddMv(Scalar alpha, const MV& A,
                                        const Matrix& B, Scalar beta, MV& mv) {
                mv.gemm(alpha,A,B,beta);
            }


            static void MvAddMv(Scalar alpha, const MV& A,
                                Scalar beta,  const MV& B, MV& mv) {

                mv.set(A);
                mv.scale(alpha);
                mv.axpy(beta,B);
            }

            /** \brief Scale every element in the MultiVector by alpha
            */
            static void MvScale(MV& mv, Scalar alpha) {
                mv.scale(alpha);
            }

            /** \brief Scale the columns of the MultiVector by the values in alpha
            */
            static void MvScale(MV& mv, const std::vector<Scalar>& alphas) {
                mv.scale(alphas);
            }

            /** \brief Compute \f$\alpha A^\top B$\f$ and store the result in
                       the matrix\f$C\f$
            */
            static void MvTransMv(const Scalar alpha, const MV& A,
                                  const MV& B, Matrix& C) {
                B.innerProducts(alpha,A,C);
            }

            /** \brief Compute the columnwise dot products of A and B and store the
                       result in dots
            */
            static void MvDot(const MV& A, const MV& B, std::vector<Scalar> &dots) {
                A.dots(B,dots);
            }

            /** \brief Compute the columnwise norms of mv and store the
                       result in dots
            */
            static void MvNorm(const MV& mv,
                               std::vector<magnitudeType>& normvec,
                               NormType type=TwoNorm) {

                ROL_TEST_FOR_EXCEPTION( (type !=TwoNorm) , std::invalid_argument,
                    "Belos::MultiVecTraits<Scalar,ROL::MultiVector<Scalar> >::MvNorm()\n"
                    "ROL::MultiVector supports only Euclidean norm");

                mv.norms(normvec);

            }

            /** \brief Copy the specified columns of A into mv
            */
            static void SetBlock(const MV& A, const std::vector<int>& index, MV& mv) {
                mv.set(A,index);
            }

            /** \brief Copy the contiguous block of columns of A into mv
            */
            static void SetBlock(const MV& A, Teuchos::Range1D& index, MV& mv) {
                int n = index.size();
                int l = index.lbound();

                std::vector<int> stdex(n);
                for( int i=0; i<n; ++i ) {
                    stdex[i] = l+i;
                }
                mv.set(A,stdex);
            }

            /** \brief Copy the columns of A into mv
            */
            static void Assign(const MV& A, MV& mv) {
                mv.set(A);
            }

            /** \brief Randomize the vectors in mv (not implemented)
            */
            static void MvRandom(MV& mv) {
                ROL_TEST_FOR_EXCEPTION(true,std::logic_error,
                    "Belos::MultiVecTraits<Scalar,ROL::MultiVector<Scalar> >::Random()\n"
                    "Random initialization not implemented for ROL::MultiVector");
            }


            static void MvInit(MV& mv,
                               const Scalar alpha = ROL::ScalarTraits<Scalar>::zero()) {
                ROL_TEST_FOR_EXCEPTION(alpha != 0,std::invalid_argument,
                    "Belos::MultiVecTraits<Scalar,ROL::MultiVector<Scalar> >::MvInit()\n");
                mv.zero();
            }


            static void MvPrint(const MV& mv, std::ostream& os) {
                ROL_TEST_FOR_EXCEPTION(true,std::logic_error,
                    "Belos::MultiVecTraits<Scalar,ROL::MultiVector<Scalar> >::MvPrint()\n"
                    "Print not implemented for ROL::MultiVector");
            }

#ifdef HAVE_BELOS_TSQR
      /// \typedef tsqr_adaptor_type
      /// \brief TsqrAdaptor specialization for ROL::MultiVector
      typedef RolStubTsqrAdapter<Scalar> tsqr_adaptor_type;
#endif // HAVE_BELOS_TSQR

    };
}


#endif
