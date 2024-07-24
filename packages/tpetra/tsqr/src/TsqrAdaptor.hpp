// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef __TSQR_Trilinos_TsqrAdaptor_hpp
#define __TSQR_Trilinos_TsqrAdaptor_hpp

/// \file TsqrAdaptor.hpp
/// \brief Abstract interface between TSQR and multivector type

#include "Tsqr_ConfigDefs.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "TsqrTypeAdaptor.hpp"
#include "TsqrCommFactory.hpp"
#include "Tsqr_GlobalVerify.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include <stdexcept>
#include <sstream>

namespace TSQR {
  /// \namespace Trilinos
  /// \brief Interface between TSQR implementation and "the rest of Trilinos."
  ///
  /// "The rest of Trilinos" in this case means the specific linear
  /// algebra library: Epetra, Tpetra, or Thyra (which is more of an
  /// interface to other linear algebra libraries, but requires its
  /// own special TSQR adaptor).
  namespace Trilinos {
    /// \class TsqrAdaptor
    /// \brief Abstract interface between TSQR and multivector type
    ///
    /// Child classes of TsqrAdaptor tell TSQR how to compute a
    /// factorization of a specific Trilinos multivector class MV.
    /// Currently, \c Tpetra::MultiVector<S, LO, GO, NodeType> for any
    /// NodeType is supported.  At the moment, the latter will only be
    /// efficient if NodeType is not a GPU node.  Support for \c
    /// Epetra_MultiVector and Thyra multivectors may be added on
    /// request.
    ///
    /// TsqrAdaptor uses the appropriate specialization of
    /// TsqrTypeAdaptor to figure out which variant of TSQR to use on
    /// the given multivector type.  The caller is responsible for
    /// constructing the intranode and internode TSQR objects.
    ///
    /// \tparam S Scalar type
    /// \tparam LO Local ordinal type
    /// \tparam GO Global ordinal type: TSQR doesn't use it, but MV does.
    /// \tparam MV Multivector type
    ///
    /// Implementers who want to support TSQR with a new MultiVector
    /// (MV) type must create a subclass of that type, using e.g., \c
    /// TsqrTpetraAdaptor as a model.  They must then create a new \c
    /// TsqrTypeAdaptor specialization (with the appropriate
    /// typedefs), and a new \c TsqrCommFactory subclass.  The
    /// TsqrCommFactory subclass gets the underlying communicator
    /// object (e.g., \c Teuchos::Comm<int>) from a "prototype"
    /// multivector and wraps it into \c TSQR::MessengerBase<S> and \c
    /// TSQR::MessengerBase<LO> objects for TSQR.
    ///
    /// Implementers who wish to change which TSQR implementation is
    /// used for a particular MultiVector type (for which a
    /// TsqrAdaptor child class exists) should change the
    /// corresponding (possibly partial) specialization of \c
    /// TsqrTypeAdaptor.  Certainly the node_tsqr_type (and perhaps
    /// also the dist_tsqr_type) typedef(s) in the \c TsqrTypeAdaptor
    /// specialization must be changed.  If no corresponding \c
    /// TsqrFactory subclass exists for that combination of
    /// node_tsqr_type and dist_tsqr_type, a new \c TsqrFactory
    /// subclass may also have to be created, to tell us how to
    /// instantiate those node_tsqr_type and dist_tsqr_type objects.
    ///
    /// Implementers who wish to add a new TSQR factorization must
    /// create a new \c TsqrFactory subclass.
    template<class S, class LO, class GO, class MV>
    class TsqrAdaptor {
    public:
      typedef S   scalar_type;
      typedef LO  local_ordinal_type;
      typedef GO  global_ordinal_type;
      typedef MV  multivector_type;

      typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

      typedef TsqrTypeAdaptor<S, LO, GO, MV>        type_adaptor;
      typedef typename type_adaptor::factory_type   factory_type;

      typedef typename type_adaptor::node_tsqr_type node_tsqr_type;
      typedef typename type_adaptor::node_tsqr_ptr  node_tsqr_ptr;

      typedef typename type_adaptor::comm_type      comm_type;
      typedef typename type_adaptor::comm_ptr       comm_ptr;

      typedef typename type_adaptor::dist_tsqr_type dist_tsqr_type;
      typedef typename type_adaptor::dist_tsqr_ptr  dist_tsqr_ptr;

      typedef typename type_adaptor::tsqr_type      tsqr_type;
      typedef typename type_adaptor::tsqr_ptr       tsqr_ptr;

      typedef typename tsqr_type::FactorOutput      factor_output_type;
      typedef Teuchos::SerialDenseMatrix<LO, S>     dense_matrix_type;
      typedef Teuchos::RCP< MessengerBase<S> >      scalar_messenger_ptr;
      typedef Teuchos::RCP< MessengerBase<LO> >     ordinal_messenger_ptr;

      //! Virtual destructor ensures memory safety for derived classes.
      virtual ~TsqrAdaptor() = default;

      /// \brief Compute explicit "thin" QR factorization of A.
      ///
      /// \param A [in/out] On input, the multivector to factor.
      ///   Overwritten with nonuseful data on output.
      ///
      /// \param Q [out] On output, the explicit "thin" Q factor of A.
      ///
      /// \param R [out] On output, the square upper triangular R
      ///   factor in the QR factorization of A.
      ///
      /// \param contiguousCacheBlocks [in] Whether the data in A (and
      ///   Q) has been reorganized so that the elements of each cache
      ///   block are stored contiguously (i.e., via the output of
      ///   cacheBlock()).  The default is false, which means that
      ///   each process' row block of A (and Q) is stored as a matrix
      ///   in column-major order, with leading dimension >= the
      ///   number of rows in the row block.
      void
      factorExplicit (multivector_type& A,
                      multivector_type& Q,
                      dense_matrix_type& R,
                      const bool contiguousCacheBlocks = false)
      {
        factor_output_type output = factor (A, R, contiguousCacheBlocks);
        explicitQ (A, output, Q, contiguousCacheBlocks);
      }

      /// \brief Compute QR factorization of the multivector A.
      ///
      /// Compute the QR factorization in place of the multivector A.
      /// The Q factor is represented implicitly; part of that is
      /// stored in place in A (overwriting the input), and the other
      /// part is returned.  The returned object as well as the
      /// representation in A are both inputs of \c explicitQ().  The R
      /// factor is copied into R.
      ///
      /// \param A [in/out] On input, the multivector whose QR
      ///   factorization is to be computed.  Overwritten on output
      ///   with part of the implicit representation of the Q factor.
      ///
      /// \param R [out] On output, the R factor from the QR
      ///   factorization of A.  Represented as a square dense matrix
      ///   (not in packed form) with the same number of columns as A.
      ///   The lower triangle of R is overwritten with zeros on
      ///   output.
      ///
      /// \param contiguousCacheBlocks [in] Whether the data in A has
      ///   been reorganized so that the elements of each cache block
      ///   are stored contiguously (i.e., via the output of
      ///   cacheBlock()).  The default is false, which means that
      ///   each process' row block of A is stored as a matrix in
      ///   column-major order, with leading dimension >= the number
      ///   of rows in the row block.
      ///
      /// \return Additional information that, together with the A
      ///   output, encodes the implicitly represented Q factor from
      ///   the QR factorization of the A input.
      ///
      /// \note Virtual but implemented, because this default
      /// implementation is correct for all multivector_type types,
      /// but not necessarily efficient.  It should be efficient if
      /// fetchNonConstView(A) does not require copying the contents
      /// of A (e.g., from GPU memory to CPU memory).
      virtual factor_output_type
      factor (multivector_type& A,
              dense_matrix_type& R,
              const bool contiguousCacheBlocks = false)
      {
        local_ordinal_type nrowsLocal, ncols, LDA;
        fetchDims (A, nrowsLocal, ncols, LDA);
        // This is guaranteed to be _correct_ for any Node type, but
        // won't necessary be efficient.  The desired model is that
        // A_local requires no copying.
        Teuchos::ArrayRCP<scalar_type> A_local = fetchNonConstView (A);

        // Reshape R if necessary.  This operation zeros out all the
        // entries of R, which is what we want anyway.
        if (R.numRows() != ncols || R.numCols() != ncols) {
          if (0 != R.shape (ncols, ncols)) {
            throw std::runtime_error ("Failed to reshape matrix R");
          }
        }
        return pTsqr_->factor (nrowsLocal, ncols, A_local.data(), LDA,
                               R.values(), R.stride(), contiguousCacheBlocks);
      }

      /// \brief Compute the explicit Q factor.
      ///
      /// Compute the explicit (multivector) "thin" (same number of
      /// columns as the input) representation of the Q factor
      /// computed by factor(), using the implicit representation
      /// returned by factor().
      ///
      /// \param Q_in [in] Same as the "A" input of factor()
      /// \param factorOutput [in] Return value of factor()
      ///   corresponding to Q_in
      /// \param Q_out [out] Explicit "thin" representation of the Q
      ///   factor.  "Explicit" means as a regular matrix (in the same
      ///   multivector storage format as the "A" input of factor()).
      ///   "Thin" (terminology used by Golub and Van Loan) means that
      ///   the dimensions of Q_out are the same as the dimensions of
      ///   the "A" input of factor().
      /// \param contiguousCacheBlocks [in] See the epinonymous
      ///   argument of factor().  In this case, it applies to both
      ///   Q_in and Q_out, which must have the same data layout.
      ///
      /// \note Virtual but implemented, because this default
      /// implementation is correct for all multivector_type types,
      /// but not necessarily efficient.  It should be efficient if
      /// fetchNonConstView(Q_out) and fetchConstView(Q_in) do not
      /// require copying (e.g., from GPU memory to CPU memory) the
      /// contents of their respective multivector inputs.
      virtual void
      explicitQ (const multivector_type& Q_in,
                 const factor_output_type& factorOutput,
                 multivector_type& Q_out,
                 const bool contiguousCacheBlocks = false)
      {
        using Teuchos::ArrayRCP;

        local_ordinal_type nrowsLocal, ncols_in, LDQ_in;
        fetchDims (Q_in, nrowsLocal, ncols_in, LDQ_in);
        local_ordinal_type nrowsLocal_out, ncols_out, LDQ_out;
        fetchDims (Q_out, nrowsLocal_out, ncols_out, LDQ_out);

        if (nrowsLocal_out != nrowsLocal) {
          std::ostringstream os;
          os << "TSQR explicit Q: input Q factor\'s node-local part has a di"
            "fferent number of rows (" << nrowsLocal << ") than output Q fac"
            "tor\'s node-local part (" << nrowsLocal_out << ").";
          throw std::runtime_error (os.str());
        }
        ArrayRCP<const scalar_type> pQin = fetchConstView (Q_in);
        ArrayRCP<scalar_type> pQout = fetchNonConstView (Q_out);
        pTsqr_->explicit_Q (nrowsLocal,
                            ncols_in, pQin.data(), LDQ_in,
                            factorOutput,
                            ncols_out, pQout.data(), LDQ_out,
                            contiguousCacheBlocks);
      }

      /// \brief Rank-revealing decomposition.
      ///
      /// Using the R factor from factor() and the explicit Q factor
      /// from explicitQ(), compute the SVD of R (\f$R = U \Sigma
      /// V^*\f$).  R.  If R is full rank (with respect to the given
      /// relative tolerance), don't change Q or R.  Otherwise,
      /// compute \f$Q := Q \cdot U\f$ and \f$R := \Sigma V^*\f$ in
      /// place (the latter may be no longer upper triangular).
      ///
      /// \param Q [in/out] On input: the explicit Q factor computed
      ///   by explicitQ().  On output: unchanged if R has full
      ///   (numerical) rank, else \f$Q := Q \cdot U\f$, where \f$U\f$
      ///   is the ncols by ncols matrix of R's left singular vectors.
      ///
      /// \param R [in/out] On input: ncols by ncols upper triangular
      ///   matrix stored in column-major order.  On output: if input
      ///   has full (numerical) rank, R is unchanged on output.
      ///   Otherwise, if \f$R = U \Sigma V^*\f$ is the SVD of R, on
      ///   output R is overwritten with \f$\Sigma \cdot V^*\f$.  This
      ///   is also an ncols by ncols matrix, but may not necessarily
      ///   be upper triangular.
      ///
      /// \return Rank \f$r\f$ of R: \f$ 0 \leq r \leq ncols\f$.
      ///
      local_ordinal_type
      revealRank (multivector_type& Q,
                  dense_matrix_type& R,
                  const magnitude_type relativeTolerance,
                  const bool contiguousCacheBlocks = false) const
      {
        using Teuchos::ArrayRCP;

        local_ordinal_type nrowsLocal, ncols, ldqLocal;
        fetchDims (Q, nrowsLocal, ncols, ldqLocal);

        ArrayRCP< scalar_type > Q_ptr = fetchNonConstView (Q);
        return pTsqr_->reveal_rank (nrowsLocal, ncols,
                                    Q_ptr.data(), ldqLocal,
                                    R.values(), R.stride(),
                                    relativeTolerance,
                                    contiguousCacheBlocks);
      }

      /// \brief Cache-block A_in into A_out.
      ///
      /// Copy A_in into A_out, in a reorganized way that improves
      /// locality of cache blocks.
      ///
      /// \warning This may invalidate some invariants of A_out, such
      ///   as the mapping from index pair (i,j) to element of A_out.
      ///   Another way to say this is that the multivector object may
      ///   not be aware that its data has been reorganized underneath
      ///   it.
      virtual void
      cacheBlock (const multivector_type& A_in,
                  multivector_type& A_out)
      {
        using Teuchos::ArrayRCP;

        local_ordinal_type nrowsLocal, ncols, LDA_in;
        fetchDims (A_in, nrowsLocal, ncols, LDA_in);
        local_ordinal_type nrowsLocal_out, ncols_out, LDA_out;
        fetchDims (A_out, nrowsLocal_out, ncols_out, LDA_out);

        if (nrowsLocal_out != nrowsLocal) {
          std::ostringstream os;
          os << "TSQR cache block: the input matrix\'s node-local part has a"
            " different number of rows (" << nrowsLocal << ") than the outpu"
            "t matrix\'s node-local part (" << nrowsLocal_out << ").";
          throw std::runtime_error (os.str());
        }
        else if (ncols_out != ncols) {
          std::ostringstream os;
          os << "TSQR cache block: the input matrix\'s node-local part has a"
            " different number of columns (" << ncols << ") than the output "
            "matrix\'s node-local part (" << ncols_out << ").";
          throw std::runtime_error (os.str());
        }
        ArrayRCP<const scalar_type> pA_in = fetchConstView (A_in);
        ArrayRCP<scalar_type> pA_out = fetchNonConstView (A_out);
        pTsqr_->cache_block (nrowsLocal, ncols, pA_out.data(),
                             pA_in.data(), LDA_in);
      }

      /// \brief Un-cache-block A_in into A_out.
      ///
      /// Undo the transformation performed by \c cacheBlock(), by
      /// copying the contiguously cache blocked data in A_in into the
      /// conventionally stored A_out.
      virtual void
      unCacheBlock (const multivector_type& A_in,
                    multivector_type& A_out)
      {
        using Teuchos::ArrayRCP;

        local_ordinal_type nrowsLocal, ncols, LDA_in;
        fetchDims (A_in, nrowsLocal, ncols, LDA_in);
        local_ordinal_type nrowsLocal_out, ncols_out, LDA_out;
        fetchDims (A_out, nrowsLocal_out, ncols_out, LDA_out);

        if (nrowsLocal_out != nrowsLocal) {
          std::ostringstream os;
          os << "TSQR un-cache-block: the input matrix\'s node-local part ha"
            "s a different number of rows (" << nrowsLocal << ") than the ou"
            "tput matrix\'s node-local part (" << nrowsLocal_out << ").";
          throw std::runtime_error (os.str());
        }
        else if (ncols_out != ncols) {
          std::ostringstream os;
          os << "TSQR cache block: the input matrix\'s node-local part has a"
            " different number of columns (" << ncols << ") than the output "
            "matrix\'s node-local part (" << ncols_out << ").";
          throw std::runtime_error (os.str());
        }
        ArrayRCP<const scalar_type> pA_in = fetchConstView (A_in);
        ArrayRCP<scalar_type> pA_out = fetchNonConstView (A_out);
        pTsqr_->un_cache_block (nrowsLocal, ncols, pA_out.data(),
                                LDA_out, pA_in.data());
      }

      /// \brief Verify the result of the "thin" QR factorization \f$A = QR\f$.
      ///
      /// This method returns a list of three magnitudes:
      /// - \f$\| A - QR \|_F\f$
      /// - \f$\|I - Q^* Q\|_F\f$
      /// - \f$\|A\|_F\f$
      ///
      /// The notation $\f\| X \|\f$ denotes the Frobenius norm
      /// (square root of sum of squares) of a matrix \f$X\f$.
      /// Returning the Frobenius norm of \f$A\f$ allows you to scale
      /// or not scale the residual \f$\|A - QR\|\f$ as you prefer.
      virtual std::vector< magnitude_type >
      verify (const multivector_type& A,
              const multivector_type& Q,
              const Teuchos::SerialDenseMatrix< local_ordinal_type, scalar_type >& R)
      {
        using Teuchos::ArrayRCP;

        local_ordinal_type nrowsLocal_A, ncols_A, LDA;
        local_ordinal_type nrowsLocal_Q, ncols_Q, LDQ;
        fetchDims (A, nrowsLocal_A, ncols_A, LDA);
        fetchDims (Q, nrowsLocal_Q, ncols_Q, LDQ);
        if (nrowsLocal_A != nrowsLocal_Q)
          throw std::runtime_error ("A and Q must have same number of rows");
        else if (ncols_A != ncols_Q)
          throw std::runtime_error ("A and Q must have same number of columns");
        else if (ncols_A != R.numCols())
          throw std::runtime_error ("A and R must have same number of columns");
        else if (R.numRows() < R.numCols())
          throw std::runtime_error ("R must have no fewer rows than columns");

        // Const views suffice for verification
        ArrayRCP<const scalar_type> A_ptr = fetchConstView (A);
        ArrayRCP<const scalar_type> Q_ptr = fetchConstView (Q);
        return global_verify (nrowsLocal_A, ncols_A, A_ptr.data(), LDA,
                              Q_ptr.data(), LDQ, R.values(), R.stride(),
                              pScalarMessenger_.get());
      }

    protected:
      /// \brief A "nonconstructor constructor."
      ///
      /// This method initializes the adaptor, as a constructor would
      /// normally do.  However, we have to make this method not a
      /// constructor, since you're not supposed to call the
      /// constructor of a pure virtual class.  ("Call the constructor
      /// of" is a synonym for "instantiate," and you can't
      /// instantiate an instance of a pure virtual class.)
      void
      init (const multivector_type& mv,
            const Teuchos::RCP<Teuchos::ParameterList>& plist)
      {
        // This is done in a multivector type - dependent way.
        fetchMessengers (mv, pScalarMessenger_, pOrdinalMessenger_);

        factory_type factory;
        // plist and pScalarMessenger_ are inputs.  Construct *pTsqr_.
        factory.makeTsqr (plist, pScalarMessenger_, pTsqr_);
      }

    private:
      /// \brief Return dimensions of a multivector object.
      ///
      /// For a given multivector A, return the number of rows stored
      /// locally on this process, the number of columns (multivectors
      /// are stored in a block row layout, so all columns of this
      /// process' row block are stored on this process), and the
      /// leading dimension of this process' row block (>= # rows on
      /// this process).
      ///
      /// \param A [in] The multivector object
      /// \param nrowsLocal [out] Number of rows of A stored locally
      ///   on this process
      /// \param ncols [out] Number of columns of A
      /// \param LDA [out] Leading dimension of this process' row
      ///   block of A
      virtual void
      fetchDims (const multivector_type& A,
                 local_ordinal_type& nrowsLocal,
                 local_ordinal_type& ncols,
                 local_ordinal_type& LDA) const = 0;

      /// \brief Get nonconst pointer to the node-local data in A.
      ///
      /// \note Child classes should implement this in such a way as
      ///   to make the above public methods always correct (though
      ///   not necessarily efficient) for all multivector types.  (It
      ///   may not be efficient if the \c Teuchos::ArrayRCP copies
      ///   between different memory spaces.)
      virtual Teuchos::ArrayRCP<scalar_type>
      fetchNonConstView (multivector_type& A) const = 0;

      /// \brief Get const pointer to the node-local data in A.
      ///
      /// \note Child classes should implement this in such a way as
      ///   to make the above public methods always correct (though
      ///   not necessarily efficient) for all multivector types.  (It
      ///   may not be efficient if the \c Teuchos::ArrayRCP copies
      ///   between different memory spaces.)
      virtual Teuchos::ArrayRCP<const scalar_type>
      fetchConstView (const multivector_type& A) const = 0;

      /// \brief Get "messenger" objects associated with the given multivector.
      ///
      /// Trilinos multivectors typically store a (distributed-memory)
      /// communicator of some kind.  For example, \c
      /// Epetra_MultiVector instances store an \c Epetra_Comm, and \c
      /// Tpetra::MultiVector instances store a \c Teuchos::Comm<int>.
      /// This method wraps the communicator in two "messenger"
      /// objects, one for communicating scalars and the other for
      /// communicating ordinals.
      ///
      /// \note The messenger objects may or may not be valid once the
      ///   given multivector falls out of scope, depending on the
      ///   multivector type.
      virtual void
      fetchMessengers (const multivector_type& mv,
                       scalar_messenger_ptr& pScalarMessenger,
                       ordinal_messenger_ptr& pOrdinalMessenger) const = 0;

      /// Object that knows how to communicate scalar_type objects.
      scalar_messenger_ptr pScalarMessenger_;

      /// Object that knows how to communicate local_ordinal_type
      /// objects.
      ordinal_messenger_ptr pOrdinalMessenger_;

      /// The \c Tsqr object that implements the Tall Skinny QR (TSQR)
      /// factorization.
      tsqr_ptr pTsqr_;
    };
  } // namespace Trilinos
} // namespace TSQR

#endif // __TSQR_Trilinos_TsqrAdaptor_hpp
