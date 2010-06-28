#ifndef __TSQR_Trilinos_TsqrAdaptor_hpp
#define __TSQR_Trilinos_TsqrAdaptor_hpp

#include "AnasaziConfigDefs.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "TsqrTrilinosMessenger.hpp"
#include "TsqrTypeAdaptor.hpp"

#include <TSQR/Tsqr.hpp>
#include <string>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Trilinos {

    /// \class TsqrAdaptor
    /// \brief Interface between a Trilinos multivector class and TSQR
    ///
    /// TsqrAdaptor tells TSQR how to compute a factorization of a
    /// Trilinos multivector class MV.  Currently, Epetra_MultiVector
    /// and Tpetra::MultiVector< S, LO, GO, NodeType > for any
    /// NodeType are supported.  At the moment, the latter will only
    /// be efficient if NodeType is not a GPU node.  TsqrAdaptor uses
    /// TsqrTypeAdaptor to figure out which variant of TSQR to use on
    /// the given multivector type.  For example, with
    /// Tpetra::MultiVector< S, LO, GO, NodeType >, if NodeType is
    /// Kokkos::TBBNode, the TBB-parallel intranode variant of TSQR
    /// will be used.  The caller is responsible for constructing
    /// the intranode and internode TSQR objects.
    ///
    /// S: scalar type
    /// LO: local ordinal type
    /// GO: global ordinal type: TSQR doesn't use it, but MV does.
    /// MV: multivector type
    ///
    /// \note This is the interface to TSQR that Trilinos sees.
    ///   Implementers who want to add a new MultiVector (MV) type
    ///   must create three different (partial) instantiations (in
    ///   this order): a TsqrTypeAdaptor, a TsqrFactory, and a
    ///   TsqrAdaptor.  Implementers who wish to change which TSQR
    ///   implementation is used for a particular MultiVector type
    ///   (for which a (partial) TsqrAdaptor instantiation exists)
    ///   should change the corresponding (partial) instantiation of
    ///   TsqrTypeAdaptor (which maps the MultiVector type to the TSQR
    ///   implementation type).
    template< class S, class LO, class GO, class MV >
    class TsqrAdaptor {
    public:
      typedef S   scalar_type;
      typedef LO  local_ordinal_type;
      typedef GO  global_ordinal_type;
      typedef MV  multivector_type;

      typedef TsqrTypeAdaptor< S, LO, GO, MV >::node_tsqr_type node_tsqr_type;
      typedef TsqrTypeAdaptor< S, LO, GO, MV >::tsqr_type      type_type;

      typedef Teuchos::RCP< node_tsqr_type >             node_tsqr_ptr;
      typedef Teuchos::RCP< tsqr_type >                  tsqr_ptr;
      typedef Teuchos::RCP< TrilinosMessenger< S > >     messenger_ptr;
      typedef typename tsqr_type::FactorOutput           factor_output_type;
      typedef Teuchos::SerialDenseMatrix< LO, S >        dense_matrix_type;

      /// \brief Constructor
      ///
      /// TsqrAdaptor constructor.  The reference-counted pointer
      /// inputs should be to already constructed objects.
      /// pComm.get() and *pNodeTsqr should have been used to
      /// initialize *pTsqr.  (We do this because different TSQR
      /// intranode and internode types have different constructors.)
      ///
      /// \note The caller is responsible for ensuring that the types
      /// of node_tsqr_type and tsqr_type are correct for the given
      /// multivector_type.  TsqrTypeAdaptor may be useful for that
      /// purpose.
      ///
      /// \param pComm [in] Shared pointer to a wrapper around
      ///   Teuchos::Comm.  We need to keep it in this class because
      ///   the already constructed Tsqr object pTsqr keeps a raw
      ///   pointer to it internally; we don't want the smart pointer
      ///   to fall out of scope.
      ///
      /// \param pNodeTsqr [in] Shared pointer to the already
      ///   constructed node TSQR object.  It was used to construct
      ///   *pTsqr, but we keep it here because pTsqr used a raw
      ///   pointer to node_tsqr_type; it doesn't have a
      ///   reference-counted pointer.
      ///
      /// \param pTsqr [in] Shared pointer to the already constructed
      ///   object that performs the QR factorization.  pComm.get()
      ///   and *pNodeTsqr should have been used to construct it.
      TsqrAdaptor (const Teuchos::ParameterList& plist)
      {
	typedef TsqrFactory< LO, S, node_tsqr_type, tsqr_type > factory_type;
	factory_type::makeTsqr (plist, pComm_, pNodeTsqr_, pTsqr_);
      }

      /// \brief Compute QR factorization of the multivector A
      ///
      /// Compute the QR factorization in place of the multivector A.
      /// The Q factor is represented implicitly; part of that is
      /// stored in place in A (overwriting the input), and the other
      /// part is returned.  The returned object as well as the
      /// representation in A are both inputs of explicitQ().  The R
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
      factor_output_type
      factor (multivector_type& A, 
	      dense_matrix_type& R,
	      const bool contiguousCacheBlocks = false);

      /// \brief Compute the explicit Q factor
      ///
      /// Compute the explicit (multivector) "thin" (same number of
      /// columns as the input) representation of the Q factor
      /// computed by factor(), using the implicit representation
      /// returned by factor().
      ///
      /// \param Q_in [in] Same as the "A" input of factor()
      /// \param factorOutput [in] Return value of factor() 
      ///  corresponding to Q_in
      /// \param Q_out [out] Explicit "thin" representation of the Q
      ///   factor.  "Explicit" means as a regular matrix (in the same
      ///   multivector storage format as the "A" input of factor()).
      ///   "Thin" (terminology used by Golub and Van Loan) means that
      ///   the dimensions of Q_out are the same as the dimensions of
      ///   the "A" input of factor().
      /// \param contiguousCacheBlocks [in] See the epinonymous
      ///   argument of factor().  In this case, it applies to both
      ///   Q_in and Q_out, which must have the same data layout.
      void 
      explicitQ (const multivector_type& Q_in, 
		 const factor_output_type& factorOutput,
		 multivector_type& Q_out, 
		 const bool contiguousCacheBlocks = false);

      /// \brief Cache-block A_in into A_out
      ///
      /// Copy A_in into A_out, in a reorganized way that improves
      /// locality of cache blocks.
      ///
      /// \warning This may invalidate some invariants of A_out, such
      ///   as the mapping from index pair (i,j) to element of A_out.
      ///   Another way to say this is that the multivector object may
      ///   not be aware that its data has been reorganized underneath
      ///   it.
      void 
      cacheBlock (const multivector_type& A_in, 
		  multivector_type& A_out) const;

      /// \brief Un-cache-block A_in into A_out
      ///
      /// Undo the transformation performed by cacheBlock(), by
      /// copying the contiguously cache blocked data in A_in into the
      /// conventionally stored A_out.
      void 
      unCacheBlock (const multivector_type& A_in, 
		    multivector_type& A_out) const;

    private:
      /// \brief Return dimensions of a multivector object
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
      ///
      /// \note (mfh 15 June 2010) The reason why we can't just use
      /// fetch_MV_dims() as the adaptor class for TSQR, is because of
      /// Tpetra::MultiVector's use of smart pointers.
      /// fetch_MV_dims() can't just return raw pointers to the
      /// internal storage of Tpetra::MultiVector, since that would
      /// require letting the raw pointers escape the scope of their
      /// parent smart pointer objects.  If we return the smart
      /// pointer, then the adaptor class depends on the type of smart
      /// pointer (or raw pointer, in the case of Epetra_MultiVector,
      /// whose Values() method returns just that).  If I could use
      /// boost::shared_ptr in TSQR then I would just use it
      /// throughout, which would remove the need for the TsqrAdaptor
      /// class; however, Boost is only an optional dependency in
      /// Trilinos (thanks to one particular compiler vendor *ahem*).
      /// I don't want TSQR to depend on Teuchos::RCP because I want
      /// to make TSQR fully portable.  So I'm stuck writing adaptors.
      void 
      fetch_MV_dims (const multivector_type& A, 
		     local_ordinal_type& nrowsLocal, 
		     local_ordinal_type& ncols, 
		     local_ordinal_type& LDA);

      /// Shared pointer to a wrapper around Teuchos::Comm.  We need
      /// to keep it in this class because *pTsqr keeps a raw pointer
      /// to it internally; we don't want *pComm_ to fall out of
      /// scope.
      messenger_ptr pComm_;

      /// Shared pointer to the "Node TSQR" (node_tsqr_type) object.
      /// It was used to construct *pTsqr, but we keep it here because
      /// pTsqr used a raw pointer to node_tsqr_type.  This ensures
      /// that *pNodeTsqr_ stays in scope until we are done with pTsqr_.
      node_tsqr_ptr pNodeTsqr_;

      /// Shared pointer to the Tsqr object that performs the QR
      /// factorization.  pComm_.get() and *pNodeTsqr_ should have been
      /// used to construct it.
      tsqr_ptr pTsqr_;
    };

  } // namespace Trilinos
} // namespace TSQR

#ifdef HAVE_ANASAZI_EPETRA
#  include "TsqrAdaptor_Epetra_MultiVector.hpp"
#endif // HAVE_ANASAZI_EPETRA

#ifdef HAVE_ANASAZI_TPETRA
#  include "TsqrAdaptor_Tpetra_MultiVector_SerialNode.hpp"
#  include "TsqrAdaptor_Tpetra_MultiVector_TBBNode.hpp"
#endif // HAVE_ANASAZI_TPETRA

#endif // __TSQR_Trilinos_TsqrAdaptor_hpp
