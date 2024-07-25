// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Tsqr_CombineNodeTsqr.hpp
/// \brief Declaration and definition of an implementation of NodeTsqr
///   (intranode TSQR) that just uses Combine for all the operations
///   on an MPI process.

#ifndef TSQR_COMBINENODETSQR_HPP
#define TSQR_COMBINENODETSQR_HPP

#include "Tsqr_NodeTsqr.hpp"
#include "Tsqr_Impl_CombineUser.hpp"
#include "Tsqr_Impl_SystemBlas.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include <memory>

namespace TSQR {
  namespace Impl {
    template<class T>
    using span = Kokkos::View<T*, Kokkos::HostSpace,
      Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    template<class Ordinal, class Scalar>
    class CombineNodeFactorOutput :
      public NodeFactorOutput<Ordinal, Scalar> {
    public:
      CombineNodeFactorOutput (std::vector<Scalar>&& tau) :
        tau_ (tau)
      {}
      ~CombineNodeFactorOutput () override = default;
      span<const Scalar> tau () const {
        return span<const Scalar> (tau_.data (), tau_.size ());
      }
    private:
      std::vector<Scalar> tau_;
    };
  } // namespace Impl

  /// \class CombineNodeTsqr
  /// \brief Implementation of NodeTsqr (intranode TSQR) that just
  ///   uses Combine for all the operations on an MPI process.
  template<class Ordinal, class Scalar>
  class CombineNodeTsqr :
    public NodeTsqr<Ordinal, Scalar>,
    private Impl::CombineUser<Ordinal, Scalar> {
  private:
    using base_type = NodeTsqr<Ordinal, Scalar>;
    using my_factor_output_type =
      Impl::CombineNodeFactorOutput<Ordinal, Scalar>;

  public:
    using ordinal_type = typename base_type::ordinal_type;
    using scalar_type = typename base_type::scalar_type;
    using mat_view_type = typename base_type::mat_view_type;
    using const_mat_view_type =
      typename base_type::const_mat_view_type;
    using magnitude_type = typename base_type::magnitude_type;
    using factor_output_type = typename base_type::factor_output_type;

    ~CombineNodeTsqr () override = default;

    Teuchos::RCP<const Teuchos::ParameterList>
    getValidParameters () const override {
      return Teuchos::parameterList ("CombineNodeTsqr");
    }

    void
    setParameterList (const Teuchos::RCP<Teuchos::ParameterList>&) override
    {}

    bool ready() const override {
      return true;
    }

    size_t cache_size_hint() const override {
      return size_t (0);
    }

    std::string description () const override {
      using Teuchos::TypeNameTraits;
      std::ostringstream os;
      os << "CombineNodeTsqr<Ordinal="
         << TypeNameTraits<Ordinal>::name() << ", Scalar="
         << TypeNameTraits<Scalar>::name() << ">: Intranode "
        "Intraprocess TSQR based on TSQR::Combine";
      return os.str();
    }

  private:
    void
    factorImpl (const mat_view_type& R,
                const mat_view_type& A,
                std::vector<Scalar>& tau) const
    {
      const ordinal_type ncols = A.extent (1);
      TEUCHOS_ASSERT( R.extent (0) == ncols &&
                      R.extent (1) == ncols );
      auto& combine = this->getCombine (ncols);
      const ordinal_type lwork =
        combine.work_size (A.extent (0), ncols, ncols);
      std::vector<Scalar> work (lwork);
      combine.factor_first (A, tau.data (), work.data (), lwork);

      // Copy the R factor resulting from the factorization out of the
      // topmost block of A) into the R output argument.
      deep_copy (R, Scalar {});
      copy_upper_triangle (R, A);
    }

  public:
    Teuchos::RCP<factor_output_type>
    factor (const ordinal_type nrows,
            const ordinal_type ncols,
            Scalar A[],
            const ordinal_type lda,
            Scalar R[],
            const ordinal_type ldr,
            const bool /* contiguousCacheBlocks */) const override
    {
      // The "contiguous cache blocks" option does nothing here, since
      // we just defer to an internal library that expects
      // column-major matrices.
      mat_view_type A_view (nrows, ncols, A, lda);
      mat_view_type R_view (ncols, ncols, R, ldr);
      std::vector<Scalar> tau (ncols);
      factorImpl (R_view, A_view, tau);
      using Teuchos::rcp;
      return rcp (new my_factor_output_type (std::move (tau)));
    }

    void
    apply (const ApplyType& applyType,
           const ordinal_type nrows,
           const ordinal_type ncols_Q,
           const Scalar Q[],
           const ordinal_type ldq,
           const factor_output_type& factorOutput,
           const ordinal_type ncols_C,
           Scalar C[],
           const ordinal_type ldc,
           const bool /* contiguousCacheBlocks */) const override
    {
      const char prefix[] = "TSQR::CombineNodeTsqr::apply: ";

      // Quick exit and error tests
      if (ncols_Q == 0 || ncols_C == 0 || nrows == 0) {
        return;
      }
      else if (ldc < nrows) {
        std::ostringstream os;
        os << prefix << "ldc (= " << ldc << ") < nrows (= "
           << nrows << ")";
        throw std::invalid_argument (os.str());
      }
      else if (ldq < nrows) {
        std::ostringstream os;
        os << prefix << "ldq (= " << ldq << ") < nrows (= "
           << nrows << ")";
        throw std::invalid_argument (os.str());
      }

      const my_factor_output_type& output = [&] () {
        const my_factor_output_type* output_ptr =
          dynamic_cast<const my_factor_output_type*> (&factorOutput);
        if (output_ptr == nullptr) {
          using Teuchos::demangleName;
          using Teuchos::TypeNameTraits;
          using Teuchos::typeName;
          std::ostringstream os;
          os << prefix << "Input factor_output_type object was not "
            "created by the same type of NodeTsqr object as this "
            "one.  This object has type " << typeName (*this) <<
            " and its subclass of factor_output_type has type " <<
            TypeNameTraits<my_factor_output_type>::name () << ", but "
            "the input factor_output_type object has dynamic type "
            << demangleName (typeid (factorOutput).name ());
          throw std::invalid_argument (os.str ());
        }
        return *output_ptr;
      } ();

      auto& combine = this->getCombine (std::max (ncols_Q, ncols_C));
      const ordinal_type lwork =
        combine.work_size (nrows, ncols_C, ncols_C);
      std::vector<Scalar> work (lwork);

      const_mat_view_type Q_view (nrows, ncols_Q, Q, ldq);
      mat_view_type C_view (nrows, ncols_C, C, ldc);
      const auto tau = output.tau ();
      combine.apply_first (applyType, Q_view, tau.data (),
                           C_view, work.data (), lwork);
    }

    void
    explicit_Q (const ordinal_type nrows,
                const ordinal_type ncols_Q,
                const Scalar Q[],
                const ordinal_type ldq,
                const factor_output_type& factorOutput,
                const ordinal_type ncols_C,
                Scalar C[],
                const ordinal_type ldc,
                const bool contiguousCacheBlocks) const override
    {
      mat_view_type C_view (nrows, ncols_C, C, ldc);
      deep_copy (C_view, Scalar {});
      this->set_diagonal_entries_to_one (C_view);
      // Apply the Q factor to C, to extract the first ncols_C columns
      // of Q in explicit form.
      apply (ApplyType::NoTranspose,
             nrows, ncols_Q, Q, ldq, factorOutput,
             ncols_C, C, ldc, contiguousCacheBlocks);
    }

    void
    cache_block (const ordinal_type /* nrows */,
                 const ordinal_type /* ncols */,
                 Scalar /* A_out */ [],
                 const Scalar /* A_in */ [],
                 const ordinal_type /* lda_in */) const override
    {}

    void
    un_cache_block (const ordinal_type /* nrows */,
                    const ordinal_type /* ncols */,
                    Scalar /* A_out */ [],
                    const ordinal_type /* lda_out */,
                    const Scalar /* A_in */ []) const override
    {}

    void
    Q_times_B (const ordinal_type nrows,
               const ordinal_type ncols,
               Scalar Q[],
               const ordinal_type ldq,
               const Scalar B[],
               const ordinal_type ldb,
               const bool /* contiguousCacheBlocks */) const override
    {
      using Teuchos::NO_TRANS;

      // We don't do any other error checking here (e.g., matrix
      // dimensions), though it would be a good idea to do so.

      // Take the easy exit if available.
      if (ncols == 0 || nrows == 0) {
        return;
      }

      Impl::SystemBlas<Scalar> blas;
      mat_view_type Q_view (nrows, ncols, Q, ldq);
      // GEMM doesn't like its input and output arguments to alias
      // each other, so we use a (deep) copy.
      Matrix<ordinal_type, Scalar> Q_copy (Q_view);

      // Q_view := Q_copy * B.
      blas.GEMM (NO_TRANS, NO_TRANS,
                 nrows, ncols, ncols,
                 Scalar (1.0), Q_copy.data (), Q_copy.stride (1),
                 B, ldb,
                 Scalar {}, Q_view.data (), Q_view.stride (1));
    }

    void
    fill_with_zeros (const ordinal_type nrows,
                     const ordinal_type ncols,
                     Scalar A[],
                     const ordinal_type lda,
                     const bool /* contiguousCacheBlocks */) const override
    {
      mat_view_type A_view (nrows, ncols, A, lda);
      deep_copy (A_view, Scalar {});
    }

  protected:
    const_mat_view_type
    const_top_block (const const_mat_view_type& C,
                     const bool /* contiguousCacheBlocks */) const override
    {
      return C; // For this class, "cache blocking" does nothing.
    }

  public:
    bool
    QR_produces_R_factor_with_nonnegative_diagonal () const override
    {
      // FIXME (19 Dec 2019) If the combine type is dynamic, we can't
      // answer this question without knowing the number of columns.
      // Just guess for now.
      constexpr ordinal_type fakeNumCols = 10;
      auto& c = this->getCombine (fakeNumCols);
      return c.QR_produces_R_factor_with_nonnegative_diagonal ();
    }
  };
} // namespace TSQR

#endif // TSQR_COMBINENODETSQR_HPP
