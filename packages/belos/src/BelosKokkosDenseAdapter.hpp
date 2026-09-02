// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//
#ifndef BELOS_KOKKOS_DENSE_MAT_TRAITS_HPP
#define BELOS_KOKKOS_DENSE_MAT_TRAITS_HPP

/*! \file BelosKokkosDenseAdapter.hpp
  \brief Full specialization of Belos::DenseMatTraits for Kokkos::DualView
  with arbitrary scalarType. All views are expected to
  have 2 dimensions.
*/

#include "KokkosBatched_Iamax.hpp"
#include "KokkosBlas1_dot.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_LAPACK.hpp"
#include "BelosDenseMatTraits.hpp"

#include "Kokkos_DualView.hpp"
#include "Kokkos_Random.hpp"
#include "KokkosKernels_ArithTraits.hpp"
#include "KokkosBlas1_scal.hpp"
#include "KokkosBlas1_axpby.hpp"
#include "KokkosBlas1_nrm2.hpp"
#include "KokkosBlas3_gemm.hpp"
#include "KokkosBlas3_trsm.hpp"
#include "KokkosLapack_geqrf.hpp"
#include "KokkosLapack_gegqr.hpp"
#include "KokkosLapack_potrf.hpp"
#include "KokkosLapack_potrs.hpp"
#include "KokkosBatched_Axpy.hpp"
#include "KokkosBatched_Dot.hpp"
#include "KokkosBatched_Rot.hpp"
#include "KokkosBatched_Rotg.hpp"

#include <vector>

namespace Belos {

  //! Helper function for copying Kokkos::DualView into conjugate Kokkos::DualView
  template<typename V>
  void kokkos_transpose(const V& dst, const V& src)
  {
    Kokkos::parallel_for(Kokkos::MDRangePolicy<typename V::execution_space, Kokkos::Rank<2>>({0, 0}, {dst.extent(0), dst.extent(1)}),
    KOKKOS_LAMBDA(int i, int j)
    {
      dst(i, j) = KokkosKernels::ArithTraits<typename V::non_const_value_type>::conj( src(j, i) );
    });
  }

  //! Full specialization of Belos::DenseMatSolver for Kokkos::DualView.
  template<class Scalar, class DM>
  class KokkosDenseSolver : public DenseSolver<Scalar, DM>
  {
  public:
    typedef typename KokkosKernels::ArithTraits<Scalar>::val_type IST; //Impl Scalar Type, as used in Tpetra
    typedef DenseMatTraits<Scalar, DM> DMT;

    //! @name Constructor/Destructor Methods
    //@{
    //! Default constructor; matrix should be set using setMatrix(), LHS and RHS set with setVectors().
    KokkosDenseSolver() {}

    //! KokkosDenseSolver destructor.
    virtual ~KokkosDenseSolver() {}
    //@}

    //! @name Set Methods
    //@{
    //! Sets the pointers for coefficient matrix
    using DenseSolver<Scalar, DM>::setMatrix;

    //! Sets the pointers for left and right hand side vector(s).
    using DenseSolver<Scalar, DM>::setVectors;
    //@}

    //! @name Strategy Modifying Methods
    //@{

    //! Set if dense matrix is symmetric positive definite
    //! \note This method must be called before the factorization is performed, otherwise it will have no effect.
    using DenseSolver<Scalar, DM>::setSPD;

    //! Causes equilibration to be called just before the matrix factorization as part of the call to \c factor.
    /*! \note This method must be called before the factorization is performed, otherwise it will have no effect.
    */
    using DenseSolver<Scalar, DM>::factorWithEquilibration;

    //! All subsequent function calls will work with the transpose-type set by this method (\c Teuchos::NO_TRANS, \c Teuchos::TRANS, and \c Teuchos::CONJ_TRANS).
    /*! \note This interface will set correct transpose flag for matrix, including complex-valued linear systems.
    */
    using DenseSolver<Scalar, DM>::solveWithTransposeFlag;
    //@}

    //! @name Factor/Solve Methods
    //@{

    //! Computes the in-place LU factorization of the matrix.
    /*!
      \return Integer error code, set to 0 if successful.
    */
    int factor()
    {
      int INFO = 0;

      // Only factor matrix if it is new, otherwise factors have been computed
      if (newMatrix_)
      {
        Teuchos::LAPACK<int,Scalar> lapack;

        if (equilibrate_)
        {
          int M = DMT::GetNumRows(*A_);
          int N = DMT::GetNumCols(*A_);

          {
            int Min_MN = TEUCHOS_MIN(M,N);
            int LDA = DMT::GetStride(*A_);

            IPIV_.resize( Min_MN );

            DMT::SyncDeviceToHost(*A_);
            const Scalar * Aptr = DMT::GetConstRawHostPtr(*A_);

            // Compute equilibration scalings
            R_ = MDMT::Create(M, 1, false);
            if (!spd_)
              C_ = MDMT::Create(N, 1, false);

            MagnitudeType ROWCND, COLCND, AMAX;
            if (spd_)
              lapack.POEQU (M, Aptr, LDA, R_->view_host().data(), &ROWCND, &AMAX, &INFO);
            else
              lapack.GEEQU (M, N, Aptr, LDA, R_->view_host().data(), C_->view_host().data(), &ROWCND, &COLCND, &AMAX, &INFO);

            if (INFO)
              return INFO;

            MDMT::SyncHostToDevice(*R_);
            if (!spd_)
              MDMT::SyncHostToDevice(*C_);
          }

          // Apply equilibration to matrix
          TEUCHOS_ASSERT(!A_->need_sync_device());
          auto A_d = A_->view_device();
          auto R_d = R_->view_device();
          if (spd_) {
            Kokkos::parallel_for(Kokkos::MDRangePolicy<typename DM::t_dev::execution_space, Kokkos::Rank<2>>({0, 0}, {A_d.extent(0), A_d.extent(1)}),
                                 KOKKOS_LAMBDA(int i, int j) {
              A_d(i, j) *= R_d(i, 0)*R_d(j, 0);
            });
          } else {
            auto C_d = C_->view_device();
            Kokkos::parallel_for(Kokkos::MDRangePolicy<typename DM::t_dev::execution_space, Kokkos::Rank<2>>({0, 0}, {A_d.extent(0), A_d.extent(1)}),
                                 KOKKOS_LAMBDA(int i, int j) {
              A_d(i, j) *= R_d(i, 0)*C_d(j, 0);
            });
          }
        }

        // Compute LU factor
        if (spd_) {
          DMT::potrf("U", *A_);
        }
        else {
          int M = DMT::GetNumRows(*A_);
          int N = DMT::GetNumCols(*A_);
          int LDA = DMT::GetStride(*A_);
          DMT::SyncDeviceToHost(*A_);
          Scalar * Aptr = DMT::GetRawHostPtr(*A_);
          lapack.GETRF(M, N, Aptr, LDA, &IPIV_[0], &INFO);
          DMT::SyncHostToDevice(*A_);
        }

      }

      return INFO;
    }

    //! Computes the solution X to AX = B for the \e this matrix and the B provided.
    /*!
      \return Integer error code, set to 0 if successful.
    */
    int solve()
    {
      bool transpose = (TRANS_ != Teuchos::NO_TRANS) ? true : false;

      // LAPACK overwrites RHS vector with solution vector, so copy if necessary
      if (B_ != X_)
        DMT::Assign(*X_, *B_); // Copy B to X if needed

      if (equilibrate_)
      {
        // Since B_ = X_, perform operations on X_.

        // Apply equilibration scalings to RHS vector
        typename MDM::t_dev R_tmp;
        if (transpose && !spd_)
          R_tmp = C_->view_device();
        else
          R_tmp = R_->view_device();
        TEUCHOS_ASSERT(!X_->need_sync_device());
        auto X_d = X_->view_device();
        Kokkos::parallel_for(Kokkos::MDRangePolicy<typename DM::t_dev::execution_space, Kokkos::Rank<2>>({0, 0}, {X_d.extent(0), X_d.extent(1)}),
                             KOKKOS_LAMBDA(int i, int j) {
          X_d(i, j) *= R_tmp(i, 0);
                             });
      }

      int INFO = 0;
      if (spd_) {
        DMT::potrs("U", *A_, *X_);
      }
      else {
        DMT::SyncDeviceToHost(*X_);

        Scalar * X = DMT::GetRawHostPtr(*X_);
        int LDX = DMT::GetStride(*X_);
        int M = DMT::GetNumRows(*X_);
        int NRHS = DMT::GetNumCols(*X_);

        int LDA = DMT::GetStride(*A_);
        Scalar * Aptr = DMT::GetRawHostPtr(*A_);
        Teuchos::LAPACK<int,Scalar> lapack;
        lapack.GETRS(Teuchos::ETranspChar[TRANS_], M, NRHS,
                     Aptr, LDA, &IPIV_[0], X, LDX, &INFO);
        DMT::SyncHostToDevice(*X_);
      }

      if (equilibrate_)
      {
        // Apply equilibration scalings to X vector
        typename MDM::t_dev C_tmp;
        if (spd_ || transpose)
          C_tmp = R_->view_device();
        else
          C_tmp = C_->view_device();
        TEUCHOS_ASSERT(!X_->need_sync_device());
        auto X_d = X_->view_device();
        Kokkos::parallel_for(Kokkos::MDRangePolicy<typename DM::t_dev::execution_space, Kokkos::Rank<2>>({0, 0}, {X_d.extent(0), X_d.extent(1)}),
                             KOKKOS_LAMBDA(int i, int j) {
          X_d(i, j) *= C_tmp(i, 0);
                             });
      }

      return INFO;
    }
    //@}

  private:

    using ST = KokkosKernels::ArithTraits<Scalar>;
    using MagnitudeType = typename KokkosKernels::ArithTraits<IST>::mag_type;
    using MDM = Kokkos::DualView<MagnitudeType **, Kokkos::LayoutLeft>;
    using MDMT = DenseMatTraits<MagnitudeType, MDM>;

    std::vector<int> IPIV_;
    Teuchos::RCP<MDM> R_, C_;

    using DenseSolver<Scalar, DM>::newMatrix_;
    using DenseSolver<Scalar, DM>::equilibrate_;
    using DenseSolver<Scalar, DM>::TRANS_;
    using DenseSolver<Scalar, DM>::spd_;

    using DenseSolver<Scalar, DM>::A_;
    using DenseSolver<Scalar, DM>::X_;
    using DenseSolver<Scalar, DM>::B_;

  };

  //! Full specialization of Belos::DenseMatTraits for Kokkos::DualView.
  //
  template<class Scalar, class... Properties>
  class DenseMatTraits<Scalar, Kokkos::DualView<typename KokkosKernels::ArithTraits<Scalar>::val_type **,Properties...>>{

  public:
    using ST = KokkosKernels::ArithTraits<Scalar>;
    using IST = typename ST::val_type;
    using DM = Kokkos::DualView<IST**,Properties...>;
    using MagnitudeType = typename KokkosKernels::ArithTraits<IST>::mag_type;
    using MDM = Kokkos::DualView<MagnitudeType**, Properties...>;
    using MDMT = DenseMatTraits<MagnitudeType, MDM>;

    //@{ \name Creation methods

    /*! \brief Creates a new empty \c Kokkos::DualView with no dimension.

    \return Reference-counted pointer to a new dense matrix of type \c Kokkos::DualView.
    */
    static Teuchos::RCP<DM> Create() {
      return Teuchos::rcp(new DM("BelosDenseCreate",0,0));
    }

    /*! \brief Creates a new empty \c Kokkos::DualView containing
     *     and \c numRows rows and \c numCols columns.
     *         Will be initialized to zeros if last parameter is true.

    \return Reference-counted pointer to a new dense matrix of type \c Kokkos::DualView.
    */
    static Teuchos::RCP<DM>
           Create( const int numRows, const int numCols, bool initZero = true) {
      if(initZero){
        return Teuchos::rcp(new DM("BelosDenseCreate2",numRows,numCols));
      }
      else {
        return Teuchos::rcp(new DM(Kokkos::view_alloc(Kokkos::WithoutInitializing,"BelosDenseCreate2"),numRows,numCols));
      }
    }

    /*! \brief Create a new copy \c DM, possibly transposed.

       \return Reference-counted pointer to a new dense matrix of type \c DM.
    */
    static Teuchos::RCP<DM>
           CreateCopy(const DM & dm, bool transpose=false)
    {
      Teuchos::RCP<DM> tmpCopyRCP = Teuchos::null;

      if (transpose) {
        // want tmpCopyRCP to end up as dm^H in the end. Prefer doing transpose on device
        tmpCopyRCP = Teuchos::rcp(new DM
                                 (Kokkos::view_alloc(Kokkos::WithoutInitializing,"BelosDenseCreateCopy"),dm.extent_int(1),dm.extent_int(0)));
        if(tmpCopyRCP->need_sync_device()) {
          // tmpCopyRCP is only up to date on the host
          kokkos_transpose(tmpCopyRCP->view_host(), dm.view_host());
          tmpCopyRCP->clear_sync_state();
          tmpCopyRCP->modify_host();
        }
        else {
          kokkos_transpose(tmpCopyRCP->view_device(), dm.view_device());
          tmpCopyRCP->clear_sync_state();
          tmpCopyRCP->modify_device();
        }
      }
      else {
        tmpCopyRCP = Teuchos::rcp(new DM
                                 (Kokkos::view_alloc(Kokkos::WithoutInitializing,"BelosDenseCreateCopy"),dm.extent_int(0),dm.extent_int(1)));
        Kokkos::deep_copy(*tmpCopyRCP, dm);
      }
      return tmpCopyRCP;
    }

    //! \brief Returns a raw pointer to the (non-const) data on the host.
    /// \note We assume that return data in in a column-major format
    /// because this is what is expected by LAPACK.
    /// \note This raw pointer is intended only for passing data to LAPACK
    /// functions. Other operations on the raw data may result in undefined behavior!
    /// \note We assume that getting this pointer to non-const data means that data
    /// will be modified. If that is not true, use of this function could result in
    /// extra data syncs.
    /// The user must NOT hold onto this pointer after any data syncs. If you modify the
    /// data again after a sync, the view will not be marked as modified the second time.
    static Scalar* GetRawHostPtr(DM & dm ) {
      dm.sync_host();
      dm.modify_host();
      return reinterpret_cast<Scalar*>(dm.view_host().data());
      //TODO: Is there any way that the user could hold on to this pointer...
      // and everything works fine the first time they pass to LAPACK.
      // But then... they call MvTimesMatAddMv which syncs to device.
      // But then they keep the same pointer and pass to LAPACK again.
      // Then they call MvTimesMatAddMv... but since they didn't call this
      // function again, we miss the sync... See thread with Heidi on this.
    }

    //! \brief Returns a raw pointer to const data on the host.
    static Scalar const * GetConstRawHostPtr(const DM & dm ) {
      // CAG: This is a bit naughty.
      const_cast<DM*>(&dm)->sync_host();
      return reinterpret_cast<Scalar const *>(dm.view_host().data());
    }

    //! \brief Returns an RCP to a Kokkos::DualView which has a subview of the given Kokkos::DualView.
    //        Row and column indexing is zero-based.
    static Teuchos::RCP<DM>
    Subview( DM & source, int numRows, int numCols, int startRow=0, int startCol=0){
      return Teuchos::rcp(new DM(source,
          Kokkos::pair<int,int>(startRow,startRow+numRows), Kokkos::pair<int,int>(startCol,startCol+numCols)));
    }

    static Teuchos::RCP<const DM>
    SubviewConst( const DM& source, int numRows, int numCols, int startRow=0, int startCol=0){
      return Teuchos::rcp(new DM(source,
          Kokkos::pair<int,int>(startRow,startRow+numRows), Kokkos::pair<int,int>(startCol,startCol+numCols)));
    }

    //! \brief Returns a deep copy of the requested subview.
    static Teuchos::RCP<DM>
    SubviewCopy( const DM& source, int numRows, int numCols, int startRow=0, int startCol=0){
      //Maaybe we could get away with just a host copy here??
      //Hmmm... but it says we return a dual view.
      //Maybe it should return a dual view with only host data copied in. Require sync to work on device.
      //This is related to the functionality of the Assign function.
      auto tmpViewRCP = Teuchos::rcp(new DM
                                  (Kokkos::view_alloc(Kokkos::WithoutInitializing,"BelosDenseSubViewCopy"),numRows,numCols));
      // I am keeping this where it copies the whole view on host and device because:
      // a) I feel like we might be inviting some weird bugs later if we don't.
      // But TODO Clarify to developer that this function needs to work on both host and device.
      // b) It's not a host-device copy or vice versa. Its a copy from device to same device and from host to host.
      // So shouldn't add much extra overhead.
      Kokkos::deep_copy(*tmpViewRCP, Kokkos::subview(source,
          Kokkos::pair<int,int>(startRow,startRow+numRows), Kokkos::pair<int,int>(startCol,startCol+numCols)));
      return tmpViewRCP;
    }
    //@}

    //@{ \name Attribute methods

    //! \brief Obtain the number of rows of \c dm.
    static int GetNumRows( const DM& dm ) {
      return dm.extent_int(0);
    }

    //! \brief Obtain the number of columns of \c dm.
    static int GetNumCols( const DM& dm ) {
      return dm.extent_int(1);
    }

    //! \brief Obtain the stride between the columns of \c dm.
    static int GetStride( const DM& dm ) {
      // Note: We force LayoutLeft, which is column major, so the stride_0 is always 1.
      // (This is the distance between two elts in same col, different rows.)
      // Lapack wants stride_1, the distance from one col to the next col if we stay
      // in the same row.
      int strides[8]; // There are 8 possible strides and all will be returned
      dm.stride(strides);
      return strides[1];
      //return dm.stride_1();  //This shortcut doesn't work for dualView.
    }

    //@}

    //@{ \name Shaping methods

    /* \brief Reshaping method for changing the size of \c dm to have \c numRows rows and \c numCols columns.
     *        All values will be initialized to zero if the final argument is true.
     *        If the final argument is fale, the previous entries in
     *        the matrix will be maintained. For new entries that did not exist in the previous matrix, values will
     *        contain noise from memory.
    */
    static void Reshape( DM& dm, const int numRows, const int numCols, bool initZero = false) {
      if(initZero){
        dm.realloc(numRows,numCols);
        Kokkos::deep_copy(dm.view_device(), 0.0);
        dm.modify_device();
      } else{
        dm.resize(numRows,numCols); //keeps values in old array.
      }
    }

    //@}

    //@{ \name Data access methods

    //! \brief Access a reference to the (i,j) entry of \c dm, \c e_i^T dm e_j.
    static Scalar & Value( DM& dm, const int i, const int j )
    {
      // Mark as modified on host, since we don't know if it will be.
      dm.sync_host();
      dm.modify_host();
      return reinterpret_cast<Scalar&>((dm.view_host())(i,j));
    }

    //! \brief Access a const reference to the (i,j) entry of \c dm, \c e_i^T dm e_j.
    static const Scalar & ValueConst( const DM& dm, const int i, const int j ) {
      // CAG: This is a bit naughty.
      const_cast<DM*>(&dm)->sync_host();
      return reinterpret_cast<Scalar const &>((dm.view_host())(i,j));
    }

    //! \brief If an accelorator is in use, sync it to device on this call.
    //
    //  \note The only Belos function that results in a need to sync to
    //  host is MvTransMv. You MUST call SyncDeviceToHost before calling
    //  any other DenseMatTraits functions after a call to MvTransMv.
    //  All DenseMatTraits functions assume the necessary data is on host
    //  and perform computations only on the host.
    //
    static void SyncDeviceToHost(DM & dm) {
      if(dm.need_sync_host()){
        if(dm.view_host().span_is_contiguous() && dm.view_device().span_is_contiguous()){
        dm.sync_host();}
        else{
          DM compat_view("compat view",dm.extent_int(0),dm.extent_int(1));
          Kokkos::deep_copy(compat_view,dm);
          compat_view.sync_host();
          Kokkos::deep_copy(dm,compat_view);
          dm.clear_sync_state();
        }
      }
    }

    static void SyncHostToDevice(DM & dm) {
      if(dm.need_sync_device()){
        if(dm.view_host().span_is_contiguous() && dm.view_device().span_is_contiguous()){
          dm.sync_device();
        }
        else{
          DM compat_view("compat view",dm.extent_int(0),dm.extent_int(1));
          Kokkos::deep_copy(compat_view,dm);
          compat_view.sync_device();
          Kokkos::deep_copy(dm,compat_view);
          dm.clear_sync_state();
        }
      }
    }
    //@}
    //@{ \name Operator methods

    //!  \brief Adds sourceDM to thisDM and returns answer in thisDM.
    static void Add( DM& thisDM, const DM& sourceDM) {
      thisDM.sync_device();
      // CAG: This is a bit naughty.
      const_cast<DM*>(&sourceDM)->sync_device();
      KokkosBlas::axpy(1.0,sourceDM.view_device(), thisDM.view_device()); //axpy(alpha,x,y), y = y + alpha*x
      thisDM.modify_device();
    }

    //!  \brief Add source to diagonal entries of dest
    static void AddDiag( DM& dest, Scalar source) {
      TEUCHOS_ASSERT(GetNumRows(dest) == GetNumCols(dest));
      TEUCHOS_ASSERT(!dest.need_sync_device());

      auto dest_d = dest.view_device();
      Kokkos::parallel_for(Kokkos::RangePolicy<typename DM::t_dev::execution_space>(0, GetNumRows(dest)),
          KOKKOS_LAMBDA(size_t i)
          {
            dest_d(i, i) += source;
          });
      dest.modify_device();
    }

    //!  \brief Fill all entries with \c value. Value is zero if not specified.
    static void PutScalar( DM& dm, Scalar value = Teuchos::ScalarTraits<Scalar>::zero()){
      dm.clear_sync_state();
      Kokkos::deep_copy( dm.view_device(), value);
      dm.modify_device();
    }

    //!  \brief Multiply all entries by a scalar. DM = value.*DM
    static void Scale( DM& dm, Scalar value) {
      dm.sync_device();
      KokkosBlas::scal( dm.view_device(), value, dm.view_device());
      dm.modify_device();
    }

    //!  \brief Multiply two dense matrices. C = beta*C + alpha*op(A)*op(B)
    static void Multiply( bool transposeA, bool transposeB, Scalar alpha, const DM &A, const DM &B, Scalar beta, DM& C)
    {
      TEUCHOS_ASSERT(!A.need_sync_device());
      TEUCHOS_ASSERT(!B.need_sync_device());
      TEUCHOS_ASSERT(!C.need_sync_device());
      KokkosBlas::gemm(transposeA ? "T" : "N", transposeB ? "T" : "N", alpha, A.view_device(), B.view_device(), beta, C.view_device());
      C.modify_device();
    }

    //!  \brief Fill the Kokkos::DualView with random entries.
    //!   Entries are assumed to be the same on each MPI rank (each matrix copy).
    static void Randomize( DM& dm) {
      int rand_seed = std::rand();
      Kokkos::Random_XorShift64_Pool<> pool(rand_seed);
      dm.clear_sync_state();
      Kokkos::fill_random(dm.view_device(), pool, -1,1);
      dm.modify_device();
    }

    //!  \brief Copies entries of source to dest (deep copy).
    static void Assign( DM& dest, const DM& source) {
      Kokkos::deep_copy(dest,source);
    }

    //!  \brief Assign source to diagonal entries of dest
    static void AssignDiag( DM& dest, Scalar source) {
      TEUCHOS_ASSERT(GetNumRows(dest) == GetNumCols(dest));
      TEUCHOS_ASSERT(!dest.need_sync_device());

      auto dest_d = dest.view_device();
      Kokkos::parallel_for(Kokkos::RangePolicy<typename DM::t_dev::execution_space>(0, GetNumRows(dest)),
          KOKKOS_LAMBDA(size_t i)
          {
            dest_d(i, i) = source;
          });
      dest.modify_device();
    }

    //!  \brief Assign source to diagonal entries of dest
    static void AssignDiag( DM& dest, const DM& source) {
      TEUCHOS_ASSERT(GetNumRows(dest) == GetNumRows(source));
      TEUCHOS_ASSERT(GetNumCols(dest) == GetNumRows(source));
      TEUCHOS_ASSERT(GetNumCols(source) == 1);
      TEUCHOS_ASSERT(!source.need_sync_device());
      TEUCHOS_ASSERT(!dest.need_sync_device());

      auto dest_d = dest.view_device();
      auto source_d = source.view_device();
      Kokkos::parallel_for(Kokkos::RangePolicy<typename DM::t_dev::execution_space>(0, GetNumRows(dest)),
          KOKKOS_LAMBDA(size_t i)
          {
            dest_d(i, i) = source_d(i, 0);
          });
      dest.modify_device();
    }

    //!  \brief Assign upper triangular entries of source to dest
    static void AssignUpperTri( DM& dest, const DM& source) {
      TEUCHOS_ASSERT(GetNumRows(dest) == GetNumRows(source));
      TEUCHOS_ASSERT(GetNumCols(dest) == GetNumCols(source));
      TEUCHOS_ASSERT(!source.need_sync_device());
      TEUCHOS_ASSERT(!dest.need_sync_device());

      auto dest_d = dest.view_device();
      auto source_d = source.view_device();
      Kokkos::parallel_for(Kokkos::MDRangePolicy<typename DM::t_dev::execution_space, Kokkos::Rank<2>>({0, 0}, {GetNumRows(dest), GetNumCols(dest)}),
          KOKKOS_LAMBDA(size_t i, size_t j)
          {
            if (i <= j)
              dest_d(i, j) = source_d(i, j);
          });
      dest.modify_device();
    }

    //!  \brief Returns the Frobenius norm of the dense matrix.
    static typename Teuchos::ScalarTraits<Scalar>::magnitudeType NormFrobenius(const DM& dm) {
      using KAT = KokkosKernels::ArithTraits<IST>;
      using mag_t = typename KAT::mag_type;
      // CAG: This is a bit naughty.
      const_cast<DM*>(&dm)->sync_device();
      mag_t frobNorm;
      Kokkos::parallel_reduce(Kokkos::MDRangePolicy<typename DM::t_dev::execution_space, Kokkos::Rank<2>>({0, 0}, {dm.extent(0), dm.extent(1)}),
          KOKKOS_LAMBDA(size_t i, size_t j, mag_t& lfrobNorm)
          {
          mag_t absVal = KAT::abs((dm.view_device())(i, j));
          lfrobNorm += absVal * absVal;
          }, frobNorm);
      return Kokkos::sqrt(frobNorm);
    }

    //!  \brief Returns the one-norm of the dense matrix.
    static typename Teuchos::ScalarTraits<Scalar>::magnitudeType NormOne(const DM& dm) {
      using KAT = KokkosKernels::ArithTraits<IST>;
      using mag_t = typename KAT::mag_type;
      // CAG: This is a bit naughty.
      const_cast<DM*>(&dm)->sync_device();
      mag_t max_sum = 0;

      Kokkos::parallel_reduce(dm.extent(1), KOKKOS_LAMBDA(const int j, mag_t& norm) {
        mag_t sum = 0;
        for(int i = 0; i < dm.extent_int(0); i++){  //rows
          sum += KAT::abs((dm.view_device())(i,j));
        }
        norm = Kokkos::max(norm, sum);
      }, Kokkos::Max<mag_t>(max_sum));
      return KAT::abs(max_sum);
    }
    //@}

    //@{ \name Solver methods

    //!  \brief Returns a dense solver object for the dense matrix.
    static Teuchos::RCP<DenseSolver<Scalar, DM>>
      createDenseSolver() {

      Teuchos::RCP<DenseSolver<Scalar, DM>> newSolver
        = Teuchos::rcp( new KokkosDenseSolver<Scalar, DM>() );
      return newSolver;
    }
    //@}

  static void trsm(const char side[], const char uplo[], const char trans[], const char diag[],
                   const Scalar& alpha, const DM& A, DM& B) {
    TEUCHOS_ASSERT(!A.need_sync_device());
    TEUCHOS_ASSERT(!B.need_sync_device());
    KokkosBlas::trsm(side, uplo, trans, diag, alpha, A.view_device(), B.view_device());
    B.modify_device();
  }

  static void nrm2(std::vector<MagnitudeType>&R, const DM &X) {
    TEUCHOS_ASSERT(!X.need_sync_device());
    Kokkos::View<MagnitudeType*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged >> R_view(R.data(), R.size());
    KokkosBlas::nrm2(R_view, X.view_device());
  }

  static void geqrf(DM &A, DM &tau) {
    TEUCHOS_ASSERT(!A.need_sync_device());
    TEUCHOS_ASSERT(GetNumCols(tau) == 1);
    Kokkos::View<int*, Kokkos::LayoutLeft, typename DM::t_dev::memory_space> info(Kokkos::ViewAllocateWithoutInitializing("geqrf_info"), 1);
    auto tau_view = Kokkos::subview(tau.view_device(), Kokkos::ALL(), 0);
    KokkosLapack::geqrf(A.view_device(),
                        tau_view,
                        info);
    auto info_h = Kokkos::create_mirror_view_and_copy(info);
    TEUCHOS_ASSERT(info_h(0) == 0);
    A.modify_device();
    tau.modify_device();
  }

  static void ungqr(const int k, DM& A, const DM& tau) {
    TEUCHOS_ASSERT(!A.need_sync_device());
    TEUCHOS_ASSERT(GetNumCols(tau) == 1);
    Kokkos::View<int*, Kokkos::LayoutLeft, typename DM::t_dev::memory_space> info(Kokkos::ViewAllocateWithoutInitializing("geqrf_info"), 1);
    KokkosLapack::gegqr(k, A.view_device(), Kokkos::subview(tau.view_device(), Kokkos::ALL(), 0), info);
    auto info_h = Kokkos::create_mirror_view_and_copy(info);
    TEUCHOS_ASSERT(info_h(0) == 0);
    A.modify_device();
  }

  static void potrf(const char uplo[], DM& A) {
    TEUCHOS_ASSERT(!A.need_sync_device());
    KokkosLapack::potrf(uplo, A.view_device());
    A.modify_device();
  }

  static void potrs(const char uplo[], const DM& A, DM &X) {
    TEUCHOS_ASSERT(!A.need_sync_device());
    TEUCHOS_ASSERT(!X.need_sync_device());
    KokkosLapack::potrs(uplo, A.view_device(), X.view_device());
    X.modify_device();
  }

  static void updateLSQR(DM &H, DM &z, Teuchos::RCP<MDM> cs, Teuchos::RCP<DM> sn, Teuchos::RCP<DM> beta, int dim, int blockSize) {
    const Scalar zero = ST::zero();

    if (blockSize == 1) {
      TEUCHOS_ASSERT(!H.need_sync_device())
      TEUCHOS_ASSERT(!z.need_sync_device())
      TEUCHOS_ASSERT(!cs.is_null())
      TEUCHOS_ASSERT(!sn.is_null())
      TEUCHOS_ASSERT(!cs->need_sync_device())
      TEUCHOS_ASSERT(!sn->need_sync_device())

      auto H_dv = H.view_device();
      auto z_dv = z.view_device();
      auto cs_dv = cs->view_device();
      auto sn_dv = sn->view_device();

      Kokkos::parallel_for("Belos::updateLSQR", Kokkos::RangePolicy<>(0, 1), KOKKOS_LAMBDA(const int k) {
        for (int i = 0; i<dim; ++i) {
          KokkosBatched::SerialRot<true>::invoke(Kokkos::subview(H_dv, Kokkos::make_pair(i, i+1), dim),
                                                 Kokkos::subview(H_dv, Kokkos::make_pair(i+1, i+2), dim),
                                                 cs_dv(i, 0),
                                                 sn_dv(i, 0));
        }
        KokkosBatched::Rotg::invoke(Kokkos::subview(H_dv, dim, dim),
                                    Kokkos::subview(H_dv, dim+1, dim),
                                    Kokkos::subview(cs_dv, dim, 0),
                                    Kokkos::subview(sn_dv, dim, 0));
        H_dv(dim+1, dim) = zero;
        KokkosBatched::SerialRot<true>::invoke(Kokkos::subview(z_dv, Kokkos::make_pair(dim, dim+1), 0),
                                               Kokkos::subview(z_dv, Kokkos::make_pair(dim+1, dim+2), 0),
                                               cs_dv(dim, 0),
                                               sn_dv(dim, 0));
      });
      H.modify_device();
      z.modify_device();
      cs->modify_device();
      sn->modify_device();
    } else {
      TEUCHOS_ASSERT(!H.need_sync_device())
      TEUCHOS_ASSERT(!z.need_sync_device())
      TEUCHOS_ASSERT(!beta.is_null())
      TEUCHOS_ASSERT(!beta->need_sync_device())

      auto H_dv = H.view_device();
      auto z_dv = z.view_device();
      auto beta_dv = beta->view_device();

      Kokkos::parallel_for("Belos::updateLSQR", Kokkos::RangePolicy<>(0, 1), KOKKOS_LAMBDA(const int k) {
        IST sigma, mu, vscale;
        Kokkos::View<IST, Kokkos::MemoryTraits<Kokkos::Unmanaged >> sigma_view(&sigma);
        //
        // QR factorization of Least-Squares system with Householder reflectors
        //
        for (int j=0; j<blockSize; j++) {
          //
          // Apply previous Householder reflectors to new block of Hessenberg matrix
          //
          for (int i=0; i<dim+j; i++) {
            auto X = Kokkos::subview(H_dv, Kokkos::make_pair(i+1, i+1+blockSize), i);
            auto Y = Kokkos::subview(H_dv, Kokkos::make_pair(i+1, i+1+blockSize), dim+j);
            KokkosBatched::SerialDot<KokkosBatched::Trans::NoTranspose, 0>::invoke(X, Y, sigma_view);
            sigma += H_dv(i, dim+j);
            sigma *= ST::conjugate(beta_dv(i, 0));
            KokkosBatched::SerialAxpy::invoke(-sigma, X, Y);
            H_dv(i, dim+j) -= sigma;
          }
          //
          // Compute new Householder reflector
          //
          auto XX = Kokkos::subview(H_dv, Kokkos::make_pair(dim+j, dim+j+blockSize+1), dim+j);
          auto maxidx = KokkosBatched::SerialIamax::invoke(XX);
          auto maxelem = ST::magnitude(XX(maxidx));
          for (int i=0; i<blockSize+1; i++)
            H_dv(dim+j+i,dim+j) /= maxelem;
          auto Z = Kokkos::subview(H_dv, Kokkos::make_pair(dim+j+1, dim+j+1+blockSize), dim+j);
          KokkosBatched::SerialDot<KokkosBatched::Trans::NoTranspose, 0>::invoke(Z, Z, sigma_view);
          MagnitudeType sign_Rjj = -ST::real(H_dv(dim+j,dim+j)) /
                                   ST::magnitude(ST::real((H_dv(dim+j,dim+j))));
          if (sigma == zero) {
            beta_dv(dim + j, 0) = zero;
          } else {
            mu = ST::squareroot(ST::conjugate(H_dv(dim+j,dim+j))*H_dv(dim+j,dim+j)+sigma);
            vscale = H_dv(dim+j,dim+j) - Teuchos::as<Scalar>(sign_Rjj)*mu;
            beta_dv(dim + j, 0) = -Teuchos::as<Scalar>(sign_Rjj) * vscale / mu;
            H_dv(dim+j,dim+j) = Teuchos::as<Scalar>(sign_Rjj)*maxelem*mu;
            for (int i=0; i<blockSize; i++)
              H_dv(dim+j+1+i,dim+j) /= vscale;
          }
          //
          // Apply new Householder reflector to rhs
          //
          for (int i=0; i<blockSize; i++) {
            auto X = Kokkos::subview(H_dv, Kokkos::make_pair(dim+j+1, dim+j+1+blockSize), dim+j);
            auto Y = Kokkos::subview(z_dv, Kokkos::make_pair(dim+j+1, dim+j+1+blockSize), i);
            KokkosBatched::SerialDot<KokkosBatched::Trans::NoTranspose, 0>::invoke(X, Y, sigma_view);
            sigma += z_dv(dim+j,i);
            sigma *= ST::conjugate(beta_dv(dim+j, 0));
            KokkosBatched::SerialAxpy::invoke(-sigma, X, Y);
            z_dv(dim+j,i) -= sigma;
          }
        }
      });

      H.modify_device();
      z.modify_device();
      beta->modify_device();
    }
  }
};


} // namespace Belos

#endif // end file BELOS_KOKKOS_DENSE_MAT_TRAITS_HPP
