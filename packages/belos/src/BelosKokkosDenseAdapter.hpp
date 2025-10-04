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
  have 2 dimensions. We hard-code LayoutLeft because LAPACK expects
  column-major matrices. 
*/

#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_LAPACK.hpp"
#include "BelosDenseMatTraits.hpp"

#include "Kokkos_Random.hpp"
#include "Kokkos_ArithTraits.hpp"
//Kokkos BLAS files:
#include "KokkosBlas1_scal.hpp"
#include "KokkosBlas1_axpby.hpp"

#include <vector>

namespace Belos {

  //! Helper function for copying Kokkos::DualView into conjugate Kokkos::DualView
  template<typename V>
  void kokkos_transpose(V& dst, const V& src)
  {
    Kokkos::parallel_for(Kokkos::MDRangePolicy<typename V::execution_space, Kokkos::Rank<2>>({0, 0}, {dst.extent(0), dst.extent(1)}),
    KOKKOS_LAMBDA(int i, int j)
    {
      dst(i, j) = Kokkos::ArithTraits<typename V::non_const_value_type>::conj( src(j, i) );
    });
  }

  //! Full specialization of Belos::DenseMatSolver for Kokkos::DualView.
  template<class Scalar>
  class KokkosDenseSolver : public DenseSolver<Scalar, Kokkos::DualView<typename Kokkos::Details::ArithTraits<Scalar>::val_type **,Kokkos::LayoutLeft>>
  {
  public:
    typedef typename Kokkos::Details::ArithTraits<Scalar>::val_type IST; //Impl Scalar Type, as used in Tpetra
    typedef DenseMatTraits<Scalar, Kokkos::DualView<IST**,Kokkos::LayoutLeft>>     DMT;
 
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
    using DenseSolver<Scalar, Kokkos::DualView<IST**,Kokkos::LayoutLeft>>::setMatrix;

    //! Sets the pointers for left and right hand side vector(s).
    using DenseSolver<Scalar, Kokkos::DualView<IST**,Kokkos::LayoutLeft>>::setVectors;
    //@}

    //! @name Strategy Modifying Methods
    //@{

    //! Set if dense matrix is symmetric positive definite
    //! \note This method must be called before the factorization is performed, otherwise it will have no effect.
    using DenseSolver<Scalar, Kokkos::DualView<IST**,Kokkos::LayoutLeft>>::setSPD;

    //! Causes equilibration to be called just before the matrix factorization as part of the call to \c factor.
    /*! \note This method must be called before the factorization is performed, otherwise it will have no effect.
    */
    using DenseSolver<Scalar, Kokkos::DualView<IST**,Kokkos::LayoutLeft>>::factorWithEquilibration;

    //! All subsequent function calls will work with the transpose-type set by this method (\c Teuchos::NO_TRANS, \c Teuchos::TRANS, and \c Teuchos::CONJ_TRANS).
    /*! \note This interface will set correct transpose flag for matrix, including complex-valued linear systems.
    */
    using DenseSolver<Scalar, Kokkos::DualView<IST**,Kokkos::LayoutLeft>>::solveWithTransposeFlag;
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

        int M = DMT::GetNumRows(*A_);
        int N = DMT::GetNumCols(*A_);
        int Min_MN = TEUCHOS_MIN(M,N);
        int LDA = DMT::GetStride(*A_);

        IPIV_.resize( Min_MN );

        DMT::SyncDeviceToHost(*A_);
        Scalar * Aptr = DMT::GetRawHostPtr(*A_);

        if (equilibrate_) 
        {
          // Compute equilibration scalings
          R_.resize(M);
          if (!spd_)
            C_.resize(N);

          MagnitudeType ROWCND, COLCND, AMAX;
          if (spd_)
            lapack.POEQU (M, Aptr, LDA, &R_[0], &ROWCND, &AMAX, &INFO);
          else
            lapack.GEEQU (M, N, Aptr, LDA, &R_[0], &C_[0], &ROWCND, &COLCND, &AMAX, &INFO); 

          if (INFO)
            return INFO;

          // Apply equilibration to matrix
          if (spd_) {
            Scalar * ptr = 0;
            for (int j=0; j<N; j++) {
              ptr = Aptr + j*LDA;
              Scalar s1 = R_[j];
              for (int i=0; i<=j; i++) {
                *ptr = *ptr*s1*R_[i];
                ptr++;
              }
            }
          }
          else {
            Scalar * ptr = 0;
            for (int j=0; j<N; j++) {
              ptr = Aptr + j*LDA;
              Scalar s1 = C_[j];
              for (int i=0; i<M; i++) {
                *ptr = *ptr*s1*R_[i];
                ptr++;
              }
            }
          }
        }

        // Compute LU factor
        if (spd_) {
          lapack.POTRF('U', M, Aptr, LDA, &INFO);
        }
        else {
          lapack.GETRF(M, N, Aptr, LDA, &IPIV_[0], &INFO);
        }

        DMT::SyncHostToDevice(*A_);
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

      DMT::SyncDeviceToHost(*X_);

      // LAPACK overwrites RHS vector with solution vector, so copy if necessary
      if (B_ != X_)
        DMT::Assign(*X_, *B_); // Copy B to X if needed

      // Since B_ = X_, perform operations on X_.

      int M = DMT::GetNumRows(*X_); 
      int NRHS = DMT::GetNumCols(*X_);
      int LDX = DMT::GetStride(*X_); 
      Scalar * X = DMT::GetRawHostPtr(*X_);

      if (equilibrate_) 
      {
        // Apply equilibration scalings to RHS vector
        MagnitudeType * R_tmp = &R_[0];
        if (transpose && !spd_) R_tmp = &C_[0];
  
        Scalar * ptr = 0;              
        for (int j=0; j<NRHS; j++) {
          ptr = X + j*LDX;
          for (int i=0; i<M; i++) {
            *ptr = *ptr*R_tmp[i];
            ptr++;
          }
        } 
      }

      int INFO = 0;

      int LDA = DMT::GetStride(*A_);
      Scalar * Aptr = DMT::GetRawHostPtr(*A_);
      Teuchos::LAPACK<int,Scalar> lapack;

      if (spd_) {
        lapack.POTRS('U', M, NRHS, Aptr, LDA, X, LDX, &INFO);
      }
      else {
        lapack.GETRS(Teuchos::ETranspChar[TRANS_], M, NRHS, 
                     Aptr, LDA, &IPIV_[0], X, LDX, &INFO);
      }

      if (equilibrate_)
      {
        Scalar * ptr = 0;
        for (int j=0; j<NRHS; j++) {
          ptr = X + j*LDX; 
            for (int i=0; i<M; i++) {
            *ptr = *ptr*R_[i];
            ptr++;
          }
        }
      }

      // Synchronize solution vector to the device
      DMT::SyncHostToDevice(*X_);

      return INFO;
    }
    //@}

  private:

    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType MagnitudeType;
  
    std::vector<int> IPIV_;
    std::vector<MagnitudeType> R_, C_;

    using DenseSolver<Scalar, Kokkos::DualView<IST**,Kokkos::LayoutLeft>>::newMatrix_;
    using DenseSolver<Scalar, Kokkos::DualView<IST**,Kokkos::LayoutLeft>>::equilibrate_;
    using DenseSolver<Scalar, Kokkos::DualView<IST**,Kokkos::LayoutLeft>>::TRANS_;
    using DenseSolver<Scalar, Kokkos::DualView<IST**,Kokkos::LayoutLeft>>::spd_;

    using DenseSolver<Scalar, Kokkos::DualView<IST**,Kokkos::LayoutLeft>>::A_;
    using DenseSolver<Scalar, Kokkos::DualView<IST**,Kokkos::LayoutLeft>>::X_;
    using DenseSolver<Scalar, Kokkos::DualView<IST**,Kokkos::LayoutLeft>>::B_;

  };

  //! Full specialization of Belos::DenseMatTraits for Kokkos::DualView.
  //
  //TODO: It seems like all of Tpetra Details returning views 
  // (e.g. getStatic2dDualView) return a LayoutLeft view.  Should we 
  // hard-code that parameter he
  template<class Scalar>
  class DenseMatTraits<Scalar, Kokkos::DualView<typename Kokkos::Details::ArithTraits<Scalar>::val_type **,Kokkos::LayoutLeft>>{

  public:
    typedef typename Kokkos::Details::ArithTraits<Scalar>::val_type IST; //Impl Scalar Type, as used in Tpetra
    
    //@{ \name Creation methods

    /*! \brief Creates a new empty \c Kokkos::DualView<IST**,Kokkos::LayoutLeft> with no dimension.

    \return Reference-counted pointer to a new dense matrix of type \c Kokkos::DualView<IST**,Kokkos::LayoutLeft>.
    */
    static Teuchos::RCP<Kokkos::DualView<IST**,Kokkos::LayoutLeft>> Create() { 
      return Teuchos::rcp(new Kokkos::DualView<IST**,Kokkos::LayoutLeft>("BelosDenseCreate",0,0));
    }     

    /*! \brief Creates a new empty \c Kokkos::DualView<IST**,Kokkos::LayoutLeft> containing 
     *     and \c numRows rows and \c numCols columns.
     *         Will be initialized to zeros if last parameter is true.

    \return Reference-counted pointer to a new dense matrix of type \c Kokkos::DualView<IST**,Kokkos::LayoutLeft>.
    */
    static Teuchos::RCP<Kokkos::DualView<IST**,Kokkos::LayoutLeft>> 
           Create( const int numRows, const int numCols, bool initZero = true) { 
      if(initZero){
        return Teuchos::rcp(new Kokkos::DualView<IST**,Kokkos::LayoutLeft>("BelosDenseCreate2",numRows,numCols));
      }
      else {
        return Teuchos::rcp(new Kokkos::DualView<IST**,Kokkos::LayoutLeft>(Kokkos::view_alloc(Kokkos::WithoutInitializing,"BelosDenseCreate2"),numRows,numCols));
      }
    }     
   
    /*! \brief Create a new copy \c DM, possibly transposed.
   
       \return Reference-counted pointer to a new dense matrix of type \c DM.
    */
    static Teuchos::RCP<Kokkos::DualView<IST**,Kokkos::LayoutLeft>> 
           CreateCopy(const Kokkos::DualView<IST**,Kokkos::LayoutLeft> & dm, bool transpose=false)
    {
      Teuchos::RCP<Kokkos::DualView<IST**,Kokkos::LayoutLeft>> tmpCopyRCP = Teuchos::null;

      if (transpose) {
        // want tmpCopyRCP to end up as dm^H in the end. Prefer doing transpose on device
        tmpCopyRCP = Teuchos::rcp(new Kokkos::DualView<IST**,Kokkos::LayoutLeft>
                                 (Kokkos::view_alloc(Kokkos::WithoutInitializing,"BelosDenseCreateCopy"),dm.extent_int(1),dm.extent_int(0)));
        if(tmpCopyRCP->need_sync_device()) {
          // tmpCopyRCP is only up to date on the host
          kokkos_transpose(tmpCopyRCP->h_view, dm.h_view);
          tmpCopyRCP->clear_sync_state();
          tmpCopyRCP->modify_host();
        }
        else {
          kokkos_transpose(tmpCopyRCP->d_view, dm.d_view);
          tmpCopyRCP->clear_sync_state();
          tmpCopyRCP->modify_device();
        }
      }
      else {
        tmpCopyRCP = Teuchos::rcp(new Kokkos::DualView<IST**,Kokkos::LayoutLeft>
                                 (Kokkos::view_alloc(Kokkos::WithoutInitializing,"BelosDenseCreateCopy"),dm.extent_int(0),dm.extent_int(1)));
        Kokkos::deep_copy(*tmpCopyRCP, dm);
      }
      return tmpCopyRCP;
    }
 
    //TODO make const and non-const function for raw pointer. 
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
    static Scalar* GetRawHostPtr(Kokkos::DualView<IST**,Kokkos::LayoutLeft> & dm ) { 
    //LAPACK could use the host ptr to modify entries, so mark as modified.
    //TODO: Is there a better way to handle this?
        std::cout << "GetRawHostPtr: Modifying on host." << std::endl;
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
    static Scalar const * GetConstRawHostPtr(const Kokkos::DualView<IST**,Kokkos::LayoutLeft> & dm ) { 
        return reinterpret_cast<Scalar const *>(dm.view_host().data());
    }

    //! \brief Marks host data modified to avoid device sync errors. 
    /// \note Belos developers must call this function after EVERY
    ///   call to LAPACK that modifies dense matrix data accessed via raw pointer. 
    //static void RawPtrDataModified(Kokkos::DualView<IST**,Kokkos::LayoutLeft> & dm ) {
      //dm.modify_host();
    //}

    //! \brief Returns an RCP to a Kokkos::DualView<IST**,Kokkos::LayoutLeft> which has a subview of the given Kokkos::DualView<IST**,Kokkos::LayoutLeft>.
    //        Row and column indexing is zero-based.
    static Teuchos::RCP<Kokkos::DualView<IST**,Kokkos::LayoutLeft>> 
    Subview( Kokkos::DualView<IST**,Kokkos::LayoutLeft> & source, int numRows, int numCols, int startRow=0, int startCol=0){
      return Teuchos::rcp(new Kokkos::DualView<IST**,Kokkos::LayoutLeft>(source, 
          Kokkos::pair<int,int>(startRow,startRow+numRows), Kokkos::pair<int,int>(startCol,startCol+numCols))); 
      //Does subview work on dual views? Brian says yes. 
    }     

    static Teuchos::RCP<const Kokkos::DualView<IST**,Kokkos::LayoutLeft>> 
    SubviewConst( const Kokkos::DualView<IST**,Kokkos::LayoutLeft>& source, int numRows, int numCols, int startRow=0, int startCol=0){
      return Teuchos::rcp(new Kokkos::DualView<IST**,Kokkos::LayoutLeft>(source, 
          Kokkos::pair<int,int>(startRow,startRow+numRows), Kokkos::pair<int,int>(startCol,startCol+numCols))); 
          //TODO: Check const-ness semantics here?
    }     

    //! \brief Returns a deep copy of the requested subview.
    static Teuchos::RCP<Kokkos::DualView<IST**,Kokkos::LayoutLeft>> 
    SubviewCopy( const Kokkos::DualView<IST**,Kokkos::LayoutLeft>& source, int numRows, int numCols, int startRow=0, int startCol=0){ 
      //Maaybe we could get away with just a host copy here??
      //Hmmm... but it says we return a dual view. 
      //Maybe it should return a dual view with only host data copied in. Require sync to work on device. 
      //This is related to the functionality of the Assign function. 
      auto tmpViewRCP = Teuchos::rcp(new Kokkos::DualView<IST**,Kokkos::LayoutLeft>
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
    static int GetNumRows( const Kokkos::DualView<IST**,Kokkos::LayoutLeft>& dm ) { 
      return dm.extent_int(0); 
    }     

    //! \brief Obtain the number of columns of \c dm.
    static int GetNumCols( const Kokkos::DualView<IST**,Kokkos::LayoutLeft>& dm ) { 
      return dm.extent_int(1); 
    }     

    //! \brief Obtain the stride between the columns of \c dm.
    static int GetStride( const Kokkos::DualView<IST**,Kokkos::LayoutLeft>& dm ) { 
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
    static void Reshape( Kokkos::DualView<IST**,Kokkos::LayoutLeft>& dm, const int numRows, const int numCols, bool initZero = false) {
      std::cout << "Calling reshape." << std::endl;
      if(initZero){
        std::cout << "in reshape initZero." << std::endl;
        dm.realloc(numRows,numCols); //changes size of both host and device view.
        Kokkos::deep_copy(dm.view_device(), 0.0);
        dm.modify_device();
        std::cout << "Reshape: Modified on device." << std::endl;
      }
      else{
      std::cout << "in reshape keep vals." << std::endl;
        dm.resize(numRows,numCols); //keeps values in old array.
      }
    }     

    //@}

    //@{ \name Data access methods

    //! \brief Access a reference to the (i,j) entry of \c dm, \c e_i^T dm e_j.
    static Scalar & Value( Kokkos::DualView<IST**,Kokkos::LayoutLeft>& dm, const int i, const int j )
    { 
    //Mark as modified on host, since we don't know if it will be. 
        std::cout << "Value: Modifying on host." << std::endl;
      dm.modify_host();
      return reinterpret_cast<Scalar&>((dm.view_host())(i,j));
      // TODO Will this result in extra syncs? Is always marking modified the best way?
    }

    //! \brief Access a const reference to the (i,j) entry of \c dm, \c e_i^T dm e_j.
    static const Scalar & ValueConst( const Kokkos::DualView<IST**,Kokkos::LayoutLeft>& dm, const int i, const int j ) { 
      return reinterpret_cast<Scalar const &>((dm.view_host())(i,j));
      //TODO check const semantics here?
    }

    //TODO: Check formatting of ALL doxygen comments!
    //
    //! \brief If an accelorator is in use, sync it to device on this call.
    //  
    //  \note The only Belos function that results in a need to sync to 
    //  host is MvTransMv. You MUST call SyncDeviceToHost before calling
    //  any other DenseMatTraits functions after a call to MvTransMv. 
    //  All DenseMatTraits functions assume the necessary data is on host
    //  and perform computations only on the host. 
    //
    static void SyncDeviceToHost(Kokkos::DualView<IST**,Kokkos::LayoutLeft> & dm) { 
      std::cout << "Called Sync Device to Host. " << std::endl;
        std::cout << "Matrix extents are: " << dm.extent_int(0) << " , " << dm.extent_int(1) << std::endl;
      if(dm.need_sync_host()){
        if(dm.view_host().span_is_contiguous() && dm.view_device().span_is_contiguous()){
        std::cout << "Syncing d2h the easy way... " << dm.extent_int(0) << " , " << dm.extent_int(1) << std::endl;
        //Stupidness to print type info: int int1 = dm.d_view;
        dm.sync_host();}
        else{
          std::cout << "d2h sync the hard way in progress..." << std::endl;
          Kokkos::DualView<IST**,Kokkos::LayoutLeft> compat_view("compat view",dm.extent_int(0),dm.extent_int(1));
          Kokkos::deep_copy(compat_view,dm);
          compat_view.sync_host();
          Kokkos::deep_copy(dm,compat_view);
          dm.clear_sync_state();
        }
      }    
    }

    static void SyncHostToDevice(Kokkos::DualView<IST**,Kokkos::LayoutLeft> & dm) { 
      std::cout << "Called Sync Host to Device. " << std::endl;
      if(dm.need_sync_device()){
        if(dm.view_host().span_is_contiguous() && dm.view_device().span_is_contiguous()){
          std::cout << "h2d sync easy way..." << std::endl;
          dm.sync_device();
        }
        else{
          std::cout << "h2d sync the hard way in progress..." << std::endl;
          Kokkos::DualView<IST**,Kokkos::LayoutLeft> compat_view("compat view",dm.extent_int(0),dm.extent_int(1));
          Kokkos::deep_copy(compat_view,dm);
          compat_view.sync_host();
          Kokkos::deep_copy(dm,compat_view);
          dm.clear_sync_state();
        }
      }    
    }
    //@}
    //@{ \name Operator methods
    
    //!  \brief Adds sourceDM to thisDM and returns answer in thisDM.
    static void Add( Kokkos::DualView<IST**,Kokkos::LayoutLeft>& thisDM, const Kokkos::DualView<IST**,Kokkos::LayoutLeft>& sourceDM) {
      KokkosBlas::axpy(1.0,sourceDM.view_device(), thisDM.view_device()); //axpy(alpha,x,y), y = y + alpha*x
      thisDM.modify_device();
        std::cout << "Add: Modified on device." << std::endl;
    }

    //!  \brief Fill all entries with \c value. Value is zero if not specified.
    static void PutScalar( Kokkos::DualView<IST**,Kokkos::LayoutLeft>& dm, Scalar value = Teuchos::ScalarTraits<Scalar>::zero()){ 
      Kokkos::deep_copy( dm.view_device(), value);
      dm.modify_device();
        std::cout << "PutScalar: Modified on device." << std::endl;
    }

    //!  \brief Multiply all entries by a scalar. DM = value.*DM
    static void Scale( Kokkos::DualView<IST**,Kokkos::LayoutLeft>& dm, Scalar value) { 
      KokkosBlas::scal( dm.view_device(), value, dm.view_device());
      dm.modify_device();
        std::cout << "Scale: Modified on device." << std::endl;
    }

    //!  \brief Fill the Kokkos::DualView with random entries.
    //!   Entries are assumed to be the same on each MPI rank (each matrix copy). 
    static void Randomize( Kokkos::DualView<IST**,Kokkos::LayoutLeft>& dm) { 
      int rand_seed = std::rand();
      Kokkos::Random_XorShift64_Pool<> pool(rand_seed); 
      Kokkos::fill_random(dm.view_device(), pool, -1,1);
      dm.modify_device();
      std::cout << "Randomize: Modified on device." << std::endl;
    }

    //!  \brief Copies entries of source to dest (deep copy). 
    static void Assign( Kokkos::DualView<IST**,Kokkos::LayoutLeft>& dest, const Kokkos::DualView<IST**,Kokkos::LayoutLeft>& source) { 
      Kokkos::deep_copy(dest,source); //Brian Kelley says this works on dual views.
      //TODO: But do we really want to do this? What if we only need the host pieces copied right now?
      // Note: going with first solution. See discussion on SubViewCopy. 
    }

    //!  \brief Returns the Frobenius norm of the dense matrix.
    static typename Teuchos::ScalarTraits<Scalar>::magnitudeType NormFrobenius(const Kokkos::DualView<const IST**,Kokkos::LayoutLeft>& dm) { 
      using KAT = Kokkos::ArithTraits<IST>;
      using mag_t = typename KAT::mag_type;
      mag_t frobNorm;
      Kokkos::parallel_reduce(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {dm.extent(0), dm.extent(1)}),
          KOKKOS_LAMBDA(size_t i, size_t j, mag_t& lfrobNorm)
          {
          mag_t absVal = KAT::abs((dm.view_device())(i, j));
          lfrobNorm += absVal * absVal;
          }, frobNorm);
      return Kokkos::sqrt(frobNorm);
    }

    //!  \brief Returns the one-norm of the dense matrix.
    static typename Teuchos::ScalarTraits<Scalar>::magnitudeType NormOne( Kokkos::DualView<IST**,Kokkos::LayoutLeft>& dm) {  
      using KAT = Kokkos::ArithTraits<IST>;
      IST sum = 0, max_sum = 0; 
      SyncDeviceToHost(dm); //TODO Kokkos-ify this
      //Kokkos::parallel_for(dm.extent_int(1), 
          //KOKKOS_LAMBDA(
      for(int j = 0; j < dm.extent_int(1); j++){ //cols
        for(int i = 0; i < dm.extent_int(0); i++){  //rows
          sum += KAT::abs((dm.view_host())(i,j));
        }
        if(KAT::abs(sum) > KAT::abs(max_sum)){
          max_sum = sum;
        }
        sum = 0;
      }
      return KAT::abs(max_sum); //TODO: Check this impl is correct
    }
    //@}

    //@{ \name Solver methods

    //!  \brief Returns a dense solver object for the dense matrix.
    static Teuchos::RCP<DenseSolver<Scalar, Kokkos::DualView<IST**,Kokkos::LayoutLeft>>>
      createDenseSolver() {

      Teuchos::RCP<DenseSolver<Scalar, Kokkos::DualView<IST**,Kokkos::LayoutLeft>>> newSolver
        = Teuchos::rcp( new KokkosDenseSolver<Scalar>() );
      return newSolver;
    }
    //@}

  };

} // namespace Belos

#endif // end file BELOS_KOKKOS_DENSE_MAT_TRAITS_HPP
