// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file BelosKokkosAdapter.hpp
    \brief Implementation of the interface between Belos virtual classes and Kokkos concrete classes.
    Allows the user to use a Kokkos::view as a Belos::MultiVector in the Belos solvers.
    \warning Note: This class does not currently allow the user to run Belos solvers which
    require accessing non-contiguous columns of data in memory.
*/

#include<Kokkos_Core.hpp>
#include<Kokkos_Random.hpp>
#include<KokkosBlas.hpp>
#include<KokkosSparse_spmv.hpp>

#include<BelosTypes.hpp>
#include<BelosMultiVec.hpp>
#include<BelosOperator.hpp>

#include <Teuchos_SerialDenseMatrix.hpp>


#ifndef BELOS_KOKKOS_ADAPTER_HPP
#define BELOS_KOKKOS_ADAPTER_HPP
namespace Belos {

//Forward class declaration of KokkosCrsOperator:
template<class ScalarType, class OrdinalType, class Device>
class KokkosCrsOperator;

/// \class KokkosMultiVec
/// \brief Implementation of Belos::MultiVec using Kokkos::View.
/// \author Jennifer Loe
///
/// \tparam ScalarType The type of entries of the multivector.
/// \tparam Device The Kokkos::ExecutionSpace where the multivector
/// should reside.
///
/// Belos::MultiVec offers a simple abstract interface for
/// multivector operations in Belos solver algorithms.  This class
/// implements Belos::MultiVec using Kokkos::View.
template<class ScalarType, class Device = Kokkos::DefaultExecutionSpace >
class KokkosMultiVec : public MultiVec<ScalarType> {

public:

  using ViewVectorType = Kokkos::View<ScalarType*,Kokkos::LayoutLeft, Device>;
  using ConstViewVectorType = Kokkos::View<const ScalarType*,Kokkos::LayoutLeft, Device>;
  using ViewMatrixType = Kokkos::View<ScalarType**,Kokkos::LayoutLeft, Device>;
  using ConstViewMatrixType = Kokkos::View<const ScalarType**,Kokkos::LayoutLeft, Device>;

  //Unmanaged view types:
  using UMHostViewVectorType =
        Kokkos::View<ScalarType*,Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using UMHostConstViewVectorType =
        Kokkos::View<const ScalarType*,Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using UMHostViewMatrixType =
        Kokkos::View<ScalarType**,Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  using UMHostConstViewMatrixType =
        Kokkos::View<const ScalarType**,Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

protected:
  ViewMatrixType myView;

public:
  //! @name Constructors/Destructor
  //@{

  /// \brief Returns a multivector with `numrows` rows and `numvecs` columns.
  ///
  /// The `label` string indicates the label for the internal `Kokkos::view`.
  /// If `zeroOut` is set to `true`, the multivector will be initialized to zeros.
  KokkosMultiVec<ScalarType, Device> (const std::string label, const int numrows, const int numvecs, const bool zeroOut = true) :
    myView (Kokkos::view_alloc(Kokkos::WithoutInitializing,label),numrows,numvecs)
    { if (zeroOut) { Kokkos::deep_copy(myView,0); } }

  /// \brief Returns a multivector with `numrows` rows and `numvecs` columns.
  ///
  /// If `zeroOut` is set to `true`, the multivector will be initialized to zeros.
  KokkosMultiVec<ScalarType, Device> (const int numrows, const int numvecs, const bool zeroOut = true) :
    myView (Kokkos::view_alloc(Kokkos::WithoutInitializing,"MV"),numrows,numvecs)
    { if (zeroOut) { Kokkos::deep_copy(myView,0); } }

  /// \brief Returns a single column multivector with `numrows` rows.
  ///
  /// If `zeroOut` is set to `true`, the multivector will be initialized to zeros.
  KokkosMultiVec<ScalarType, Device> (const int numrows, const bool zeroOut = true) :
    myView(Kokkos::view_alloc(Kokkos::WithoutInitializing,"MV"),numrows,1)
    { if (zeroOut) { Kokkos::deep_copy(myView,0); } }

  //! Copy constructor (performs deep copy).

  /// This copy constructor returns a new KokksMultiVec containing a
  /// deep copy of the multivector given by the user.
  KokkosMultiVec<ScalarType, Device> (const KokkosMultiVec<ScalarType, Device> &sourceVec) :
    myView(Kokkos::view_alloc(Kokkos::WithoutInitializing,"MV"),(int)sourceVec.GetGlobalLength(),sourceVec.GetNumberVecs())
  { Kokkos::deep_copy(myView,sourceVec.GetInternalViewConst()); }

  //! Copy constructor for type conversion. (Performs deep copy.)

  /// This copy constructor returns a new KokksMultiVec containing a
  /// deep copy of the multivector given by the user.
  /// The internal data of the multivector is converted from
  /// `ScalarType` to `ScalarType2`.
  template < class ScalarType2 >
    KokkosMultiVec<ScalarType, Device> (const KokkosMultiVec<ScalarType2, Device> &sourceVec) :
    myView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "MV"),(int)sourceVec.GetGlobalLength(),sourceVec.GetNumberVecs())
  { Kokkos::deep_copy(myView,sourceVec.GetInternalViewConst()); }

  //! Assignment operator (performs deep copy).

  /// This `=` operator performs a deep copy of
  /// the right-hand side KokkosMultiVec to the
  /// left-hand side KokkosMultiVec. The left-hand
  /// side MultiVec will be resized if necessary.
  KokkosMultiVec<ScalarType, Device> & operator=(const KokkosMultiVec<ScalarType, Device> & sourceVec) {
    int len = sourceVec.GetGlobalLength();
    int cols = sourceVec.GetNumberVecs();
    if( len != (int)myView.extent(0) || cols != (int)myView.extent(1) ){
      Kokkos::resize(myView, len, cols);
    }
    Kokkos::deep_copy(myView,sourceVec.GetInternalViewConst());
    return *this;
  }

  //! Assignment operator for type conversion. (Performs deep copy.)

  /// This `=` operator performs a deep copy of
  /// the right-hand side KokkosMultiVec to the
  /// left-hand side KokkosMultiVec. The left-hand
  /// side MultiVec will be resized if necessary.
  /// The internal data of the right-hand side multivec
  /// is converted from `ScalarType` to `ScalarType2`.
  template < class ScalarType2 >
  KokkosMultiVec<ScalarType, Device> & operator=(const KokkosMultiVec<ScalarType2, Device> & sourceVec) {
    int len = sourceVec.GetGlobalLength();
    int cols = sourceVec.GetNumberVecs();
    if( len != (int)myView.extent(0) || cols != (int)myView.extent(1) ){
      Kokkos::resize(myView, len, cols);
    }
    Kokkos::deep_copy(myView,sourceVec.GetInternalViewConst());
    return *this;
  }

  //! Create a KokkosMultiVec from a given `Kokkos::view`.

  /// Returns a KokkosMultiVec that internally stores the
  /// data given in `sourceView`. If `makeCopy` has value
  /// `true`, then this function makes a deep copy of the
  /// `Kokkos::view`.  If `false`, then the KokkosMultiVec stores a
  /// shallow copy of the given `Kokkos::view`.  (This option assumes that
  /// the user will make no changes to that view outside of
  /// the KokkosMultiVec interface.)
  KokkosMultiVec<ScalarType, Device> (const ViewMatrixType & sourceView, bool makeCopy = true) {
    if( makeCopy ){
      if( sourceView.extent(0) != myView.extent(0) || sourceView.extent(1) != myView.extent(1) ){
        Kokkos::resize(myView, sourceView.extent(0), sourceView.extent(1));
      }
      Kokkos::deep_copy(myView, sourceView);
    }
    else{
      myView = sourceView;
    }
  }

  //! Create a KokkosMultiVec with ScalarType2 from a given `Kokkos::view` with ScalarType1.

  /// Returns a KokkosMultiVec that internally stores the
  /// data given in `sourceView`. This function always makes
  /// a deep copy of the sourceView in order to change scalar
  /// types.
  template < class ScalarType2 > //TODO: Fix this so that passing in a view without device specified actually compiles...
  KokkosMultiVec<ScalarType, Device> (const Kokkos::View<ScalarType2**,Kokkos::LayoutLeft, Device> & sourceView) :
    myView(Kokkos::view_alloc(Kokkos::WithoutInitializing, "MV"),sourceView.extent(0),sourceView.extent(1))
  { Kokkos::deep_copy(myView,sourceView); }

  //! Destructor (default)
  ~KokkosMultiVec<ScalarType, Device>(){}

  //@}

  //! @name Data Access Methods
  //@{

  //! Reader (const) method for internal view.

  /// Returns the (const) view that stores internal data
  /// for the KokkosMultiVec.
  ConstViewMatrixType GetInternalViewConst() const {
    return myView;
  }

  //! Reader/Writer method for internal view.

  /// Returns the view that stores internal data
  /// for the KokkosMultiVec.
  /// \warning ** Be careful!! ** The KokkosMultiVec does
  /// NOT expect users to use this function to
  /// make changes to its internal data. Most
  /// necessary changes can be made via other functions.
  /// Make sure you know what you are doing!
  ViewMatrixType GetInternalViewNonConst(){
    return myView;
  }
  //@}

  //! @name Copy/View functions inherited from Belos::MultiVec
  //@{

  //! Creates a MultiVec with same data distribution as `this`.

  /// A virtual "copy constructor" that returns a pointer to a new
  /// Belos::MultiVec.  (Underneath, the vector is a KokkosMultiVec,
  /// and its properties can be accessed via `dynamic_cast`.)
  /// This vector's entries are
  /// not copied; instead, a new MultiVec is created with the same
  /// data distribution, but with numvecs columns (numvecs > 0).
  /// Multivector entries are not initialized.
  ///
  /// \param numvecs [in] The number of columns in the output
  ///   multivector.  Must be positive.
  MultiVec<ScalarType> * Clone ( const int numvecs ) const{
    KokkosMultiVec<ScalarType, Device> * ptr = new KokkosMultiVec<ScalarType, Device>(myView.extent(0),numvecs, false);
    return ptr;
  }

  //! \brief Creates a MultiVec which is a (deep) copy of `this`.
  ///
  /// A virtual "copy constructor" returning a pointer to a new
  /// object of Belos::MultiVec. (KokkosMultiVec underneath.)
  /// All of this vector's entries are
  /// copied and a new stand-alone multivector is created.  (deep
  /// copy).
  MultiVec<ScalarType> * CloneCopy () const{
    KokkosMultiVec<ScalarType, Device> * ptr = new KokkosMultiVec<ScalarType, Device>(myView.extent(0),myView.extent(1), false);
    Kokkos::deep_copy(ptr->GetInternalViewNonConst(),myView);
    return ptr;
  }

  //! \brief Creates a MultiVec which is a (deep) copy of selected columns of `this`.
  ///
  /// A virtual "copy constructor" returning a pointer to a new Belos::MultiVec.
  /// This vector's entries are copied and a new
  /// stand-alone MultiVector is created where only selected columns
  /// are chosen.  (deep copy).
  /// Indexing is zero-based, so an `std::vector index` with values 0, 3, 5 would
  /// indicate copying the first, 4th, and 6th columns of the original multivector.
  /// Indices need not be contiguous or ordered.
  /// Result is `output[:,i] = (*this)[:,index[i]]`.
  MultiVec<ScalarType> * CloneCopy ( const std::vector<int>& index ) const{
    // JAL- If debug options needed, could add validity checks of index.
    // See debug code in belos/src/tpetra/BelosMultiVecTraits_Tpetra.hpp.
    int numvecs = index.size();
    KokkosMultiVec<ScalarType, Device> * B = new KokkosMultiVec<ScalarType, Device>("B",myView.extent(0),numvecs, false);
    bool isContigAscending = true;

    //Check whether the given indices are contiguous and ascending.
    for(unsigned int i=0; i< (index.size()-1); i++){
      if( index[i+1] != index[i]+1 ){
        isContigAscending = false;
      }
    }

    //Copy the vectors: (Case depends on indices.)
    if(isContigAscending && index.size()==(unsigned)this->GetNumberVecs()){ //Copy entire multivec.
      Kokkos::deep_copy(B->GetInternalViewNonConst(),myView);
    }
    else if (isContigAscending){ //Copy contiguous subset
      ViewMatrixType ThisSub = Kokkos::subview(myView, Kokkos::ALL, std::make_pair(index.front(), index.back()+1));
      Kokkos::deep_copy(B->GetInternalViewNonConst(),ThisSub);
    }
    else{ //Copy columns one by one
      for(unsigned int i=0; i<index.size(); i++){
        auto Bsub = Kokkos::subview(B->GetInternalViewNonConst(), Kokkos::ALL, i);
        auto ThisSub = Kokkos::subview(myView, Kokkos::ALL, index[i]);
        Kokkos::deep_copy(Bsub, ThisSub);
      }
    }
    return B;
  }

  //! \brief Creates a (read-only) MultiVec which is a shallow copy of selected columns of `this`.
  ///
  /// A virtual view constructor returning a pointer to a new Belos::MultiVec with
  /// selected columns. (Column indexing is zero-based.) The view is read-only.
  /// This vector's entries are shared and hence no
  /// memory is allocated for the columns. (Internally, we create
  /// a `Kokkos::subview`.)
  /// \warning At this time, the Kokkos-Belos adapter only supports
  /// viewing column indices that form a contiguous subset in memory.
  /// Thus, the values in `index` must be contiguous and ascending (e.g. 0,1,2,3).
  const MultiVec<ScalarType> * CloneView ( const std::vector<int>& index ) const { //TODO This isn't const!!
    bool isContigAscending = true;
    //Check whether the given indices are contiguous and ascending.
    for(unsigned int i=0; i< (index.size()-1); i++){
      if( index[i+1] != index[i]+1 ){
        isContigAscending = false;
      }
    }
    if(isContigAscending ){
    const KokkosMultiVec<ScalarType, Device> * B =
        new KokkosMultiVec<ScalarType, Device>(Kokkos::subview(myView, Kokkos::ALL, std::make_pair(index.front(), index.back()+1)),false);
      return B;
    }
    else{
      throw std::runtime_error("CloneView asked for non-contiguous subset. \n This feature is not yet supported in Belos for Kokkos.");
    }
  }


  //! \brief Creates a nonconst MultiVec which is a shallow copy of selected columns of `this`.
  ///
  /// A virtual view constructor returning a pointer to a new Belos::MultiVec with
  /// selected columns. (Column indexing is zero-based.)
  /// This vector's entries are shared and hence no
  /// memory is allocated for the columns. (Internally, we create
  /// a `Kokkos::subview`.) Writing to this view will change the
  /// entries of the original multivector.
  /// \warning At this time, the Kokkos-Belos adapter only supports
  /// viewing column indices that form a contiguous subset in memory.
  /// Thus, the values in `index` must be contiguous and ascending (e.g. 0,1,2,3).
  MultiVec<ScalarType> * CloneViewNonConst ( const std::vector<int>& index ){
    bool isContigAscending = true;
    //Check whether the given indices are contiguous and ascending.
    for(unsigned int i=0; i< (index.size()-1); i++){
      if( index[i+1] != index[i]+1 ){
        isContigAscending = false;
      }
    }
    if(isContigAscending ){
    KokkosMultiVec<ScalarType, Device> * B =
        new KokkosMultiVec<ScalarType, Device>(Kokkos::subview(myView, Kokkos::ALL, std::make_pair(index.front(), index.back()+1)),false);
      return B;
    }
    else{
      throw std::runtime_error("CloneViewNonConst asked for non-contiguous subset. \n This feature is not yet supported in Belos for Kokkos.");
    }
  }

  //! Copy the vectors in A to the vectors in (*this) specified by index. (Deep copy.)

  /// Copies vectors of A to a sub-block of vectors of this multivector. (Deep copy.)
  /// The sub-block to be overwritten is given by the indices and need not be contiguous.
  /// Result is (*this)[:,index[i]] = A[:,i].
  /// Column indexing is zero-based.
  ///
  void SetBlock ( const MultiVec<ScalarType>& A, const std::vector<int>& index ){
    KokkosMultiVec<ScalarType, Device> *A_vec = dynamic_cast<KokkosMultiVec *>(&const_cast<MultiVec<ScalarType> &>(A));

    if( index.size() > myView.extent(1) ){
      throw std::runtime_error("Error in KokkosMultiVec::SetBlock. A cannot have more vectors than (*this).");
    }
    bool isContigAscending = true;
    //Check whether the given indices are contiguous and ascending.
    for(unsigned int i=0; i< (index.size()-1); i++){
      if( index[i+1] != index[i]+1 ){
        isContigAscending = false;
      }
    }

    //Perform deep copy of sub block:
    if(isContigAscending && index.size()==(unsigned)this->GetNumberVecs()){ //Copy entire multivec.
      Kokkos::deep_copy(myView,A_vec->GetInternalViewConst());
    }
    else if (isContigAscending){ //Copy contiguous subset
      ConstViewMatrixType Asub = Kokkos::subview(A_vec->GetInternalViewConst(), Kokkos::ALL, std::make_pair(0,(int)index.size()));
      ViewMatrixType ThisSub = Kokkos::subview(myView, Kokkos::ALL, std::make_pair(index.front(), index.back()+1));
      Kokkos::deep_copy(ThisSub, Asub);
    }
    else{ //Copy columns one by one
      for(unsigned int i=0; i<index.size(); i++){
        ConstViewVectorType Asub = Kokkos::subview(A_vec->GetInternalViewConst(), Kokkos::ALL, i);
        ViewVectorType ThisSub = Kokkos::subview(myView, Kokkos::ALL, index[i]);
        Kokkos::deep_copy(ThisSub, Asub);
      }
    }
  }
  //@}

  //! @name Attribue functions inherited from Belos::MultiVec
  //@{
  //! Returns the number of rows in the multivector.
  ptrdiff_t GetGlobalLength () const {
    return static_cast<ptrdiff_t>(myView.extent(0));
  }

  //! Returns the number of columns in the multivector.
  int GetNumberVecs () const { return myView.extent(1); }
  //@}

  //! @name Mathematical functions inherited from Belos::MultiVec
  //@{
  //! \brief `*this <- alpha * A * B + beta * (*this)`
  ///
  /// where alpha and beta are scalars and the dimensions of `A*B` match
  /// the dimensions of `(*this)`.
  void MvTimesMatAddMv ( const ScalarType alpha, const MultiVec<ScalarType>& A,
                         const Teuchos::SerialDenseMatrix<int,ScalarType>& B, const ScalarType beta ){
    KokkosMultiVec<ScalarType, Device> *A_vec = dynamic_cast<KokkosMultiVec *>(&const_cast<MultiVec<ScalarType> &>(A));
    if( myView.extent(1) == 1 && A_vec->GetInternalViewConst().extent(1) == 1){ //B is a scalar.
      ScalarType scal1 = alpha*B(0,0);
      ViewVectorType mysub = Kokkos::subview(myView, Kokkos::ALL, 0);
      ConstViewVectorType Asub = Kokkos::subview(A_vec->GetInternalViewConst(), Kokkos::ALL, 0);
      KokkosBlas::axpby(scal1, Asub, beta, mysub);
    }
    else{
      UMHostConstViewMatrixType mat_h(B.values(), A_vec->GetInternalViewConst().extent(1), myView.extent(1));
      ViewMatrixType mat_d(Kokkos::view_alloc(Kokkos::WithoutInitializing,"mat"), A_vec->GetInternalViewConst().extent(1), myView.extent(1));
      Kokkos::deep_copy(mat_d, mat_h);
      if( myView.extent(1) == 1 ){ // B has only 1 col
          ConstViewVectorType Bsub = Kokkos::subview(mat_d, Kokkos::ALL, 0);
          ViewVectorType mysub = Kokkos::subview(myView, Kokkos::ALL, 0);
          KokkosBlas::gemv("N", alpha, A_vec->GetInternalViewConst(), Bsub, beta, mysub);
      }
      else{
        KokkosBlas::gemm("N", "N", alpha, A_vec->GetInternalViewConst(), mat_d, beta, myView);
      }
    }
  }

  //! `*this <- alpha * A + beta * B`
  ///
  /// Scale and add two vectors.  Store the result in `*this`.
  void MvAddMv ( const ScalarType alpha, const MultiVec<ScalarType>& A, const ScalarType beta,
                 const MultiVec<ScalarType>& B){
    KokkosMultiVec<ScalarType, Device> *A_vec = dynamic_cast<KokkosMultiVec *>(&const_cast<MultiVec<ScalarType> &>(A));
    KokkosMultiVec<ScalarType, Device> *B_vec = dynamic_cast<KokkosMultiVec *>(&const_cast<MultiVec<ScalarType> &>(B));

    KokkosBlas::update(alpha, A_vec->GetInternalViewConst(), beta, B_vec->GetInternalViewConst(), (ScalarType) 0.0, myView);
  }

  /// `*this <- alpha * this`
  ///
  //! Scale (multiply) each element of the vectors in \c *this with \c alpha.
  void MvScale ( const ScalarType alpha ) {
    KokkosBlas::scal(myView, alpha, myView);
  }

  /// `*this[:,i] <- alpha[i] * this[:,i]`
  ///
  //! Scale (multiply) each element of the \c i-th vector in \c *this with \c alpha[i].
  void MvScale ( const std::vector<ScalarType>& alpha ){

    //Move scalar values to a Kokkos View:
    UMHostConstViewVectorType scalars_h(alpha.data(), alpha.size());
    ViewVectorType scalars_d(Kokkos::view_alloc(Kokkos::WithoutInitializing,"scalars_d"), alpha.size());
    Kokkos::deep_copy(scalars_d, scalars_h);

    KokkosBlas::scal(myView, scalars_d, myView);
  }

  //! B <- alpha * A^T * (*this)
  ///
  /// Computes matrix product with transpose.  Result is a dense matrix.
  /// Conjugate transpose is used as appropriate.
  void MvTransMv ( const ScalarType alpha, const MultiVec<ScalarType>& A, Teuchos::SerialDenseMatrix<int,ScalarType>& B ) const{
    KokkosMultiVec<ScalarType, Device> *A_vec = dynamic_cast<KokkosMultiVec *>(&const_cast<MultiVec<ScalarType> &>(A));
    if(A_vec->myView.extent(1) == 1 && myView.extent(1) == 1){
      ConstViewVectorType Asub = Kokkos::subview(A_vec->GetInternalViewConst(), Kokkos::ALL, 0);
      ViewVectorType mysub = Kokkos::subview(myView, Kokkos::ALL, 0);
      ScalarType soln = KokkosBlas::dot(Asub, mysub);
      soln = alpha*soln;
      B(0,0) = soln;
    }
   // ***
   // For MvTransMv, this option runs slower than GEMM on NVIDIA V100.
   // Do not enable for now.
   // ****
   // else if( myView.extent(1) == 1 ){ // Only 1 col in soln vec
   //   ViewVectorType soln(Kokkos::view_alloc(Kokkos::WithoutInitializing,"soln"), A_vec->GetInternalViewConst().extent(1));
   //   ViewVectorType mysub = Kokkos::subview(myView, Kokkos::ALL, 0);
   //   KokkosBlas::gemv("C", alpha, A_vec->GetInternalViewConst(), mysub, ScalarType(0.0), soln);
   //   for( unsigned int i = 0; i < soln.extent(0); i++){
   //     B(i,0) = soln(i);
   //   }
   // }
    else{
      UMHostViewMatrixType soln_h(B.values(), A_vec->GetInternalViewConst().extent(1), myView.extent(1));
      ViewMatrixType soln_d(Kokkos::view_alloc(Kokkos::WithoutInitializing,"mat"), A_vec->GetInternalViewConst().extent(1), myView.extent(1));
      KokkosBlas::gemm("C", "N", alpha, A_vec->GetInternalViewConst(), myView, ScalarType(0.0), soln_d);
      Kokkos::deep_copy(soln_h, soln_d);
    }
  }


  //! b[i] = A[i]^T * this[i]

  /// Performs a dot product between A and (*this).
  /// Uses conjugate transpose when appropriate.
  /// Output is a vector.
  void MvDot ( const MultiVec<ScalarType>& A, std::vector<ScalarType>& b ) const{
    //Put output vector in unmanaged Kokkos view:
    UMHostViewVectorType dotView_h(b.data(),myView.extent(1));
    ViewVectorType dotView_d(Kokkos::view_alloc(Kokkos::WithoutInitializing,"Dot"),myView.extent(1));

    KokkosMultiVec<ScalarType, Device> *A_vec = dynamic_cast<KokkosMultiVec *>(&const_cast<MultiVec<ScalarType> &>(A));

    KokkosBlas::dot(dotView_d, A_vec->GetInternalViewConst(), myView);
    Kokkos::deep_copy(dotView_h, dotView_d);
  }

  //! alpha[i] = norm of i-th column of (*this)

  /// Valid norm types are Belos::TwoNorm, Belos::OneNorm,
  /// and Belos::InfNorm.
  void MvNorm ( std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>& normvec, NormType norm_type = TwoNorm ) const{

    //Put output vector in unmanaged Kokkos view:
    using magnitudeType = typename Teuchos::ScalarTraits<ScalarType>::magnitudeType;
    Kokkos::View<magnitudeType*,Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> normView_h(normvec.data(),myView.extent(1));
    Kokkos::View<magnitudeType*,Kokkos::LayoutLeft, Device> normView_d(Kokkos::view_alloc(Kokkos::WithoutInitializing,"Norm"),myView.extent(1));

    switch( norm_type ) {
      case ( OneNorm ) :
        KokkosBlas::nrm1(normView_d, myView);
        break;
      case ( TwoNorm ) :
        KokkosBlas::nrm2(normView_d, myView);
        break;
      case ( InfNorm ) :
        KokkosBlas::nrminf(normView_d, myView);
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
            "Belos::KokkosMultiVec::MvNorm: Invalid norm_type "
            << norm_type << ".  The current list of valid norm "
            "types is {OneNorm, TwoNorm, InfNorm}.");
    }
    Kokkos::deep_copy(normView_h, normView_d);
  }
  //@}

  //! @name Initialization functions inherited from Belos::MultiVec
  //@{
  //! Fill all columns of *this with random values.
  void MvRandom() {
    int rand_seed = std::rand();
    Kokkos::Random_XorShift64_Pool<> pool(rand_seed);
    Kokkos::fill_random(myView, pool, -1,1);
  }

  //! Initialize each element of (*this) to the scalar value alpha.
  void MvInit ( const ScalarType alpha ) {
     Kokkos::deep_copy(myView,alpha);
  }
  //@}

  //! @name Print function inherited from Belos::MultiVec
  //@{
  //! Print (*this) to the given output stream.

  /// (This function will first copy the multivector to host space
  /// if needed.)
  void MvPrint( std::ostream& os ) const {
    typename ViewMatrixType::HostMirror hostView("myViewMirror", myView.extent(0), myView.extent(1));
    Kokkos::deep_copy(hostView, myView);
    for(unsigned int i = 0; i < (hostView.extent(0)); i++){
      for (unsigned int j = 0; j < (hostView.extent(1)); j++){
        os << hostView(i , j) << "  ";
      }
      os << std::endl;
    }
    os << std::endl;
  }
  //@}

};

/// \class KokkosCrsOperator
/// \brief Implementation of Belos::Operator using KokkosSparse::CrsMatrix.
template<class ScalarType, class OrdinalType=int, class Device=Kokkos::DefaultExecutionSpace>
class KokkosCrsOperator : public Operator<ScalarType> {

private:
  // Shallow copy of the CrsMatrix used for SpMV.
  KokkosSparse::CrsMatrix<ScalarType, OrdinalType, Device> myMatrix;

public:
  //! @name Constructor/Destructor
  //@{
  //! Constructor obtains a shallow copy of the given CrsMatrix.
  KokkosCrsOperator<ScalarType, OrdinalType, Device> (const KokkosSparse::CrsMatrix<ScalarType, OrdinalType, Device> mat)
   : myMatrix(mat) {}

  //! Destructor.
  ~KokkosCrsOperator<ScalarType, OrdinalType, Device>(){}
  //@}

  //! @name Methods relating to applying the operator
  //@{

  /// \brief Apply the operator to x, putting the result in y.
  ///
  /// Take the KokkosMultiVec \c x and apply the operator (or its
  /// transpose or Hermitian transpose) to it, writing the result
  /// into the KokkosMultiVec \c y. (This performs a sparse matrix-
  /// vector product.)
  ///
  /// \param x [in] The input multivector.
  ///
  /// \param y [out] The output multivector.  x and y may not alias
  ///   (i.e., be views of) one another.
  ///
  /// \param trans [in] Whether to apply the operator (Belos::NOTRANS), its
  ///   transpose (Belos::TRANS), or its Hermitian transpose (Belos::CONJTRANS).
  ///   The default is Belos::NOTRANS. (Defined in BelosTypes.hpp.)
  ///
  void Apply (const MultiVec<ScalarType>& x,  MultiVec<ScalarType>& y,  ETrans trans=NOTRANS) const{

    // Determine transpose mode:
    char mode[] = "X";
    switch(trans){
      case NOTRANS:
        mode[0]='N';
        break;
      case TRANS:
        mode[0]='T';
        break;
      case CONJTRANS:
        mode[0]='C';
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
            "Belos::KokkosCrsOperator::Apply: Invalid ETrans type ");
    }

    //Use dynamic_cast to tell the compiler these are Kokkos Multivecs.
    KokkosMultiVec<ScalarType, Device> *x_vec =
            dynamic_cast<KokkosMultiVec<ScalarType, Device> *>(&const_cast<MultiVec<ScalarType> &>(x));
    KokkosMultiVec<ScalarType, Device> *y_vec = dynamic_cast<KokkosMultiVec<ScalarType, Device> *>(&y);

    //KokkosSparse::spmv computes y = beta*y + alpha*Op(A)*x
    ScalarType alpha = 1.0;
    ScalarType beta = 0;
    KokkosSparse::spmv(mode, alpha, myMatrix, x_vec->GetInternalViewConst(), beta, y_vec->GetInternalViewNonConst());
  }

  /// \brief Whether this operator implements applying the transpose.
  ///
  /// This function returns true since we can always apply the transpose
  /// of a Kokkos::CrsMatrix.
  ///
  bool HasApplyTranspose () const {
    return true;
  }
  //@}
};
}// end namespace Belos
#endif
