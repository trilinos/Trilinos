/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Details_Behavior.hpp"
//#include "Teuchos_SerialDenseMatrix.hpp"
#include "KokkosBlas.hpp"
#include <type_traits>

namespace { // (anonymous)

template<class MapType>
Teuchos::RCP<const MapType>
makeLocalMap (const MapType& inputMap,
	      const typename MapType::local_ordinal_type lclNumRows)
{
  Teuchos::RCP<MapType> map (new MapType (lclNumRows,
					  inputMap.getIndexBase (),
					  inputMap.getComm (),
					  Tpetra::LocallyReplicated));
  return Teuchos::rcp_const_cast<const MapType> (map);
}

template<class DualViewType, class SizeType>
DualViewType
makeUninitializedDualView (const SizeType numRows, const SizeType numCols)
{
  using Kokkos::view_alloc;
  using Kokkos::WithoutInitializing;
  using dual_view_type = DualViewType;
  using dev_view_type = typename dual_view_type::t_dev;

  // This needs to be a string and not a char*, if given as an
  // argument to Kokkos::view_alloc.  This is because view_alloc also
  // allows a raw pointer as its first argument.  See
  // https://github.com/kokkos/kokkos/issues/434.
  const std::string label ("MV::DualView");

  // NOTE (mfh 18 Feb 2015, 12 Apr 2015, 22 Sep 2016, 12 Mar 2019)
  // Separate creation of the DualView's Views works around
  // Kokkos::DualView's current inability to accept an
  // AllocationProperties initial argument (as Kokkos::View does).
  // However, the work-around is harmless, since it does what the
  // (currently nonexistent) equivalent DualView constructor would
  // have done anyway.
  dev_view_type d_view (view_alloc (label, WithoutInitializing),
			numRows, numCols);

  if (::Tpetra::Details::Behavior::debug ()) {
    // Filling with NaN is a cheap and effective way to tell if
    // downstream code is trying to use a MultiVector's data without
    // them having been initialized.  ArithTraits lets us call nan()
    // even if the scalar type doesn't define it; it just returns some
    // undefined value in the latter case.
    using IST = typename DualViewType::non_const_value_type;
    const IST nan = Kokkos::ArithTraits<IST>::nan ();
    KokkosBlas::fill (d_view, nan);
  }
  auto h_view = Kokkos::create_mirror_view (d_view);

  dual_view_type dv (d_view, h_view);
  // Whether or not the user cares about the initial contents of the
  // MultiVector, the device and host views are out of sync.  We
  // prefer to work in device memory.  The way to ensure this happens
  // is to mark the device view as modified.
  dv.modify_device ();
  return dv;
}

// NOTE (mfh 14 Mar 2019) Complication avoids "you let a View persist
// after Kokkos::finalize" run-time warnings.  Note also that
// assigning an empty DualView to an existing DualView won't actually
// deallocate the modified flags View.
template<class DualViewType, class SizeType>
class DualViewPool {
private:
  DualViewType* dualViewPtr_ = nullptr;

  static DualViewType
  makeDualView (const SizeType numRows, const SizeType numCols)
  {
    // Motivating use cases for these initial sizes:
    //
    // 1. GMRES (need restart length number of rows)
    // 2. Single reduce CG (need 2 x 2)
    constexpr SizeType initNumRows = 30;
    constexpr SizeType initNumCols = 2;
    return makeUninitializedDualView<DualViewType, SizeType>
      (std::max (initNumRows, numRows),
       std::max (initNumCols, numCols));
  }

public:
  DualViewPool () : dualViewPtr_ (new DualViewType ()) {
    DualViewType* dvp = dualViewPtr_; // Avoid 'this' capture concerns.
    // Don't let DualView instances persist after Kokkos::finalize.
    Kokkos::push_finalize_hook ([=] () { delete dvp; });
  }

  DualViewPool (const DualViewPool&) = delete;
  DualViewPool& operator= (const DualViewPool&) = delete;
  ~DualViewPool () = default;

  DualViewType get (const SizeType numRows, const SizeType numCols)
  {
    TEUCHOS_ASSERT( dualViewPtr_ != nullptr );

    if (SizeType (dualViewPtr_->extent (0)) < numRows ||
	SizeType (dualViewPtr_->extent (1)) < numCols) {
      *dualViewPtr_ = DualViewType ();
      *dualViewPtr_ = makeDualView (numRows, numCols);
    }
    TEUCHOS_ASSERT( SizeType (dualViewPtr_->extent (0)) >= numRows );
    TEUCHOS_ASSERT( SizeType (dualViewPtr_->extent (1)) >= numCols );

    using pair_type = std::pair<SizeType, SizeType>;
    return Kokkos::subview (*dualViewPtr_, pair_type (0, numRows),
			    pair_type (0, numCols));
  }
};

template<class MultiVectorType>
MultiVectorType
makeLocalMultiVector (const MultiVectorType& gblMv,
		      const typename MultiVectorType::local_ordinal_type lclNumRows,
		      const typename MultiVectorType::local_ordinal_type numCols)
{
  using dual_view_type = typename MultiVectorType::dual_view_type;
  using LO = typename MultiVectorType::local_ordinal_type;

  static DualViewPool<dual_view_type, LO> dualViewPool;
  dual_view_type dv = dualViewPool.get (lclNumRows, numCols);

  TEUCHOS_ASSERT( LO (dv.extent (0)) == lclNumRows );
  TEUCHOS_ASSERT( LO (dv.extent (1)) == numCols );
  TEUCHOS_ASSERT( LO (dv.d_view.extent (0)) == lclNumRows );
  TEUCHOS_ASSERT( LO (dv.d_view.extent (1)) == numCols );
  TEUCHOS_ASSERT( LO (dv.h_view.extent (0)) == lclNumRows );
  TEUCHOS_ASSERT( LO (dv.h_view.extent (1)) == numCols );

  auto lclMap = makeLocalMap (* (gblMv.getMap ()), lclNumRows);
  TEUCHOS_ASSERT( LO (lclMap->getNodeNumElements ()) == lclNumRows );
  MultiVectorType X (lclMap, dv);
  TEUCHOS_ASSERT( LO (X.getLocalLength ()) == lclNumRows );
  TEUCHOS_ASSERT( LO (X.getNumVectors ()) == numCols );
  TEUCHOS_ASSERT( X.isConstantStride () );

  if (LO (X.getStride ()) != LO (dv.d_view.stride (1)) ||
      LO (X.getStride ()) != LO (dv.h_view.stride (1))) {
    using std::endl;
    std::ostringstream os;
    os << "*** makeLocalMultiVector:" << endl
       << "X: {lclNumRows: " << X.getLocalLength ()
       << ", numCols: " << X.getNumVectors ()
       << ", stride: " << X.getStride ()
       << "}" << endl
       << "dv.d_view.stride(0): " << dv.d_view.stride (0) << endl
       << "dv.d_view.stride(1): " << dv.d_view.stride (1) << endl
       << "dv.h_view.stride(0): " << dv.h_view.stride (0) << endl
       << "dv.h_view.stride(1): " << dv.h_view.stride (1) << endl;
    std::cerr << os.str ();
    MPI_Abort (MPI_COMM_WORLD, -1);
  }
  TEUCHOS_ASSERT( LO (X.getStride ()) == LO (dv.d_view.stride (1)) );
  TEUCHOS_ASSERT( LO (X.getStride ()) == LO (dv.h_view.stride (1)) );
  return X;
}

#if 0

template<class MultiVectorType>
void
copyMultiVectorToSerialDenseMatrix
(Teuchos::SerialDenseMatrix<int, typename MultiVectorType::scalar_type>& X_out,
 const MultiVectorType& X_in)
{
  using IST = typename MultiVectorType::impl_scalar_type;
  using output_view_type =
    Kokkos::View<IST**, Kokkos::LayoutLeft, Kokkos::HostSpace,
		 Kokkos::MemoryUnmanaged>;
  using LO = typename MultiVectorType::local_ordinal_type;
  using pair_type = std::pair<LO, LO>;

  const LO numRows = static_cast<LO> (X_out.numRows ());
  const LO numCols = static_cast<LO> (X_out.numCols ());
  TEUCHOS_ASSERT( numRows == LO (X_in.getLocalLength ()) );
  TEUCHOS_ASSERT( numCols == LO (X_in.getNumVectors ()) );

  output_view_type X_out_orig (reinterpret_cast<IST*> (X_out.values ()),
			       X_out.stride (), numCols);
  auto X_out_lcl =
    Kokkos::subview (X_out_orig, pair_type (0, numRows), Kokkos::ALL ());
  if (X_in.need_sync_device ()) {
    Kokkos::deep_copy (X_out_lcl, X_in.getLocalViewHost ());
  }
  else {
    Kokkos::deep_copy (X_out_lcl, X_in.getLocalViewDevice ());
  }
}

template<class MultiVectorType>
void
copySerialDenseMatrixToMultiVector
(MultiVectorType& X_out,
 const Teuchos::SerialDenseMatrix<int, typename MultiVectorType::scalar_type>& X_in)
{
  using IST = typename MultiVectorType::impl_scalar_type;
  using input_view_type =
    Kokkos::View<const IST**, Kokkos::LayoutLeft, Kokkos::HostSpace,
		 Kokkos::MemoryUnmanaged>;
  using LO = typename MultiVectorType::local_ordinal_type;
  using pair_type = std::pair<LO, LO>;

  const LO numRows = static_cast<LO> (X_in.numRows ());
  const LO numCols = static_cast<LO> (X_in.numCols ());
  TEUCHOS_ASSERT( numRows == LO (X_out.getLocalLength ()) );
  TEUCHOS_ASSERT( numCols == LO (X_out.getNumVectors ()) );

  input_view_type X_in_orig
    (reinterpret_cast<const IST*> (X_in.values ()),
     X_in.stride (), numCols);
  auto X_in_lcl =
    Kokkos::subview (X_in_orig, pair_type (0, numRows),
		     pair_type (0, numCols));
  if (X_out.need_sync_device ()) {
    X_out.modify_host ();
    Kokkos::deep_copy (X_out.getLocalViewHost (), X_in_lcl);
    // Make sure that the output MV is up-to-date on device.
    X_out.sync_device ();
  }
  else {
    X_out.modify_device ();
    Kokkos::deep_copy (X_out.getLocalViewDevice (), X_in_lcl);
  }
}

#endif // 0

template<class MV>
bool multiVectorsLocallyEqual (const MV& X, const MV& Y)
{
  const size_t lclNumRows = X.getLocalLength ();
  if (Y.getLocalLength () != lclNumRows) {
    return false;
  }
  const size_t numCols = X.getNumVectors ();
  if (Y.getNumVectors () != numCols) {
    return false;
  }

  // We don't want to change the user's sync state, so if we find
  // ourselves needing to sync to host, we make a deep copy first.
  MV X_copy;
  MV Y_copy;
  if (X.need_sync_host ()) {
    X_copy = MV (X, Teuchos::Copy);
    X_copy.sync_host ();
  }
  else {
    X_copy = X; // harmless shallow copy
  }
  if (Y.need_sync_host ()) {
    Y_copy = MV (Y, Teuchos::Copy);
    Y_copy.sync_host ();
  }
  else {
    Y_copy = Y; // harmless shallow copy
  }

  // Comparing a vector at a time avoids issues with noncontiguous MVs.
  for (size_t j = 0; j < numCols; ++j) {
    auto X_j = X_copy.getVector (j);
    auto Y_j = Y_copy.getVector (j);
    auto X_j_lcl_2d = X_j->getLocalViewHost ();
    auto Y_j_lcl_2d = Y_j->getLocalViewHost ();
    auto X_j_lcl = Kokkos::subview (X_j_lcl_2d, Kokkos::ALL (), 0);
    auto Y_j_lcl = Kokkos::subview (Y_j_lcl_2d, Kokkos::ALL (), 0);

    for (size_t i = 0; i < lclNumRows; ++i) {
      if (X_j_lcl(i) != Y_j_lcl(i)) {
	return false;
      }
    }
  }

  return true;
}

template<class MV>
bool multiVectorsEqual (const MV& X, const MV& Y)
{
  const auto& X_map = * (X.getMap ());
  const auto& Y_map = * (Y.getMap ());
  if (! X_map.isSameAs (Y_map)) {
    return false;
  }

  const bool lclEqual = multiVectorsLocallyEqual (X, Y);
  const int lclEq = lclEqual ? 1 : 0;
  int gblEq = 0;

  using Teuchos::outArg;  
  using Teuchos::reduceAll;
  using Teuchos::REDUCE_MIN;
  const auto& comm = * (X_map.getComm ());
  reduceAll<int, int> (comm, REDUCE_MIN, lclEq, outArg (gblEq));

  return gblEq == 1;
}

template<class MV>
void
XH_times_Y (const typename MV::scalar_type& alpha, const MV& X, const MV& Y,
	    MV& Z)
{
  using LO = typename MV::local_ordinal_type;
  using IST = typename MV::impl_scalar_type;
  using KAT = Kokkos::ArithTraits<IST>;
  using ST = typename MV::scalar_type;  
  using STS = Teuchos::ScalarTraits<ST>;

  // beta == 0 always here

  if (alpha == STS::zero ()) {
    //if (beta == STS::zero ()) {
      Z.putScalar (STS::zero ());
    // }
    // else {
    //   Z.scale (beta);
    // }
  }
  else {
    const IST alpha_IST = static_cast<IST> (alpha);
    // const IST beta_IST = static_cast<IST> (beta);    

    // Make everything constant stride.
    MV X2 (X, Teuchos::Copy);
    MV Y2 (Y, Teuchos::Copy);
    MV Z2 (Z, Teuchos::Copy);
    
    X2.sync_host ();
    Y2.sync_host ();
    Z2.sync_host ();
    Z2.modify_host ();

    auto X_lcl = X2.getLocalViewHost ();
    auto Y_lcl = Y2.getLocalViewHost ();
    auto Z_lcl = Z2.getLocalViewHost ();
    const LO numTerms = X_lcl.extent (0);

    for (LO j = 0; j < Z_lcl.extent (1); ++j) {
      for (LO i = 0; i < Z_lcl.extent (0); ++i) {
	// IST tmp = (beta_IST == KAT::zero ()) ? KAT::zero () :
	//   beta_IST * Z_lcl(i,j);
	IST tmp = KAT::zero ();
	for (LO k = 0; k < numTerms; ++k) {
	  tmp += KAT::conj (X_lcl(k,i)) * Y_lcl(k,j);
	}
	Z_lcl(i,j) = tmp;	
      }
    }      
    
    Z2.sync_device ();
    Z.sync_device ();
    Z.modify_device ();
    Tpetra::deep_copy (Z, Z2);
    Z.reduce ();
  }
}

template<class MV>
void
XH_times_Y_2 (const typename MV::scalar_type& alpha, const MV& X, const MV& Y,
	      MV& Z)
{
  using LO = typename MV::local_ordinal_type;
  using IST = typename MV::impl_scalar_type;
  using KAT = Kokkos::ArithTraits<IST>;
  using ST = typename MV::scalar_type;  
  using STS = Teuchos::ScalarTraits<ST>;

  // beta == 0 always here

  if (alpha == STS::zero ()) {
    //if (beta == STS::zero ()) {
      Z.putScalar (STS::zero ());
    // }
    // else {
    //   Z.scale (beta);
    // }
  }
  else {
    const IST alpha_IST = static_cast<IST> (alpha);
    // const IST beta_IST = static_cast<IST> (beta);    

    // Make everything constant stride.
    MV X2 (X, Teuchos::Copy);
    MV Y2 (Y, Teuchos::Copy);
    MV Z2 (Z, Teuchos::Copy);
    
    X2.sync_host ();
    Y2.sync_host ();
    Z2.sync_host ();
    Z2.modify_host ();

    auto X_lcl = X2.getLocalViewHost ();
    auto Y_lcl = Y2.getLocalViewHost ();
    auto Z_lcl = Z2.getLocalViewHost ();
    const LO numTerms = X_lcl.extent (0);

    for (LO j = 0; j < Z_lcl.extent (1); ++j) {
      for (LO i = 0; i < Z_lcl.extent (0); ++i) {
	// IST tmp = (beta_IST == KAT::zero ()) ? KAT::zero () :
	//   beta_IST * Z_lcl(i,j);
	IST tmp = KAT::zero ();
	for (LO k = 0; k < numTerms; ++k) {
	  tmp += KAT::conj (X_lcl(k,i)) * Y_lcl(k,j);
	}
	Z_lcl(i,j) = tmp;	
      }
    }      
    
    Z2.sync_device ();
    Z2.modify_device ();
    Z2.reduce ();
    
    Z.sync_device ();
    Z.modify_device ();
    Tpetra::deep_copy (Z, Z2);
  }
}

//
// UNIT TESTS
//

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( MultiVector, DualViewPool, Scalar, LocalOrdinal, GlobalOrdinal, Node )
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;
  using SC = Scalar;
  using LO = LocalOrdinal;
  using GO = GlobalOrdinal;
  using NT = Node;
  using MV = Tpetra::MultiVector<SC, LO, GO, NT>;
  using map_type = typename MV::map_type;
  using STS = Teuchos::ScalarTraits<SC>;
  using dual_view_type = typename MV::dual_view_type;
  using pair_type = std::pair<size_t, size_t>;

  out << "Test Tpetra::MultiVector::multiply with DualViewPool" << endl;
  Teuchos::OSTab tab1 (out);

  const auto comm = Tpetra::TestingUtilities::getDefaultComm ();
  constexpr LO lclNumRows = 5;
  constexpr LO numCols = 2;
  constexpr GO indexBase = 0;
  const GO gblNumRows = GO (comm->getSize ()) * GO (lclNumRows);

  RCP<const map_type> map =
    rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));

  MV X (map, numCols);
  MV Y (map, numCols);

  const size_t X_numCols = X.getNumVectors ();
  const size_t Y_numCols = Y.getNumVectors ();

  // Fill each column of X and Y with a different number.  Start with
  // 1 instead of 0, because 0 is the default fill value.
  {
    SC curVal = STS::one ();
    X.sync_host ();
    X.modify_host ();
    for (size_t j = 0; j < X_numCols; ++j) {
      X.getVectorNonConst (j)->putScalar (curVal);
      curVal += STS::one ();
    }
    X.sync_device ();
  }
  {
    SC curVal = STS::one ();
    Y.sync_host ();
    Y.modify_host ();
    for (size_t j = 0; j < Y_numCols; ++j) {
      Y.getVectorNonConst (j)->putScalar (SC (10.0) * curVal);
      curVal += STS::one ();
    }
    Y.sync_device ();
  }
  
  const auto lclMap = makeLocalMap (*map, X_numCols);

  MV Z1 (lclMap, Y_numCols);
  Z1.multiply (Teuchos::CONJ_TRANS, Teuchos::NO_TRANS,
	       STS::one (), X, Y, STS::zero ());
  {
    MV Z1a (lclMap, Y_numCols);
    XH_times_Y (STS::one (), X, Y, Z1a);
    
    const bool Z1_Z1a_equal = multiVectorsEqual (Z1, Z1a);
    TEST_ASSERT( Z1_Z1a_equal );

    Z1a.putScalar (-STS::one ()); // flag value
    XH_times_Y (STS::one (), X, Y, Z1a);    
    const bool Z1_Z1a_again_equal = multiVectorsEqual (Z1, Z1a);
    TEST_ASSERT( Z1_Z1a_again_equal );

    Z1a.putScalar (-STS::one ()); // flag value
    XH_times_Y_2 (STS::one (), X, Y, Z1a);    
    const bool Z1_Z1a_again_equal2 = multiVectorsEqual (Z1, Z1a);
    TEST_ASSERT( Z1_Z1a_again_equal2 );
    
    MV Z2 = makeLocalMultiVector (X, X_numCols, Y_numCols);
    Z2.multiply (Teuchos::CONJ_TRANS, Teuchos::NO_TRANS,
		 STS::one (), X, Y, STS::zero ());

    const bool Z1_Z2_equal = multiVectorsEqual (Z1, Z2);
    TEST_ASSERT( Z1_Z2_equal );

    Z2.putScalar (-STS::one ()); // flag value
    XH_times_Y (STS::one (), X, Y, Z2);    
    const bool Z1_Z2_again_equal = multiVectorsEqual (Z1, Z2);
    TEST_ASSERT( Z1_Z2_again_equal );

    Z2.putScalar (-STS::one ()); // flag value
    XH_times_Y_2 (STS::one (), X, Y, Z2);    
    const bool Z1_Z2_again_equal2 = multiVectorsEqual (Z1, Z2);
    TEST_ASSERT( Z1_Z2_again_equal2 );
  }

  {
    MV Z3 = makeLocalMultiVector (X, X_numCols, Y_numCols);
    TEST_ASSERT( Z3.isConstantStride () );

    auto Z3_lcl_d = Z3.getLocalViewDevice ();
    auto Z3_lcl_h = Z3.getLocalViewHost ();
    
    TEST_ASSERT( size_t (Z3_lcl_d.stride (1)) >= X_numCols );
    TEST_ASSERT( size_t (Z3_lcl_h.stride (1)) >= X_numCols );
    TEST_ASSERT( size_t (Z3_lcl_d.extent (0)) == X_numCols );
    TEST_ASSERT( size_t (Z3_lcl_h.extent (0)) == X_numCols );  
    TEST_ASSERT( size_t (Z3_lcl_d.extent (1)) == Y_numCols );
    TEST_ASSERT( size_t (Z3_lcl_h.extent (1)) == Y_numCols );  
    Z3.multiply (Teuchos::CONJ_TRANS, Teuchos::NO_TRANS,
		 STS::one (), X, Y, STS::zero ());

    const bool Z1_Z3_equal = multiVectorsEqual (Z1, Z3);
    TEST_ASSERT( Z1_Z3_equal );
  }

  {
    dual_view_type Z4_dv ("Z4", X_numCols, Y_numCols);
    MV Z4 (lclMap, Z4_dv);
    Z4.multiply (Teuchos::CONJ_TRANS, Teuchos::NO_TRANS,
		 STS::one (), X, Y, STS::zero ());
    const bool Z1_Z4_equal = multiVectorsEqual (Z1, Z4);
    TEST_ASSERT( Z1_Z4_equal );
  }

  {
    dual_view_type Z5_dv ("Z5", X_numCols, Y_numCols);
    MV Z5 (lclMap, Z5_dv);
    Z5.multiply (Teuchos::CONJ_TRANS, Teuchos::NO_TRANS,
		 STS::one (), X, Y, STS::zero ());
    const bool Z1_Z5_equal = multiVectorsEqual (Z1, Z5);
    TEST_ASSERT( Z1_Z5_equal );
  }
  
  {
    // Make sure Z6 has a stride greater than its number of rows.
    const size_t Z6_stride = X_numCols + 31;
    dual_view_type Z6_dv_extra ("Z6", Z6_stride, Y_numCols);
    auto Z6_dv = Kokkos::subview (Z6_dv_extra,
				  pair_type (0, X_numCols),
				  pair_type (0, Y_numCols));
    TEST_ASSERT( size_t (Z6_dv.extent (0) == X_numCols ) );
    TEST_ASSERT( size_t (Z6_dv.extent (1) == Y_numCols ) );		 
    // Kokkos could in theory insert padding in the row dimension.
    TEST_ASSERT( size_t (Z6_dv.d_view.stride (1)) >= Z6_stride );
    TEST_ASSERT( size_t (Z6_dv.h_view.stride (1)) >= Z6_stride );  
  
    MV Z6 (lclMap, Z6_dv);
    Z6.multiply (Teuchos::CONJ_TRANS, Teuchos::NO_TRANS,
		 STS::one (), X, Y, STS::zero ());
    const bool Z1_Z6_equal = multiVectorsEqual (Z1, Z6);
    TEST_ASSERT( Z1_Z6_equal );

    Z6.putScalar (-STS::one ()); // flag value
    XH_times_Y (STS::one (), X, Y, Z6);
    const bool Z1_Z6_again_equal = multiVectorsEqual (Z1, Z6);
    TEST_ASSERT( Z1_Z6_again_equal );

    Z6.putScalar (-STS::one ()); // flag value
    XH_times_Y_2 (STS::one (), X, Y, Z6);
    const bool Z1_Z6_again_equal2 = multiVectorsEqual (Z1, Z6);
    TEST_ASSERT( Z1_Z6_again_equal2 );
  }
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( MultiVector, DualViewPool, SCALAR, LO, GO, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )

} // namespace (anonymous)
