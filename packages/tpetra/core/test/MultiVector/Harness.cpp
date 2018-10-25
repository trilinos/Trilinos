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

#include <Tpetra_TestingUtilities.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <Teuchos_DefaultSerialComm.hpp>

namespace Harness {

  enum AccessMode {
    ReadOnly,
    WriteOnly,
    ReadWrite
  };

  namespace Impl {
    // Given a global object, get its default memory space (both the
    // type and the default instance thereof).
    template<class GlobalObjectType>
    struct DefaultMemorySpace {
      using type = typename GlobalObjectType::device_type::memory_space;

      // Given a global object, get its (default) memory space instance.
      static type space (const GlobalObjectType& /* G */) {
	// This stub just assumes that 'type' is default constructible.
	// In Kokkos, default-constructing a memory space instance just
	// gives the default memory space.
	return type ();
      }
    };

    // Struct that tells withLocalAccess how to access a global object's
    // local data.  Do not use this directly; start with readOnly,
    // writeOnly, or readWrite.
    template<class GlobalObjectType,
	     class MemorySpace,
	     const AccessMode am>
    class LocalAccess; // forward declaration

    // Mapping from LocalAccess to the "master" local object type.  The
    // latter gets the local data from a global object, and holds on to
    // it until after the user's function (input to withLocalAccess)
    // returns.
    template<class LocalAccessType>
    struct GetMasterLocalObject {};

    // Given a LocalAccess instance (which has a reference to a global
    // object), get an instance of its master local object.  This may be
    // a heavyweight operation.
    //
    // If you need to specialize this, just specialize get() in
    // GetMasterLocalObject above.
    template<class LocalAccessType>
    typename GetMasterLocalObject<LocalAccessType>::master_local_object_type
    getMasterLocalObject (LocalAccessType LA) {
      return GetMasterLocalObject<LocalAccessType>::get (LA);
    }

    // Mapping from "master" local object type to the nonowning "local
    // view" type that users see (as arguments to the function that they
    // give to withLocalAccess).  The master local object may encode the
    // memory space and access mode, but the mapping to local view type
    // may also need run-time information.
    template<class MasterLocalObjectType>
    struct GetNonowningLocalObject {};

    // Given a master local object, get an instance of a nonowning local
    // object.  Users only ever see the nonowning local object, and
    // subviews (slices) thereof.  This is supposed to be a lightweight
    // operation.
    //
    // If you need to specialize this, just specialize get() in
    // GetNonowningLocalObject above.
    template<class MasterLocalObjectType>
    typename GetNonowningLocalObject<MasterLocalObjectType>::nonowning_local_object_type
    getNonowningLocalObject (const MasterLocalObjectType& master) {
      return GetNonowningLocalObject<MasterLocalObjectType>::get (master);
    }

    // Use the LocalAccess type as the template parameter to determine
    // the type of the nonowning local view to the global object's data.
    // This only works if GetMasterLocalObject has been specialized for
    // these template parameters, and if GetNonowningLocalObject has
    // been specialized for the resulting "master" local object type.
    template<class LocalAccessType>
    class LocalAccessFunctionArgument {
    private:
      using gmlo = GetMasterLocalObject<LocalAccessType>;
      using master_local_object_type = typename gmlo::master_local_object_type;
      using gnlo = GetNonowningLocalObject<master_local_object_type>;
    public:
      using type = typename gnlo::nonowning_local_object_type;
    };
  } // namespace Impl

  //////////////////////////////////////////////////////////////////////
  // Users call readOnly, writeOnly, and readWrite, in order to declare
  // how they intend to access a global object's local data.
  //////////////////////////////////////////////////////////////////////

  // Declare that you want to access the given global object's local
  // data in read-only mode.
  template<class GlobalObjectType>
  Impl::LocalAccess<GlobalObjectType,
		    typename Impl::DefaultMemorySpace<GlobalObjectType>::type,
		    ReadOnly>
  readOnly (GlobalObjectType&);

  template<class GlobalObjectType>
  Impl::LocalAccess<GlobalObjectType,
		    typename Impl::DefaultMemorySpace<GlobalObjectType>::type,
		    ReadOnly>
  readOnly (const GlobalObjectType&);

  // Declare that you want to access the given global object's local
  // data in write-only mode.
  template<class GlobalObjectType>
  Impl::LocalAccess<GlobalObjectType,
		    typename Impl::DefaultMemorySpace<GlobalObjectType>::type,
		    WriteOnly>
  writeOnly (GlobalObjectType&);

  // Declare that you want to access the given global object's local
  // data in read-and-write mode.
  template<class GlobalObjectType>
  Impl::LocalAccess<GlobalObjectType,
		    typename Impl::DefaultMemorySpace<GlobalObjectType>::type,
		    ReadWrite>
  readWrite (GlobalObjectType&);

  namespace Impl {
    // Declaration of access intent for a global object.
    //
    // Users aren't supposed to make instances of this class.  They
    // should use readOnly, writeOnly, or readWrite instead, then call
    // methods like on() and valid() on the resulting LocalAccess
    // instance.
    template<class GlobalObjectType,
	     class MemorySpace,
	     const AccessMode am>
    class LocalAccess {
    public:
      using global_object_type = GlobalObjectType;
      using memory_space = typename MemorySpace::memory_space;    
      static constexpr AccessMode access_mode = am;

      // Users must NOT call the LocalAccess constructor directly.  They
      // should instead start by calling readOnly, writeOnly, or
      // readWrite above.  They may then use methods like on() or
      // valid() (see below).
      //
      // G is a reference, because we only access it in a delimited
      // scope.  G is nonconst, because even read-only local access may
      // modify G.  For example, G may give access to its local data via
      // lazy allocation of a data structure that differs from its
      // normal internal storage format.
      //
      // Memory spaces should behave like Kokkos memory spaces.  Default
      // construction should work and should get the default instance of
      // the space.  Otherwise, it may make sense to get the default
      // memory space from G.
      LocalAccess (global_object_type& G,
		   memory_space space = memory_space (),
		   const bool isValid = true) :
	G_ (G),
	space_ (space),
	valid_ (isValid)
      {}

      // Type users see, that's an argument to the function that they give
      // to withLocalAccess.
      using function_argument_type =
	typename LocalAccessFunctionArgument<LocalAccess<global_object_type, memory_space, access_mode> >::type;

    public:
      // Give users run-time control over whether they actually want to
      // access the object.  If isValid is false, implementations should
      // not spend any effort getting the master local object.  This may
      // save time on allocating temporary space, copying from device to
      // host, etc.  This implies that implementations must be able to
      // construct "null" / empty master local objects.
      LocalAccess<GlobalObjectType, MemorySpace, am>
      valid (const bool isValid) const {
	std::cout << "  .valid(" << (isValid ? "true" : "false") << ")" << std::endl;
	return {this->G_, this->space_, isValid};
      }

      // Let users access this object in a different memory space.
      //
      // NOTE: This also works for PGAS.  'space' in that case could be
      // something like an MPI_Win, or a UPC-style "shared pointer."
      //
      // NewMemorySpace should behave like a Kokkos memory space
      // instance.
      template<class NewMemorySpace>
      LocalAccess<GlobalObjectType, NewMemorySpace, am>
      on (NewMemorySpace space) const {
	std::cout << "  .on(" << space.name () << ")" << std::endl;
	return {this->G_, space, this->valid_};
      }

      // Is access supposed to be valid?  (See valid() above.)
      bool isValid () const { return this->valid_; }

      memory_space getSpace () const { return space_; }
  
    public:
      // Keep by reference, because this struct is only valid in a
      // delimited scope.
      global_object_type& G_;
      memory_space space_; // assume shallow-copy semantics
      bool valid_; // will I actually need to access this object?

    private:
      // Nonmember "constructors"; see above for declarations.  This are
      // friends, because they are the only ways that users are supposed
      // to construct LocalAccess instances.
      template<class GOT> friend
      LocalAccess<GOT, typename Impl::DefaultMemorySpace<GOT>::type, ReadOnly> readOnly (GOT&);
      template<class GOT> friend
      LocalAccess<GOT, typename Impl::DefaultMemorySpace<GOT>::type, ReadOnly> readOnly (const GOT&);
      template<class GOT> friend
      LocalAccess<GOT, typename Impl::DefaultMemorySpace<GOT>::type, WriteOnly> writeOnly (GOT&);
      template<class GOT> friend
      LocalAccess<GOT, typename Impl::DefaultMemorySpace<GOT>::type, ReadWrite> readWrite (GOT&);
    };
  } // namespace Impl

  template<class GOT>
  Impl::LocalAccess<GOT, typename Impl::DefaultMemorySpace<GOT>::type, ReadOnly>
  readOnly (GOT& G)
  {
    std::cout << "readOnly" << std::endl;
    return {G, Impl::DefaultMemorySpace<GOT>::space (G), true};
  }

  template<class GOT>
  Impl::LocalAccess<GOT, typename Impl::DefaultMemorySpace<GOT>::type, ReadOnly>
  readOnly (const GOT& G)
  {
    std::cout << "readOnly" << std::endl;
    GOT& G_nc = const_cast<GOT&> (G);
    return {G_nc, Impl::DefaultMemorySpace<GOT>::space (G_nc), true};
  }
  
  template<class GOT>
  Impl::LocalAccess<GOT, typename Impl::DefaultMemorySpace<GOT>::type, WriteOnly>
  writeOnly (GOT& G)
  {
    std::cout << "writeOnly" << std::endl;
    return {G, Impl::DefaultMemorySpace<GOT>::space (G), true};
  }

  template<class GOT>
  Impl::LocalAccess<GOT, typename Impl::DefaultMemorySpace<GOT>::type, ReadWrite>
  readWrite (GOT& G)
  {
    std::cout << "readWrite" << std::endl;  
    return {G, Impl::DefaultMemorySpace<GOT>::space (G), true};
  }

  namespace Impl {
    template<class SC, class LO, class GO, class NT,
	     class MemorySpace>
    using multivector_nonconst_nonowning_local_object_type =
      typename std::conditional<
        std::is_same<typename MemorySpace::memory_space, Kokkos::HostSpace>::value,
        typename Tpetra::MultiVector<SC, LO, GO, NT>::dual_view_type::t_host,
        typename Tpetra::MultiVector<SC, LO, GO, NT>::dual_view_type::t_dev>::type;

    template<class SC, class LO, class GO, class NT,
	     class MemorySpace>
    using multivector_const_nonowning_local_object_type =
      typename multivector_nonconst_nonowning_local_object_type<SC, LO, GO, NT, MemorySpace>::const_type;

    template<class SC, class LO, class GO, class NT,
	     class MemorySpace, 
	     const AccessMode am>
    using multivector_nonowning_local_object_type =
      typename std::conditional<
        am == ReadOnly, 
	multivector_const_nonowning_local_object_type<SC, LO, GO, NT, MemorySpace>,
	multivector_nonconst_nonowning_local_object_type<SC, LO, GO, NT, MemorySpace> >::type;
  } // namespace Impl

  template<class SC, class LO, class GO, class NT,
	   class MemorySpace, 
	   const AccessMode access_mode>
  Impl::multivector_nonowning_local_object_type<SC, LO, GO, NT, MemorySpace, access_mode>
  getMultiVector (Impl::LocalAccess<Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace, access_mode> LA)
  {
    using memory_space = typename MemorySpace::memory_space;
    using ret_type =
      Impl::multivector_nonowning_local_object_type<SC, LO, GO, NT,
						    memory_space, access_mode>;
    std::cout << "getMultiVector" << std::endl;

    // if (access_mode == WriteOnly) { ...} // FIXME (mfh 22 Oct 2018)
    if (LA.G_.template need_sync<memory_space> ()) {
      LA.G_.template sync<memory_space> ();
    }
    if (access_mode != ReadWrite) {
      LA.G_.template modify<memory_space> ();
    }
    // FIXME (mfh 22 Oct 2018) This might break if we need copy-back
    // semantics, e.g., for a memory space for which the
    // Tpetra::MultiVector does not store data.  In that case, we
    // would need some kind of object whose destructor copies back,
    // and it would need to have the whole DualView, not just the
    // View on one side.  Watch out for copy elision.  The object
    // could just be std::shared_ptr and could handle copy-back via
    // custom deleter.
    if (LA.isValid ()) {
      // this converts to const if applicable
      return ret_type (LA.G_.template getLocalView<memory_space> ()); 
    }
    else {
      return ret_type (); // "null" Kokkos::View
    }
  }

#if 0
  template<class SC, class LO, class GO, class NT,
	   class MemorySpace, 
	   const AccessMode access_mode>
  Impl::multivector_nonowning_local_object_type<SC, LO, GO, NT, MemorySpace, access_mode>
  getMultiVector (Impl::LocalAccess<const Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace, access_mode> LA)
  {
    using MV = Tpetra::MultiVector<SC, LO, GO, NT>;
    MV& G = const_cast<MV&> (LA.G_);
    Impl::LocalAccess<const MV, MemorySpace, access_mode> LA2 (G, LA.space_, LA.valid_);
    return getMultiVector (LA2);
  }
#endif // 0

  template<class SC, class LO, class GO, class NT,
	   class MemorySpace, 
	   const AccessMode access_mode>
  auto
  getVector (Impl::LocalAccess<Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace, access_mode> LA)
    -> decltype (Kokkos::subview (getMultiVector (LA), Kokkos::ALL (), 0))
  {
    return Kokkos::subview (getMultiVector (LA), Kokkos::ALL (), 0);
  }

  template<class SC, class LO, class GO, class NT,
	   class MemorySpace, 
	   const AccessMode access_mode>
  auto
  getVector (Impl::LocalAccess<Tpetra::Vector<SC, LO, GO, NT>, MemorySpace, access_mode> LA)
    -> decltype (Kokkos::subview (getMultiVector (Impl::LocalAccess<Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace, access_mode> (LA.G_, LA.space_,  LA.valid_)), Kokkos::ALL (), 0))
  {
    return Kokkos::subview (getMultiVector (Impl::LocalAccess<Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace, access_mode> (LA.G_, LA.space_,  LA.valid_)), Kokkos::ALL (), 0);
  }
  
  template<class SC, class LO, class GO, class NT,
	   class UnaryFunction,
	   class ExecutionSpace,
	   class MemorySpace = typename ExecutionSpace::memory_space>
  void
  transform (Tpetra::MultiVector<SC, LO, GO, NT>& X,
	     UnaryFunction f,
	     ExecutionSpace execSpace = ExecutionSpace (),
	     MemorySpace memSpace = MemorySpace ())
  {
    using execution_space =
      typename Tpetra::MultiVector<SC, LO, GO, NT>::device_type::execution_space;
    using range_type = Kokkos::RangePolicy<execution_space, LO>;
    if (X.getNumVectors () == size_t (1)) {
      auto X_lcl = getVector (readWrite (X).on (memSpace));
      Kokkos::parallel_for ("transform",
			    range_type (execSpace, 0, X_lcl.extent (0)),
			    KOKKOS_LAMBDA (const LO i) {
			      X_lcl(i) = f (X_lcl(i), i, 0);
			    });
    }
    else {
      auto X_lcl = getMultiVector (readWrite (X).on (memSpace));
      Kokkos::parallel_for ("transform",
			    range_type (execSpace, 0, X_lcl.extent (0)),
			    KOKKOS_LAMBDA (const LO i) {
			      const LO numVecs = X_lcl.extent (1);
			      for (LO j = 0; j < numVecs; ++j) {
				X_lcl(i,j) = f (X_lcl(i,j), i, j);
			      }
			    });
    }
  }

  template<class SC, class LO, class GO, class NT,
	   class UnaryFunction>
  void
  transform (Tpetra::MultiVector<SC, LO, GO, NT>& X,
	     UnaryFunction f)
  {
    using execution_space = typename Tpetra::MultiVector<SC, LO, GO, NT>::device_type::execution_space;
    using memory_space = typename Tpetra::MultiVector<SC, LO, GO, NT>::device_type::memory_space;    
    transform (X, f, execution_space (), memory_space ());
  }

  namespace Impl {
    // Specialization of GetMasterLocalObject for Tpetra::MultiVector.
    template<class SC, class LO, class GO, class NT,
	     class MemorySpace,
	     const AccessMode am>
    struct GetMasterLocalObject<LocalAccess<Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace, am> > {
    public:
      using local_access_type = LocalAccess<Tpetra::MultiVector<SC, LO, GO, NT>, MemorySpace, am>;        
    private:
      using global_object_type = typename local_access_type::global_object_type;
      using memory_space = typename local_access_type::memory_space;
      static constexpr AccessMode access_mode = local_access_type::access_mode;
      using non_const_value_type = typename Tpetra::MultiVector<SC, LO, GO, NT>::impl_scalar_type;
      using value_type = typename std::conditional<
        access_mode == ReadOnly,
        const non_const_value_type,
        non_const_value_type
	>::type;

    public:
      // FIXME (mfh 22 Oct 2018) See FIXME below.
      using master_local_object_type =
	Kokkos::View<value_type**, 
		     typename global_object_type::dual_view_type::t_dev::array_layout,
		     MemorySpace>; // FIXME (mfh 22 Oct 2018) need to make sure execution_space matches

      static master_local_object_type get (local_access_type LA) {
	std::cout << "Get master local object" << std::endl;

	// if (access_mode == WriteOnly) { ...} // FIXME (mfh 22 Oct 2018)
	if (LA.G_.template need_sync<memory_space> ()) {
	  LA.G_.template sync<memory_space> ();
	}
	if (access_mode != ReadWrite) {
	  LA.G_.template modify<memory_space> ();
	}
	// FIXME (mfh 22 Oct 2018) This might break if we need copy-back
	// semantics, e.g., for a memory space for which the
	// Tpetra::MultiVector does not store data.  In that case, we
	// would need some kind of object whose destructor copies back,
	// and it would need to have the whole DualView, not just the
	// View on one side.  Watch out for copy elision.  The object
	// could just be std::shared_ptr and could handle copy-back via
	// custom deleter.
	if (LA.isValid ()) {
	  // this converts to const if applicable
	  return master_local_object_type (LA.G_.template getLocalView<memory_space> ()); 
	}
	else {
	  return master_local_object_type (); // "null" Kokkos::View
	}
      }
    };

    // Specialization of GetMasterLocalObject for Tpetra::Vector.
    template<class SC, class LO, class GO, class NT,
	     class MemorySpace,
	     const AccessMode am>
    struct GetMasterLocalObject<LocalAccess<Tpetra::Vector<SC, LO, GO, NT>, MemorySpace, am> > {
    public:
      using local_access_type = LocalAccess<Tpetra::Vector<SC, LO, GO, NT>, MemorySpace, am>;        
    private:
      using global_object_type = typename local_access_type::global_object_type;
      using memory_space = typename local_access_type::memory_space;
      static constexpr AccessMode access_mode = local_access_type::access_mode;
      using non_const_value_type = typename global_object_type::impl_scalar_type;
      using value_type = typename std::conditional<
        access_mode == ReadOnly,
        const non_const_value_type,
        non_const_value_type
	>::type;

    public:
      // FIXME (mfh 22 Oct 2018) See FIXME below.
      using master_local_object_type =
	Kokkos::View<value_type**, 
		     typename global_object_type::dual_view_type::t_dev::array_layout,
		     MemorySpace>; // FIXME (mfh 22 Oct 2018) need to make sure execution_space matches

      static master_local_object_type get (local_access_type LA) {
	std::cout << "Get master local object" << std::endl;

	// if (access_mode == WriteOnly) { ...} // FIXME (mfh 22 Oct 2018)
	if (LA.G_.template need_sync<memory_space> ()) {
	  LA.G_.template sync<memory_space> ();
	}
	if (access_mode != ReadWrite) {
	  LA.G_.template modify<memory_space> ();
	}
	// FIXME (mfh 22 Oct 2018) This might break if we need copy-back
	// semantics, e.g., for a memory space for which the
	// Tpetra::MultiVector does not store data.  In that case, we
	// would need some kind of object whose destructor copies back,
	// and it would need to have the whole DualView, not just the
	// View on one side.  Watch out for copy elision.  The object
	// could just be std::shared_ptr and could handle copy-back via
	// custom deleter.
	if (LA.isValid ()) {
	  // this converts to const if applicable
	  return master_local_object_type (LA.G_.template getLocalView<memory_space> ()); 
	}
	else {
	  return master_local_object_type (); // "null" Kokkos::View
	}
      }
    };

    // Specialization of GetNonowningLocalObject for StubMasterLocalObject.
    template<class ValueType,
	     class LayoutType,
	     class MemorySpace>
    struct GetNonowningLocalObject<Kokkos::View<ValueType**, LayoutType, MemorySpace> > {
      using master_local_object_type = Kokkos::View<ValueType**, LayoutType, MemorySpace>;
      using nonowning_local_object_type = Kokkos::View<ValueType**, LayoutType, MemorySpace, Kokkos::MemoryUnmanaged>;
      static nonowning_local_object_type get (const master_local_object_type& M) {
	return nonowning_local_object_type (M); // standard Kokkos::View assignment
      }
    };
  
  } // namespace Impl
  
} // namespace Harness

namespace { // (anonymous)

  using Tpetra::TestingUtilities::getDefaultComm;
  using Teuchos::outArg;
  using Teuchos::RCP;  
  using Teuchos::rcp;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::endl;
  using GST = Tpetra::global_size_t;
  using map_type = Tpetra::Map<>;
  using multivec_type = Tpetra::MultiVector<>;
  using vec_type = Tpetra::Vector<>;
  using GO = map_type::global_ordinal_type;

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST( VectorHarness, GetVector )
  {
    constexpr bool debug = true;

    RCP<Teuchos::FancyOStream> outPtr = debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Test Harness:{getMultiVector, getVector}" << endl;
    Teuchos::OSTab tab0 (myOut);

    myOut << "Create a Map" << endl;
    auto comm = getDefaultComm ();
    const auto INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const size_t numLocal = 13;
    const size_t numVecs  = 3;
    const GO indexBase = 0;
    auto map = rcp (new map_type (INVALID, numLocal, indexBase, comm));

    myOut << "Create a MultiVector, and make sure that it has "
      "the right number of vectors (columns)" << endl;
    multivec_type mvec (map, numVecs);
    TEST_EQUALITY( mvec.getNumVectors (), numVecs );

    myOut << "Create a Vector, and make sure that "
      "it has exactly one vector (column)" << endl;
    vec_type vec (map);
    TEST_EQUALITY_CONST(vec.getNumVectors (), size_t (1));

    // asMultiVector converts LocalAccess struct to nonowning rank-2 Kokkos::View.
    //
    // auto X_lcl_2d = asMultiVector (readOnly (mvec)); // rank 2
    // auto X_lcl_2d_host = asMultiVector (readOnly (mvec).on (Kokkos::HostSpace)); // rank 2
    {
      auto X_lcl_ro = Harness::getMultiVector (Harness::readOnly (mvec));
      // Make sure X_lcl_ro can be assigned to the type we expect it to
      // be.  It doesn't have to be that type, it just has to be
      // assignable to that type.
      Kokkos::View<const double**, Kokkos::LayoutLeft, multivec_type::device_type, Kokkos::MemoryUnmanaged> X_lcl_ro2 = X_lcl_ro;
      static_assert (decltype (X_lcl_ro)::Rank == 2, "Rank is not 2");
      TEST_ASSERT( size_t (X_lcl_ro.extent (0)) == numLocal );
      TEST_ASSERT( size_t (X_lcl_ro.extent (1)) == numVecs );
    }

    // Test whether read-only access works with a const MultiVector&.
    {
      auto X_lcl_ro = Harness::getMultiVector (Harness::readOnly (const_cast<const multivec_type&> (mvec)));
      // Make sure X_lcl_ro can be assigned to the type we expect it to
      // be.  It doesn't have to be that type, it just has to be
      // assignable to that type.
      Kokkos::View<const double**, Kokkos::LayoutLeft, multivec_type::device_type, Kokkos::MemoryUnmanaged> X_lcl_ro2 = X_lcl_ro;
      static_assert (decltype (X_lcl_ro)::Rank == 2, "Rank is not 2");
      TEST_ASSERT( size_t (X_lcl_ro.extent (0)) == numLocal );
      TEST_ASSERT( size_t (X_lcl_ro.extent (1)) == numVecs );
    }

    {
      auto X_lcl_wo = Harness::getMultiVector (Harness::writeOnly (mvec));
      Kokkos::View<double**, Kokkos::LayoutLeft, multivec_type::device_type, Kokkos::MemoryUnmanaged> X_lcl_wo2 = X_lcl_wo;
      static_assert (decltype (X_lcl_wo)::Rank == 2, "Rank is not 2");
      TEST_ASSERT( size_t (X_lcl_wo.extent (0)) == numLocal );
      TEST_ASSERT( size_t (X_lcl_wo.extent (1)) == numVecs );
    }

    {
      auto X_lcl_rw = Harness::getMultiVector (Harness::readWrite (mvec));
      Kokkos::View<double**, Kokkos::LayoutLeft, multivec_type::device_type, Kokkos::MemoryUnmanaged> X_lcl_rw2 = X_lcl_rw;
      static_assert (decltype (X_lcl_rw)::Rank == 2, "Rank is not 2");
      TEST_ASSERT( size_t (X_lcl_rw.extent (0)) == numLocal );
      TEST_ASSERT( size_t (X_lcl_rw.extent (1)) == numVecs );
    }

    // NOTE (mfh 23 Oct 2018) Don't use the commented-out interface
    // below.  It's really only meant to be an implementation detail
    // of withLocalAccess, which we're not using for VectorHarness.
#if 0
    {
      auto X_lcl_ro = Harness::Impl::getNonowningLocalObject (Harness::Impl::getMasterLocalObject (Harness::readOnly (mvec)));

      // Make sure X_lcl_ro can be assigned to the type we expect it to
      // be.  It doesn't have to be that type, it just has to be
      // assignable to that type.
      Kokkos::View<const double**, Kokkos::LayoutLeft, multivec_type::device_type, Kokkos::MemoryUnmanaged> X_lcl_ro2 = X_lcl_ro;
      static_assert (decltype (X_lcl_ro)::Rank == 2, "Rank is not 2");
      TEST_ASSERT( size_t (X_lcl_ro.extent (0)) == numLocal );
      TEST_ASSERT( size_t (X_lcl_ro.extent (1)) == numVecs );
    }
#endif // 0

    // asVector converts LocalAccess struct to nonowning rank-1 Kokkos::View.
    //
    // auto X_lcl_1d = asVector (readOnly (mvec)); // rank 1
    // auto X_lcl_1d_host = asVector (readOnly (mvec).on (Kokkos::HostSpace)); // rank 1
    {
      auto X_lcl_1d_ro = Harness::getVector (Harness::readOnly (vec));
      Kokkos::View<const double*, Kokkos::LayoutLeft, multivec_type::device_type, Kokkos::MemoryUnmanaged> X_lcl_1d_ro2 = X_lcl_1d_ro;
      static_assert (decltype (X_lcl_1d_ro)::Rank == 1, "Rank is not 1");
      TEST_ASSERT( size_t (X_lcl_1d_ro.extent (0)) == numLocal );
    }

    {
      auto X_lcl_1d_wo = Harness::getVector (Harness::writeOnly (vec));
      Kokkos::View<double*, Kokkos::LayoutLeft, multivec_type::device_type, Kokkos::MemoryUnmanaged> X_lcl_1d_wo2 = X_lcl_1d_wo;
      static_assert (decltype (X_lcl_1d_wo)::Rank == 1, "Rank is not 1");
      TEST_ASSERT( size_t (X_lcl_1d_wo.extent (0)) == numLocal );
    }

    {
      auto X_lcl_1d_wr = Harness::getVector (Harness::readWrite (vec));
      Kokkos::View<double*, Kokkos::LayoutLeft, multivec_type::device_type, Kokkos::MemoryUnmanaged> X_lcl_1d_wr2 = X_lcl_1d_wr;
      static_assert (decltype (X_lcl_1d_wr)::Rank == 1, "Rank is not 1");
      TEST_ASSERT( size_t (X_lcl_1d_wr.extent (0)) == numLocal );
    }

    //
    // Examples of using the result of getVector in Kokkos::parallel_for kernels.
    //

    {
      using execution_space = vec_type::device_type::execution_space;
      using memory_space = vec_type::device_type::memory_space;
      using LO = vec_type::local_ordinal_type;
      using range_type = Kokkos::RangePolicy<execution_space, LO>;

      auto X_lcl_1d_wo = Harness::getVector (Harness::writeOnly (vec).on (memory_space ()));
      static_assert (std::is_same<
		       decltype (X_lcl_1d_wo)::device_type::execution_space,
		       vec_type::dual_view_type::t_host::execution_space>::value,
		     "Not host execution space");
      Kokkos::parallel_for (
	"Device kernel for write-only getVector",
	range_type (0, LO (numLocal)),
	KOKKOS_LAMBDA (const LO lclRow) { X_lcl_1d_wo(lclRow) = 42.0; });
    }

    {
      using host_execution_space = vec_type::dual_view_type::t_host::execution_space;
      using LO = vec_type::local_ordinal_type;
      using range_type = Kokkos::RangePolicy<host_execution_space, LO>;

      auto X_lcl_1d_wo = Harness::getVector (Harness::writeOnly (vec).on (Kokkos::HostSpace ()));
      static_assert (std::is_same<
		       decltype (X_lcl_1d_wo)::device_type::execution_space,
		       vec_type::dual_view_type::t_host::execution_space>::value,
		     "Not host execution space");
      Kokkos::parallel_for (
	"Host kernel for write-only getVector",
	range_type (0, LO (numLocal)),
	[=] (const LO lclRow) {
	  std::pair<double, double> p {3.0, 4.0}; // some not-device function
	  X_lcl_1d_wo(lclRow) = p.first * p.second;
	});

      // Just plain modify some entries, in some sequential order.
      // Just in case LO is unsigned, start at +1.
      for (LO lclRowPlusOne = LO (numLocal); lclRowPlusOne > LO (0); --lclRowPlusOne) {
	const LO lclRow = lclRowPlusOne - LO (1);
	// operator[] for 1-D Kokkos::View does the same thing as
	// operator().
	X_lcl_1d_wo[lclRow] = double (lclRow) + 42.0;
      }
    }

    // Make sure that the test passed on all processes, not just Proc 0.
    int lclSuccess = success ? 1 : 0;
    int gblSuccess = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
  }

  TEUCHOS_UNIT_TEST( VectorHarness, Transform )
  {
    constexpr bool debug = true;
    using LO = typename multivec_type::local_ordinal_type;

    RCP<Teuchos::FancyOStream> outPtr = debug ?
      Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
      Teuchos::rcpFromRef (out);
    Teuchos::FancyOStream& myOut = *outPtr;

    myOut << "Test Harness::transform" << endl;
    Teuchos::OSTab tab0 (myOut);

    myOut << "Create a Map" << endl;
    auto comm = getDefaultComm ();
    const auto INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const size_t numLocal = 13;
    const size_t numVecs  = 3;
    const GO indexBase = 0;
    auto map = rcp (new map_type (INVALID, numLocal, indexBase, comm));

    myOut << "Create a MultiVector, and make sure that it has "
      "the right number of vectors (columns)" << endl;
    multivec_type X (map, numVecs);
    TEST_EQUALITY( X.getNumVectors (), numVecs );

    using Harness::transform;
    transform (X, [] (const double X_ij, const LO i, const LO j) {
	return X_ij + double (i+1.0) + double (j+1.0);
      });
    transform (X, [] (const double X_ij, const LO i, const LO j) {
	return X_ij + double (i+1.0) + double (j+1.0);
      }, Kokkos::DefaultHostExecutionSpace ());
    transform (X, [] (const double X_ij, const LO i, const LO j) {
	return X_ij + double (i+1.0) + double (j+1.0);
      });
    transform (X, [] (const double X_ij, const LO i, const LO j) {
	return X_ij + double (i+1.0) + double (j+1.0);
      }, Kokkos::DefaultHostExecutionSpace (), Kokkos::HostSpace ());


    X.sync<Kokkos::HostSpace> ();
    auto X_lcl = X.getLocalView<Kokkos::HostSpace> ();
    bool ok = true;    
    for (LO j = 0; j < LO (X.getNumVectors ()); ++j) {
      for (LO i = 0; i < LO (X.getLocalLength ()); ++i) {
	if (X_lcl(i,j) != 4.0 * (double (i+1.0) + double(j+1.0))) {
	  ok = false;
	}
      }
    }
    TEST_ASSERT( ok );

    myOut << "Create a Vector, and make sure that "
      "it has exactly one vector (column)" << endl;
    vec_type vec (map);
    TEST_EQUALITY_CONST(vec.getNumVectors (), size_t (1));

    transform (vec, [] (const double X_ij, const LO i, const LO j) {
	return X_ij + double (i+1.0) + double (j+1.0);
      });
  }  
} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  const int errCode =
    Teuchos::UnitTestRepository::runUnitTestsFromMain (argc, argv);
  return errCode;
}

