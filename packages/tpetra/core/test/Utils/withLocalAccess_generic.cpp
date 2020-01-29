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
// ************************************************************************
// @HEADER
*/

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_withLocalAccess.hpp"
#include "Kokkos_Core.hpp"
#include <memory>

namespace { // (anonymous)

// In this test, we provide specializations of GetMasterLocalObject
// and GetNonowningLocalObject for this "stub global object" class, as
// a way to test that we can make withLocalAccess work for types other
// than Tpetra::MultiVector and Tpetra::Vector.
class StubGlobalObject {
private:
  using memory_space = Kokkos::HostSpace;
  using execution_space = Kokkos::DefaultHostExecutionSpace;

public:
  using device_type = Kokkos::Device<execution_space, memory_space>;

  StubGlobalObject() = delete;
  StubGlobalObject(const int input) :
    localData_(std::make_shared<int> (input))
  {}

  int value() const {
    return *localData_;
  }

  static constexpr int flagValue = -1;

private:
  std::shared_ptr<int> localData_;

  // Notice how I'm not letting users get the pointer directly.  If
  // your class must expose a public method for getting a pointer to
  // the data, try your best to make that a nonowning pointer.
  friend std::shared_ptr<int>
  getStubGlobalObjectOwningLocalData(StubGlobalObject& G)
  {
    return G.localData_;
  }
};

} // namespace (anonymous)

namespace Tpetra {
namespace Details {

template<class AccessMode>
struct StubGlobalObjectLocalObjectImpl {};

template<>
struct StubGlobalObjectLocalObjectImpl<Access<EAccess::ReadOnly> >
{
  // In practice, prefer things that behave like std::unique_ptr,
  // since you don't actually need to share state.
  using master_local_object_type = std::shared_ptr<const int>;

  static master_local_object_type
  getOwning(StubGlobalObject& G)
  {
    std::shared_ptr<int> p_nc =
      getStubGlobalObjectOwningLocalData(G);
    return std::shared_ptr<const int>(p_nc);
  }

  using nonowning_local_object_type = const int*;
  static nonowning_local_object_type
  getNonowning(StubGlobalObject& /* G */,
               const master_local_object_type& p)
  {
    return p.get ();
  }
};

// Example of the "copy-back" model.
template<>
struct StubGlobalObjectLocalObjectImpl<Access<EAccess::WriteOnly> >
{
private:
  struct Deleter {
    // Capturing the std::shared_ptr means that we increment the
    // reference count, and thus keep it alive.
    Deleter(int* raw, std::shared_ptr<int> sp) : raw_ (raw), sp_ (sp)
    {}
    void operator() (int* p) {
      if (p != nullptr) {
        *sp_ = *raw_; // "copy back" happens here
        delete p;
      }
    }
    int* raw_ = nullptr;
    std::shared_ptr<int> sp_;
  };

public:
  using master_local_object_type = std::unique_ptr<int, Deleter>;

  static master_local_object_type
  getOwning(StubGlobalObject& G)
  {
    std::shared_ptr<int> p_nc = getStubGlobalObjectOwningLocalData(G);
    // Test write-only-ness, by deliberately ignoring any existing
    // value.  In practice, your class could do this in debug mode.
    int* raw = new int(StubGlobalObject::flagValue);
    Deleter deleter(raw, p_nc);
    return {raw, deleter};
  }

  using nonowning_local_object_type = int*;

  static nonowning_local_object_type
  getNonowning(StubGlobalObject& /* G */,
               const master_local_object_type& p)
  {
    return p.get ();
  }
};

template<>
struct StubGlobalObjectLocalObjectImpl<Access<EAccess::ReadWrite> >
{
  // In practice, prefer things that behave like std::unique_ptr,
  // since you don't actually need to share state.
  using master_local_object_type = std::shared_ptr<int>;

  static master_local_object_type
  getOwning(StubGlobalObject& G)
  {
    return getStubGlobalObjectOwningLocalData(G);
  }

  using nonowning_local_object_type = int*;

  static nonowning_local_object_type
  getNonowning(StubGlobalObject& /* G */,
               const master_local_object_type& p)
  {
    return p.get ();
  }
};

//////////////////////////////////////////////////////////////////////

template<class ... Args>
struct GetMasterLocalObject<LocalAccess<StubGlobalObject, Args...>>
{
private:
  using local_access_type = LocalAccess<StubGlobalObject, Args...>;
  using access_mode = typename local_access_type::access_mode;
  using impl_type =
    StubGlobalObjectLocalObjectImpl<access_mode>;

public:
  // In practice, prefer things that behave like std::unique_ptr,
  // since you don't actually need to share state.
  using master_local_object_type =
    typename impl_type::master_local_object_type;

  static master_local_object_type get(local_access_type LA) {
    return impl_type::getOwning(LA.G_);
  }
};

//////////////////////////////////////////////////////////////////////

template<class ... Args>
struct GetNonowningLocalObject<LocalAccess<StubGlobalObject, Args...>>
{
private:
  using local_access_type = LocalAccess<StubGlobalObject, Args...>;
  using access_mode = typename local_access_type::access_mode;
  using impl_type = StubGlobalObjectLocalObjectImpl<access_mode>;

public:
  using master_local_object_type =
    typename impl_type::master_local_object_type;
  using nonowning_local_object_type =
    typename impl_type::nonowning_local_object_type;

  static nonowning_local_object_type
  get(local_access_type LA,
      const master_local_object_type& p)
  {
    return impl_type::getNonowning(LA.G_, p);
  }
};

} // namespace Details
} // namespace Tpetra

//////////////////////////////////////////////////////////////////////

namespace { // (anonymous)

TEUCHOS_UNIT_TEST(WithLocalAccess, generic)
{
  using ::Tpetra::withLocalAccess;
  using ::Tpetra::readOnly;
  using ::Tpetra::writeOnly;
  using ::Tpetra::readWrite;

  StubGlobalObject g1 (1);
  StubGlobalObject g2 (2);
  StubGlobalObject g3 (3);

  TEUCHOS_ASSERT( g1.value () == 1 );
  TEUCHOS_ASSERT( g2.value () == 2 );
  TEUCHOS_ASSERT( g3.value () == 3 );

  withLocalAccess
    ([] (int* g1_wo) {
       TEUCHOS_ASSERT( g1_wo != nullptr );
       // Write-only access means that users may not count on the
       // value read.  We enforce this in this generic test by
       // checking for a "flag value."
       const int flagValue = StubGlobalObject::flagValue;
       TEUCHOS_ASSERT( *g1_wo == flagValue );
       *g1_wo = 111;
     },
     writeOnly (g1));
  withLocalAccess
    ([] (const int* g1_ro) {
      TEUCHOS_ASSERT( g1_ro != nullptr );
      TEUCHOS_ASSERT( *g1_ro == 111 );
    },
    readOnly (g1));

  withLocalAccess
    ([] (int* g1_rw, int* g2_wo) {
      TEUCHOS_ASSERT( g1_rw != nullptr );
      TEUCHOS_ASSERT( g2_wo != nullptr );

      TEUCHOS_ASSERT( *g1_rw == 111 );
      *g1_rw = 666;

      const int flagValue = StubGlobalObject::flagValue;
      TEUCHOS_ASSERT( *g2_wo == flagValue );
      *g2_wo = 777;
    },
    readWrite (g1), writeOnly (g2));

  withLocalAccess
    ([] (const int* g1_ro, const int* g2_ro, const int* g3_ro) {
      TEUCHOS_ASSERT( g1_ro != nullptr );
      TEUCHOS_ASSERT( g2_ro != nullptr );
      TEUCHOS_ASSERT( g3_ro != nullptr );

      TEUCHOS_ASSERT( *g1_ro == 666 );
      TEUCHOS_ASSERT( *g2_ro == 777 );
      TEUCHOS_ASSERT( *g3_ro == 3 );
    },
    readOnly (g1), readOnly (g2), readOnly (g3));
}

} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  // Initialize MPI (if enabled) before initializing Kokkos.  This
  // lets MPI control things like pinning processes to sockets.
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  const int errCode =
    Teuchos::UnitTestRepository::runUnitTestsFromMain (argc, argv);
  return errCode;
}
