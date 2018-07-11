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

#include <Tpetra_ConfigDefs.hpp>

// mfh 29 Jun 2014: It makes life easier for Sierra developers if
// Trilinos just builds all .cpp files unconditionally, rather than
// making the decision whether to build them dependent on CMake
// options.  Thus, we just exclude all the content of this file if MPI
// is not enabled.
#ifdef HAVE_TPETRA_MPI
#  include <Tpetra_MpiPlatform.hpp>

namespace Tpetra {

  MpiPlatform<Tpetra::Details::DefaultTypes::node_type>::
  MpiPlatform () :
    comm_ (Teuchos::createMpiComm<int> (Teuchos::opaqueWrapper<MPI_Comm> (MPI_COMM_WORLD)))
  {}

  MpiPlatform<Tpetra::Details::DefaultTypes::node_type>::
  MpiPlatform (int* argc, char*** argv) :
    comm_ (Teuchos::null)
  {
    initialize (argc, argv);
    comm_ = getDefaultComm ();
  }

  MpiPlatform<Tpetra::Details::DefaultTypes::node_type>::
  MpiPlatform (const Teuchos::RCP<NodeType>& /* node */) :
    comm_ (Teuchos::createMpiComm<int> (Teuchos::opaqueWrapper<MPI_Comm> (MPI_COMM_WORLD)))
  {
  }

  MpiPlatform<Tpetra::Details::DefaultTypes::node_type>::
  MpiPlatform (int* argc, char*** argv, const Teuchos::RCP<NodeType>& /* node */) :
    comm_ (Teuchos::null)
  {
    initialize (argc, argv);
    comm_ = getDefaultComm ();
  }

  MpiPlatform<Tpetra::Details::DefaultTypes::node_type>::
  MpiPlatform (const Teuchos::RCP<NodeType>& /* node */,
               const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> >& rawMpiComm)
    : comm_ (Teuchos::null)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      rawMpiComm.is_null (), std::invalid_argument, "Tpetra::MpiPlatform "
      "constructor: The input RCP<OpaqueWrapper<MPI_Comm> > is null.  That "
      "means something different than MPI_COMM_NULL.  If you want to give "
      "MPI_COMM_NULL to this constructor, please wrap MPI_COMM_NULL in a "
      "nonnull Teuchos::OpaqueWrapper by using the "
      "Teuchos::opaqueWrapper<MPI_Comm>() nonmember constructor.");
    comm_ = Teuchos::createMpiComm<int> (rawMpiComm);
  }

  MpiPlatform<Tpetra::Details::DefaultTypes::node_type>::
  MpiPlatform (int* argc,
               char*** argv,
               const Teuchos::RCP<NodeType>& /* node */,
               const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> >& rawMpiComm)
    : comm_ (Teuchos::null)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      rawMpiComm.is_null (), std::invalid_argument, "Tpetra::MpiPlatform "
      "constructor: The input RCP<OpaqueWrapper<MPI_Comm> > is null.  That "
      "means something different than MPI_COMM_NULL.  If you want to give "
      "MPI_COMM_NULL to this constructor, please wrap MPI_COMM_NULL in a "
      "nonnull Teuchos::OpaqueWrapper by using the "
      "Teuchos::opaqueWrapper<MPI_Comm>() nonmember constructor.");
    comm_ = Teuchos::createMpiComm<int> (rawMpiComm);

    // NOTE (mfh 29 Jun 2014): The OpaqueWrapper might wrap the
    // MPI_Comm in something that calls MPI_Comm_free.  Thus, we
    // can't just ignore it; we have to give it to the comm_ so that
    // its destructor (which might call MPI_Comm_free) will be
    // called at the right time.  This is why we don't set comm
    // using getDefaultComm().  This is also why we pass comm_
    // directly to initialize(): that way there aren't two
    // references to the raw MPI_Comm floating around, and comm_'s
    // destructor will get to do the right thing.
    initialize (argc, argv, comm_);
  }

  MpiPlatform<Tpetra::Details::DefaultTypes::node_type>::
  MpiPlatform (const Teuchos::RCP<NodeType>& /* node */, MPI_Comm rawMpiComm)
    : comm_ (Teuchos::createMpiComm<int> (Teuchos::opaqueWrapper<MPI_Comm> (rawMpiComm)))
  {
  }

  MpiPlatform<Tpetra::Details::DefaultTypes::node_type>::
  MpiPlatform (int* argc,
               char*** argv,
               const Teuchos::RCP<NodeType>& /* node */,
               MPI_Comm rawMpiComm)
    : comm_ (Teuchos::null)
  {
    initialize (argc, argv, rawMpiComm);
    comm_ = getDefaultComm ();
  }

  MpiPlatform<Tpetra::Details::DefaultTypes::node_type>::
  ~MpiPlatform () {}

  Teuchos::RCP<const Teuchos::Comm<int> >
  MpiPlatform<Tpetra::Details::DefaultTypes::node_type>::
  getComm () const {
    TEUCHOS_TEST_FOR_EXCEPTION(
      comm_.is_null (), std::logic_error, "Tpetra::MpiPlatform::getComm: "
      "The default communicator is null.  This should never happen.  "
      "Please report this bug to the Tpetra developers.");
    return comm_;
  }

  Teuchos::RCP<MpiPlatform< ::Tpetra::Details::DefaultTypes::node_type>::NodeType>
  MpiPlatform<Tpetra::Details::DefaultTypes::node_type>::
  getNode () const
  {
    return Teuchos::rcp (new Tpetra::Details::DefaultTypes::node_type);
  }

} // namespace Tpetra

#endif // HAVE_TPETRA_MPI
