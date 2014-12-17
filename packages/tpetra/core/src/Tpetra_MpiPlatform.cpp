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

  MpiPlatform<KokkosClassic::DefaultNode::DefaultNodeType>::
  MpiPlatform () :
    comm_ (Teuchos::createMpiComm<int> (Teuchos::opaqueWrapper<MPI_Comm> (MPI_COMM_WORLD))),
    node_ (Teuchos::null)
  {
    // mfh 29 Jun 2014: Don't initialize the Node yet.  This ensures
    // that (new) Kokkos won't get initialized with the wrong
    // command-line arguments, at least not until getNode() is called.
    // Initializing Kokkos with the wrong command-line arguments may
    // result in poor performance due to the wrong assignment of
    // software threads to hardware execution units.
    //
    // if (node_.is_null ()) {
    //   node_ = KokkosClassic::DefaultNode::getDefaultNode ();
    // }
  }

  MpiPlatform<KokkosClassic::DefaultNode::DefaultNodeType>::
  MpiPlatform (int* argc, char*** argv) :
    comm_ (Teuchos::null),
    node_ (Teuchos::null)
  {
    initialize (argc, argv);
    comm_ = getDefaultComm ();

    // mfh 29 Jun 2014: Don't initialize the Node yet.  See above note.
    //
    // if (node_.is_null ()) {
    //   node_ = KokkosClassic::DefaultNode::getDefaultNode ();
    // }
  }

  MpiPlatform<KokkosClassic::DefaultNode::DefaultNodeType>::
  MpiPlatform (const Teuchos::RCP<NodeType>& node) :
    comm_ (Teuchos::createMpiComm<int> (Teuchos::opaqueWrapper<MPI_Comm> (MPI_COMM_WORLD))),
    node_ (node)
  {
    // mfh 29 Jun 2014: Don't initialize the Node yet.  This ensures
    // that (new) Kokkos won't get initialized with the wrong
    // command-line arguments, at least not until getNode() is
    // called.  Initializing Kokkos with the wrong command-line
    // arguments may result in poor performance due to the wrong
    // assignment of software threads to hardware execution units.
    //
    // if (node_.is_null ()) {
    //   node_ = KokkosClassic::DefaultNode::getDefaultNode ();
    // }
  }

  MpiPlatform<KokkosClassic::DefaultNode::DefaultNodeType>::
  MpiPlatform (int* argc, char*** argv, const Teuchos::RCP<NodeType>& node) :
    comm_ (Teuchos::null),
    node_ (node)
  {
    initialize (argc, argv);
    comm_ = getDefaultComm ();

    // mfh 29 Jun 2014: Don't initialize the Node yet.  See above note.
    //
    // if (node_.is_null ()) {
    //   node_ = KokkosClassic::DefaultNode::getDefaultNode ();
    // }
  }

  MpiPlatform<KokkosClassic::DefaultNode::DefaultNodeType>::
  MpiPlatform (const Teuchos::RCP<NodeType>& node,
               const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> >& rawMpiComm)
    : comm_ (Teuchos::null),
      node_ (node)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      rawMpiComm.is_null (), std::invalid_argument, "Tpetra::MpiPlatform "
      "constructor: The input RCP<OpaqueWrapper<MPI_Comm> > is null.  That "
      "means something different than MPI_COMM_NULL.  If you want to give "
      "MPI_COMM_NULL to this constructor, please wrap MPI_COMM_NULL in a "
      "nonnull Teuchos::OpaqueWrapper by using the "
      "Teuchos::opaqueWrapper<MPI_Comm>() nonmember constructor.");
    comm_ = Teuchos::createMpiComm<int> (rawMpiComm);

    // mfh 29 Jun 2014: Don't initialize the Node yet.  See above note.
    //
    // if (node_.is_null ()) {
    //   node_ = KokkosClassic::DefaultNode::getDefaultNode ();
    // }
  }

  MpiPlatform<KokkosClassic::DefaultNode::DefaultNodeType>::
  MpiPlatform (int* argc,
               char*** argv,
               const Teuchos::RCP<NodeType>& node,
               const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> >& rawMpiComm)
    : comm_ (Teuchos::null),
      node_ (node)
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

    // mfh 29 Jun 2014: Don't initialize the Node yet.  See above note.
    //
    // if (node_.is_null ()) {
    //   node_ = KokkosClassic::DefaultNode::getDefaultNode ();
    // }
  }

  MpiPlatform<KokkosClassic::DefaultNode::DefaultNodeType>::
  MpiPlatform (const Teuchos::RCP<NodeType>& node, MPI_Comm rawMpiComm)
    : comm_ (Teuchos::createMpiComm<int> (Teuchos::opaqueWrapper<MPI_Comm> (rawMpiComm))),
      node_ (node)
  {
    // mfh 29 Jun 2014: Don't initialize the Node yet.  See above note.
    //
    // if (node_.is_null ()) {
    //   node_ = KokkosClassic::DefaultNode::getDefaultNode ();
    // }
  }

  MpiPlatform<KokkosClassic::DefaultNode::DefaultNodeType>::
  MpiPlatform (int* argc,
               char*** argv,
               const Teuchos::RCP<NodeType>& node,
               MPI_Comm rawMpiComm)
    : comm_ (Teuchos::null),
      node_ (node)
  {
    initialize (argc, argv, rawMpiComm);
    comm_ = getDefaultComm ();

    // mfh 29 Jun 2014: Don't initialize the Node yet.  See above note.
    //
    // if (node_.is_null ()) {
    //   node_ = KokkosClassic::Details::getNode<Node> ();
    // }
  }

  MpiPlatform<KokkosClassic::DefaultNode::DefaultNodeType>::
  ~MpiPlatform () {}

  Teuchos::RCP<const Teuchos::Comm<int> >
  MpiPlatform<KokkosClassic::DefaultNode::DefaultNodeType>::
  getComm () const {
    TEUCHOS_TEST_FOR_EXCEPTION(
      comm_.is_null (), std::logic_error, "Tpetra::MpiPlatform::getComm: "
      "The default communicator is null.  This should never happen.  "
      "Please report this bug to the Tpetra developers.");
    return comm_;
  }

  Teuchos::RCP<MpiPlatform<KokkosClassic::DefaultNode::DefaultNodeType>::NodeType>
  MpiPlatform<KokkosClassic::DefaultNode::DefaultNodeType>::
  getNode () const
  {
    typedef MpiPlatform<NodeType> this_type;
    if (node_.is_null ()) {
      // NOTE (mfh 29 Jun 2014): Creating an instance of one of the
      // new Kokkos wrapper Nodes _must_ call Kokkos::initialize.
      // If Kokkos has not been initialized yet, this may result in
      // Kokkos being initialized correctly, since we have no way to
      // pass it the command-line arguments at this point.  This is
      // why we should prefer the *Platform constructors that take
      // argc and argv, since they can (and do) call
      // Kokkos::initialize (by calling Tpetra::initialize).
      //
      // mfh 29 Jun 2014: We're only keeping the *Platform classes
      // for backwards compatibility anyway, so I don't feel bad
      // about the const_cast here.
      const_cast<this_type*> (this)->node_ =
        KokkosClassic::DefaultNode::getDefaultNode ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        node_.is_null (), std::logic_error, "Tpetra::MpiPlatform::getNode: "
        "KokkosClassic::DefaultNode::getDefaultNode() returned null.  "
        "This should never happen.  "
        "Please report this bug to the Tpetra developers.");
    }
    return node_;
  }

} // namespace Tpetra

#endif // HAVE_TPETRA_MPI
