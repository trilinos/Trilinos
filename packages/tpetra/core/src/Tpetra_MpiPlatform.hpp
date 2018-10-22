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

#ifndef TPETRA_MPIPLATFORM_HPP
#define TPETRA_MPIPLATFORM_HPP

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_Core.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_Describable.hpp>

namespace Tpetra {

  /// \brief Implementation of the Platform concept for MPI-based
  ///   platforms.
  ///
  ///  \warning This class is DEPRECATED and will be REMOVED SOON.  Do
  ///    not use <tt>*Platform</tt> classes any more.  To initialize
  ///    Tpetra, include <tt>Tpetra_Core.hpp</tt> and use
  ///    Tpetra::ScopeGuard, or Tpetra::initialize and
  ///    Tpetra::finalize.  To get Tpetra's default Comm instance,
  ///    include <tt>Tpetra_Core.hpp</tt> and call
  ///    <tt>Tpetra::getDefaultComm()</tt>.  For the default Node
  ///    type, use <tt>Tpetra::Map<>::node_type</tt>.  Do not create
  ///    Node instances yourself.  It is OK for Node instances to be
  ///    null.
  template <class Node>
  class TPETRA_DEPRECATED MpiPlatform : public Teuchos::Describable {
  public:
    //! @name Typedefs
    //@{

    //! Kokkos Node type; the template parameter of this class.
    typedef Node NodeType;

    //@}
    //! @name Constructors and destructor
    //@{

    /// \brief Constructor that accepts a Kokkos Node.
    ///
    /// This version of the constructor uses MPI_COMM_WORLD as the
    /// communicator.  It is declared "explicit" to forbid silent
    /// conversions via assignment from the Node instance to an
    /// MpiPlatform.
    ///
    /// \param node [in/out] The Kokkos Node instance.
    ///
    /// Node will be deprecated and removed in favor of
    /// Kokkos::Device, so the Node can and should always be
    /// Teuchos::null.
    explicit MpiPlatform (const Teuchos::RCP<NodeType>& node) :
      comm_ (new Teuchos::MpiComm<int> (MPI_COMM_WORLD))
    {}

    /// \brief Constructor that accepts the same arguments as
    ///   MPI_Init(), plus a Kokkos Node.
    ///
    /// This version of the constructor uses MPI_COMM_WORLD as the
    /// default communicator.
    ///
    /// \param argc [in/out] First argument of MPI_Init().
    /// \param argv [in/out] Second argument of MPI_Init().
    /// \param node [in/out] The Kokkos Node instance.
    ///
    /// Node will be deprecated and removed in favor of
    /// Kokkos::Device, so the Node can and should always be
    /// Teuchos::null.
    MpiPlatform (int* argc,
                 char*** argv,
                 const Teuchos::RCP<NodeType>& /* node */) :
      comm_ (Teuchos::null)
    {
      initialize (argc, argv);
      comm_ = getDefaultComm ();
    }

    /// \brief Constructor that accepts a Kokkos Node and a wrapped
    ///   MPI communicator.
    ///
    /// This version of the constructor accepts an arbitrary MPI
    /// communicator.  It requires that you first wrap the MPI
    /// communicator in a Teuchos::OpaqueWrapper.  This is helpful if
    /// you want to "free" the communicator automatically after use.
    ///
    /// If you just have a raw \c MPI_Comm and you want to be
    /// responsible for freeing it after use, use the constructor
    /// version that takes a raw \c MPI_Comm.  Otherwise, see the
    /// documentation of Teuchos::OpaqueWrapper to learn how to wrap
    /// an \c MPI_Comm and how to set the wrapper to free the
    /// communicator automatically after use.
    ///
    /// \param node [in/out] The Kokkos Node instance.
    /// \param rawMpiComm [in] The MPI communicator, wrapped in a
    ///   <tt>Teuchos::OpaqueWrapper</tt>.
    ///
    /// Node will be deprecated and removed in favor of
    /// Kokkos::Device, so the Node can and should always be
    /// Teuchos::null.
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

    /// \brief Constructor that accepts the same arguments as
    ///   MPI_Init(), plus a Kokkos Node and a wrapped MPI
    ///   communicator.
    ///
    /// This version of the constructor accepts an arbitrary MPI
    /// communicator.  It requires that you first wrap the MPI
    /// communicator in a Teuchos::OpaqueWrapper.  This is helpful if
    /// you want to "free" the communicator automatically after use.
    ///
    /// If you just have a raw \c MPI_Comm and you want to be
    /// responsible for freeing it after use, use the constructor
    /// version that takes a raw \c MPI_Comm.  Otherwise, see the
    /// documentation of Teuchos::OpaqueWrapper to learn how to wrap
    /// an \c MPI_Comm and how to set the wrapper to free the
    /// communicator automatically after use.
    ///
    /// \param argc [in/out] First argument of MPI_Init().
    /// \param argv [in/out] Second argument of MPI_Init().
    /// \param node [in/out] The Kokkos Node instance.
    /// \param rawMpiComm [in] The MPI communicator, wrapped in a
    ///   <tt>Teuchos::OpaqueWrapper</tt>.
    ///
    /// Node will be deprecated and removed in favor of
    /// Kokkos::Device, so the Node can and should always be
    /// Teuchos::null.
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

    /// \brief Constructor that accepts a Kokkos Node and a raw MPI
    ///   communicator.
    ///
    /// This version of the constructor accepts an arbitrary "raw"
    /// (not wrapped) MPI communicator.  You are responsible for
    /// freeing the MPI communicator after use, if necessary.
    ///
    /// \param node [in/out] The Kokkos Node instance.
    /// \param rawMpiComm [in] The "raw" (not wrapped) MPI
    ///   communicator.
    ///
    /// Node will be deprecated and removed in favor of
    /// Kokkos::Device, so the Node can and should always be
    /// Teuchos::null.
    MpiPlatform (const Teuchos::RCP<NodeType>& node, MPI_Comm rawMpiComm)
      : comm_ (new Teuchos::MpiComm<int> (rawMpiComm))
    {}

    /// \brief Constructor that accepts the same arguments as
    ///   MPI_Init(), plus a Kokkos Node and a raw MPI communicator.
    ///
    /// This version of the constructor accepts an arbitrary "raw"
    /// (not wrapped) MPI communicator.  You are responsible for
    /// freeing the MPI communicator after use, if necessary.
    ///
    /// \param argc [in/out] First argument of MPI_Init().
    /// \param argv [in/out] Second argument of MPI_Init().
    /// \param node [in/out] The Kokkos Node instance.
    /// \param rawMpiComm [in] The "raw" (not wrapped) MPI
    ///   communicator.
    ///
    /// Node will be deprecated and removed in favor of
    /// Kokkos::Device, so the Node can and should always be
    /// Teuchos::null.
    MpiPlatform (int* argc,
                 char*** argv,
                 const Teuchos::RCP<NodeType>& /* node */,
                 MPI_Comm rawMpiComm)
      : comm_ (Teuchos::null)
    {
      initialize (argc, argv, rawMpiComm);
      comm_ = getDefaultComm ();
    }

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~MpiPlatform () {}

    //@}
    //! @name Methods to access the communicator and Kokkos Node.
    //@{

    //! The Teuchos::Comm instance with which this object was created.
    Teuchos::RCP<const Teuchos::Comm<int> > getComm () const {
      TEUCHOS_TEST_FOR_EXCEPTION(
        comm_.is_null (), std::logic_error, "Tpetra::MpiPlatform::getComm: "
        "The default communicator is null.  This should never happen.  "
        "Please report this bug to the Tpetra developers.");
      return comm_;
    }

    /// \brief The default Kokkos Node instance.
    ///
    /// Node will be deprecated and removed in favor of
    /// Kokkos::Device.  Thus, the Node may be Teuchos::null, and
    /// users should not depend on \c getNode.get().
    Teuchos::RCP<NodeType> getNode () const {
      return Teuchos::rcp (new NodeType);
    }

    //@}
  protected:
    //! Teuchos::Comm object instantiated for the platform.
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;

  private:
    //! Unimplemented copy constructor (syntactically forbidden).
    MpiPlatform (const MpiPlatform<NodeType>& platform);
    //! Unimplemented assignment operator (syntactically forbidden).
    MpiPlatform& operator= (const MpiPlatform<NodeType>& platform);
  };

  /// \class MpiPlatform< ::Tpetra::Details::DefaultTypes::node_type>
  /// \brief MpiPlatform specialization for the default Node type.
  ///
  ///  \warning This class is DEPRECATED and will be REMOVED SOON.  Do
  ///    not use <tt>*Platform</tt> classes any more.  To initialize
  ///    Tpetra, include <tt>Tpetra_Core.hpp</tt> and use
  ///    Tpetra::ScopeGuard, or Tpetra::initialize and
  ///    Tpetra::finalize.  To get Tpetra's default Comm instance,
  ///    include <tt>Tpetra_Core.hpp</tt> and call
  ///    <tt>Tpetra::getDefaultComm()</tt>.  For the default Node
  ///    type, use <tt>Tpetra::Map<>::node_type</tt>.  Do not create
  ///    Node instances yourself.  It is OK for Node instances to be
  ///    null.
  template <>
  class TPETRA_DEPRECATED MpiPlatform< ::Tpetra::Details::DefaultTypes::node_type> :
    public Teuchos::Describable {
  public:
    //! @name Typedefs
    //@{

    //! Kokkos Node type; the template parameter of this class.
    typedef ::Tpetra::Details::DefaultTypes::node_type NodeType;

    //@}
    //! @name Constructors and destructor
    //@{

    //! Default constructor: uses Kokkos default node and MPI_COMM_WORLD.
    MpiPlatform ();

    /// \brief Constructor that accepts the same arguments as
    ///   MPI_Init().
    ///
    /// \param argc [in/out] First argument of MPI_Init().
    /// \param argv [in/out] Second argument of MPI_Init().
    MpiPlatform (int* argc, char*** argv);

    /// \brief Constructor that accepts a Kokkos Node.
    ///
    /// This version of the constructor uses MPI_COMM_WORLD as the
    /// communicator.  It is declared "explicit" to forbid silent type
    /// conversions from the Node instance to an MpiPlatform.  (A
    /// single-argument constructor that is not declared "explicit"
    /// defines a type conversion method from the input type to the
    /// constructor's class's type.)  The "explicit" declaration does
    /// not affect typical use of this constructor.
    ///
    /// \param node [in/out] The Kokkos Node instance.
    ///
    /// Node will be deprecated and removed in favor of
    /// Kokkos::Device, so the Node can and should always be
    /// Teuchos::null.
    explicit MpiPlatform (const Teuchos::RCP<NodeType>& node);

    /// \brief Constructor that accepts the same arguments as
    ///   MPI_Init(), plus a Kokkos Node.
    ///
    /// This version of the constructor uses MPI_COMM_WORLD as the
    /// default communicator.
    ///
    /// \param argc [in/out] First argument of MPI_Init().
    /// \param argv [in/out] Second argument of MPI_Init().
    /// \param node [in/out] The Kokkos Node instance.
    ///
    /// Node will be deprecated and removed in favor of
    /// Kokkos::Device, so the Node can and should always be
    /// Teuchos::null.
    MpiPlatform (int* argc, char*** argv, const Teuchos::RCP<NodeType>& node);

    /// \brief Constructor that accepts a Kokkos Node and a wrapped
    ///   MPI communicator.
    ///
    /// This version of the constructor accepts an arbitrary MPI
    /// communicator.  It requires that you first wrap the MPI
    /// communicator in a Teuchos::OpaqueWrapper.  This is helpful if
    /// you want to "free" the communicator automatically after use.
    ///
    /// If you just have a raw \c MPI_Comm and you want to be
    /// responsible for freeing it after use, use the constructor
    /// version that takes a raw \c MPI_Comm.  Otherwise, see the
    /// documentation of Teuchos::OpaqueWrapper to learn how to wrap
    /// an \c MPI_Comm and how to set the wrapper to free the
    /// communicator automatically after use.
    ///
    /// \param node [in/out] The Kokkos Node instance.
    /// \param rawMpiComm [in] The MPI communicator, wrapped in a
    ///   <tt>Teuchos::OpaqueWrapper</tt>.
    ///
    /// Node will be deprecated and removed in favor of
    /// Kokkos::Device, so the Node can and should always be
    /// Teuchos::null.
    MpiPlatform (const Teuchos::RCP<NodeType>& node,
                 const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> >& rawMpiComm);

    /// \brief Constructor that accepts the same arguments as
    ///   MPI_Init(), plus a Kokkos Node and a wrapped MPI
    ///   communicator.
    ///
    /// This version of the constructor accepts an arbitrary MPI
    /// communicator.  It requires that you first wrap the MPI
    /// communicator in a Teuchos::OpaqueWrapper.  This is helpful if
    /// you want to "free" the communicator automatically after use.
    ///
    /// If you just have a raw \c MPI_Comm and you want to be
    /// responsible for freeing it after use, use the constructor
    /// version that takes a raw \c MPI_Comm.  Otherwise, see the
    /// documentation of Teuchos::OpaqueWrapper to learn how to wrap
    /// an \c MPI_Comm and how to set the wrapper to free the
    /// communicator automatically after use.
    ///
    /// \param argc [in/out] First argument of MPI_Init().
    /// \param argv [in/out] Second argument of MPI_Init().
    /// \param node [in/out] The Kokkos Node instance.
    /// \param rawMpiComm [in] The MPI communicator, wrapped in a
    ///   <tt>Teuchos::OpaqueWrapper</tt>.
    ///
    /// Node will be deprecated and removed in favor of
    /// Kokkos::Device, so the Node can and should always be
    /// Teuchos::null.
    MpiPlatform (int* argc,
                 char*** argv,
                 const Teuchos::RCP<NodeType>& node,
                 const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> >& rawMpiComm);

    /// \brief Constructor that accepts a Kokkos Node and a raw MPI
    ///   communicator.
    ///
    /// This version of constructor accepts an arbitrary "raw" (not
    /// wrapped) MPI communicator.  You are responsible for freeing
    /// the MPI communicator after use, if necessary.
    ///
    /// \param node [in/out] The Kokkos Node instance.
    /// \param rawMpiComm [in] The "raw" (not wrapped) MPI
    ///   communicator.
    ///
    /// Node will be deprecated and removed in favor of
    /// Kokkos::Device, so the Node can and should always be
    /// Teuchos::null.
    MpiPlatform (const Teuchos::RCP<NodeType>& node, MPI_Comm rawMpiComm);

    /// \brief Constructor that accepts the same arguments as
    ///   MPI_Init(), plus a Kokkos Node and a raw MPI communicator.
    ///
    /// This version of the constructor accepts an arbitrary "raw"
    /// (not wrapped) MPI communicator.  You are responsible for
    /// freeing the MPI communicator after use, if necessary.
    ///
    /// \param argc [in/out] First argument of MPI_Init().
    /// \param argv [in/out] Second argument of MPI_Init().
    /// \param node [in/out] The Kokkos Node instance.
    /// \param rawMpiComm [in] The "raw" (not wrapped) MPI
    ///   communicator.
    ///
    /// Node will be deprecated and removed in favor of
    /// Kokkos::Device, so the Node can and should always be
    /// Teuchos::null.
    MpiPlatform (int* argc,
                 char*** argv,
                 const Teuchos::RCP<NodeType>& node,
                 MPI_Comm rawMpiComm);

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~MpiPlatform ();

    //@}
    //! @name Methods to access the communicator and Kokkos Node.
    //@{

    //! The Teuchos::Comm instance with which this object was created.
    Teuchos::RCP<const Teuchos::Comm<int> > getComm () const;

    /// \brief The default Kokkos Node instance.
    ///
    /// Node will be deprecated and removed in favor of
    /// Kokkos::Device, so the Node can and should always be
    /// Teuchos::null.
    Teuchos::RCP<NodeType> getNode () const;

    //@}
  private:
    //! Unimplemented copy constructor (syntactically forbidden).
    MpiPlatform (const MpiPlatform<NodeType>& platform);

    //! Unimplemented assignment operator (syntactically forbidden).
    MpiPlatform& operator= (const MpiPlatform<NodeType>& platform);

  protected:
    //! Teuchos::Comm object instantiated for the platform.
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;
  };

} // namespace Tpetra

#endif // TPETRA_MPIPLATFORM_HPP
