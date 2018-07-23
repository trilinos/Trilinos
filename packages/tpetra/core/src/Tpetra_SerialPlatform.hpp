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

#ifndef TPETRA_SERIALPLATFORM_HPP
#define TPETRA_SERIALPLATFORM_HPP

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_Core.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_DefaultSerialComm.hpp>
#include <Teuchos_Describable.hpp>

namespace Tpetra {

  /// \brief Implementation of the Platform concept for non-MPI platforms.
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
  class TPETRA_DEPRECATED SerialPlatform : public Teuchos::Describable {
  public:
    //! @name Typedefs
    //@{

    //! Kokkos Node type; the template parameter of this class.
    typedef Node NodeType;

    //@}
    //! @name Constructors and destructor
    //@{

    /// \brief Constructor that accepts the same arguments as
    ///   Tpetra::initialize().
    ///
    /// \param argc [in/out] First argument of Tpetra::initialize().
    /// \param argv [in/out] Second argument of Tpetra::initialize().
    explicit SerialPlatform (int* argc, char*** argv) :
      comm_ (Teuchos::null)
    {
      initialize (argc, argv);
      comm_ = getDefaultComm ();
    }

    /// \brief Constructor that accepts a Kokkos Node instance.
    ///
    /// Node will be deprecated and removed in favor of
    /// Kokkos::Device, so the Node can and should always be
    /// Teuchos::null.
    explicit SerialPlatform (const Teuchos::RCP<NodeType>& /* node */) :
      comm_ (Teuchos::rcp (new Teuchos::SerialComm<int> ()))
    {
    }

    /// \brief Constructor that accepts the same arguments as
    ///   Tpetra::initialize(), plus a Kokkos Node.
    ///
    /// \param argc [in/out] First argument of Tpetra::initialize().
    /// \param argv [in/out] Second argument of Tpetra::initialize().
    /// \param node [in/out] The Kokkos Node instance.
    ///
    /// Node will be deprecated and removed in favor of
    /// Kokkos::Device, so the Node can and should always be
    /// Teuchos::null.
    explicit SerialPlatform (int* argc, char*** argv,
                             const Teuchos::RCP<NodeType>& /* node */) :
      comm_ (Teuchos::null)
    {
      initialize (argc, argv);
      comm_ = getDefaultComm ();
    }

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~SerialPlatform () {}

    //@}
    //! @name Methods to access the communicator and Kokkos Node.
    //@{

    //! The Teuchos::Comm instance with which this object was created.
    Teuchos::RCP<const Teuchos::Comm<int> > getComm () const {
      return comm_;
    }

    /// \brief The Kokkos Node instance.
    ///
    /// Since Node will be deprecated and removed in favor of
    /// Kokkos::Device, this method may return Teuchos::null.
    Teuchos::RCP<Node> getNode () const {
      return Teuchos::rcp (new Node);
    }

    //@}
  protected:
    //! Teuchos::Comm object instantiated for the platform.
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;

  private:
    //! Unimplemented copy constructor (syntactically forbidden).
    SerialPlatform (const SerialPlatform<NodeType>& platform);
    //! Unimplemented assignment operator (syntactically forbidden).
    SerialPlatform& operator= (const SerialPlatform<NodeType>& platform);
  };

  /// \class SerialPlatform<Tpetra::Details::DefaultTypes::node_type>
  /// \brief SerialPlatform specialization for the default Node type.
  ///
  /// \note Tpetra::Details::DefaultTypes::node_type is a typedef, and
  ///   may have a different type, depending on Trilinos' build
  ///   options.
  template <>
  class TPETRA_DEPRECATED SerialPlatform<Tpetra::Details::DefaultTypes::node_type> :
    public Teuchos::Describable {
  public:
    //! @name Typedefs
    //@{

    //! Kokkos Node type; the template parameter of this class.
    typedef Tpetra::Details::DefaultTypes::node_type NodeType;

    //@}
    //! @name Constructors and destructor
    //@{

    /// \brief Default constructor: uses Kokkos default node.
    ///
    /// The specialization of SerialPlatform for the default Node type
    /// includes a default constructor.
    SerialPlatform () :
      comm_ (Teuchos::rcp (new Teuchos::SerialComm<int> ()))
    {}

    /// \brief Constructor that accepts the same arguments as
    ///   Tpetra::initialize().
    ///
    /// \param argc [in/out] First argument of Tpetra::initialize().
    /// \param argv [in/out] Second argument of Tpetra::initialize().
    explicit SerialPlatform (int* argc, char*** argv) :
      comm_ (Teuchos::null)
    {
      initialize (argc, argv);
      comm_ = getDefaultComm ();
    }

    /// \brief Constructor that accepts a Kokkos Node.
    ///
    /// This version of the constructor is declared "explicit" to
    /// forbid silent type conversions from the Node instance to a
    /// SerialPlatform.  (A single-argument constructor that is not
    /// declared "explicit" defines a type conversion method from the
    /// input type to the constructor's class's type.)  The "explicit"
    /// declaration does not affect typical use of this constructor.
    ///
    /// Node will be deprecated and removed in favor of
    /// Kokkos::Device, so the Node can and should always be
    /// Teuchos::null.
    explicit SerialPlatform (const Teuchos::RCP<NodeType>& /* node */) :
      comm_ (Teuchos::rcp (new Teuchos::SerialComm<int> ()))
    {
    }

    /// \brief Constructor that accepts the same arguments as
    ///   Tpetra::initialize(), plus a Kokkos Node.
    ///
    /// \param argc [in/out] First argument of Tpetra::initialize().
    /// \param argv [in/out] Second argument of Tpetra::initialize().
    /// \param node [in/out] The Kokkos Node instance.
    ///
    /// Node will be deprecated and removed in favor of
    /// Kokkos::Device, so the Node can and should always be
    /// Teuchos::null.
    explicit SerialPlatform (int* argc, char*** argv,
                             const Teuchos::RCP<NodeType>& /* node */) :
      comm_ (Teuchos::null)
    {
      initialize (argc, argv);
      comm_ = getDefaultComm ();
    }

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~SerialPlatform () {}

    //@}
    //! @name Methods to access the communicator and Kokkos Node.
    //@{

    //! The Teuchos::Comm instance with which this object was created.
    Teuchos::RCP<const Teuchos::Comm<int> > getComm() const {
      return comm_;
    }

    /// \brief The Kokkos Node instance.
    ///
    /// Node will be deprecated and removed in favor of
    /// Kokkos::Device, so the Node can and should always be
    /// Teuchos::null.
    Teuchos::RCP<Tpetra::Details::DefaultTypes::node_type>
    getNode () const
    {
      return Teuchos::rcp (new Tpetra::Details::DefaultTypes::node_type);
    }

    //@}
  private:
    //! Unimplemented copy constructor (syntactically forbidden).
    SerialPlatform (const SerialPlatform<NodeType>& platform);

    //! Unimplemented assignment operator (syntactically forbidden).
    SerialPlatform& operator= (const SerialPlatform<NodeType>& platform);

  protected:
    //! Teuchos::Comm object instantiated for the platform.
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;
  };

} // namespace Tpetra

#endif // TPETRA_SERIALPLATFORM_HPP
