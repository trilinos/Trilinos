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
  /// \warning This class will be DEPRECATED, in favor of the
  ///   initialize() functions in Tpetra_Core.hpp.  Please use those
  ///   functions for safe, consistent initialization Kokkos, on which
  ///   Tpetra depends.  If you must use this class, please prefer the
  ///   constructors that take \c argc and \c argv.  Those
  ///   constructors will call initialize() for you.
  ///
  /// SerialPlatform is an implementation of Tpetra's Platform
  /// concept.  Classes implementing the Platform concept are
  /// templated on the Kokkos Node type.  They have at least the
  /// following public interface:
  /// \code
  /// // This is not a real class; it just illustrates the concept.
  /// template<class Node>
  /// class Platform {
  /// public:
  ///   typedef Node NodeType;
  ///   explicit Platform (const RCP<Node>& node);
  ///   RCP<const Comm<int> > getComm() const;
  ///   RCP<Node> getNode() const;
  /// };
  /// \endcode
  /// SerialPlatform uses a "communicator" containing one process.  It
  /// is available whether or not Trilinos was built with MPI (the
  /// Message-Passing Interface which provides a distributed-memory
  /// parallel programming model).
  template <class Node>
  class SerialPlatform : public Teuchos::Describable {
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
      comm_ (Teuchos::null),
      node_ (Teuchos::null)
    {
      initialize (argc, argv);
      comm_ = getDefaultComm ();

      // mfh 29 Jun 2014: Don't initialize the Node yet.  This ensures
      // that (new) Kokkos won't get initialized with the wrong
      // command-line arguments, at least not until getNode() is
      // called.  Initializing Kokkos with the wrong command-line
      // arguments may result in poor performance due to the wrong
      // assignment of software threads to hardware execution units.
      //
      // if (node_.is_null ()) {
      //   node_ = KokkosClassic::Details::getNode<NodeType> ();
      // }
    }

    /// Constructor that accepts a Kokkos Node.
    ///
    /// \param node [in/out] The Kokkos Node instance.  If null, this
    ///   class will create a Node with default parameters, at some
    ///   time no later than during the first call to getNode().
    explicit SerialPlatform (const Teuchos::RCP<NodeType>& node) :
      comm_ (Teuchos::rcp (new Teuchos::SerialComm<int> ())),
      node_ (node)
    {
      // mfh 29 Jun 2014: Don't initialize the Node yet.  See above note.
    }

    /// \brief Constructor that accepts the same arguments as
    ///   Tpetra::initialize(), plus a Kokkos Node.
    ///
    /// \param argc [in/out] First argument of Tpetra::initialize().
    /// \param argv [in/out] Second argument of Tpetra::initialize().
    /// \param node [in/out] The Kokkos Node instance.  If null, this
    ///   class will create a Node with default parameters, at some
    ///   time no later than during the first call to getNode().
    explicit SerialPlatform (int* argc, char*** argv,
                             const Teuchos::RCP<NodeType>& node) :
      comm_ (Teuchos::null),
      node_ (node)
    {
      initialize (argc, argv);
      comm_ = getDefaultComm ();
      // mfh 29 Jun 2014: Don't initialize the Node yet.  See above note.
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

    //! The Kokkos Node instance with which this object was created.
    Teuchos::RCP<Node> getNode () const {
      typedef SerialPlatform<NodeType> this_type;
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
          KokkosClassic::Details::getNode<NodeType> ();
        TEUCHOS_TEST_FOR_EXCEPTION(
          node_.is_null (), std::logic_error, "Tpetra::MpiPlatform::getNode: "
          "KokkosClassic::Details::getNode<NodeType>() returned null.  "
          "This should never happen.  "
          "Please report this bug to the Tpetra developers.");
      }
      return node_;
    }

    //@}
  protected:
    //! Teuchos::Comm object instantiated for the platform.
    Teuchos::RCP<const Teuchos::Comm<int> > comm_;
    //! Kokkos Node object instantiated for the platform.
    Teuchos::RCP<NodeType> node_;

  private:
    //! Unimplemented copy constructor (syntactically forbidden).
    SerialPlatform (const SerialPlatform<NodeType>& platform);
    //! Unimplemented assignment operator (syntactically forbidden).
    SerialPlatform& operator= (const SerialPlatform<NodeType>& platform);
  };

  /// \class SerialPlatform<KokkosClassic::DefaultNode::DefaultNodeType>
  /// \brief SerialPlatform specialization for KokkosClassic::DefaultNode::DefaultNodeType.
  ///
  /// \warning KokkosClassic::DefaultNode::DefaultNodeType is a typedef, and
  ///   may have a different type, depending on Trilinos' build
  ///   options.  For example, it may be KokkosClassic::SerialNode if
  ///   Trilinos was built without a threading library, or
  ///   KokkosClassic::TPINode if Trilinos was built with Pthreads.
  ///
  /// \note In the past (up to and including the 10.8 Trilinos
  ///   release), the specialization of SerialPlatform for the default
  ///   Node type delayed instantiation of the default Node instance
  ///   until getNode() was called.  We have changed this behavior to
  ///   simplify the code and make the specialization of
  ///   SerialPlatform conform more closely to the generic version of
  ///   SerialPlatform.
  template <>
  class SerialPlatform<KokkosClassic::DefaultNode::DefaultNodeType> :
    public Teuchos::Describable {
  public:
    //! @name Typedefs
    //@{

    //! Kokkos Node type; the template parameter of this class.
    typedef KokkosClassic::DefaultNode::DefaultNodeType NodeType;

    //@}
    //! @name Constructors and destructor
    //@{

    /// \brief Default constructor: uses Kokkos default node.
    ///
    /// The specialization of SerialPlatform for the default Node type
    /// includes a default constructor.  At some point before the
    /// first call to getNode() returns, this class will create a Node
    /// with default parameters.
    SerialPlatform () :
      comm_ (Teuchos::rcp (new Teuchos::SerialComm<int> ())),
      node_ (Teuchos::null)
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

    /// \brief Constructor that accepts the same arguments as
    ///   Tpetra::initialize().
    ///
    /// \param argc [in/out] First argument of Tpetra::initialize().
    /// \param argv [in/out] Second argument of Tpetra::initialize().
    explicit SerialPlatform (int* argc, char*** argv) :
      comm_ (Teuchos::null),
      node_ (Teuchos::null)
    {
      initialize (argc, argv);
      comm_ = getDefaultComm ();
      // mfh 29 Jun 2014: Don't initialize the Node yet.  See above note.
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
    /// \param node [in/out] The Kokkos Node instance.  If null, this
    ///   class will create a Node with default parameters, at some
    ///   time no later than during the first call to getNode().
    explicit SerialPlatform (const Teuchos::RCP<NodeType>& node) :
      comm_ (Teuchos::rcp (new Teuchos::SerialComm<int> ())),
      node_ (node)
    {
      // mfh 29 Jun 2014: Don't initialize the Node yet.  See above note.
    }

    /// \brief Constructor that accepts the same arguments as
    ///   Tpetra::initialize(), plus a Kokkos Node.
    ///
    /// \param argc [in/out] First argument of Tpetra::initialize().
    /// \param argv [in/out] Second argument of Tpetra::initialize().
    /// \param node [in/out] The Kokkos Node instance.  If null, this
    ///   class will create a Node with default parameters, at some
    ///   time no later than during the first call to getNode().
    explicit SerialPlatform (int* argc, char*** argv,
                             const Teuchos::RCP<NodeType>& node) :
      comm_ (Teuchos::null),
      node_ (node)
    {
      initialize (argc, argv);
      comm_ = getDefaultComm ();
      // mfh 29 Jun 2014: Don't initialize the Node yet.  See above note.
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

    //! The Kokkos Node instance with which this object was created.
    Teuchos::RCP<KokkosClassic::DefaultNode::DefaultNodeType> getNode () const {
      typedef SerialPlatform<NodeType> this_type;
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
          KokkosClassic::Details::getNode<NodeType> ();
        TEUCHOS_TEST_FOR_EXCEPTION(
          node_.is_null (), std::logic_error, "Tpetra::MpiPlatform::getNode: "
          "KokkosClassic::Details::getNode<NodeType>() returned null.  "
          "This should never happen.  "
          "Please report this bug to the Tpetra developers.");
      }
      return node_;
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

    //! Node object instantiated for the platform.
    Teuchos::RCP<NodeType> node_;
  };

} // namespace Tpetra

#endif // TPETRA_SERIALPLATFORM_HPP
