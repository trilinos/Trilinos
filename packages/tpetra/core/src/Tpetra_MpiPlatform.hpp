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
  /// \warning This class will be DEPRECATED, in favor of the
  ///   initialize() functions in Tpetra_Core.hpp.  Please use those
  ///   functions for safe, consistent initialization of MPI and
  ///   Kokkos, on which Tpetra depends.  If you must use this class,
  ///   please prefer the constructors that take \c argc and \c argv.
  ///   Those constructors will call initialize() for you.
  ///
  /// MpiPlatform is an implementation of Tpetra's Platform concept.
  /// Classes implementing the Platform concept are templated on the
  /// Kokkos Node type.  They have at least the following public
  /// interface:
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
  /// MpiPlatform also has a constructor that accepts an MPI
  /// communicator, over which the application using the platform
  /// should perform communication.  The default communicator is
  /// MPI_COMM_WORLD.  MpiPlatform is only available if Trilinos was
  /// built with MPI.
  template <class Node>
  class MpiPlatform : public Teuchos::Describable {
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
    /// \param node [in/out] The Kokkos Node instance.  If null, this
    ///   class will create a Node with default parameters, at some
    ///   time no later than during the first call to getNode().
    explicit MpiPlatform (const Teuchos::RCP<NodeType>& node) :
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
      //   node_ = KokkosClassic::Details::getNode<NodeType> ();
      // }
    }

    /// \brief Constructor that accepts the same arguments as
    ///   MPI_Init(), plus a Kokkos Node.
    ///
    /// This version of the constructor uses MPI_COMM_WORLD as the
    /// default communicator.
    ///
    /// \param argc [in/out] First argument of MPI_Init().
    /// \param argv [in/out] Second argument of MPI_Init().
    /// \param node [in/out] The Kokkos Node instance.  If null, this
    ///   class will create a Node with default parameters, at some
    ///   time no later than during the first call to getNode().
    MpiPlatform (int* argc,
                 char*** argv,
                 const Teuchos::RCP<NodeType>& node) :
      comm_ (Teuchos::null),
      node_ (node)
    {
      initialize (argc, argv);
      comm_ = getDefaultComm ();

      // mfh 29 Jun 2014: Don't initialize the Node yet.  See above note.
      //
      // if (node_.is_null ()) {
      //   node_ = KokkosClassic::Details::getNode<Node> ();
      // }
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
    /// \param node [in/out] The Kokkos Node instance.  If null, this
    ///   class will create a Node with default parameters, at some
    ///   time no later than during the first call to getNode().
    /// \param rawMpiComm [in] The MPI communicator, wrapped in a
    ///   <tt>Teuchos::OpaqueWrapper</tt>.
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
      //   node_ = KokkosClassic::Details::getNode<NodeType> ();
      // }
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
    /// \param node [in/out] The Kokkos Node instance.  If null, this
    ///   class will create a Node with default parameters, at some
    ///   time no later than during the first call to getNode().
    /// \param rawMpiComm [in] The MPI communicator, wrapped in a
    ///   <tt>Teuchos::OpaqueWrapper</tt>.
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
      //   node_ = KokkosClassic::Details::getNode<NodeType> ();
      // }
    }

    /// \brief Constructor that accepts a Kokkos Node and a raw MPI
    ///   communicator.
    ///
    /// This version of the constructor accepts an arbitrary "raw"
    /// (not wrapped) MPI communicator.  You are responsible for
    /// freeing the MPI communicator after use, if necessary.
    ///
    /// \param node [in/out] The Kokkos Node instance.  If null, this
    ///   class will create a Node with default parameters, at some
    ///   time no later than during the first call to getNode().
    /// \param rawMpiComm [in] The "raw" (not wrapped) MPI
    ///   communicator.
    MpiPlatform (const Teuchos::RCP<NodeType>& node, MPI_Comm rawMpiComm)
      : comm_ (Teuchos::createMpiComm<int> (Teuchos::opaqueWrapper<MPI_Comm> (rawMpiComm))),
        node_ (node)
    {
      // mfh 29 Jun 2014: Don't initialize the Node yet.  See above note.
      //
      // if (node_.is_null ()) {
      //   node_ = KokkosClassic::Details::getNode<NodeType> ();
      // }
    }

    /// \brief Constructor that accepts the same arguments as
    ///   MPI_Init(), plus a Kokkos Node and a raw MPI communicator.
    ///
    /// This version of the constructor accepts an arbitrary "raw"
    /// (not wrapped) MPI communicator.  You are responsible for
    /// freeing the MPI communicator after use, if necessary.
    ///
    /// \param argc [in/out] First argument of MPI_Init().
    /// \param argv [in/out] Second argument of MPI_Init().
    /// \param node [in/out] The Kokkos Node instance.  If null, this
    ///   class will create a Node with default parameters, at some
    ///   time no later than during the first call to getNode().
    /// \param rawMpiComm [in] The "raw" (not wrapped) MPI
    ///   communicator.
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
      //   node_ = KokkosClassic::Details::getNode<NodeType> ();
      // }
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
    /// If a constructor was called with a nonnull Node instance, this
    ///   just returns that instance.  Otherwise, this method returns
    ///   a Node set up with default parameters.  It only initializes
    ///   that Node once, no matter how many times this method is
    ///   called.
    Teuchos::RCP<NodeType> getNode () const {
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
    MpiPlatform (const MpiPlatform<NodeType>& platform);
    //! Unimplemented assignment operator (syntactically forbidden).
    MpiPlatform& operator= (const MpiPlatform<NodeType>& platform);
  };

  /// \class MpiPlatform<Tpetra::Details::DefaultTypes::node_type>
  /// \brief MpiPlatform specialization for the default Node type.
  ///
  /// \note Tpetra::Details::DefaultTypes::node_type is a typedef,
  ///   and may have a different type, depending on Trilinos' build
  ///   options.
  template <>
  class MpiPlatform<Tpetra::Details::DefaultTypes::node_type> :
    public Teuchos::Describable {
  public:
    //! @name Typedefs
    //@{

    //! Kokkos Node type; the template parameter of this class.
    typedef Tpetra::Details::DefaultTypes::node_type NodeType;

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
    /// This specialization of MpiPlatform for the default node type
    /// will instantiate a default Node if node.is_null().
    ///
    /// \param node [in/out] The Kokkos Node instance.  If null, this
    ///   class will create a Node with default parameters, at some
    ///   time no later than during the first call to getNode().
    explicit MpiPlatform (const Teuchos::RCP<NodeType>& node);

    /// \brief Constructor that accepts the same arguments as
    ///   MPI_Init(), plus a Kokkos Node.
    ///
    /// This version of the constructor uses MPI_COMM_WORLD as the
    /// default communicator.
    ///
    /// \param argc [in/out] First argument of MPI_Init().
    /// \param argv [in/out] Second argument of MPI_Init().
    /// \param node [in/out] The Kokkos Node instance.  If null, this
    ///   class will create a Node with default parameters, at some
    ///   time no later than during the first call to getNode().
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
    /// \param node [in/out] The Kokkos Node instance.  If null, this
    ///   class will create a Node with default parameters, at some
    ///   time no later than during the first call to getNode().
    /// \param rawMpiComm [in] The MPI communicator, wrapped in a
    ///   <tt>Teuchos::OpaqueWrapper</tt>.
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
    /// \param node [in/out] The Kokkos Node instance.  If null, this
    ///   class will create a Node with default parameters, at some
    ///   time no later than during the first call to getNode().
    /// \param rawMpiComm [in] The MPI communicator, wrapped in a
    ///   <tt>Teuchos::OpaqueWrapper</tt>.
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
    /// \param node [in/out] The Kokkos Node instance.  If null, this
    ///   class will create a Node with default parameters, at some
    ///   time no later than during the first call to getNode().
    /// \param rawMpiComm [in] The "raw" (not wrapped) MPI
    ///   communicator.
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
    /// \param node [in/out] The Kokkos Node instance.  If null, this
    ///   class will create a Node with default parameters, at some
    ///   time no later than during the first call to getNode().
    /// \param rawMpiComm [in] The "raw" (not wrapped) MPI
    ///   communicator.
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
    /// If a constructor was called with a nonnull Node instance, this
    ///   just returns that instance.  Otherwise, this method returns
    ///   a Node set up with default parameters.  It only initializes
    ///   that Node once, no matter how many times this method is
    ///   called.
    Teuchos::RCP<NodeType> getNode () const;

    //@}
  private:
    //! Unimplemented copy constructor (syntactically forbidden).
    MpiPlatform (const MpiPlatform<NodeType>& platform);

    //! Unimplemented assignment operator (syntactically forbidden).
    MpiPlatform& operator= (const MpiPlatform<NodeType>& platform);

  protected:
    //! Teuchos::Comm object instantiated for the platform.
    RCP<const Teuchos::Comm<int> > comm_;

    //! Kokkos Node object instantiated for the platform.
    RCP<NodeType> node_;
  };

} // namespace Tpetra

#endif // TPETRA_MPIPLATFORM_HPP
