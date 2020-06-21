// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_InputTraits.hpp
    \brief Traits for application input objects.
*/

#ifndef ZOLTAN2_INPUTTRAITS_HPP
#define ZOLTAN2_INPUTTRAITS_HPP

#include <Zoltan2_Standards.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_RowMatrix.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_RowGraph.hpp>

#ifdef HAVE_ZOLTAN2_EPETRA
#include <Epetra_CrsMatrix.h>
#include <Epetra_CrsGraph.h>
#endif

#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_RowMatrix.hpp>
#include <Xpetra_TpetraRowMatrix.hpp>
#include <Xpetra_CrsGraph.hpp>

namespace Zoltan2{

// Default local ordinal
typedef ::Tpetra::Details::DefaultTypes::local_ordinal_type default_lno_t;

// Default global ordinal
typedef ::Tpetra::Details::DefaultTypes::global_ordinal_type default_gno_t;

// Default offset type
typedef size_t default_offset_t;

// Default scalar type (for weights, coordinates)
typedef ::Tpetra::Details::DefaultTypes::scalar_type default_scalar_t;

// Until Kokkos node types are supported, use default
typedef ::Tpetra::Details::DefaultTypes::node_type default_node_t;

// Default part number type.
typedef int default_part_t;  // Restrictions in MPI interface will make it
                             // somewhat difficult to change default_part_t to
                             // long long, since we use part_t for ranks
                             // and we sometimes broadcast arrays whose
                             // size has type part_t.
                             // part_t must be a signed data type.

/*! \brief A simple class that can be the User template argument
 *             for an InputAdapter.
 *
 *  BasicUserTypes is a convenience class that provides a simple way
 *  to supply the template argument when constructing an InputAdapter.
 *
 *  Typically, a C++ user will have a templated or non-templated class that
 *  represents his or her input data.  He or she will create an InputTraits
 *  specialization for this class, and then supply that as the template
 *  argument of the InputAdapter class.  (Commonly used Trilinos classes
 *  that represent vectors, graphs and matrices already have InputTraits
 *  specializations.)
 *
 *  This makes sense if you are a C++ programmer who uses templates, but
 *  is a great deal of trouble if you are not.  In this case you can
 *  construct your InputAdapter as follows.
 *
 *  Suppose you want to construct a Zoltan2::BasicVectorAdapter object and
 *  you use \c float for vector values in your application, \c long for
 *  global identifiers, and \c int for local indices.
 *
 *  You need to determine an integral data type that Zoltan2 can use internally
 *  for global identifiers. Often this is the same data type that you use for
 *  this purpose, but if you have non-integral global identifiers (such as
 *  pair<int, int>) then Zoltan2 needs you to supply an integral data
 *  type that is large enough to globally enumerate your objects.
 *  In this example let's say that you want Zoltan2 to use \c unsigned \c long
 *  for its global Identifiers.  Then you would do the following:
 *
\code
   typedef BasicUserTypes<float, long, int, unsigned long> myTypes;
   Zoltan2::BasicVectorAdapter<myTypes> myInput({constructor arguments})
\endcode
 *
 * In particular, the BasicUserTypes template parameters are:

    \li \c scalar is the data type for element values, weights and coordinates.
    \li \c lno (local number) is the integral data type used by the application and by Zoltan2 for local indices and local counts.
    \li \c gno (global number) is the integral data type used by the application and Zoltan2 to represent global identifiers and global counts.
 */

template <typename scalar=default_scalar_t,
          typename lno=default_lno_t,
          typename gno=default_gno_t,
          typename node=default_node_t>
class BasicUserTypes{
};

/*! \brief The traits required of User input classes or structures.
 *
 *  We use the symbol \em User as an abtraction of the user's input.
 *  If the user's data structure is an instance of a C++ class,
 *  or is a C-language struct,
 *  then it should have an InputTraits specialization.  The InputAdapter
 *  is templated on the User object, and Zoltan2 uses the User object
 *  traits to get the user's data types.
 *
 *  If the User object is \em not a C++ class or C-language struct, and
 *  particularly if you do not plan to pass your User object as an
 *  argument to the InputAdapter constructor, you can use the
 *  BasicUserTypes class as your User class.
 *
 *  InputTraits already exist for:
 *
 *  \li Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>
 *  \li Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>
 *  \li Epetra_CrsMatrix
 *  \li Xpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>
 *  \li Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>
 *  \li Epetra_CrsGraph
 *  \li Xpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>
 *  \li Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>
 *  \li Tpetra::RowGraphx<LocalOrdinal,GlobalOrdinal,Node>
 *  \li Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>
 *  \li Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>
 *  \li Epetra_Vector
 *  \li Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>
 *  \li Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>
 *  \li Epetra_MultiVector
 *  \li BasicUserTypes
 *
 */
template <typename User>
struct InputTraits {

  /*! \brief The data type for weights and coordinates.
   */
  typedef default_scalar_t scalar_t;

  /*! \brief The ordinal type (e.g., int, long, int64_t) that represents
   *              local counts and local indices.
   */
  typedef default_lno_t lno_t;

  /*! \brief The ordinal type (e.g., int, long, int64_t) that can represent
   *             global counts and identifiers.
   */
  typedef default_gno_t gno_t;

  /*! \brief The data type to represent offsets.
   */
  typedef default_offset_t offset_t;

  /*! \brief The data type to represent part numbers.
   */
  typedef default_part_t part_t;

  /*! \brief  The Kokkos node type.  This is only meaningful for users
   *          of Tpetra objects.
   */
  typedef default_node_t node_t;

  /*! \brief  The name of the user's input object.
   */
  static inline std::string name() {return "InputAdapter";}
};

#ifndef DOXYGEN_SHOULD_SKIP_THIS

// This combination of macros is used to define a single line
// Z2_STATIC_ASSERT_TYPES for each InputTraits with custom template types
#define Z2_ISSAME(s,type) (std::is_same< s, type >::value)

#define Z2_STYPES(s) ( Z2_ISSAME(s,float) || \
  Z2_ISSAME(s,double) || Z2_ISSAME(s,int) )

#define Z2_LTYPES(l) ( Z2_ISSAME(l,int) ||  \
  Z2_ISSAME(l,long) || Z2_ISSAME(l,long long) || Z2_ISSAME(l,ssize_t) )

#define Z2_GTYPES(g) ( Z2_ISSAME(g,int) || Z2_ISSAME(g,long) || \
  Z2_ISSAME(g,long long) || Z2_ISSAME(g,ssize_t) || \
  Z2_ISSAME(g,unsigned int) || Z2_ISSAME(g,unsigned long) || \
  Z2_ISSAME(g,unsigned long long) || Z2_ISSAME(g,size_t) )

#define Z2_SERROR "Invalid scalar type. It must be float, double, or int."

#define Z2_LERROR "Invalid local ordinal type. It must be int, long, " \
  "long long, or ssize_t."

#define Z2_GERROR "Invalid global ordinal type. It must be int, long, " \
  "long long, ssize_t, unsigned int, unsigned long long, size_t."

#ifdef Z2_INVERT_STATIC_ASSERT_FOR_UNIT_TESTING
  #define Z2_STATIC_ASSERT_TYPES static_assert( ( !Z2_STYPES(scalar_t) || \
  !Z2_LTYPES(lno_t) || !Z2_GTYPES(gno_t) ), \
  "Inverted unit test for InputTraits was supposed to fail but did not." );
#else
  #define Z2_STATIC_ASSERT_TYPES static_assert( Z2_STYPES(scalar_t), \
  Z2_SERROR ); static_assert( Z2_LTYPES(lno_t), Z2_LERROR ); \
  static_assert( Z2_GTYPES(gno_t), Z2_GERROR );
#endif

template <typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<BasicUserTypes<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
{
  typedef Scalar        scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef LocalOrdinal  offset_t;
  typedef Zoltan2::default_part_t  part_t;
  typedef Node          node_t;
  static inline std::string name() {return "BasicUserTypes";}

  Z2_STATIC_ASSERT_TYPES // validate the types
};

template <typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
{
  typedef Scalar        scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef size_t offset_t;
  typedef Zoltan2::default_part_t  part_t;
  typedef Node          node_t;
  static inline std::string name() {return "Xpetra::CrsMatrix";}

  Z2_STATIC_ASSERT_TYPES // validate the types
};

template <typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
{
  typedef Scalar        scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef size_t offset_t;
  typedef Zoltan2::default_part_t  part_t;
  typedef Node          node_t;
  static inline std::string name() {return "Tpetra::CrsMatrix";}

  Z2_STATIC_ASSERT_TYPES // validate the types
};

#ifdef HAVE_ZOLTAN2_EPETRA
template < >
struct InputTraits<Epetra_CrsMatrix>
{
  typedef double scalar_t;
  typedef int lno_t;
  typedef int gno_t;
  typedef size_t offset_t;
  typedef Zoltan2::default_part_t  part_t;
  typedef Zoltan2::default_node_t node_t;
  static inline std::string name() {return "Epetra_CrsMatrix";}
};
#endif

template <typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<Xpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
{
  typedef Scalar        scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef size_t offset_t;
  typedef Zoltan2::default_part_t  part_t;
  typedef Node          node_t;
  static inline std::string name() {return "Xpetra::RowMatrix";}

  Z2_STATIC_ASSERT_TYPES // validate the types
};

template <typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
{
  typedef Scalar        scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef size_t offset_t;
  typedef Zoltan2::default_part_t  part_t;
  typedef Node          node_t;
  static inline std::string name() {return "Tpetra::RowMatrix";}

  Z2_STATIC_ASSERT_TYPES // validate the types
};

template <typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<Tpetra::RowGraph<LocalOrdinal,GlobalOrdinal,Node> >
{
  typedef default_scalar_t scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef size_t offset_t;
  typedef Zoltan2::default_part_t  part_t;
  typedef Node          node_t;
  static inline std::string name() {return "Tpetra::RowGraph";}
};

template <typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<Xpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >
{
  typedef default_scalar_t scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef size_t offset_t;
  typedef Zoltan2::default_part_t  part_t;
  typedef Node          node_t;
  static inline std::string name() {return "Xpetra::CrsGraph";}

  Z2_STATIC_ASSERT_TYPES // validate the types
};

template <typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >
{
  typedef default_scalar_t scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef size_t offset_t;
  typedef Zoltan2::default_part_t  part_t;
  typedef Node          node_t;
  static inline std::string name() {return "Tpetra::CrsGraph";}

  Z2_STATIC_ASSERT_TYPES // validate the types
};

#ifdef HAVE_ZOLTAN2_EPETRA
template < >
struct InputTraits<Epetra_CrsGraph>
{
  typedef double scalar_t;
  typedef int   lno_t;
  typedef int   gno_t;
  typedef size_t offset_t;
  typedef Zoltan2::default_part_t  part_t;
  typedef Zoltan2::default_node_t node_t;
  static inline std::string name() {return "Epetra_CrsGraph";}
};
#endif

template <typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
{
  typedef Scalar        scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef size_t offset_t;
  typedef Zoltan2::default_part_t  part_t;
  typedef Node          node_t;
  static inline std::string name() {return "Xpetra::Vector";}

  Z2_STATIC_ASSERT_TYPES // validate the types
};

 /*! \todo A Tpetra::Vector is a Tpetra::MultiVector - can we just
 *     define MultiVector traits only?  Ditto with Xpetra.  Test this.
 */
template <typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
{
  typedef Scalar        scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef size_t offset_t;
  typedef Zoltan2::default_part_t  part_t;
  typedef Node          node_t;
  static inline std::string name() {return "Tpetra::Vector";}

  Z2_STATIC_ASSERT_TYPES  // validate the types
};

#ifdef HAVE_ZOLTAN2_EPETRA
template < >
struct InputTraits<Epetra_Vector>
{
  typedef double scalar_t;
  typedef int   lno_t;
  typedef int   gno_t;
  typedef size_t offset_t;
  typedef Zoltan2::default_part_t  part_t;
  typedef Zoltan2::default_node_t node_t;
  static inline std::string name() {return "Epetra_Vector";}
};
#endif

template <typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
{
  typedef Scalar        scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef size_t offset_t;
  typedef Zoltan2::default_part_t  part_t;
  typedef Node          node_t;
  static inline std::string name() {return "Xpetra::MultiVector";}

  Z2_STATIC_ASSERT_TYPES // validate the types
};

template <typename Scalar,
          typename LocalOrdinal,
          typename GlobalOrdinal,
          typename Node>
struct InputTraits<Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
{
  typedef Scalar        scalar_t;
  typedef LocalOrdinal  lno_t;
  typedef GlobalOrdinal gno_t;
  typedef size_t offset_t;
  typedef Zoltan2::default_part_t  part_t;
  typedef Node          node_t;
  static inline std::string name() {return "Tpetra::MultiVector";}

  Z2_STATIC_ASSERT_TYPES // validate the types
};

#ifdef HAVE_ZOLTAN2_EPETRA
template < >
struct InputTraits<Epetra_MultiVector>
{
  typedef double scalar_t;
  typedef int   lno_t;
  typedef int   gno_t;
  typedef size_t offset_t;
  typedef Zoltan2::default_part_t  part_t;
  typedef Zoltan2::default_node_t node_t;
  static inline std::string name() {return "Epetra_MultiVector";}
};
#endif

#endif // DOXYGEN_SHOULD_SKIP_THIS


}  // namespace Zoltan2
#endif // ZOLTAN2_INPUTTRAITS_HPP
