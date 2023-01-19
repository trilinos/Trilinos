// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef XPETRA_MULTIVECTORFACTORY_DECL_HPP
#define XPETRA_MULTIVECTORFACTORY_DECL_HPP

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_MultiVector_decl.hpp"

#include "Xpetra_TpetraMultiVector.hpp"

#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraMultiVector.hpp"
#include "Xpetra_EpetraIntMultiVector.hpp"
#endif

// #include "Xpetra_BlockedMap.hpp"
#include "Xpetra_Exceptions.hpp"


namespace Xpetra {

/*!
  @class MultiVectorFactory
  @brief Factory for any type of Xpetra::MultiVector and its derived classes

  Creates instances of \c Xpetra::MulitVector and \c Xpetra::BlockedMultiVector ,
  depending on the type of map, i.e. \c Xpetra::Map vs. \c Xpetra::BlockedMap .

  @tparam Scalar
  @tparam LocalOrdinal
  @tparam GlobalOrdinal
  @tparam Node

  \note Although this class can gerenate \c Xpetra::BlockedMultiVector ,
  it always returns \c Xpetra::MultiVector .
  Don't forget to cast to \c Xpetra::BlockedMultiVector , if you need the blocked layout directly.
*/
template<class Scalar,
         class LocalOrdinal,
         class GlobalOrdinal,
         class Node = KokkosClassic::DefaultNode::DefaultNodeType>
class MultiVectorFactory
{
  private:

    //! Private constructor. This is a static class.
    MultiVectorFactory() {}

  public:

    //! Constructor specifying the number of non-zeros for all rows.
    static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& map, size_t NumVectors, bool zeroOut = true);

    //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy).
    static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
          const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar>>& ArrayOfPtrs,
          size_t                                                            NumVectors);

    static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP<const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node > > &source,
        Teuchos::DataAccess copyOrView);
};



// we need the Epetra specialization only if Epetra is enabled
#if defined(HAVE_XPETRA_EPETRA)

#if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)


// Specialization for Scalar=double, LO=GO=int and Serial node
// Used both for Epetra and Tpetra
// For any other node definition the general default implementation is used which allows Tpetra only
template<>
class MultiVectorFactory<double, int, int, EpetraNode>
{

    typedef double     Scalar;
    typedef int        LocalOrdinal;
    typedef int        GlobalOrdinal;
    typedef EpetraNode Node;

  private:

    //! Private constructor. This is a static class.
    MultiVectorFactory();

  public:

    static RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& map, size_t NumVectors, bool zeroOut = true);

    //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy).
    static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
          const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar>>& ArrayOfPtrs,
          size_t                                                            NumVectors);

    static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP<const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node > > &source,
        Teuchos::DataAccess copyOrView);
};



template<>
class MultiVectorFactory<int, int, int, EpetraNode>
{

    typedef int        Scalar;
    typedef int        LocalOrdinal;
    typedef int        GlobalOrdinal;
    typedef EpetraNode Node;

  private:

    //! Private constructor. This is a static class.
    MultiVectorFactory();

  public:

    static RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& map, size_t NumVectors, bool zeroOut = true);

    //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy).
    static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
    Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
          const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar>>& ArrayOfPtrs,
          size_t                                                            NumVectors);

    static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
    Build(const Teuchos::RCP<const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node > > &source,
          Teuchos::DataAccess copyOrView);
};


// we need the Epetra specialization only if Epetra is enabled
#if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)


template <>
class MultiVectorFactory<double, int, long long, EpetraNode>
{
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef long long GlobalOrdinal;
  typedef EpetraNode Node;

private:

  //! Private constructor. This is a static class.
  MultiVectorFactory();

public:

  static RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map,
        size_t NumVectors,
        bool zeroOut = true);

  //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy).
  static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &map,
        const Teuchos::ArrayView< const Teuchos::ArrayView< const Scalar > > &ArrayOfPtrs,
        size_t NumVectors);

  static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Build(const Teuchos::RCP<const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node > > &source,
        Teuchos::DataAccess copyOrView);

};


template <>
class MultiVectorFactory<int, int, long long, EpetraNode>
{

  typedef int Scalar;
  typedef int LocalOrdinal;
  typedef long long GlobalOrdinal;
  typedef EpetraNode Node;

private:

  //! Private constructor. This is a static class.
  MultiVectorFactory();

public:

  static RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map,
        size_t NumVectors,
        bool zeroOut = true);

  //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy).
  static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Build(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &map,
        const Teuchos::ArrayView< const Teuchos::ArrayView< const Scalar > > &ArrayOfPtrs,
        size_t NumVectors);

  static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Build(const Teuchos::RCP<const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node > > &source,
        Teuchos::DataAccess copyOrView);

};


#endif      // END !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)

#endif      // END !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)

#endif      // END HAVE_XPETRA_EPETRA


}      // namespace Xpetra

#define XPETRA_MULTIVECTORFACTORY_SHORT
#endif
