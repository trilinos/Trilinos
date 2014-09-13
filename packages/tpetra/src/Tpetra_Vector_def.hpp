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

#ifndef TPETRA_VECTOR_DEF_HPP
#define TPETRA_VECTOR_DEF_HPP

#include <Kokkos_DefaultArithmetic.hpp>
#include <Kokkos_NodeTrace.hpp>
#include <Tpetra_MultiVector.hpp>

#ifdef DOXYGEN_USE_ONLY
  #include "Tpetra_Vector_decl.hpp"
#endif

namespace Tpetra {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  Vector (const Teuchos::RCP<const map_type>& map, const bool zeroOut)
    : base_type (map, 1, zeroOut)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  Vector (const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& source)
    : base_type (source)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  Vector (const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& source,
          const Teuchos::DataAccess copyOrView)
    : base_type (source, copyOrView)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  Vector (const Teuchos::RCP<const map_type>& map,
          const Teuchos::ArrayRCP<Scalar>& view,
          EPrivateHostViewConstructor /* dummy */)
    : base_type (map, view, view.size (), 1, HOST_VIEW_CONSTRUCTOR)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  Vector (const Teuchos::RCP<const map_type>& map,
          const Teuchos::ArrayRCP<Scalar>& view,
          EPrivateComputeViewConstructor /* dummy */)
    : base_type (map, view, view.size (), 1, COMPUTE_VIEW_CONSTRUCTOR)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  Vector (const Teuchos::RCP<const map_type>& map,
          const Teuchos::ArrayView<const Scalar>& values)
    : base_type (map, values, values.size (), 1)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  ~Vector ()
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceGlobalValue(GlobalOrdinal globalRow, const Scalar &value) {
    this->base_type::replaceGlobalValue(globalRow,0,value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar &value) {
    this->base_type::sumIntoGlobalValue(globalRow,0,value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceLocalValue(LocalOrdinal myRow, const Scalar &value) {
    this->base_type::replaceLocalValue(myRow,0,value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoLocalValue(LocalOrdinal myRow, const Scalar &value) {
    this->base_type::sumIntoLocalValue(myRow,0,value);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::get1dCopy(ArrayView<Scalar> A) const {
    size_t lda = this->getLocalLength();
    this->get1dCopy(A,lda);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Scalar Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dot(const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &a) const {
    using Teuchos::outArg;
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION( !this->getMap()->isCompatible(*a.getMap()), std::runtime_error,
        "Tpetra::Vector::dots(): Vectors do not have compatible Maps:" << std::endl
        << "this->getMap(): " << std::endl << *this->getMap()
        << "a.getMap(): " << std::endl << *a.getMap() << std::endl);
#else
    TEUCHOS_TEST_FOR_EXCEPTION( this->getLocalLength() != a.getLocalLength(), std::runtime_error,
        "Tpetra::Vector::dots(): Vectors do not have the same local length.");
#endif
    Scalar gbldot;
    gbldot = MVT::Dot(this->lclMV_,a.lclMV_);
    if (this->isDistributed()) {
      Scalar lcldot = gbldot;
      Teuchos::reduceAll(*this->getMap()->getComm(),Teuchos::REDUCE_SUM,lcldot,outArg(gbldot));
    }
    return gbldot;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Scalar
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::meanValue () const
  {
    Scalar mean;
    this->meanValue (Teuchos::arrayView (&mean, 1));
    return mean;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  norm1 () const
  {
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType mag_type;
    mag_type norm;
    this->norm1 (Teuchos::arrayView (&norm, 1));
    return norm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm2 () const
  {
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType mag_type;
    mag_type norm;
    this->norm2 (Teuchos::arrayView (&norm, 1));
    return norm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normInf () const
  {
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType mag_type;
    mag_type norm;
    this->normInf (Teuchos::arrayView (&norm, 1));
    return norm;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  normWeighted (const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& weights) const
  {
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType mag_type;
    mag_type norm;
    this->normWeighted (weights, Teuchos::arrayView (&norm, 1));
    return norm;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  offsetView (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& subMap,
              size_t offset) const
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> V;

    const size_t newNumRows = subMap->getNodeNumElements ();
    const bool tooManyElts = newNumRows + offset > this->lclMV_.getOrigNumRows ();
    if (tooManyElts) {
      const int myRank = this->getMap ()->getComm ()->getRank ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        newNumRows + offset > MVT::getNumRows (this->lclMV_),
        std::runtime_error,
        "Tpetra::Vector::offsetView: Invalid input Map.  Input Map owns "
        << subMap->getNodeNumElements () << " elements on process " << myRank
        << ".  offset = " << offset << ".  Yet, the Vector contains only "
        << this->lclMV_.getOrigNumRows () << " on this process.");
    }

    KokkosClassic::MultiVector<Scalar, Node> newLocalMV =
      this->lclMV_.offsetView (newNumRows, this->lclMV_.getNumCols (), offset, 0);
    return rcp (new V (subMap, newLocalMV.getValuesNonConst (), COMPUTE_VIEW_CONSTRUCTOR));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  offsetViewNonConst (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& subMap,
                      size_t offset)
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> V;

    const size_t newNumRows = subMap->getNodeNumElements ();
    const bool tooManyElts = newNumRows + offset > this->lclMV_.getOrigNumRows ();
    if (tooManyElts) {
      const int myRank = this->getMap ()->getComm ()->getRank ();
      TEUCHOS_TEST_FOR_EXCEPTION(
        newNumRows + offset > MVT::getNumRows (this->lclMV_),
        std::runtime_error,
        "Tpetra::Vector::offsetViewNonConst: Invalid input Map.  Input Map owns "
        << subMap->getNodeNumElements () << " elements on process " << myRank
        << ".  offset = " << offset << ".  Yet, the MultiVector contains only "
        << this->lclMV_.getOrigNumRows () << " on this process.");
    }

    KokkosClassic::MultiVector<Scalar, Node> newLocalMV =
      this->lclMV_.offsetViewNonConst (newNumRows, this->lclMV_.getNumCols (), offset, 0);
    return rcp (new V (subMap, newLocalMV.getValuesNonConst (), COMPUTE_VIEW_CONSTRUCTOR));
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::description() const
  {
    using Teuchos::TypeNameTraits;

    std::ostringstream oss;
    oss << "\"Tpetra::Vector\": {"
        << "Template parameters: {"
        << "Scalar: " << TypeNameTraits<Scalar>::name ()
        << ", LocalOrdinal: " << TypeNameTraits<LocalOrdinal>::name ()
        << ", GlobalOrdinal: " << TypeNameTraits<GlobalOrdinal>::name ()
        << ", Node: " << TypeNameTraits<Node>::name ()
        << "}";
    if (this->getObjectLabel () != "") {
      oss << ", Label: \"" << this->getObjectLabel () << "\", ";
    }
    oss << "Global length: " << this->getGlobalLength ()
        << "}";
    return oss.str ();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    using std::endl;
    using std::setw;
    using Teuchos::Comm;
    using Teuchos::RCP;
    using Teuchos::TypeNameTraits;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;

    const Teuchos::EVerbosityLevel vl =
      (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

    const Map<LocalOrdinal, GlobalOrdinal, Node>& map = * (this->getMap ());
    RCP<const Comm<int> > comm = map.getComm ();
    const int myImageID = comm->getRank ();
    const int numImages = comm->getSize ();
    Teuchos::OSTab tab0 (out);

    if (vl != VERB_NONE) {
      if (myImageID == 0) {
        out << "\"Tpetra::Vector\":" << endl;
      }
      Teuchos::OSTab tab1 (out);// applies to all processes
      if (myImageID == 0) {
        out << "Template parameters:";
        {
          Teuchos::OSTab tab2 (out);
          out << "Scalar: " << TypeNameTraits<Scalar>::name () << endl
              << "LocalOrdinal: " << TypeNameTraits<LocalOrdinal>::name () << endl
              << "GlobalOrdinal: " << TypeNameTraits<GlobalOrdinal>::name () << endl
              << "Node: " << TypeNameTraits<Node>::name () << endl;
        }
        out << endl;
        if (this->getObjectLabel () != "") {
          out << "Label: \"" << this->getObjectLabel () << "\"" << endl;
        }
        out << "Global length: " << this->getGlobalLength () << endl;
      }
      for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
        if (myImageID == imageCtr) {
          if (vl != VERB_LOW) {
            // VERB_MEDIUM and higher prints getLocalLength()
            out << "Process: " << myImageID << endl;
            Teuchos::OSTab tab2 (out);
            out << "Local length: " << this->getLocalLength () << endl;

            if (vl == VERB_EXTREME && this->getLocalLength () > 0) {
              // VERB_EXTREME prints values
              out << "Global indices and values:" << endl;
              Teuchos::OSTab tab3 (out);
              RCP<Node> node = this->lclMV_.getNode ();
              ArrayRCP<const Scalar> myview =
                node->template viewBuffer<Scalar> (this->getLocalLength (),
                                                   MVT::getValues (this->lclMV_));
              for (size_t i = 0; i < this->getLocalLength (); ++i) {
                out << map.getGlobalElement (i) << ": " << myview[i] << endl;
              }
            }
          }
          std::flush (out); // give output time to complete
        }
        comm->barrier (); // give output time to complete
        comm->barrier ();
        comm->barrier ();
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node >
  createCopy (const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node >& src)
  {
    if (src.getCopyOrView () == Teuchos::Copy) { // copy semantics
      return src; // Copy constructor will make a deep copy.
    } else { // view semantics
      // Create a deep copy using the two-argument copy constructor,
      // set the result 'dst' to have view semantics (which 'src' also
      // has), and return it.
      Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> dst (src, Teuchos::Copy);
      dst.setCopyOrView (Teuchos::View);
      return dst;
    }
  }
} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#if defined(TPETRA_HAVE_KOKKOS_REFACTOR)
#include "Tpetra_KokkosRefactor_Vector_def.hpp"
#endif


#define TPETRA_VECTOR_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class Vector< SCALAR , LO , GO , NODE >; \
  template Vector< SCALAR , LO , GO , NODE > createCopy( const Vector< SCALAR , LO , GO , NODE >& src); \


#endif // TPETRA_VECTOR_DEF_HPP
