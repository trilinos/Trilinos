// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_TPETRAVECTOR_DEF_HPP
#define XPETRA_TPETRAVECTOR_DEF_HPP
#include "Xpetra_TpetraVector_decl.hpp"

namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraVector(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
                 bool zeroOut)
  : TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map, 1, zeroOut) {
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraVector(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
                 const Teuchos::ArrayView<const Scalar>& A)
  : TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map, A, map->getLocalNumElements(), 1) {
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ~TpetraVector() {
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceGlobalValue(GlobalOrdinal globalRow, const Scalar& value) {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceGlobalValue");
  getTpetra_Vector()->replaceGlobalValue(globalRow, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar& value) {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoGlobalValue");
  getTpetra_Vector()->sumIntoGlobalValue(globalRow, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceLocalValue(LocalOrdinal myRow, const Scalar& value) {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceLocalValue");
  getTpetra_Vector()->replaceLocalValue(myRow, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    sumIntoLocalValue(LocalOrdinal myRow, const Scalar& value) {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoLocalValue");
  getTpetra_Vector()->sumIntoLocalValue(myRow, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    norm1() const {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm1");
  return getTpetra_Vector()->norm1();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    norm2() const {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm2");
  return getTpetra_Vector()->norm2();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    normInf() const {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normInf");
  return getTpetra_Vector()->normInf();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Scalar
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    meanValue() const {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::meanValue");
  return getTpetra_Vector()->meanValue();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    description() const {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::description");
  return getTpetra_Vector()->description();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::describe");
  getTpetra_Vector()->describe(out, verbLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Scalar
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    dot(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& a) const {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dot");
  return getTpetra_Vector()->dot(*toTpetra(a));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraVector(const Teuchos::RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& vec)
  : TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(vec) {
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getTpetra_Vector() const {
  return this->TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getTpetra_MultiVector()->getVectorNonConst(0);
}

#ifdef HAVE_XPETRA_EPETRA

#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))

// specialization of TpetraVector for GO=int and NO=SerialNode
template <class Scalar>
class TpetraVector<Scalar, int, int, EpetraNode>
  : public virtual Vector<Scalar, int, int, EpetraNode>, public TpetraMultiVector<Scalar, int, int, EpetraNode> {
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef EpetraNode Node;

#undef XPETRA_TPETRAMULTIVECTOR_SHORT
#undef XPETRA_TPETRAVECTOR_SHORT
#include "Xpetra_UseShortNames.hpp"
#define XPETRA_TPETRAMULTIVECTOR_SHORT
#define XPETRA_TPETRAVECTOR_SHORT

 public:
  using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dot;                 // overloading, not hiding
  using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::norm1;               // overloading, not hiding
  using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::norm2;               // overloading, not hiding
  using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::normInf;             // overloading, not hiding
  using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::meanValue;           // overloading, not hiding
  using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceGlobalValue;  // overloading, not hiding
  using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sumIntoGlobalValue;  // overloading, not hiding
  using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceLocalValue;   // overloading, not hiding
  using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sumIntoLocalValue;   // overloading, not hiding

  //! @name Constructor/Destructor Methods
  //@{

  //! Sets all vector entries to zero.
  TpetraVector(const Teuchos::RCP<const Map>& map, bool zeroOut = true)
    : TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map, 1, zeroOut) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                "int",
                                typeid(EpetraNode).name());
  }

  //! Set multi-vector values from an array using Teuchos memory management classes. (copy)
  TpetraVector(const Teuchos::RCP<const Map>& map, const Teuchos::ArrayView<const Scalar>& A)
    : TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map, A, map->getLocalNumElements(), 1) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                "int",
                                typeid(EpetraNode).name());
  }

  virtual ~TpetraVector() {}

  void replaceGlobalValue(GlobalOrdinal globalRow, const Scalar& value) {}

  //! Adds specified value to existing value at the specified location.
  void sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar& value) {}

  //! Replace current value at the specified location with specified values.
  void replaceLocalValue(LocalOrdinal myRow, const Scalar& value) {}

  //! Adds specified value to existing value at the specified location.
  void sumIntoLocalValue(LocalOrdinal myRow, const Scalar& value) {}

  //@}

  //! @name Mathematical methods
  //@{

  //! Return 1-norm of this Vector.
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm1() const {
    return Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero());
  }

  //! Compute 2-norm of this Vector.
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm2() const {
    return Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero());
  }

  //! Compute Inf-norm of this Vector.
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType normInf() const {
    return Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero());
  }

  //! Compute mean (average) value of this Vector.
  Scalar meanValue() const { return Teuchos::ScalarTraits<Scalar>::zero(); }

  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const { return std::string(""); }

  //! Print the object with some verbosity level to an FancyOStream object.
  void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const {}

  //@}

  //! Computes dot product of this Vector against input Vector x.
  Scalar dot(const Vector& a) const { return Teuchos::ScalarTraits<Scalar>::zero(); }

  //! @name Xpetra specific
  //@{

  //! TpetraMultiVector constructor to wrap a Tpetra::MultiVector object
  TpetraVector(const Teuchos::RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& vec)
    : TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(vec) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                "int",
                                typeid(EpetraNode).name());
  }

  //! Get the underlying Tpetra multivector
  RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  getTpetra_Vector() const {
    return Teuchos::null;
  }

  //@}

};  // TpetraVector class (specialization on GO=int, NO=EpetraNode)

#endif  // #if((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT)))
        //    || (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))

#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))))

// specialization of TpetraVector for GO=int and NO=SerialNode
template <class Scalar>
class TpetraVector<Scalar, int, long long, EpetraNode>
  : public virtual Vector<Scalar, int, long long, EpetraNode>, public TpetraMultiVector<Scalar, int, long long, EpetraNode> {
  typedef int LocalOrdinal;
  typedef long long GlobalOrdinal;
  typedef EpetraNode Node;

#undef XPETRA_TPETRAMULTIVECTOR_SHORT
#undef XPETRA_TPETRAVECTOR_SHORT
#include "Xpetra_UseShortNames.hpp"
#define XPETRA_TPETRAMULTIVECTOR_SHORT
#define XPETRA_TPETRAVECTOR_SHORT

 public:
  using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dot;                 // overloading, not hiding
  using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::norm1;               // overloading, not hiding
  using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::norm2;               // overloading, not hiding
  using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::normInf;             // overloading, not hiding
  using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::meanValue;           // overloading, not hiding
  using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceGlobalValue;  // overloading, not hiding
  using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sumIntoGlobalValue;  // overloading, not hiding
  using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::replaceLocalValue;   // overloading, not hiding
  using TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::sumIntoLocalValue;   // overloading, not hiding

  //! @name Constructor/Destructor Methods
  //@{

  //! Sets all vector entries to zero.
  TpetraVector(const Teuchos::RCP<const Map>& map, bool zeroOut = true)
    : TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map, 1, zeroOut) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                "long long",
                                typeid(EpetraNode).name());
  }

  //! Set multi-vector values from an array using Teuchos memory management classes. (copy)
  TpetraVector(const Teuchos::RCP<const Map>& map, const Teuchos::ArrayView<const Scalar>& A)
    : TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map, A, map->getLocalNumElements(), 1) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                "long long",
                                typeid(EpetraNode).name());
  }

  //! Destructor.
  virtual ~TpetraVector() {}

  //@}

  //! @name Post-construction modification routines
  //@{

  //! Replace current value at the specified location with specified value.
  void replaceGlobalValue(GlobalOrdinal globalRow, const Scalar& value) {}

  //! Adds specified value to existing value at the specified location.
  void sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar& value) {}

  //! Replace current value at the specified location with specified values.
  void replaceLocalValue(LocalOrdinal myRow, const Scalar& value) {}

  //! Adds specified value to existing value at the specified location.
  void sumIntoLocalValue(LocalOrdinal myRow, const Scalar& value) {}

  //@}

  //! @name Mathematical methods
  //@{

  //! Return 1-norm of this Vector.
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm1() const {
    return Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero());
  }

  //! Compute 2-norm of this Vector.
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm2() const {
    return Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero());
  }

  //! Compute Inf-norm of this Vector.
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType normInf() const {
    return Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero());
  }

  //! Compute mean (average) value of this Vector.
  Scalar meanValue() const { return Teuchos::ScalarTraits<Scalar>::zero(); }

  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const { return std::string(""); }

  //! Print the object with some verbosity level to an FancyOStream object.
  void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const {}

  //@}

  //! Computes dot product of this Vector against input Vector x.
  Scalar dot(const Vector& a) const { return Teuchos::ScalarTraits<Scalar>::zero(); }

  //! Compute Weighted 2-norm (RMS Norm) of this Vector.

  //! @name Xpetra specific
  //@{

  //! TpetraMultiVector constructor to wrap a Tpetra::MultiVector object
  TpetraVector(const Teuchos::RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& vec) {
    XPETRA_TPETRA_ETI_EXCEPTION(typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                typeid(TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, EpetraNode>).name(),
                                "long long",
                                typeid(EpetraNode).name());
  }

  //! Get the underlying Tpetra multivector
  RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  getTpetra_Vector() const {
    return Teuchos::null;
  }

  //@}

};  // TpetraVector class (specialization on GO=long long, NO=EpetraNode)

#endif  // #if((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG)))
        //    || (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))))

#endif  // HAVE_XPETRA_EPETRA
}  // namespace Xpetra

#endif
