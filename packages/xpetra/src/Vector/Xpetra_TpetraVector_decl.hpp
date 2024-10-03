// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_TPETRAVECTOR_DECL_HPP
#define XPETRA_TPETRAVECTOR_DECL_HPP

#include "Xpetra_TpetraConfigDefs.hpp"

#include "Xpetra_Vector.hpp"
#include "Xpetra_MultiVector_decl.hpp"
#include "Xpetra_TpetraMultiVector_decl.hpp"

#include "Xpetra_TpetraMap_decl.hpp"
#include "Xpetra_Utils.hpp"

#include "Tpetra_Vector.hpp"

namespace Xpetra {

// TODO: move that elsewhere
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> toTpetra(Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&);

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> toTpetra(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&);

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> toXpetra(RCP<const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> vec);

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> toXpetra(RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> vec);

//
//

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class TpetraVector
  : public virtual Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>,
    public TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
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
  TpetraVector(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& map, bool zeroOut = true);

  //! Set multi-vector values from an array using Teuchos memory management classes. (copy)
  TpetraVector(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& map, const Teuchos::ArrayView<const Scalar>& A);

  //! Destructor.
  virtual ~TpetraVector();

  //@}

  //! @name Post-construction modification routines
  //@{

  //! Replace current value at the specified location with specified value.
  void replaceGlobalValue(GlobalOrdinal globalRow, const Scalar& value);

  //! Adds specified value to existing value at the specified location.
  void sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar& value);

  //! Replace current value at the specified location with specified values.
  void replaceLocalValue(LocalOrdinal myRow, const Scalar& value);

  //! Adds specified value to existing value at the specified location.
  void sumIntoLocalValue(LocalOrdinal myRow, const Scalar& value);

  //@}

  //! @name Mathematical methods
  //@{

  //! Return 1-norm of this Vector.
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm1() const;

  //! Compute 2-norm of this Vector.
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm2() const;

  //! Compute Inf-norm of this Vector.
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType normInf() const;

  //! Compute mean (average) value of this Vector.
  Scalar meanValue() const;

  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  //! Return a simple one-line description of this object.
  std::string description() const;

  //! Print the object with some verbosity level to an FancyOStream object.
  void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const;

  //@}

  //! Computes dot product of this Vector against input Vector x.
  Scalar dot(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& a) const;

  //! @name Xpetra specific
  //@{

  //! TpetraMultiVector constructor to wrap a Tpetra::MultiVector object
  TpetraVector(const Teuchos::RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& vec);

  //! Get the underlying Tpetra multivector
  RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> getTpetra_Vector() const;

  //@}

};  // TpetraVector class

// TODO: move that elsewhere
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
toTpetra(Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
  typedef TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraVectorClass;
  XPETRA_DYNAMIC_CAST(TpetraVectorClass, x, tX, "toTpetra");
  return tX.getTpetra_Vector();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
toTpetra(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) {
  typedef TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraVectorClass;
  XPETRA_DYNAMIC_CAST(const TpetraVectorClass, x, tX, "toTpetra");
  return tX.getTpetra_Vector();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
toXpetra(RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> vec) {
  if (!vec.is_null())
    return rcp(new TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(vec));

  return Teuchos::null;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
toXpetra(RCP<const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> vec) {
  // We cast away the const to wrap the Tpetra vector into an Xpetra object. But it's OK because the Xpetra vector is returned as const.
  return toXpetra(Teuchos::rcp_const_cast<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(vec));
}

}  // namespace Xpetra

#define XPETRA_TPETRAVECTOR_SHORT
#endif  // XPETRA_TPETRAVECTOR_DECL_HPP
