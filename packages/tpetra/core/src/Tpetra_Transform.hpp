// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_TRANSFORM_HPP
#define TPETRA_TRANSFORM_HPP

/// \file Tpetra_Tranform.hpp
/// \brief Declaration of the Tpetra::Transform class

#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>

namespace Tpetra {

///
/** Base Class for all Tpetra Transforms.
 *
 * This is the abstract definition for all Tpetra Transforms.
 * Depending on the type of Transform, several specializations are
 * available: Structural, SameType, InPlace, View.
 */

template <typename T, typename U>
class Transform {
 public:
  /** @name Typedefs for templated classes */
  //@{

  using OriginalType      = Teuchos::RCP<T>;
  using OriginalConstType = Teuchos::RCP<const T>;
  using NewType           = Teuchos::RCP<U>;
  using NewConstType      = Teuchos::RCP<const U>;

  //@}

  ///
  virtual ~Transform() = default;

  /** @name Pure Virtual Methods which must be implemented by subclasses */
  //@{

  /// \brief Analysis of transform operation on original object and
  /// construction of new object.
  ///
  /// @return Returns an RCP to the newly created object of type NewType.
  virtual NewType operator()(const OriginalType& orig) = 0;

  /// \brief Forward transfer of data.
  ///
  /// Forward transfer of data from <tt>orig</tt> object input in the
  /// <tt>operator()</tt> method call to the new object created in this
  /// same call.
  virtual void fwd() = 0;

  /// \brief Reverse transfer of data.
  ///
  /// Reverse transfer of data from new object created in the
  /// <tt>operator()</tt> method call to the <tt>orig</tt> object input
  /// to this same method.
  virtual void rvs() = 0;

  //@}

  /** @name Virtual functions with default implements allowing for optional
   * implementation by the Transform developer */
  //@{

  /// \brief Initial analysis phase of transform.
  ///
  /// Initial analysis phase of transform to confirm the transform
  /// is possible allowing methods <tt>construct()</tt>,
  /// <tt>fwd()</tt> and <tt>rvs()</tt> to be successfully utilized.
  ///
  /// The default implementation calls method <tt>operator()</tt>
  /// and stores the resulting object in an internal attribute
  /// <tt>newObj_</tt>.
  virtual void analyze(const OriginalType& orig);

  /// \brief Construction of new object as a result of the transform.
  ///
  /// The default implementation returns internal attribute <tt>newObj_</tt>.
  virtual NewType construct();

  /// \brief Check if transformed object has been constructed.
  ///
  /// The default implementation returns <tt>true</tt> if
  /// <tt>newObj_</tt> != 0.
  virtual bool isConstructed();

  //@}

 protected:
  /// \brief Default constructor
  ///
  /// Protected to allow only derived classes to use.
  ///
  /// Initializes attributes <tt>origObj_</tt> and <tt>newObj_</tt>
  /// to Teuchos::null.
  Transform()
    : origObj_(Teuchos::null)
    , newObj_(Teuchos::null) {}

  OriginalType origObj_;
  NewType newObj_;

 private:
  Transform(const Transform<T, U>& src)
    : origObj_(src.origObj_)
    , newObj_(src.newObj_) {}

  Transform<T, U>& operator=(const Transform<T, U>& src) {
    // Not currently supported.
    Teuchos::GlobalMPISession::abort();
    return (*this);
  }

};  // end class Transform

template <typename T, typename U>
void Transform<T, U>::
    analyze(const OriginalType& orig) {
  origObj_ = orig;
  newObj_  = Teuchos::rcp_dynamic_cast<U>(origObj_);
  return;
}

template <typename T, typename U>
typename Transform<T, U>::NewType
Transform<T, U>::
    construct() {
  return newObj_;
}

template <typename T, typename U>
bool Transform<T, U>::
    isConstructed() {
  return (newObj_ != Teuchos::null);
}

template <typename T, typename U>
class StructuralTransform : public Transform<T, U> {
 public:
  void fwd() { return; }
  void rvs() { return; }

  virtual ~StructuralTransform() = default;
};

template <typename T>
class SameTypeTransform : public Transform<T, T> {
 public:
  using TransformType = Teuchos::RCP<T>;

  virtual ~SameTypeTransform() = default;
};

template <typename T>
class StructuralSameTypeTransform : public SameTypeTransform<T> {
 public:
  void fwd() { return; }
  void rvs() { return; }

  virtual ~StructuralSameTypeTransform() = default;
};

//  template<typename T>
//  class InPlaceTransform : public SameTypeTransform<T>
//  {
//   public:
//    typename Transform<T,T>::NewTypeRef
//    operator()
//    ( typename Transform<T,T>::OriginalTypeRef orig )
//    { this->origObj_ = &orig;
//      this->newObj_ = &orig;
//      return orig;
//    }
//
//    virtual ~InPlaceTransform() = default;
//  };

template <typename T>
class ViewTransform : public SameTypeTransform<T> {
 public:
  void fwd() { return; }
  void rvs() { return; }

  virtual ~ViewTransform() = default;
};

}  // namespace Tpetra

#endif  // TPETRA_TRANSFORM_HPP
