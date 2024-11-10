// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef BELOS_OPERATORT_HPP
#define BELOS_OPERATORT_HPP

#include "BelosOperatorTraits.hpp"

namespace Belos {
  // Base class for the Belos OP template type. Similar to Belos::Operator<> but this one deals with any kind of vector (not only Belos::MultiVec as the Belos::Operator interface)
  template <class MV>
  class OperatorT {

  public:

    //! @name Constructor/Destructor
    //@{

    //! Default constructor
    OperatorT() {};

    //! Destructor.
    virtual ~OperatorT() {};
    //@}

    //! @name Operator application method
    //@{

    /*! \brief This routine takes the Belos::MultiVec \c x and applies the operator
      to it resulting in the Belos::MultiVec \c y, which is returned.
      \note It is expected that any problem with applying this operator to \c x will be
      indicated by an std::exception being thrown.
    */
    virtual void Apply ( const MV & x, MV & y, ETrans trans=NOTRANS ) const = 0;
  };

  /// \brief Specialization of OperatorTraits for OperatorT.
  ///
  /// This is a partial template specialization of
  /// Belos::OperatorTraits class using the Belos::OperatorT
  /// abstract interface. Any class that inherits
  /// from Belos::OperatorT will be accepted by the Belos templated
  /// solvers, due to this specialization of Belos::OperatorTraits.
  template <class ScalarType, class MV>
  class OperatorTraits<ScalarType, MV, OperatorT<MV> >
  {

  public:
    //! Specialization of Apply() for OperatorT.
    static void Apply (const OperatorT<MV>& Op,
                       const MV& x,
                       MV& y, ETrans trans=NOTRANS) {
      Op.Apply (x, y, trans);
    }
  };

} // namespace Belos

#endif // BELOS_OPERATORT_HPP
