// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_SCALARPARAMETERENTRY_HPP
#define SACADO_SCALARPARAMETERENTRY_HPP

#include "Sacado_Traits.hpp"
#include "Sacado_mpl_apply.hpp"

namespace Sacado {

  //! Base traits definition mapping evaluation types to value types
  /*!
   * Users should provide evaluation type tags for various evalutions (e.g.,
   * Residual, Jacobian, Tangent, ...) and specialize SPE_ValueType on
   * those tags.  Each specialization should have a public typedef called type
   * specifying the ValueType for that evaluation type, e.g.,
   * \code
   *
   * struct ResidualType {};
   * template <> struct SPE_ValueType<ResidualType> {
   *  typedef double type;
   * };
   *
   * struct JacobianType {};
   * template <> struct SPE_ValueType<JacobianType> {
   *  typedef Sacado::Fad::DFad<double> type;
   * };
   *
   * \endcode
   */
  struct DefaultEvalTypeTraits {
    template <class EvalType> struct apply {
      typedef EvalType type; };
  };

  //! Abstract interface for all entries in Sacado::ScalarParameterFamily
  class AbstractScalarParameterEntry {
  public:

    //! Default contructor
    AbstractScalarParameterEntry() {}

    //! Destructor
    virtual ~AbstractScalarParameterEntry() {}

    //! Set real parameter value
    virtual void setRealValue(double value) = 0;

    //! Get real parameter value
    virtual double getRealValue() const = 0;

    //! Print entry
    virtual void print(std::ostream& os) const = 0;
  };

  /*!
   * \brief A base class for scalar parameter values
   */
  template <typename EvalType, typename EvalTypeTraits = DefaultEvalTypeTraits>
  class ScalarParameterEntry : public AbstractScalarParameterEntry {

  public:

    typedef typename EvalTypeTraits::template apply<EvalType>::type ScalarT;

    //! Default constructor
    ScalarParameterEntry() {}

    //! Destructor
    virtual ~ScalarParameterEntry() {}

    //! Set parameter this object represents to \em value
    /*!
     * Treat the set parameter as an independent for derivative computations
     * (use setRealValue() otherwise).
     */
    virtual void setValue(const ScalarT& value) = 0;

    //! Get parameter value this object represents
    virtual const ScalarT& getValue() const = 0;

    //! Get real parameter value
    /*!
     * Default implementation should work in most cases.
     */
    virtual double getRealValue() const {
      return Sacado::ScalarValue<ScalarT>::eval(this->getValue());
    }

    //! Print entry
    /*!
     * Default implementation should work in most cases.
     */
    virtual void print(std::ostream& os) const {
      os << getValue();
    }

  };
}

#endif
