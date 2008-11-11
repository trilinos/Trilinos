// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
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

  };
}

#endif
