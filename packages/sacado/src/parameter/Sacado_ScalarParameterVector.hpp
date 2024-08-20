// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_SCALARPARAMETERVECTOR_HPP
#define SACADO_SCALARPARAMETERVECTOR_HPP

#include "Sacado_ParameterVectorBase.hpp"
#include "Sacado_ScalarParameterFamily.hpp"

namespace Sacado {

  /*!
   * \brief Specialization of Sacado::ParameterVectorBase for scalar parameters
   */
  template <typename EvalTypeTraits = DefaultEvalTypeTraits>
  class ScalarParameterVector :
    public ParameterVectorBase<ScalarParameterFamily<EvalTypeTraits>, double> {

  public:

    //! Default constructor
    ScalarParameterVector() {}

    //! Copy constructor
    ScalarParameterVector(const ScalarParameterVector& source) :
      ParameterVectorBase<ScalarParameterFamily<EvalTypeTraits>, double>(source) {}

    //! Destructor
    virtual ~ScalarParameterVector() {}

    //! Assignment operator
    ScalarParameterVector& operator = (const ScalarParameterVector& source) {
      ParameterVectorBase<ScalarParameterFamily<EvalTypeTraits>, double>::operator=(source);
      return *this;
    }

  };

}

#endif
