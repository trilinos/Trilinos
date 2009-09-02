// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#ifndef SACADO_PARAMETERACCESSOR_HPP
#define SACADO_PARAMETERACCESSOR_HPP

#include "Sacado_ScalarParameterEntry.hpp"

  /*!
   * \brief Abstract class that provides access to a parameter
   * value in a code for the parameter library. An object of this
   * type is required to construct a ParameterRegistration object.
   */

namespace Sacado {
  template<typename EvalType, typename EvalTypeTraits = DefaultEvalTypeTraits> class ParameterAccessor {
  private:
    typedef typename EvalTypeTraits::template apply<EvalType>::type ScalarT;

  public:

    virtual ~ParameterAccessor() {};

    //! Method that returns a reference to the parameter value given the name
    //! The ParameterLibrary call this method when a parameter value cahnges
    virtual ScalarT& getValue(const std::string &n) = 0;
  };
}

#endif
