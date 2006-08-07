// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
// (etphipp@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#ifndef LOCA_MULTIPREDICTOR_FACTORY_H
#define LOCA_MULTIPREDICTOR_FACTORY_H

#include "Teuchos_RefCountPtr.hpp"

// forward declarations
namespace Teuchos {
  class ParameterList;
}
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
  namespace MultiPredictor {
    class AbstractStrategy;
  }
}

namespace LOCA {

  namespace MultiPredictor {

    //! Factory for creating %Predictor strategy objects
    /*!
     * The parameters passed to the create() through the \c predictorParams
     * argument method should specify the
     * "Method" as described below, as well as any additional parameters
     * for the particular strategy.
     * <ul>
     * <li> "Method" - Name of the predictor method. Valid choices are
     *   <ul>
     *   <li> "Constant" (LOCA::MultiPredictor::Constant)
     *   <li> "Tangent" (LOCA::MultiPredictor::Tangent) 
     *   <li> "Secant" (LOCA::MultiPredictor::Secant) [Default]
     *   <li> "Random" (LOCA::MultiPredictor::Random)
     *   <li> "Restart" (LOCA::MultiPredictor::Restart)
     *   </ul>
     * </ul>
     */
    class Factory {

    public:

      //! Constructor
      Factory(const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data);

      //! Destructor
      virtual ~Factory();

      //! Create predictor strategy
      /*!
       * \param topParams [in] Parsed top-level parameter list.  
       * \param predictorParams [in] %Predictor parameters as described above
       */
      Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy>
      create(
	const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RefCountPtr<Teuchos::ParameterList>& predictorParams);

      //! Return strategy name given by \c predictorParams
      const string& strategyName(Teuchos::ParameterList& predictorParams) const;

    private:

      //! Private to prohibit copying
      Factory(const Factory&);

      //! Private to prohibit copying
      Factory& operator = (const Factory&);

    protected:

      //! Global data
      Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

    }; // Class Factory

  } // Namespace MultiPredictor

} // Namespace LOCA

#endif // LOCA_MULTIPREDICTOR_FACTORY_H
