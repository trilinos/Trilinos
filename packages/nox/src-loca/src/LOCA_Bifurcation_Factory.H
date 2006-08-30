// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef LOCA_BIFURCATION_FACTORY_H
#define LOCA_BIFURCATION_FACTORY_H

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
  namespace MultiContinuation {
    class AbstractGroup;
  }
}

namespace LOCA {

  namespace Bifurcation {

    //! Factory for creating bifurcation strategy objects
    /*!
     * The parameters passed to the create() through the \c bifurcationParams
     * argument method should specify the "Type" and "Formulation" as 
     * described below, as well as any additional parameters for the 
     * particular strategy.
     * <ul>
     * <li> "Type" - Name of the bifurcation type. Valid choices are
     *      <ul>
     *      <li> "None" - No bifurcation [Default]
     *      <li> "Turning Point" Turning point bifurations
     *      <li> "Pitchfork" Pitchfork bifurcations
     *      <li> "Hopf" Hopf bifurcations
     *      <li> "User-Defined" User defined bifurcation.  Set 
     *           "User-Defined Name" to be the parameter list name containing
     *           user-defined strategy, which must be of type
     *           Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>.
     *      </ul>
     * <li> "Formulation" - Name of the bifurcation formulation. Valid 
     *      choices are
     *      <ul> 
     *      <li> For turning point bifurcations:
     *          <ul>
     *          <li> "Moore-Spence" [Default]
     *               (LOCA::TurningPoint::MooreSpence::ExtendedGroup)
     *          <li> "Minimally Augmented" 
     *               (LOCA::TurningPoint::MinimallyAugmented::ExtendedGroup)
     *          </ul>
     *      <li> For pitchfork bifurcations:
     *           <ul>
     *           <li> "Moore-Spence" [Default]
     *                (LOCA::Pitchfork::MooreSpence::ExtendedGroup)
     *           </ul>
     *      <li> For Hopf bifurcations:
     *           <ul>
     *           <li> "Moore-Spence" [Default]
     *                (LOCA::Hopf::MooreSpence::ExtendedGroup)
     *           </ul>
     *      </ul>
     * </ul>
     */
    class Factory {

    public:

      //! Constructor
      Factory(const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data);

      //! Destructor
      virtual ~Factory();

      //! Create bifurcation strategy
      /*!
       * \param topParams [in] Parsed top-level parameter list.  
       * \param bifurcationParams [in] Bifurcation parameters as described above
       * \param grp [in] Underlying group
       */
      Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>
      create(
      const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RefCountPtr<Teuchos::ParameterList>& bifurcationParams,
      const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& grp);

      //! Return strategy name given by \c bifurcationParams
      string 
      strategyName(Teuchos::ParameterList& bifurcationParams) const;

    private:

      //! Private to prohibit copying
      Factory(const Factory&);

      //! Private to prohibit copying
      Factory& operator = (const Factory&);

    protected:

      //! Global data
      Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

    }; // Class Factory

  } // Namespace Bifurcation

} // Namespace LOCA

#endif // LOCA_BIFURCATION_FACTORY_H
