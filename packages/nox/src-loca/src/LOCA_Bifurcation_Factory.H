// $Id: LOCA_Bifurcation_Factory.H,v 1.11 2007/06/21 16:22:52 rhoope Exp $
// $Source: /space/CVS/Trilinos/packages/nox/src-loca/src/LOCA_Bifurcation_Factory.H,v $

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#include "Teuchos_RCP.hpp"

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
     *      <li> "Phase Transition" Phase Transition points
     *      <li> "User-Defined" User defined bifurcation.  Set 
     *           "User-Defined Name" to be the parameter list name containing
     *           user-defined strategy, which must be of type
     *           Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>.
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
     *           <li> "Minimally Augmented" 
     *               (LOCA::Pitchfork::MinimallyAugmented::ExtendedGroup)
     *           </ul>
     *      <li> For Hopf bifurcations:
     *           <ul>
     *           <li> "Moore-Spence" [Default]
     *                (LOCA::Hopf::MooreSpence::ExtendedGroup)
     *           <li> "Minimally Augmented" 
     *               (LOCA::Hopf::MinimallyAugmented::ExtendedGroup)
     *           </ul>
     *      </ul>
     * </ul>
     */
    class Factory {

    public:

      //! Constructor
      Factory(const Teuchos::RCP<LOCA::GlobalData>& global_data);

      //! Destructor
      virtual ~Factory();

      //! Create bifurcation strategy
      /*!
       * \param topParams [in] Parsed top-level parameter list.  
       * \param bifurcationParams [in] Bifurcation parameters as described above
       * \param grp [in] Underlying group
       */
      Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
      create(
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RCP<Teuchos::ParameterList>& bifurcationParams,
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp);

      //! Return strategy name given by \c bifurcationParams
      std::string 
      strategyName(Teuchos::ParameterList& bifurcationParams) const;

    private:

      //! Private to prohibit copying
      Factory(const Factory&);

      //! Private to prohibit copying
      Factory& operator = (const Factory&);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

    }; // Class Factory

  } // Namespace Bifurcation

} // Namespace LOCA

#endif // LOCA_BIFURCATION_FACTORY_H
