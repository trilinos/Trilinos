// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
