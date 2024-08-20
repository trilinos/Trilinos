// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_MULTICONTINUATION_NATURALGROUP_H
#define LOCA_MULTICONTINUATION_NATURALGROUP_H

#include "LOCA_MultiContinuation_ExtendedGroup.H"       // base class

namespace LOCA {

  namespace MultiContinuation {

    /*!
     * \brief Specialization of LOCA::MultiContinuation::ExtendedGroup
     * to natural continuation.
     */
    /*!
     * Natural continuation corresponds to a continuation equation
     * \f$g(x,p,x_0,p_0,x^\ast,p^\ast,v,\Delta s)=0\f$ with \f$g\f$ given by
     * \f[
     *     g(x,p,x_0,p_0,x^\ast,p^\ast,v,\Delta s)=
     *         p-p_0-v_p \Delta s
     * \f]
     * where \f$v_p\f$ is the parameter component of the predictor direction
     * \f$v\f$. This corresponds geometrically to constraining the nonlinear
     * solver steps used in calculating \f$F(x,p)=0\f$ to be orthogonal to the
     * parameter axis.  The natural constraint \f$g\f$ is
     * represented by a LOCA::MultiContinuation::NaturalConstraint object.
     */
    class NaturalGroup :
      public virtual LOCA::MultiContinuation::ExtendedGroup {

    public:

      //! Constructor
      /*!
       * \param global_data [in] Global data object
       * \param topParams [in] Parsed top-level parameter list.
       * \param continuationParams [in] Continuation parameters.
       * \param grp [in] Group representing \f$F\f$.
       * \param pred [in] Predictor strategy.
       * \param paramIDs [in] Parameter IDs of continuation parameters.
       */
      NaturalGroup(
      const Teuchos::RCP<LOCA::GlobalData>& global_data,
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RCP<Teuchos::ParameterList>& continuationParams,
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const std::vector<int>& paramIDs);

      //! Copy constructor
      NaturalGroup(const NaturalGroup& source,
             NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~NaturalGroup();

      /*!
       * @name Implementation of NOX::Abstract::Group virtual methods
       */
      //@{

      //! Assignment operator
      virtual NOX::Abstract::Group&
      operator=(const NOX::Abstract::Group& source);

      //! Clone function
      virtual Teuchos::RCP<NOX::Abstract::Group>
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //@}

      /*!
       * @name Implementation of LOCA::MultiContinuation::AbstractStrategy virtual methods
       */
      //@{

      //! Copy
      virtual void copy(const NOX::Abstract::Group& source);

      //@}

    private:

      //! Prohibit generation and use of operator=()
      NaturalGroup& operator=(const NaturalGroup& source);

    };

  } // namespace MultiContinuation
} // namespace LOCA

#endif
