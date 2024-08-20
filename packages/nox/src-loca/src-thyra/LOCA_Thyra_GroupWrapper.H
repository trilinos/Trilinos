// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_THYRA_GROUPWRAPPER_H
#define LOCA_THYRA_GROUPWRAPPER_H

#include "LOCA_Thyra_Group.H"                   // base class

// Forward declares
namespace Teuchos {
  class ParameterList;
}

namespace LOCA {

  //! %LOCA %Thyra support
  namespace Thyra {

    //! Extension of the LOCA::Thyra::Group to provide access to solution data
    class GroupWrapper : public LOCA::Thyra::Group {

    public:

      //! Constructor
      GroupWrapper(const Teuchos::RCP<LOCA::GlobalData>& global_data,
        const NOX::Thyra::Vector& initial_guess,
        const Teuchos::RCP< ::Thyra::ModelEvaluator<double> >& model,
        const LOCA::ParameterVector& p,
        int p_index,
        bool implement_dfdp = false);

      /*!
       * \brief Copy constructor. If type is DeepCopy, takes ownership of
       * valid shared Jacobian.
       */
      GroupWrapper(const GroupWrapper& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~GroupWrapper();

      //! Assignment operator
      virtual GroupWrapper& operator=(const GroupWrapper& source);

      /*!
       * @name Overloaded NOX::Thyra::Group  methods.
       */
      //@{

      //! Assignment operator.
      virtual NOX::Abstract::Group&
      operator=(const NOX::Abstract::Group& source);

      //! Assignment operator.
      virtual NOX::Abstract::Group&
      operator=(const NOX::Thyra::Group& source);

      //! Assignment operator.
      virtual LOCA::Thyra::Group&
      operator=(const LOCA::Thyra::Group& source);

      //! Cloning function
      virtual Teuchos::RCP<NOX::Abstract::Group>
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Access the solution
      const Teuchos::RCP<const NOX::Thyra::Vector> getNOXThyraVecRCPX(){ return x_vec_; }

    };

  } // namespace Thyra
} // namespace LOCA


#endif
