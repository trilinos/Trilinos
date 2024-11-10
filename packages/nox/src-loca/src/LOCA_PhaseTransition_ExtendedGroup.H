// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_PHASETRANSITION_EXTENDEDGROUP_H
#define LOCA_PHASETRANSITION_EXTENDEDGROUP_H

#include "LOCA_Abstract_Group.H"  // base class
#include "LOCA_Parameter_Vector.H"
#include "LOCA_GlobalData.H"

#include "NOX_Common.H"             // class data element (std::string)
#include "LOCA_PhaseTransition_ExtendedVector.H"
#include "LOCA_MultiContinuation_AbstractGroup.H"         // base class
#include "LOCA_Extended_MultiAbstractGroup.H"             // base class
#include "Teuchos_ParameterList.hpp"             // base class
//#include "LOCA_Abstract_TransposeSolveGroup.H"            // base class


//extern "C" {void fill_resid_and_matrix_control(double **, int, int);}

// Forward declares
namespace NOX {
  namespace Parameter {
    class List;
  }
}

namespace LOCA {
  namespace PhaseTransition {

    //class Group : public virtual NOX::Abstract::Group {
    class ExtendedGroup :
            public virtual LOCA::Extended::MultiAbstractGroup,
            public virtual LOCA::MultiContinuation::AbstractGroup {

    public:

      //! Constructor
      ExtendedGroup(
              const Teuchos::RCP<LOCA::GlobalData> gD,
              const Teuchos::RCP<Teuchos::ParameterList>& bifurcationparams_,
              const Teuchos::RCP<LOCA::PhaseTransition::AbstractGroup>& grp_);

      //! Copy constructor
      ExtendedGroup(const LOCA::PhaseTransition::ExtendedGroup& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      ~ExtendedGroup();

      //Start methods for NOX::ABstract::Group

      NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);
      //! See above.
      NOX::Abstract::Group& operator=(const LOCA::PhaseTransition::ExtendedGroup& source);

      /** @name "Compute" functions. */
      //@{

      void setX(const NOX::Abstract::Vector& y);
      //! See above
      void setX(const LOCA::PhaseTransition::ExtendedVector& y);

      void computeX(const NOX::Abstract::Group& grp, const NOX::Abstract::Vector& d, double step);
      //! See above.
      void computeX(const LOCA::PhaseTransition::ExtendedGroup& grp,
          const LOCA::PhaseTransition::ExtendedVector& d, double step);

      NOX::Abstract::Group::ReturnType computeF();

      NOX::Abstract::Group::ReturnType computeJacobian();

      NOX::Abstract::Group::ReturnType computeNewton(Teuchos::ParameterList& params);

      //@}

      /** @name Jacobian operations.
       *
       * Operations using the Jacobian matrix. These may not be defined in
       * matrix-free scenarios. */

      //@{

      NOX::Abstract::Group::ReturnType applyJacobian(const LOCA::PhaseTransition::ExtendedVector& input,
                             LOCA::PhaseTransition::ExtendedVector& result) const;

      //! See above
      NOX::Abstract::Group::ReturnType applyJacobian(const NOX::Abstract::Vector& input,
                             NOX::Abstract::Vector& result) const;

      NOX::Abstract::Group::ReturnType applyJacobianInverse(Teuchos::ParameterList& params,
                        const LOCA::PhaseTransition::ExtendedVector& input,
                                  LOCA::PhaseTransition::ExtendedVector& result) const;

      NOX::Abstract::Group::ReturnType applyJacobianInverse(Teuchos::ParameterList& params,
                                const NOX::Abstract::Vector& input,
                                NOX::Abstract::Vector& result) const;

      //@}

      /** @name "Is" functions
       *
       * Checks to see if various objects have been computed. Returns true
       * if the corresponding "compute" function has been called since the
       * last update to the solution vector (via instantiation or
       * computeX). */
      //@{

      bool isF() const;
      bool isJacobian() const;
      bool isNewton() const;

      //@}

      /** @name "Get" functions
       *
       * Note that these function do not check whether or not the vectors
       * are valid. Must use the "Is" functions for that purpose. */
      //@{

      const NOX::Abstract::Vector& getX() const;

      const NOX::Abstract::Vector& getF() const;

      double getNormF() const;

      const NOX::Abstract::Vector& getNewton() const;

      const NOX::Abstract::Vector& getGradient() const;

      Teuchos::RCP< const NOX::Abstract::Vector > getXPtr() const;

      Teuchos::RCP< const NOX::Abstract::Vector > getFPtr() const;

      Teuchos::RCP< const NOX::Abstract::Vector > getNewtonPtr() const;

      Teuchos::RCP< const NOX::Abstract::Vector > getGradientPtr() const;

      //@}

      virtual Teuchos::RCP<NOX::Abstract::Group>
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Print out the group
      void print() const;

      //! Start methods for LOCA::Abstract::Group

      //@{
      void copy(const NOX::Abstract::Group& source) { *this = source; }

      //! Set the parameter vector in the group to p (pVector = p).
      void setParams(const LOCA::ParameterVector& p);

       //! Set parameter indexed by (integer) paramID
      void setParam(int paramID, double val);

      //! Set parameter indexed by (std::string) paramID
      void setParam(std::string paramID, double val);

      //! Return a const reference to the ParameterVector owned by the group.
      const LOCA::ParameterVector& getParams() const;

      //! Return copy of parameter indexed by (integer) paramID
      double getParam(int paramID) const;

      //! Return copy of parameter indexed by (std::string) paramID
      double getParam(std::string paramID) const;

      //! Set parameter indexed by (std::string) paramID
      void printSolution(const NOX::Abstract::Vector& solution, const double param) const;
      void printSolution(const double param) const;

      //@}

      //! Return underlying group
      virtual
      Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup>
      getUnderlyingGroup() const { return grp;};

      //! Return underlying group
      virtual
      Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
      getUnderlyingGroup() { return grp;};

     virtual void
     setParamsMulti( const std::vector<int>& paramIDs,
                     const NOX::Abstract::MultiVector::DenseMatrix& vals);

      virtual NOX::Abstract::Group::ReturnType
      computeDfDpMulti( const std::vector<int>& paramIDs,
                        NOX::Abstract::MultiVector& dfdp,
                                          bool isValid_F);

    protected:

      //! resets the isValid flags to false
      void resetIsValid();

    protected:

      //! Solution vector.
      Teuchos::RCP<LOCA::PhaseTransition::AbstractGroup> grp;

      /** @name Vectors */
      //@{
      //! Solution vector.
      Teuchos::RCP<LOCA::PhaseTransition::ExtendedVector> xVector;
      //! Right-hand-side vector (function evaluation).
      Teuchos::RCP<LOCA::PhaseTransition::ExtendedVector> fVector;
      //! Newton direction vector.
      Teuchos::RCP<LOCA::PhaseTransition::ExtendedVector> newtonVector;
      //@}

      /** @name IsValid flags
       *
       * True if the current solution is up-to-date with respect to the
       * currect xVector. */
      //@{
      bool isValidF;
      bool isValidJacobian;
      bool isValidNewton;
      //@}

      double normF;
      const Teuchos::RCP<LOCA::GlobalData> globalData;
      int  bifParamID;

    };

  } // namespace
} // namespace NOX


#endif
