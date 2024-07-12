// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_THYRA_GROUP_H
#define LOCA_THYRA_GROUP_H

#include "LOCA_Abstract_Group.H"               // base class
#include "NOX_Thyra_Group.H"                   // base class
#include "LOCA_Parameter_Vector.H"             // class data element
#include "Thyra_VectorBase.hpp"                // class data element
#include "LOCA_Thyra_SaveDataStrategy.H"       // class data element

// Forward declares
namespace Teuchos {
  class ParameterList;
}

namespace LOCA {

  //! %LOCA %Thyra support
  namespace Thyra {

    //! Extension of the NOX::Thyra::Group to %LOCA.
    class Group : public virtual NOX::Thyra::Group,
          public virtual LOCA::Abstract::Group {

    public:

      //! Constructor (internally constructs nox group)
      Group(const Teuchos::RCP<LOCA::GlobalData>& global_data,
            const NOX::Thyra::Vector& initial_guess,
            const Teuchos::RCP< ::Thyra::ModelEvaluator<double> >& model,
            const LOCA::ParameterVector& p,
            int p_index,
            bool implement_dfdp = false,
            const Teuchos::RCP<const ::Thyra::VectorBase<double> >& weight_vector = Teuchos::null,
            const bool set_transient_in_args = true);

      //! Constructor for Householder continuation where parameters are in separate param vectors.
      Group(const Teuchos::RCP<LOCA::GlobalData>& global_data,
            const NOX::Thyra::Group& nox_group,
            const LOCA::ParameterVector& p,
            const std::vector<int>& p_index,
            bool implement_dfdp = false,
            const Teuchos::RCP<const ::Thyra::VectorBase<double> >& weight_vector = Teuchos::null,
            const bool set_transient_in_args = true);

      /*!
       * \brief Copy constructor. If type is DeepCopy, takes ownership of
       * valid shared Jacobian.
       */
      Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~Group();

      //! Assignment operator
      virtual Group& operator=(const Group& source);

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

      //! Cloning function
      virtual Teuchos::RCP<NOX::Abstract::Group>
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Overloaded computeF()
      virtual NOX::Abstract::Group::ReturnType
      computeF();

      //! Overloaded computeJacobian()
      virtual NOX::Abstract::Group::ReturnType
      computeJacobian();

      //@}

      /*!
       * @name Implementation of LOCA::MultiContinuation::AbstractGroup virtual methods.
       */
      //@{

      //! Copy
      virtual void copy(const NOX::Abstract::Group& source);

      //! Set the parameters
      virtual void setParams(const ParameterVector& p);

      //! Set parameter indexed by paramID
      virtual void setParam(int paramID, double val);

      //! Set parameter indexed by paramID
      virtual void setParam(std::string paramID, double val);

      //! Return a const reference to the ParameterVector owned by the group.
      const LOCA::ParameterVector& getParams() const;

      //! Return copy of parameter indexed by paramID
      virtual double getParam(int paramID) const;

      //! Return copy of parameter indexed by paramID
      virtual double getParam(std::string paramID) const;

      /*!
       * Compute \f$\partial F/\partial p\f$ for each parameter \f$p\f$
       * indexed by paramIDs.  The first column of \em dfdp holds F,
       * which is valid if \em isValidF is true.  Otherwise F must be
       * computed.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeDfDpMulti(const std::vector<int>& paramIDs,
               NOX::Abstract::MultiVector& dfdp,
               bool isValidF);

      //! Perform any preprocessing before a continuation step starts.
      virtual void
      preProcessContinuationStep(
                 LOCA::Abstract::Iterator::StepStatus stepStatus);

      //! Perform any postprocessing after a continuation step finishes.
      virtual void
      postProcessContinuationStep(
                 LOCA::Abstract::Iterator::StepStatus stepStatus);

      //! Projects solution to a few scalars for multiparameter continuation
      virtual void projectToDraw(const NOX::Abstract::Vector& x,
                 double *px) const;

      //! Returns the dimension of the project to draw array
      virtual int projectToDrawDimension() const;

      //! Compute a scaled dot product
      /*!
       * The implementation here uses the standard dot product divided
       * by the length of the vectors.
       */
      virtual double
      computeScaledDotProduct(const NOX::Abstract::Vector& a,
                  const NOX::Abstract::Vector& b) const;

      //! Call the user interface print() routine, solution vector
      virtual void printSolution(const double conParam) const;

      //! Call the user interface print() routine, any vector
      virtual void printSolution(const NOX::Abstract::Vector& x,
                                 const double conParam) const;

      //! Scales a vector using scaling vector
      /*!
       * The implementation here rescales the vector by the square root of
       * its length.
       */
      virtual void
      scaleVector(NOX::Abstract::Vector& x) const;

      //@}

      /*!
       * @name Implementation of LOCA::TimeDependent::AbstractGroup virtual methods.
       */
      //@{

      //! Compute the shifted matrix
      virtual NOX::Abstract::Group::ReturnType
      computeShiftedMatrix(double alpha, double beta);

      //! Multiply the shifted matrix by a vector.
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrix(const NOX::Abstract::Vector& input,
             NOX::Abstract::Vector& result) const;

      //! Multiply the shifted matrix by a multi-vector.
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrixMultiVector(
                const NOX::Abstract::MultiVector& input,
                NOX::Abstract::MultiVector& result) const;

      /*!
       * \brief Apply the inverse of the shifted matrix by a multi-vector, as
       * needed by the shift-and-invert and generalized Cayley transformations.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrixInverseMultiVector(
                    Teuchos::ParameterList& params,
                const NOX::Abstract::MultiVector& input,
                NOX::Abstract::MultiVector& result) const;
      //@}

      //! Set strategy object for saving continuation data
      void setSaveDataStrategy(
             const Teuchos::RCP<LOCA::Thyra::SaveDataStrategy>& s);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Parameter vector
      LOCA::ParameterVector params;

      //! Parameter index in Thyra::ModelEvaluator
      std::vector<int> param_index;

      //! Pointer to Thyra vector storing parameters
      std::vector<Teuchos::RCP< const ::Thyra::VectorBase<double> > > param_thyra_vec;

      //! Pointer to x_dot vector of zeros
      Teuchos::RCP< const ::Thyra::VectorBase<double> > x_dot_vec;

      //! Strategy object for saving continuation data
      Teuchos::RCP<LOCA::Thyra::SaveDataStrategy> saveDataStrategy;

      //! Whether we should implement df/dp, or use DerivUtils
      bool implement_dfdp;

      /** \brief Optional weighting vector for function scaling

           This is pulled out of the initial guess vector
      */
      Teuchos::RCP<const ::Thyra::VectorBase<double> > weight_vec_;

      // If set to true, each parameter is in a separate model evaluator parameter vector
      bool paramsInSeparatePVecs;

      // If set to true, when evaluating the Jacobian, if transient
      // in_args are supported, set the values for a steady state
      // Jacobian. LOCA uses transient in_args for stability analysis
      // and otherwise normally assumes steady-state. It default sets
      // transients to make sure the Jacobian is steady-state
      // (alpha=0, beta=1, x_dot=0). This override was added for using
      // LOCA bordering algorithms within a Tempus transient loop
      // where Tempus will set the transient parameters.
      bool set_transient_in_args_;

    private:

      void updateThyraParamView();
      void updateThyraXDot();
    };

  } // namespace Thyra
} // namespace LOCA


#endif
