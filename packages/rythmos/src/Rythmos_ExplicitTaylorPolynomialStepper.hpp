//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef RYTHMOS_EXPLICIT_TAYLOR_POLYNOMIAL_STEPPER_H
#define RYTHMOS_EXPLICIT_TAYLOR_POLYNOMIAL_STEPPER_H

#include "Rythmos_Stepper.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_ModelEvaluatorHelpers.hpp"
#include "RTOpPack_RTOpTHelpers.hpp"

namespace Rythmos {

  /*!
   * \brief Implementation of Rythmos::Stepper for explicit Taylor polynomial
   * time integration of ODEs.
   */
  /*!
   * Let 
   * \f[
   *     \frac{dx}{dt} = f(x,t), \quad x(t_0) = a
   * \f]
   * be an ODE initial-value problem.  This class implements a single time
   * step of an explicit Taylor polynomial time integration method for
   * computing numerical solutions to the IVP.  The method consists of 
   * computing a local truncated Taylor series solution to the ODE (section
   * \ref Rythmos_ETI_local_TS), estimating a step size within the radius
   * of convergence of the Taylor series (section \ref Rythmos_ETI_stepsize)
   * and then summing the polynomial at that step to compute the next
   * point in the numerical integration (section \ref Rythmos_ETI_sum).
   * The algorithmic parameters to the method are controlled through the
   * <tt> params </tt> argument to the constructor which are described in
   * section \ref Rythmos_ETI_params.
   *
   * \section Rythmos_ETI_local_TS Computing the Taylor Polynomial
   * 
   * Let
   * \f[
   *      x(t) = \sum_{k=0}^\infty x_k (t-t_0)^k 
   * \f]
   * be a power series solution to the IVP above.  Then \f$f(x(t))\f$ can
   * be expaned in a power series along the solution curve \f$x(t)\f$:
   * \f[
   *      f(x(t),t) = \sum_{k=0}^\infty f_k (t-t_0)^k
   * \f]
   * where 
   * \f[
   *      f_k = \left.\frac{1}{k!}\frac{d^k}{dt^k} f(x(t),t)\right|_{t=t_0}.
   * \f]
   * By differentiating the power series for \f$x(t)\f$ to compute a power
   * series for \f$dx/dt\f$ and then comparing coefficients in the 
   * equation \f$dx/dt=f(x(t),t)\f$ we find the following recurrence
   * relationship for the Taylor coefficients \f$\{x_k\}\f$:
   * \f[
   *     x_{k+1} = \frac{1}{k+1} f_k, \quad k=0,1,\dots
   * \f]
   * where each coefficient \f$f_k\f$ is a function only of 
   * \f$x_0,\dots,x_k\f$ and can be efficiently computed using the Taylor
   * polynomial mode of automatic differentation.  This allows the Taylor
   * coefficients to be iteratively computed starting with the initial point
   * \f$x_0\f$ up to some fixed degree \f$d\f$ to yield a local truncated 
   * Taylor polynomial approximating the solution to the IVP.
   *
   * \section Rythmos_ETI_stepsize Computing a Step Size
   *
   * With the truncated Taylor polynomial solution 
   * \f$\tilde{x}(t) = \sum_{k=0}^d x_k (t-t_0)^k\f$ in hand, a step size
   * is chosen by estimating the truncation error in the polynomial solution
   * and forcing this error to be less than some prescribed tolerance.  Let
   * \f[
   *     \rho = \max_{d/2\leq k\leq d} (1+\|x_k\|_\infty)^{1/k}
   * \f]
   * so \f$\|x_k\|_\infty\leq\rho^k\f$ for \f$d/2\leq k \leq d\f$.  Assume 
   * \f$\|x_k\|\leq\rho^k\f$ for \f$k>d\f$ as well, then for any \f$h<1/\rho\f$
   * it can be shown that the truncation error is bounded by
   * \f[
   *      \frac{(\rho h)^{d+1}}{1-\rho h}.
   * \f]
   * A step size \f$h\f$ is then given by
   * \f[
   *      h = \exp\left(\frac{1}{d+1}\log\varepsilon-\log\rho\right)
   * \f]
   * for some error tolerance \f$\varepsilon\f$ given an error of approximatly
   * \f$\varepsilon\f$.
   *
   * \section Rythmos_ETI_sum Summing the Polynomial
   *
   * With a step size \f$h\f$ computed, 
   * \f[
   *     x^\ast = \sum_{k=0}^d x_k h^k
   * \f]
   * is used as the next integration point where a new Taylor series is
   * calculated.  Local error per step can also be controlled by computing
   * \f$\|dx^\ast/dt - f(x^\ast)\|_\infty\f$.  If this error is too large,
   * the step size can be reduced to an appropriate size.
   *
   * \section Rythmos_ETI_params Parameters
   *
   * This method recognizes the following algorithmic parameters that can
   * be set in the <tt> params </tt> argument to the constructor:
   * <ul>
   * <li> "Initial Time" (Scalar) [Default = 0] Initial integration time
   * <li> "Final Time"   (Scalar) [Default = 1] Final integration time
   * <li> "Local Error Tolerance" (Magnitude) [Default = 1.0e-10] Error tolerance on \f$\|dx^\ast/dt - f(x^\ast)\|_\infty\f$ as described above.
   * <li> "Minimum Step Size" (Scalar) [Default = 1.0e-10] Minimum step size
   * <li> "Maximum Step Size" (Scalar) [Default = 1.0] Maximum step size
   * <li> "Taylor Polynomial Degree" (int) [Default = 40] Degree of local Taylor polynomial approximation.
   * </ul>
   */

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType ScalarMag;

  template<class Scalar>
  class ExplicitTaylorPolynomialStepper : public Stepper<Scalar>
  {
    public:

    //! Typename of magnitude of scalars
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    
    /*! 
     * Constructor.  See description above for the list of parameters this
     * class uses in \c params.
     */
    ExplicitTaylorPolynomialStepper(
        const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model, 
        Teuchos::ParameterList& params
        );

    /** \brief . */
    void setModel(const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model);
    
    //! Destructor
    ~ExplicitTaylorPolynomialStepper();
    
    //! Take a time step of magnitude \c dt
    Scalar TakeStep(Scalar dt, StepSizeType flag);

    /** \brief . */
    const StepStatus<Scalar> getStepStatus();

    /** \brief . */
    void reInitialize();

    /** \brief . */
    std::string description() const;

    /** \brief . */
    std::ostream& describe(
      std::ostream                &out
      ,const Teuchos::EVerbosityLevel      verbLevel
      ,const std::string          leadingIndent
      ,const std::string          indentSpacer
      ) const;

    /// Redefined from InterpolationBufferBase 
    /// Add points to buffer
    bool SetPoints(
      const std::vector<Scalar>& time_list
      ,const std::vector<Thyra::VectorBase<Scalar> >& x_list
      ,const std::vector<Thyra::VectorBase<Scalar> >& xdot_list
      );

    /// Get values from buffer
    bool GetPoints(
      const std::vector<Scalar>& time_list
      ,std::vector<Thyra::VectorBase<Scalar> >* x_list
      ,std::vector<Thyra::VectorBase<Scalar> >* xdot_list
      ,std::vector<ScalarMag>* accuracy_list
      ) const;

    /// Fill data in from another interpolation buffer
    bool SetRange(
      const Scalar& time_lower
      ,const Scalar& time_upper
      ,const InterpolationBufferBase<Scalar> & IB
      );

    /// Get interpolation nodes
    bool GetNodes(std::vector<Scalar>* time_list) const;

    /// Remove interpolation nodes
    bool RemoveNodes(std::vector<Scalar>* time_list) const;

    /// Get order of interpolation
    int GetOrder() const;

    private:

    //! Computes a local Taylor series solution to the ODE
    void computeTaylorSeriesSolution_();

    /*! 
     * \brief Computes of log of the estimated radius of convergence of the 
     * Taylor series.
     */
    magnitude_type estimateLogRadius_();

    //! Underlying model
    Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > model_;

    //! Current solution vector
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x_vector_;

    //! Vector store approximation to \f$dx/dt\f$
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > x_dot_vector_;

    //! Vector store ODE residual
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > f_vector_;

    //! Polynomial for x
    Teuchos::RefCountPtr<Teuchos::Polynomial<Thyra::VectorBase<Scalar> > > x_poly_;

    //! Polynomial for f
    Teuchos::RefCountPtr<Teuchos::Polynomial<Thyra::VectorBase<Scalar> > > f_poly_;

    //! Current time
    Scalar t_;

    //! Current step size
    Scalar dt_;

    //! Initial integration time
    Scalar t_initial_;

    //! Final integration time
    Scalar t_final_;

    //! Local error tolerance for each time step
    magnitude_type local_error_tolerance_;

    //! Smallest acceptable time step size
    Scalar min_step_size_;

    //! Largest acceptable time step size
    Scalar max_step_size_;

    //! Degree of local Taylor series expansion
    unsigned int degree_;

    //! Used in time step size computation
    Scalar linc_;
  };

  //! Reduction operator for a logarithmic infinity norm
  /*!
   * This class implements a reduction operator for computing the 
   * logarithmic infinity norm of a vector:
   * \f[
   *      \|1 + log(x)\|_\infty.
   * \f]
   */
  template <typename Scalar>
    class ROpLogNormInf : public RTOpPack::ROpScalarReductionBase<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> {
      public:

        /** \brief . */
        ROpLogNormInf() : RTOpPack::RTOpT<Scalar>("ROpLogInfNorm"){}

        /** \brief . */
        typename Teuchos::ScalarTraits<Scalar>::magnitudeType 
          operator() (const RTOpPack::ReductTarget& reduct_obj) const { 
            return this->getRawVal(reduct_obj); 
          }

        /** @name Overridden from RTOpT */
        //@{

        /** \brief . */
        void reduce_reduct_objs(const RTOpPack::ReductTarget& _in_reduct_obj, 
            RTOpPack::ReductTarget* _inout_reduct_obj) const
        {
          const RTOpPack::ReductTargetScalar<Scalar>& in_reduct_obj = 
            Teuchos::dyn_cast<const RTOpPack::ReductTargetScalar<Scalar> >(_in_reduct_obj);
          RTOpPack::ReductTargetScalar<Scalar>& inout_reduct_obj = 
            Teuchos::dyn_cast<RTOpPack::ReductTargetScalar<Scalar> >(*_inout_reduct_obj);
          if(in_reduct_obj.get() > inout_reduct_obj.get()) 
            inout_reduct_obj.set(in_reduct_obj.get());
        }

        /** \brief . */
        void apply_op(const int num_vecs, 
            const RTOpPack::SubVectorT<Scalar> sub_vecs[],
            const int num_targ_vecs, 
            const RTOpPack::MutableSubVectorT<Scalar> targ_sub_vecs[], 
            RTOpPack::ReductTarget *_reduct_obj) const
        {
          RTOpPack::ReductTargetScalar<Scalar>& reduct_obj = 
            Teuchos::dyn_cast<RTOpPack::ReductTargetScalar<Scalar> >(*_reduct_obj);

          RTOP_APPLY_OP_1_0(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);

          typename Teuchos::ScalarTraits<Scalar>::magnitudeType norm_inf = 
            reduct_obj.get();

          // unit stride
          if(v0_s == 1) {
            for(RTOp_index_type i=0; i<subDim; i++) {
              typename Teuchos::ScalarTraits<Scalar>::magnitudeType
                mag = Teuchos::ScalarTraits<Scalar>::magnitude(*v0_val++);
              mag = std::log(Teuchos::ScalarTraits<Scalar>::one() + mag);
              norm_inf = mag > norm_inf ? mag : norm_inf;
            }
          } else {
            for(RTOp_index_type i=0; i<subDim; i++, v0_val+=v0_s) {
              typename Teuchos::ScalarTraits<Scalar>::magnitudeType
                mag = Teuchos::ScalarTraits<Scalar>::magnitude(*v0_val);
              mag = std::log(Teuchos::ScalarTraits<Scalar>::one() + mag);
              norm_inf = mag > norm_inf ? mag : norm_inf;
            }
          }
          reduct_obj.set(norm_inf);
        }
        //@}

    }; // class ROpLogNormInf

  //! Computs logarithmic infinity norm of a vector using ROpLogNormInf.
  template <typename Scalar>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  log_norm_inf(const Thyra::VectorBase<Scalar>& x)
  {
    ROpLogNormInf<Scalar> log_norm_inf_op;
    Teuchos::RefCountPtr<RTOpPack::ReductTarget> log_norm_inf_targ = 
      log_norm_inf_op.reduct_obj_create();
    const Thyra::VectorBase<Scalar>* vecs[] = { &x };
    Thyra::applyOp<Scalar>(log_norm_inf_op,1,vecs,0,
			   (Thyra::VectorBase<Scalar>**)NULL,
			   log_norm_inf_targ.get());
    
    return log_norm_inf_op(*log_norm_inf_targ);
  }

  template<class Scalar>
  ExplicitTaylorPolynomialStepper<Scalar>::ExplicitTaylorPolynomialStepper(
	const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &m,
	Teuchos::ParameterList& params)
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;

    // Get initial time
    t_initial_ = params.get("Initial Time", ST::zero());

    // Get final time
    t_final_ = params.get("Final Time", ST::one());

    // Get local error tolerance
    local_error_tolerance_ = 
      params.get("Local Error Tolerance", magnitude_type(1.0e-10));

    // Get minimum step size
    min_step_size_ = params.get("Minimum Step Size", Scalar(1.0e-10));

    // Get maximum step size
    max_step_size_ = params.get("Maximum Step Size", Scalar(1.0));

    // Get degree_ of Taylor polynomial expansion
    degree_ = params.get("Taylor Polynomial Degree", 40);

    linc_ = Scalar(-16.0*std::log(10.0)/degree_);
  
    model_ = m;
    t_ = t_initial_;
    x_vector_ = model_->getNominalValues().get_x()->clone_v();
    x_dot_vector_ = Thyra::createMember(model_->get_x_space());
    f_vector_ = Thyra::createMember(model_->get_f_space());
    x_poly_ = 
      Teuchos::rcp(new Teuchos::Polynomial<Thyra::VectorBase<Scalar> >(0, 
								     *x_vector_,
								     degree_));
    f_poly_ = 
      Teuchos::rcp(new Teuchos::Polynomial<Thyra::VectorBase<Scalar> >(0, 
								     *f_vector_,
								     degree_));
  }

  template<class Scalar>
  ExplicitTaylorPolynomialStepper<Scalar>::~ExplicitTaylorPolynomialStepper()
  {
  }

  template<class Scalar>
  Scalar 
  ExplicitTaylorPolynomialStepper<Scalar>::TakeStep(Scalar dt, StepSizeType flag)
  {
    if (flag == VARIABLE_STEP) {
      // If t_ > t_final_, we're done
      if (t_ > t_final_) {
        dt_ = Teuchos::ScalarTraits<Scalar>::zero();
        return dt_;
      }

      // Compute a local truncated Taylor series solution to system
      computeTaylorSeriesSolution_();

      // Estimate log of radius of convergence of Taylor series
      Scalar rho = estimateLogRadius_();

      // Set step size
      Scalar dt = std::exp(linc_ - rho);

      // If step size is too big, reduce
      if (dt > max_step_size_) {
        dt = max_step_size_;
      }

      // If step goes past t_final_, reduce
      if (t_+dt > t_final_) {
        dt = t_final_-t_;
      }

      magnitude_type local_error;

      do {

        // compute x(t_+dt), xdot(t_+dt)
        x_poly_->evaluate(dt, x_vector_.get(), x_dot_vector_.get());

        // compute f( x(t_+dt), t_+dt )
        Thyra::eval_f(*model_, *x_vector_, t_+dt, f_vector_.get());

        // compute || xdot(t_+dt) - f( x(t_+dt), t_+dt ) ||
        Thyra::Vp_StV(x_dot_vector_.get(), -Teuchos::ScalarTraits<Scalar>::one(),
          *f_vector_);
        local_error = norm_inf(*x_dot_vector_);

        if (local_error > local_error_tolerance_) {
          dt *= 0.7;
        }

      } while (local_error > local_error_tolerance_ && dt > min_step_size_);

      // Check if minimum step size was reached
      TEST_FOR_EXCEPTION(dt < min_step_size_, 
            std::runtime_error,
            "ExplicitTaylorPolynomialStepper<Scalar>::takeStep(): " 
            << "Step size reached minimum step size " 
            << min_step_size_ << ".  Failing step." );

      // Increment t_
      t_ += dt;

      dt_ = dt;

      return dt;

    } else {

      // If t_ > t_final_, we're done
      if (t_ > t_final_) {
        dt_ = Teuchos::ScalarTraits<Scalar>::zero();
        return dt_;
      }

      // Compute a local truncated Taylor series solution to system
      computeTaylorSeriesSolution_();

      // If step size is too big, reduce
      if (dt > max_step_size_) {
        dt = max_step_size_;
      }

      // If step goes past t_final_, reduce
      if (t_+dt > t_final_) {
        dt = t_final_-t_;
      }

      // compute x(t_+dt)
      x_poly_->evaluate(dt, x_vector_.get());

      // Increment t_
      t_ += dt;

      dt_ = dt;

      return dt;
    }
  }


  template<class Scalar>
  const StepStatus<Scalar> ExplicitTaylorPolynomialStepper<Scalar>::getStepStatus()
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    StepStatus<Scalar> stepStatus;

    stepStatus.stepSize = dt_;
    stepStatus.order = this->GetOrder();
    stepStatus.time = t_;
    stepStatus.solution = x_vector_;
    stepStatus.solutionDot = x_dot_vector_;
    stepStatus.residual = f_vector_;

    return(stepStatus);
  }

  template<class Scalar>
  void ExplicitTaylorPolynomialStepper<Scalar>::reInitialize()
  {
  }

  template<class Scalar>
  void
  ExplicitTaylorPolynomialStepper<Scalar>::computeTaylorSeriesSolution_()
  {
    Teuchos::RefCountPtr<Thyra::VectorBase<Scalar> > tmp;

    // Set degree_ of polynomials to 0
    x_poly_->setDegree(0);
    f_poly_->setDegree(0);

    // Set degree_ 0 coefficient
    x_poly_->setCoefficient(0, *x_vector_);

    for (unsigned int k=1; k<=degree_; k++) {

      // compute [f] = f([x])
      Thyra::eval_f_poly(*model_, *x_poly_, t_, f_poly_.get());

      x_poly_->setDegree(k);
      f_poly_->setDegree(k);
      
      // x[k] = f[k-1] / k
      tmp = x_poly_->getCoefficient(k);
      copy(*(f_poly_->getCoefficient(k-1)), tmp.get());
      scale(Scalar(1.0)/Scalar(k), tmp.get());
    }

  }

  template<class Scalar>
  typename ExplicitTaylorPolynomialStepper<Scalar>::magnitude_type
  ExplicitTaylorPolynomialStepper<Scalar>::estimateLogRadius_()
  {
    magnitude_type rho = 0;
    magnitude_type tmp;
    for (unsigned int k=degree_/2; k<=degree_; k++) {
      tmp = log_norm_inf(*(x_poly_->getCoefficient(k))) / k;
      if (tmp > rho) {
        rho = tmp;
      }
    }
    return rho;
  }

  template<class Scalar>
  std::string ExplicitTaylorPolynomialStepper<Scalar>::description() const
  {
    std::string name = "Rythmos::ExplicitTaylorPolynomialStepper";
    return(name);
  }

  template<class Scalar>
  std::ostream& ExplicitTaylorPolynomialStepper<Scalar>::describe(
        std::ostream                &out
        ,const Teuchos::EVerbosityLevel      verbLevel
        ,const std::string          leadingIndent
        ,const std::string          indentSpacer
        ) const
  {
    if (verbLevel == Teuchos::VERB_EXTREME) {
      out << description() << "::describe" << std::endl;
      out << "model_ = " << std::endl;
      out << model_->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
      out << "x_vector_ = " << std::endl;
      out << x_vector_->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
      out << "x_dot_vector_ = " << std::endl;
      out << x_dot_vector_->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
      out << "f_vector_ = " << std::endl;
      out << f_vector_->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
      out << "x_poly_ = " << std::endl;
      out << x_poly_->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
      out << "f_poly_ = " << std::endl;
      out << f_poly_->describe(out,verbLevel,leadingIndent,indentSpacer) << std::endl;
      out << "t_ = " << t_ << std::endl;
      out << "t_initial_ = " << t_initial_ << std::endl;
      out << "t_final_ = " << t_final_ << std::endl;
      out << "local_error_tolerance_ = " << local_error_tolerance_ << std::endl;
      out << "min_step_size_ = " << min_step_size_ << std::endl;
      out << "max_step_size_ = " << max_step_size_ << std::endl;
      out << "degree_ = " << degree_ << std::endl;
      out << "linc_ = " << linc_ << std::endl;
    }
    return(out);
  }


template<class Scalar>
bool ExplicitTaylorPolynomialStepper<Scalar>::SetPoints(
    const std::vector<Scalar>& time_list
    ,const std::vector<Thyra::VectorBase<Scalar> >& x_list
    ,const std::vector<Thyra::VectorBase<Scalar> >& xdot_list)
{
  return(false);
}

template<class Scalar>
bool ExplicitTaylorPolynomialStepper<Scalar>::GetPoints(
    const std::vector<Scalar>& time_list
    ,std::vector<Thyra::VectorBase<Scalar> >* x_list
    ,std::vector<Thyra::VectorBase<Scalar> >* xdot_list
    ,std::vector<ScalarMag>* accuracy_list) const
{
  return(false);
}

template<class Scalar>
bool ExplicitTaylorPolynomialStepper<Scalar>::SetRange(
    const Scalar& time_lower
    ,const Scalar& time_upper
    ,const InterpolationBufferBase<Scalar>& IB)
{
  return(false);
}

template<class Scalar>
bool ExplicitTaylorPolynomialStepper<Scalar>::GetNodes(std::vector<Scalar>* time_list) const
{
  return(false);
}

template<class Scalar>
bool ExplicitTaylorPolynomialStepper<Scalar>::RemoveNodes(std::vector<Scalar>* time_list) const
{
  return(false);
}


template<class Scalar>
int ExplicitTaylorPolynomialStepper<Scalar>::GetOrder() const
{
  return(4);
}

template<class Scalar>
void ExplicitTaylorPolynomialStepper<Scalar>::setModel(const Teuchos::RefCountPtr<const Thyra::ModelEvaluator<Scalar> > &model)
{
  TEST_FOR_EXCEPT(model == Teuchos::null)
  model_ = model;
}


} // namespace Rythmos

#endif // RYTHMOS_EXPLICIT_TAYLOR_POLYNOMIAL_STEPPER_H
