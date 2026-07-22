// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_EVALUATION_LOGGER_MODEL_EVALUATOR_HPP
#define THYRA_DEFAULT_EVALUATION_LOGGER_MODEL_EVALUATOR_HPP

#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Teuchos_Time.hpp"

namespace Thyra {


/** \brief This class wraps any ModelEvaluator object and logs the evaluation
 * of various functions.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Nonlin_ME_support_grp
 */
template<class Scalar>
class DefaultEvaluationLoggerModelEvaluator
  : virtual public ModelEvaluatorDelegatorBase<Scalar>
{
public:

  /** \name Constructors/initializers/accessors/utilities. */
  //@{

  /** \brief . */
  DefaultEvaluationLoggerModelEvaluator();

  /** \brief . */
  DefaultEvaluationLoggerModelEvaluator(
    const RCP<ModelEvaluator<Scalar> >   &thyraModel
    ,const RCP<std::ostream>             &tableOut
    );

  /** \brief Initalize.
   *
   * \param  thyraModel     [in] Model being wrapped.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>thyraModel.get()!=NULL</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->getUnderlyingModel.get() == thyraModel.get()</tt>
   * </ul>
   */
  void initialize(
    const RCP<ModelEvaluator<Scalar> >   &thyraModel
    ,const RCP<std::ostream>             &tableOut
    );

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase */
  //@{

  /** \brief . */
  void evalModelImpl(
    const ModelEvaluatorBase::InArgs<Scalar>    &inArgs
    ,const ModelEvaluatorBase::OutArgs<Scalar>  &outArgs
    ) const;

  //@}

private:

  RCP<std::ostream> tableOut_;
  Teuchos::Time timer_;
  
  mutable bool headerPrinted_;
  mutable bool supports_f_;
  mutable bool supports_W_;
  
  static const int flt_width_;
  static const int flt_sciPrec_;
  static const int flt_prec_;
  static const char flt_line_[];
  static const int int_width_;
  static const char int_line_[];
  
  void printHeader( const ModelEvaluatorBase::OutArgs<Scalar> &outArgs ) const;
  void printLine( const ModelEvaluatorBase::OutArgs<Scalar> &outArgs ) const;
  
};

// /////////////////////////////////
// Implementations

// Constructors/initializers/accessors/utilities

template<class Scalar>
const int DefaultEvaluationLoggerModelEvaluator<Scalar>::flt_width_ = 25; 
template<class Scalar>
const int DefaultEvaluationLoggerModelEvaluator<Scalar>::flt_sciPrec_  = 16;
template<class Scalar>
const int DefaultEvaluationLoggerModelEvaluator<Scalar>::flt_prec_  = 16;
template<class Scalar>
const char DefaultEvaluationLoggerModelEvaluator<Scalar>::flt_line_[]  = "-------------------------";
template<class Scalar>
const int DefaultEvaluationLoggerModelEvaluator<Scalar>::int_width_ = 10; 
template<class Scalar>
const char DefaultEvaluationLoggerModelEvaluator<Scalar>::int_line_[]  = "----------";

template<class Scalar>
DefaultEvaluationLoggerModelEvaluator<Scalar>::DefaultEvaluationLoggerModelEvaluator()
  :timer_(""),headerPrinted_(false)
{}

template<class Scalar>
DefaultEvaluationLoggerModelEvaluator<Scalar>::DefaultEvaluationLoggerModelEvaluator(
  const RCP<ModelEvaluator<Scalar> >   &thyraModel
  ,const RCP<std::ostream>             &tableOut
  )
  :timer_(""),headerPrinted_(false), supports_f_(false), supports_W_(false)
{
  initialize(thyraModel,tableOut);
}

template<class Scalar>
void DefaultEvaluationLoggerModelEvaluator<Scalar>::initialize(
  const RCP<ModelEvaluator<Scalar> >   &thyraModel
  ,const RCP<std::ostream>             &tableOut
  )
{
  TEUCHOS_TEST_FOR_EXCEPT( tableOut.get()==NULL );
  this->ModelEvaluatorDelegatorBase<Scalar>::initialize(thyraModel);
  tableOut_ = tableOut;
  timer_.start(true);
  headerPrinted_ = false;
}


// Public functions overridden from Teuchos::Describable


template<class Scalar>
std::string DefaultEvaluationLoggerModelEvaluator<Scalar>::description() const
{
  const RCP<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  std::ostringstream oss;
  oss << "Thyra::DefaultEvaluationLoggerModelEvaluator{";
  oss << "thyraModel=";
  if(thyraModel.get())
    oss << "\'"<<thyraModel->description()<<"\'";
  else
    oss << "NULL";
  oss << "}";
  return oss.str();
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
void DefaultEvaluationLoggerModelEvaluator<Scalar>::evalModelImpl(
  const ModelEvaluatorBase::InArgs<Scalar>     &inArgs
  ,const ModelEvaluatorBase::OutArgs<Scalar>   &outArgs
  ) const
{

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_BEGIN(
    "Thyra::DefaultEvaluationLoggerModelEvaluator",inArgs,outArgs
    );

  thyraModel->evalModel(inArgs,outArgs);

  if(!headerPrinted_) {
    printHeader(outArgs);
    headerPrinted_ = true;
  }
  printLine(outArgs);

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();
  
}


// private


template<class Scalar>
void DefaultEvaluationLoggerModelEvaluator<Scalar>::printHeader(
  const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{

  using std::setw;
  using std::setprecision;
  using std::right;
  using std::left;
  typedef ModelEvaluatorBase MEB;

  supports_f_ = outArgs.supports(MEB::OUT_ARG_f);
  supports_W_ = outArgs.supports(MEB::OUT_ARG_W);

  const int Ng = outArgs.Ng();

  *tableOut_
    << "\n***"
    << "\n*** Table of function evaluations vs. CPU time"
    << "\n***\n";

  *tableOut_
    << "\nModel Evaluator Description:\n" << Teuchos::describe(*this,Teuchos::VERB_LOW);
  
  *tableOut_ << "\n";
  *tableOut_ << "  " << left << setw(flt_width_) << "time(s)";
  for( int j = 0; j < Ng; ++j ) {
    std::ostringstream oss;
    oss << "||g("<<j<<")||";
    *tableOut_ << "  " << left << setw(flt_width_) << oss.str();
  }
  if(supports_f_)
    *tableOut_ << "  " << left << setw(flt_width_) << "||f||";
  if(supports_W_)
    *tableOut_ << "  " << left << setw(int_width_) << "Calc W";
  *tableOut_ << "\n";
  
  *tableOut_ << "  " << left << setw(flt_width_) << flt_line_;   // time(s)
  for( int j = 0; j < Ng; ++j )
    *tableOut_ << "  " << left << setw(flt_width_) << flt_line_; // ||g(j)||
  if(supports_f_)
    *tableOut_ << "  " << left << setw(flt_width_) << flt_line_; // ||f||
  if(supports_W_)
    *tableOut_ << "  " << left << setw(int_width_) << int_line_; // Calc W
  *tableOut_ << "\n";

}

template<class Scalar>
void DefaultEvaluationLoggerModelEvaluator<Scalar>::printLine(
  const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{

  using std::right;
  using std::left;
  using std::setprecision;
  using std::setw;

  const int Ng = outArgs.Ng();

  RCP<const VectorBase<Scalar> > f, g_j;
  
  *tableOut_ << "  " << setprecision(flt_prec_) << right << setw(flt_width_) << timer_.totalElapsedTime(true);
  for( int j = 0; j < Ng; ++j ) {
    if((g_j=outArgs.get_g(j)).get())
      *tableOut_ << "  " << setprecision(flt_sciPrec_) << right << setw(flt_width_) << norm(*g_j);
    else
      *tableOut_ << "  " << right << setw(flt_width_) << "-";
  }
  if(supports_f_) {
    if((f=outArgs.get_f()).get())
      *tableOut_ << "  " << setprecision(flt_sciPrec_) << right << setw(flt_width_) << norm(*f);
    else
      *tableOut_ << "  " << right << setw(flt_width_) << "-";
  }
  if(supports_W_) {
    if(outArgs.get_W().get())
      *tableOut_ << "  " << right << setw(int_width_) << "1";
    else
      *tableOut_ << "  " << right << setw(int_width_) << "-";
  }
  *tableOut_ << "\n";

}

} // namespace Thyra

#endif // THYRA_DEFAULT_EVALUATION_LOGGER_MODEL_EVALUATOR_HPP
