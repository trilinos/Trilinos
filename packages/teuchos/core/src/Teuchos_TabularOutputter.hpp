// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_TABULAR_OUTPUTTER_HPP
#define TEUCHOS_TABULAR_OUTPUTTER_HPP


#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_Exceptions.hpp"


namespace Teuchos {


/** \brief Utility class that makes it easy to create formatted tables
 * of output.
 *
 */
class TEUCHOSCORE_LIB_DLL_EXPORT TabularOutputter {
public:

  /** \name Public types */
  //@{

  /** \brief . */
  enum EFieldType { DOUBLE, INT, STRING };
  enum { numFieldTypes = 3 };

  /** \brief . */
  enum EFieldJustification { LEFT, RIGHT };
  enum { numFieldJustifications = 2 };

  /** \brief . */
  enum EFloatingOutputType { SCIENTIFIC, GENERAL };
  enum { numFloatingOutputTypes = 2 };

  /** \brief .  */
  class MissingFieldsError : public ExceptionBase
  {public:MissingFieldsError(const std::string& what_arg) : ExceptionBase(what_arg) {}};

  /** \brief .  */
  class InvalidFieldSpecError : public ExceptionBase
  {public:InvalidFieldSpecError(const std::string& what_arg) : ExceptionBase(what_arg) {}};

  /** \brief .  */
  class MissingHeaderError : public ExceptionBase
  {public:MissingHeaderError(const std::string& what_arg) : ExceptionBase(what_arg) {}};

  /** \brief .  */
  class InvalidFieldOutputError : public ExceptionBase
  {public:InvalidFieldOutputError(const std::string& what_arg) : ExceptionBase(what_arg) {}};

  //@}

  /** \brief . */
  TabularOutputter(std::ostream &out);

  /** \brief . */
  TabularOutputter(const RCP<std::ostream> &out);

  /** \brief Set the ostream that all output will be sent to. */
  void setOStream( const RCP<std::ostream> &out );

  /** \brief Add a new field to be output. */
  void pushFieldSpec( const std::string &fieldName,
    const EFieldType fieldType = DOUBLE,
    const EFieldJustification fieldJustification = RIGHT,
    const EFloatingOutputType floatingOutputType = SCIENTIFIC,
    const int width = -1
    );

  /** \brief Set the precision of output for a field.
   *
   * This will also determine the width of the field.
   */
  void setFieldTypePrecision( const EFieldType fieldType, const int prec );

  /** \brief Output the headers. */
  void outputHeader();

  /** \brief Output to the next field. */
  template<typename T>
  void outputField( const T& t );

  /** \brief Finalize the row of output. */
  void nextRow(const bool allowRemainingFields = false);

private:

  // Private types

  struct FieldSpec {
    FieldSpec(std::string fieldName_in, EFieldType fieldType_in,
      EFieldJustification fieldJustification_in,
      EFloatingOutputType floatingOutputType_in,
      const int outputWidth_in
      )
      :fieldName(fieldName_in), fieldType(fieldType_in),
       fieldJustification(fieldJustification_in),
       floatingOutputType(floatingOutputType_in),
       outputWidth(outputWidth_in),
       precision(-1) // Gets set later
      {}
    std::string fieldName;
    EFieldType fieldType;
    EFieldJustification fieldJustification;
    EFloatingOutputType floatingOutputType;
    int outputWidth;
    int precision;
  };

  // Private data members

  static const std::string fieldSpacer_;

//use pragmas to disable some false-positive warnings for windows sharedlibs export
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4251)
#endif
  Array<FieldSpec> fieldSpecs_;
  RCP<FancyOStream> out_;
  Tuple<int,numFieldTypes> fieldTypePrecision_;
#ifdef _MSC_VER
#pragma warning(pop)
#endif

  int currFieldIdx_;

  Time timer_;
  int numLoops_;

  // Private member functions

  void initialize();

  double adjustTime( const double &time_in )
    {
      return ( time_in > 0.0 ? time_in : -1.0 );
    }

public: // Should be hidden

  void startTimer(const int numLoops)
    {
      timer_.reset();
      timer_.start();
      numLoops_ = numLoops;
    }

  double stopTimer()
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPT(numLoops_ == -1);
#endif
      timer_.stop();
      const double relTime =
        adjustTime(timer_.totalElapsedTime()) / numLoops_;
      numLoops_ = -1;
      return relTime;
    }

private:

  // Not defined and not to be called!
  TabularOutputter();

};


/** \brief Start a timer block using a TabularOutputter object . */
#define TEUCHOS_START_PERF_OUTPUT_TIMER(OUTPUTTER, NUMLOOPS) \
  (OUTPUTTER).startTimer(NUMLOOPS); \
  for ( int k = 0; k < (NUMLOOPS); ++k )


/** \brief Start a timer block using a TabularOutputter object . */
#define TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(OUTPUTTER, NUMLOOPS, NUMINNERLOOPS) \
  (OUTPUTTER).startTimer((NUMLOOPS)*(NUMINNERLOOPS)); \
  for ( int k = 0; k < (NUMLOOPS); ++k )


/** \brief Start a timer block using a TabularOutputter object . */
#define TEUCHOS_START_PERF_OUTPUT_TIMER_INNERLOOP(OUTPUTTER, NUMLOOPS, NUMINNERLOOPS) \
  (OUTPUTTER).startTimer((NUMLOOPS)*(NUMINNERLOOPS)); \
  for ( int k = 0; k < (NUMLOOPS); ++k )


/** \brief End a timer block, output the time field to a TabularOutputter
 * object, and set a variable with the time.
 */
#define TEUCHOS_END_PERF_OUTPUT_TIMER(OUTPUTTER, VARNAME) \
  const double VARNAME = (OUTPUTTER).stopTimer(); \
  (OUTPUTTER).outputField(VARNAME)


//
// Implementations
//


template<typename T>
void TabularOutputter::outputField( const T& t )
{

  using std::setw;

#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    currFieldIdx_ == -1,
    MissingHeaderError,
    "Error, you can not output a field until you print the header with"
    " outputHeader()."
    );
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(currFieldIdx_ < as<int>(fieldSpecs_.size())),
    InvalidFieldOutputError,
    "Error, you have already output all of the "
    << fieldSpecs_.size() << " fields for this tabular output."
    "  You must call nextRow() before outputting to the next row."
    );
#endif

  FieldSpec &fieldSpec = fieldSpecs_[currFieldIdx_];

  *out_ << fieldSpacer_ << std::setprecision(fieldSpec.precision);

  switch(fieldSpec.fieldJustification) {
    case LEFT:
      *out_ << std::left;
      break;
    case RIGHT:
      *out_ << std::right;
      break;
    default: {
      TEUCHOS_TEST_FOR_EXCEPT(true);
    }
  }

  switch(fieldSpec.floatingOutputType) {
    case SCIENTIFIC:
      *out_ << std::scientific;
      break;
    case GENERAL:
      *out_ << std::fixed;
      break;
    default: {
      TEUCHOS_TEST_FOR_EXCEPT(true);
    }
  }

  *out_ << setw(fieldSpec.outputWidth) << t;

  ++currFieldIdx_;

}



} // namespace Teuchos


#endif // TEUCHOS_TABULAR_OUTPUTTER_HPP
