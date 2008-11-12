// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_TABULAR_OUTPUTTER_HPP
#define TEUCHOS_TABULAR_OUTPUTTER_HPP


#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Time.hpp"


namespace Teuchos {


/** \brief Utility class that makes it easy to create formatted tables
 * of output.
 *
 */
class TabularOutputter {
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

  //@}

  /** \brief . */
  TabularOutputter();

  /** \brief . */
  TabularOutputter(std::ostream &out);

  /** \brief Set the ostream that all output will be sent to. */
  void setOStream( const RCP<std::ostream> &out );

  /** \brief Add a new field to be output. */
  void pushField( const std::string &fieldName,
    const EFieldType fieldType = DOUBLE,
    const EFieldJustification fieldJustification = RIGHT,
    const EFloatingOutputType floatingOutputType = SCIENTIFIC );

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
  void nextRow();

private:

  // Private types

  struct FieldSpec {
    FieldSpec(std::string fieldName_in, EFieldType fieldType_in,
      EFieldJustification fieldJustification_in,
      EFloatingOutputType floatingOutputType_in
      )
      :fieldName(fieldName_in), fieldType(fieldType_in),
       fieldJustification(fieldJustification_in),
       floatingOutputType(floatingOutputType_in),
       outputWidth(-1),
       precision(-1)
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

  Array<FieldSpec> fieldSpecs_;
  RCP<FancyOStream> out_;
  Tuple<int,numFieldTypes> fieldTypePrecision_;

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
      TEST_FOR_EXCEPT(numLoops_ == -1);
#endif      
      timer_.stop();
      const double relTime = 
        adjustTime(timer_.totalElapsedTime()) / numLoops_;
      numLoops_ = -1;
      return relTime;
    }

};


/** \brief Start a timer block using a TabularOutputter object . */
#define TEUCHOS_START_PERF_OUTPUT_TIMER(OUTPUTTER, NUMLOOPS) \
  (OUTPUTTER).startTimer(NUMLOOPS); \
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
  TEUCHOS_ASSERT(currFieldIdx_ < as<int>(fieldSpecs_.size())); 
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
      TEST_FOR_EXCEPT(true);
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
      TEST_FOR_EXCEPT(true);
    }
  }

  *out_ << setw(fieldSpec.outputWidth) << t;

  ++currFieldIdx_;
  
}



} // namespace Teuchos


#endif // TEUCHOS_TABULAR_OUTPUTTER_HPP
