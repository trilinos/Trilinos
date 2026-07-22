// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_TabularOutputter.hpp"
#include "Teuchos_as.hpp"


namespace {


int getFieldWidth(const Teuchos::TabularOutputter::EFieldType fieldType,
  const int prec)
{
  typedef Teuchos::TabularOutputter TO;
  switch(fieldType)
  {
    case TO::DOUBLE:
      return prec + 8; // Leave room sign and exponent
    case TO::INT:
      return prec+1; // leave room for sign
    case TO::STRING:
      return prec;
  }
  return -1; // Will never be called
}


const std::string getFieldLine(const int width)
{
  std::string line;
  line.append(width, '-');
  return line;
}


} // namespace


namespace Teuchos {


const std::string TabularOutputter::fieldSpacer_("  ");


TabularOutputter::TabularOutputter(std::ostream &out):
  timer_(""),
  numLoops_(0)
{
  initialize();
  setOStream(rcpFromRef(out));
}


TabularOutputter::TabularOutputter(const RCP<std::ostream> &out):
  timer_(""),
  numLoops_(0)
{
  initialize();
  setOStream(out);
}


void TabularOutputter::setOStream( const RCP<std::ostream> &out )
{
#ifdef TEUCHOS_DEBUG
  out.assert_not_null();
#endif
  out_ = fancyOStream(out);
}


void TabularOutputter::pushFieldSpec(
  const std::string &fieldName, const EFieldType fieldType,
  const EFieldJustification fieldJustification,
  const EFloatingOutputType floatingOutputType,
  const int width
  )
{
#ifdef TEUCHOS_DEBUG
  if (width > 0) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      !(as<int>(fieldName.size()) <= width),
      InvalidFieldSpecError,
      "Error, the length of the field name \""<<fieldName<<"\"\n"
      "is "<<fieldName.size()<<" which is larger than the\n"
      "specifically set field width "<<width<<"!"
      );
  }
#endif
  fieldSpecs_.push_back(
    FieldSpec(fieldName, fieldType, fieldJustification, floatingOutputType,
      TEUCHOS_MAX(as<int>(fieldName.size()), width)
      )
    );
}


void TabularOutputter::setFieldTypePrecision( const EFieldType fieldType,
  const int prec )
{
  fieldTypePrecision_[fieldType] = prec;
}


void TabularOutputter::outputHeader()
{

  using std::left;
  using std::setw;

  const int numFields = static_cast<int>(fieldSpecs_.size());

#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    numFields==0, MissingFieldsError,
    "Error, you must add at least one field spec using pushFieldSpec(...)!"
    );
#endif


  for (int i = 0; i < numFields; ++i) {
    FieldSpec &fieldSpec = fieldSpecs_[i];
    const EFieldType fieldType = fieldSpec.fieldType;
    const int fieldTypePrecision = fieldTypePrecision_[fieldType];
    fieldSpec.precision = fieldTypePrecision;
    const int fieldPrecisionWidth =
      getFieldWidth(fieldType, fieldTypePrecision);
    if (fieldSpec.outputWidth < fieldPrecisionWidth) {
      fieldSpec.outputWidth = fieldPrecisionWidth;
    }
    *out_ << fieldSpacer_ << left << setw(fieldSpec.outputWidth) << fieldSpec.fieldName;
  }
  *out_ << "\n";

  for (int i = 0; i < numFields; ++i) {
    const FieldSpec &fieldSpec = fieldSpecs_[i];
    *out_ << fieldSpacer_ << left << setw(fieldSpec.outputWidth) << getFieldLine(fieldSpec.outputWidth);
  }
  *out_ << "\n";

  currFieldIdx_ = 0;

}


void TabularOutputter::nextRow(const bool allowRemainingFields)
{
  const int numFields = static_cast<int>(fieldSpecs_.size());
  if (allowRemainingFields) {
    while (currFieldIdx_ < numFields) {
      outputField("-");
    }
  }
  else {
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      !(currFieldIdx_ == numFields),
      InvalidFieldOutputError,
      "Error, you must call outputField(...) for every field in the row\n"
      "before you call nextRow()!"
      );
#endif
  }
  *out_ << "\n";
  currFieldIdx_ = 0;
}


// Private member functions


void TabularOutputter::initialize()
{
  std::fill( fieldTypePrecision_.begin(), fieldTypePrecision_.end(), 4 );
  currFieldIdx_ = -1;
}


} // namespace Teuchos
