// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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


TabularOutputter::TabularOutputter(std::ostream &out)
  :timer_("")
{
  initialize();
  setOStream(rcpFromRef(out));
}


TabularOutputter::TabularOutputter(const RCP<std::ostream> &out)
  :timer_("")
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

  const int numFields = fieldSpecs_.size();

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
  const int numFields = fieldSpecs_.size();
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
