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


#include "Teuchos_TabularOutputter.hpp"


namespace {


const int getFieldWidth(const Teuchos::TabularOutputter::EFieldType fieldType,
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


TabularOutputter::TabularOutputter()
  :timer_("")
{
  initialize();
}


TabularOutputter::TabularOutputter(std::ostream &out)
  :timer_("")
{
  initialize();
  setOStream(rcpFromRef(out));
}


void TabularOutputter::setOStream( const RCP<std::ostream> &out )
{
  out_ = fancyOStream(out);
}


void TabularOutputter::pushField(
  const std::string &fieldName, const EFieldType fieldType,
  const EFieldJustification fieldJustification,
  const EFloatingOutputType floatingOutputType
  )
{
  fieldSpecs_.push_back(
    FieldSpec(fieldName, fieldType, fieldJustification, floatingOutputType)
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

  for (int i = 0; i < numFields; ++i) {
    FieldSpec &fieldSpec = fieldSpecs_[i];
    const EFieldType fieldType = fieldSpec.fieldType;
    const int fieldTypePrecision = fieldTypePrecision_[fieldType];
    fieldSpec.precision = fieldTypePrecision;
    fieldSpec.outputWidth = getFieldWidth(fieldType, fieldTypePrecision);
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


void TabularOutputter::nextRow()
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(currFieldIdx_, as<int>(fieldSpecs_.size()));
#endif
  *out_ << "\n";
  currFieldIdx_ = 0;
}


// Private member functions


void TabularOutputter::initialize()
{
  std::fill( fieldTypePrecision_.begin(), fieldTypePrecision_.end(), -1 );
  currFieldIdx_ = -1;
}


} // namespace Teuchos
