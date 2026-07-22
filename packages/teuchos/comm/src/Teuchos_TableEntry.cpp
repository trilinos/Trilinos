// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_TableEntry.hpp"

using namespace Teuchos;


/* --------- base class methods ---------------------------------------------- */

std::string TableEntry::toChoppedString(int maxWidth) const
{
  return toString().substr(0, maxWidth);
}



/* --------- DoubleEntry methods -------------------------------------------- */

DoubleEntry::DoubleEntry(const double& value, int precision, const std::ios_base::fmtflags& flags)
  : TableEntry(), data_(value), precision_(precision), fmtflags_(flags)
{}

std::string DoubleEntry::toString() const
{
  std::ostringstream toss;
  toss.setf(fmtflags_);
  toss << std::setprecision(precision_) << data_;
  return toss.str();
}



/* --------- IntEntry methods -------------------------------------------- */

IntEntry::IntEntry(int value, const std::ios_base::fmtflags& flags)
  : TableEntry(), data_(value), fmtflags_(flags)
{}

std::string IntEntry::toString() const
{
  std::ostringstream toss;
  toss.setf(fmtflags_);
  toss << data_;
  return toss.str();
}



/* --------- StringEntry methods -------------------------------------------- */

StringEntry::StringEntry(std::string value)
  : TableEntry(), data_(value)
{}

std::string StringEntry::toString() const
{
  return data_;
}





/* --------- CompoundEntryWithParentheses methods ------------------------- */

CompoundEntryWithParentheses
::CompoundEntryWithParentheses(const RCP<TableEntry>& first,
                                const RCP<TableEntry>& second,
                                bool spaceBeforeParens)
  : TableEntry(),
    first_(first),
    second_(second),
    spaceBeforeParens_(spaceBeforeParens)
{}

std::string CompoundEntryWithParentheses::toString() const
{
  std::ostringstream toss;

  toss << first_->toString();
  if (spaceBeforeParens_) toss << " ";
  toss << "(" << second_->toString() << ")";

  return toss.str();
}






