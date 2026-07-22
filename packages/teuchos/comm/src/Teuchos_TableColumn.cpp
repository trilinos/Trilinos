// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_TableColumn.hpp"

using namespace Teuchos;


TableColumn::TableColumn(const Array<std::string>& vals)
  : data_(vals.size())
{
  for (Array<std::string>::size_type i=0; i<vals.size(); i++)
    {
      data_[i] = rcp(new StringEntry(vals[i]));
    }
}


TableColumn::TableColumn(const Array<double>& vals,
                         int precision,
                         const std::ios_base::fmtflags& flags)
  : data_(vals.size())
{
  for (Array<double>::size_type i=0; i<vals.size(); i++)
    {
      data_[i] = rcp(new DoubleEntry(vals[i], precision, flags));
    }
}


TableColumn::TableColumn(const Array<double>& first,
                         const Array<double>& second,
                         int precision,
                         const std::ios_base::fmtflags& flags,
                         bool spaceBeforeParentheses)
  : data_(first.size())
{
  std::ios_base::fmtflags fixedflags = flags;
  fixedflags &= ~std::cout.scientific;   // unset scientific
  fixedflags &= ~std::cout.fixed;        // unset fixed
  for (Array<double>::size_type i=0; i<first.size(); i++)
    {
      RCP<DoubleEntry> x1 = rcp(new DoubleEntry(first[i],  precision, flags));
      RCP<DoubleEntry> x2 = rcp(new DoubleEntry(second[i], precision, fixedflags));
      data_[i]
        = rcp(new CompoundEntryWithParentheses(x1, x2,
                                               spaceBeforeParentheses));
    }
}

void TableColumn::addEntry(const RCP<TableEntry>& entry_in)
{
  data_.append(entry_in);
}



