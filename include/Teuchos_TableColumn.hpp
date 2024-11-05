// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_TABLECOLUMN_H
#define TEUCHOS_TABLECOLUMN_H

/** \file Teuchos_TableColumn.hpp
    \brief A column of TableEntry objects
    */

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_TableEntry.hpp"
#include "Teuchos_Array.hpp"

#include <iomanip>

namespace Teuchos
{
  /**
   *
   * KL 30 Apr 2006 -- initial design.
   */
  class TableColumn
  {
  public:
    /** \brief Empty ctor */
    TableColumn() : data_() {;}

    /** \brief  Form a column of std::string entries */
    TableColumn(const Array<std::string>& vals);

    /** \brief  Form a column of double entries */
    TableColumn(const Array<double>& vals, int precision, const std::ios_base::fmtflags& flags);

    /** \brief  Form a column of compound entries written as "first(second)" */
    TableColumn(const Array<double>& first, const Array<double>& second,
                int precision, const std::ios_base::fmtflags& flags,
                bool spaceBeforeParentheses);

    /** */
    int numRows() const {return Teuchos::as<int>(data_.size());}

    /** */
    void addEntry(const RCP<TableEntry>& entry);

    /** */
    const RCP<TableEntry>& entry(int i) const {return data_[i];}

  private:
    Array<RCP<TableEntry> > data_;
  };


}
#endif
