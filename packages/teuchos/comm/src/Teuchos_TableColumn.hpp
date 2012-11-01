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

#ifndef TEUCHOS_TABLECOLUMN_H
#define TEUCHOS_TABLECOLUMN_H

/** \file Teuchos_TableColumn.hpp
    \brief A column of TableEntry objects
    */

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_TableEntry.hpp"
#include "Teuchos_Array.hpp"

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
    TableColumn(const Array<double>& vals, int precision); 

    /** \brief  Form a column of compound entries written as "first(second)" */
    TableColumn(const Array<double>& first, const Array<double>& second,
                int precision,
                bool spaceBeforeParentheses); 

    /** */
    int numRows() const {return data_.size();}

    /** */
    void addEntry(const RCP<TableEntry>& entry);

    /** */
    const RCP<TableEntry>& entry(int i) const {return data_[i];}

  private:
    Array<RCP<TableEntry> > data_;
  };


}
#endif
