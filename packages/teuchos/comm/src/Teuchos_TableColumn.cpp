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



