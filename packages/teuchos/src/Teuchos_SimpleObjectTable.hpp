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

#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

#ifndef TEUCHOS_SIMPLEOBJECTTABLE_HPP
#define TEUCHOS_SIMPLEOBJECTTABLE_HPP

/*! \file Teuchos_SimpleObjectTable.hpp
    \brief A simple object table class for Teuchos
*/

/*! \class Teuchos::SimpleObjectTable
    \brief This class provides a central place to store objects
*/

namespace Teuchos
{

template <class T>
class SimpleObjectTable
{
  public:

    SimpleObjectTable();

    ~SimpleObjectTable();

    int storeRCP(const RCP<T> & robj);

    int storeNew(T* obj, bool owned = true);

    template <class TOld>
    int storeCastedRCP(const RCP<TOld> & robj_old);

    int removeRCP(int &index);

    const RCP<T> getRCP(int index);

    void purge();

  private:

    Array< RCP<T> > tableOfObjects;

    Array< int > freedIndices;

};

template <class T>
SimpleObjectTable<T>::SimpleObjectTable()
{

}

template <class T>
SimpleObjectTable<T>::~SimpleObjectTable()
{
  purge();
}

template <class T>
int SimpleObjectTable<T>::storeRCP(const RCP<T> & robj)
{
  robj.assert_not_null();

  int index = -1;

  if (freedIndices.size() != 0) {
    index = freedIndices.back();
    freedIndices.pop_back();
    tableOfObjects[index] = robj;
  } else {
    tableOfObjects.push_back(robj);
    index = tableOfObjects.size() - 1;
  }

  return index;
}

template <class T>
int SimpleObjectTable<T>::storeNew(T* obj, bool owned)
{
  return storeRCP(rcp(obj, owned));
}

template <class T>
template <class TOld>
int SimpleObjectTable<T>::storeCastedRCP(const RCP<TOld> & robj_old)
{
  return storeRCP(rcp_dynamic_cast<T>(robj_old, true));
}

template <class T>
int SimpleObjectTable<T>::removeRCP(int &index)
{
  if (tableOfObjects[index] == Teuchos::null) {
    throw RangeError("Item has already been deleted from SimpleObjectTable.");
  }

  int cnt = tableOfObjects[index].strong_count();

  tableOfObjects[index] = Teuchos::null;
  freedIndices.push_back(index);
  index = -1;

  return (cnt-1);
}

template <class T>
const RCP<T> SimpleObjectTable<T>::getRCP(int index)
{
  if (tableOfObjects[index] == Teuchos::null) {
    throw RangeError("Item has already been deleted from SimpleObjectTable.");
  }

  return tableOfObjects[index];
}

template <class T>
void SimpleObjectTable<T>::purge()
{
  int ocnt = tableOfObjects.size();
  for (int i=0; i<ocnt; i++) {
    tableOfObjects[i] = Teuchos::null;
  }

  if (tableOfObjects.size() > 0)
    tableOfObjects.erase(tableOfObjects.begin(), tableOfObjects.end());
  if (freedIndices.size() > 0)
    freedIndices.erase(freedIndices.begin(), freedIndices.end());
}

} // end namespace Teuchos

#endif

