// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

