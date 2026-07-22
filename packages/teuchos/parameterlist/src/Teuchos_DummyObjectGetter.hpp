// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Teuchos_DummyObjectGetter.hpp
*/

#ifndef TEUCHOS_DUMMYOBJECTGETTER_HPP
#define TEUCHOS_DUMMYOBJECTGETTER_HPP



namespace Teuchos{


/** \brief Class for retrieving a dummy object of type T.*/
template<class T>
class DummyObjectGetter{

public:

  /** \brief Retrieves a dummy object of type T. */
  static RCP<T> getDummyObject(){
    return rcp(new T);
  }

};


} // namespace Teuchos


#endif // TEUCHOS_DUMMYOBJECTGETTER_HPP
