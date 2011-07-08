//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifndef TPETRAEXT_TYPE_STACK_HPP
#define TPETRAEXT_TYPE_STACK_HPP

#include "Tpetra_ConfigDefs.hpp"
#include <Teuchos_ParameterList.hpp>

//! Macro to create a Tpetra::Ext::TypeStack of height 1
#define TPETRAEXT_TYPESTACK1(lbl,T1)              \
        typedef Tpetra::Ext::TypeStackBottom<T1>  \
                                             lbl;

//! Macro to create a Tpetra::Ext::TypeStack of height 2
#define TPETRAEXT_TYPESTACK2(lbl,T1,T2)           \
        typedef Tpetra::Ext::TypeStack<T1,        \
                                       T2>        \
                                       lbl;

//! Macro to create a Tpetra::Ext::TypeStack of height 3
#define TPETRAEXT_TYPESTACK3(lbl,T1,T2,T3)        \
        typedef Tpetra::Ext::TypeStack<T1,        \
                Tpetra::Ext::TypeStack<T2,        \
                                       T3> >      \
                                       lbl;

//! Macro to create a Tpetra::Ext::TypeStack of height 4
#define TPETRAEXT_TYPESTACK4(lbl,T1,T2,T3,T4)     \
        typedef Tpetra::Ext::TypeStack<T1,        \
                Tpetra::Ext::TypeStack<T2,        \
                Tpetra::Ext::TypeStack<T3,        \
                                       T4> > >    \
                                       lbl;

namespace Tpetra {
  namespace Ext {

    //! Implementation of a Tpetra::Ext::TypeStack, supporting the last entry.
    template <class T>
    struct TypeStackBottom {
      typedef T                  type;
      typedef TypeStackBottom<T> next; // infinite loop; better be watching for bottom
      enum { bottom = true };
      enum { height = 1 };
    };

    //! Implementation of a Tpetra::Ext::TypeStack, supporting the next to last entry.
    template <class T, class S>
    struct TypeStack {
      typedef T                  type;
      typedef TypeStackBottom<S> next;
      enum { bottom = false };
      enum { height = 2 };
    };

    //! Generic implementation of a Tpetra::Ext::TypeStack. This is the model that should be programmed to.
    template <class T, class S, class SS>
    struct TypeStack<T, TypeStack<S,SS> > {
      typedef T               type;
      typedef TypeStack<S,SS> next;
      enum { bottom = false };
      enum { height = 1 + next::height };
    };

    template <class TS,class Init>
    RCP<Teuchos::ParameterList>
    initStackDB(Teuchos::ParameterList &pl, Init &init)
    {
      using Teuchos::ParameterList;
      typedef typename TS::type T;
      RCP<ParameterList> db = init.template initDB<T>(pl);
      if (! TS::bottom) {
        RCP<ParameterList> subdb = initStackDB<typename TS::next>(pl.sublist("child"),init);
        db->set("child", *subdb);
      }
      return db;
    }

  } // namespace Tpetra::Ext
} // namespace Tpetra


/** 
  \example TypeStackTest.cpp
  An example for using the Tpetra::Ext::TypeStack class and its associated macros.
 */


#endif // TPETRAEXT_TYPE_STACK_HPP
