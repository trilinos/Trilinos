// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER

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
