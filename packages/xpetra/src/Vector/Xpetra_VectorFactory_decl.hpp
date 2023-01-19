// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef XPETRA_VECTORFACTORY_DECL_HPP
#define XPETRA_VECTORFACTORY_DECL_HPP

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Vector.hpp"

#  include "Xpetra_TpetraVector_decl.hpp"
#ifdef HAVE_XPETRA_EPETRA
#  include "Xpetra_EpetraVector.hpp"
#  include "Xpetra_EpetraIntVector.hpp"
#endif

#include "Xpetra_BlockedMap_decl.hpp"
#include "Xpetra_BlockedVector_decl.hpp"
#include "Xpetra_Exceptions.hpp"



namespace Xpetra {



  template <class Scalar        /* = Vector<>::scalar_type*/,
            class LocalOrdinal  /* = typename Vector<Scalar>::local_ordinal_type*/,
            class GlobalOrdinal /* = typename Vector<Scalar, LocalOrdinal>::local_ordinal_type*/,
            class Node          /* = typename Vector<Scalar, LocalOrdinal, GlobalOrdinal>::node_type*/>
  class VectorFactory
  {
    #undef XPETRA_VECTORFACTORY_SHORT
    #include "Xpetra_UseShortNames.hpp"

  private:

    //! Private constructor. This is a static class.
    VectorFactory() = default;

  public:

    //! Constructor specifying the number of non-zeros for all rows.
    static Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>> 
    Build(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>> &map, bool zeroOut=true)
    {
        XPETRA_MONITOR("VectorFactory::Build");

        RCP<const Xpetra::BlockedMap<LocalOrdinal,GlobalOrdinal,Node>>
            bmap = Teuchos::rcp_dynamic_cast<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>(map);

        if(!bmap.is_null())
        {
            return rcp(new Xpetra::BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bmap, zeroOut));
        }

        if(map->lib() == UseTpetra)
        {
            return rcp(new TpetraVector(map, zeroOut));
        }

        XPETRA_FACTORY_ERROR_IF_EPETRA(map->lib());
        XPETRA_FACTORY_END;
    }

  };  // class VectorFactory

#define XPETRA_VECTORFACTORY_SHORT



#if defined(HAVE_XPETRA_EPETRA)


// we need the Epetra specialization only if Epetra is enabled
#if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)

  // Specialization for Scalar=double, LO=GO=int and EpetraNode node
  // Used both for Epetra and Tpetra
  // For any other node definition the general default implementation is used which allows Tpetra only
  template <>
  class VectorFactory<double, int, int, EpetraNode> 
  {
    typedef double     Scalar;
    typedef int        LocalOrdinal;
    typedef int        GlobalOrdinal;
    typedef EpetraNode Node;

#undef XPETRA_VECTORFACTORY_SHORT
#include "Xpetra_UseShortNames.hpp"

  private:

    //! Private constructor. This is a static class.
    VectorFactory() = default;

  public:

    static RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>> 
    Build(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>>& map, 
          bool zeroOut=true);

  };
#endif  // #if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)



  // Specialization for Scalar=double, LO=int, GO=long long and EpetraNode
  // Used both for Epetra and Tpetra
  // For any other node definition the general default implementation is used which allows Tpetra only
#if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)

  template <>
  class VectorFactory<double, int, long long, EpetraNode> 
  {

    typedef double     Scalar;
    typedef int        LocalOrdinal;
    typedef long long  GlobalOrdinal;
    typedef EpetraNode Node;

#undef XPETRA_VECTORFACTORY_SHORT
#include "Xpetra_UseShortNames.hpp"

  private:

    //! Private constructor. This is a static class.
    VectorFactory() = default;

  public:

    static RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>>
    Build(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>>& map, 
          bool zeroOut=true);

  };
#endif  // #if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)
#define XPETRA_VECTORFACTORY_SHORT



// we need the Epetra specialization only if Epetra is enabled
#if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)

  // Specialization for Scalar=int, LO=GO=int and EpetraNode
  // Used both for Epetra and Tpetra
  // For any other node definition the general default implementation is used which allows Tpetra only
  template <>
  class VectorFactory<int, int, int, EpetraNode> 
  {

    typedef int        Scalar;
    typedef int        LocalOrdinal;
    typedef int        GlobalOrdinal;
    typedef EpetraNode Node;

#undef XPETRA_VECTORFACTORY_SHORT
#include "Xpetra_UseShortNames.hpp"

  private:

    //! Private constructor. This is a static class.
    VectorFactory() = default;

  public:

    static RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>> 
    Build(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>>& map, 
          bool zeroOut=true);

  };
#endif  // #if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)



// we need the Epetra specialization only if Epetra is enabled
#if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)

  // Specialization for Scalar=int, LO=int, GO=long long and Serial node
  // Used both for Epetra and Tpetra
  // For any other node definition the general default implementation is used which allows Tpetra only

  template <>
  class VectorFactory<int, int, long long, EpetraNode> 
  {
    typedef int        Scalar;
    typedef int        LocalOrdinal;
    typedef long long  GlobalOrdinal;
    typedef EpetraNode Node;

#undef XPETRA_VECTORFACTORY_SHORT
#include "Xpetra_UseShortNames.hpp"

  private:

    //! Private constructor. This is a static class.
    VectorFactory() = default;

  public:

    static RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>> 
    Build(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>>& map, 
          bool zeroOut=true);

  };
#endif    // !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)


#endif  // #if defined(HAVE_XPETRA_EPETRA)


}   // namespace Xpetra


#define XPETRA_VECTORFACTORY_SHORT
#endif    // XPETRA_VECTORFACTORY_DECL_HPP


