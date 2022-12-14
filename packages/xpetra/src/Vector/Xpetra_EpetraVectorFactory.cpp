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
#include "Xpetra_VectorFactory.hpp"
#include "Xpetra_Vector.hpp"
#include "Xpetra_BlockedVector.hpp"

namespace Xpetra {


#if defined(HAVE_XPETRA_EPETRA)


// we need the Epetra specialization only if Epetra is enabled
#if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)

RCP<Xpetra::Vector<double, int, int, EpetraNode>>
VectorFactory<double, int, int, EpetraNode>::
Build(const Teuchos::RCP<const Xpetra::Map<int, int, EpetraNode>>& map, bool zeroOut)
{
    XPETRA_MONITOR("VectorFactory::Build");

    RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> 
      bmap = Teuchos::rcp_dynamic_cast<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>(map);

    if(!bmap.is_null())
    {
        return rcp(new Xpetra::BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bmap, zeroOut));
    }


    if(map->lib() == UseTpetra)
    {
        return rcp(new TpetraVector(map, zeroOut));
    }

    if(map->lib() == UseEpetra)
    {
        return rcp(new EpetraVectorT<GlobalOrdinal, EpetraNode>(map, zeroOut));
    }

    XPETRA_FACTORY_END;
}

#endif      // #if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)



#if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)

RCP<Xpetra::Vector<double, int, long long, EpetraNode>>
VectorFactory<double, int, long long, EpetraNode>::
Build(const Teuchos::RCP<const Xpetra::Map<int, long long, EpetraNode>>& map, bool zeroOut)
{
    XPETRA_MONITOR("VectorFactory::Build");

    RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> bmap =
      Teuchos::rcp_dynamic_cast<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>(map);
    if(!bmap.is_null())
    {
        return rcp(new Xpetra::BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bmap, zeroOut));
    }

    if(map->lib() == UseTpetra)
    {
        return rcp(new TpetraVector(map, zeroOut));
    }

    if(map->lib() == UseEpetra)
    {
        return rcp(new EpetraVectorT<GlobalOrdinal, Node>(map, zeroOut));
    }

    XPETRA_FACTORY_END;
}

#endif      // #if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)



// we need the Epetra specialization only if Epetra is enabled
#if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)

RCP<Xpetra::Vector<int, int, int, EpetraNode>>
VectorFactory<int, int, int, EpetraNode>::
Build(const Teuchos::RCP<const Xpetra::Map<int, int, EpetraNode>>& map, bool zeroOut)
{
    XPETRA_MONITOR("VectorFactory::Build");

    RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> bmap =
      Teuchos::rcp_dynamic_cast<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>(map);
    if(!bmap.is_null())
    {
        return rcp(new Xpetra::BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bmap, zeroOut));
    }

    if(map->lib() == UseTpetra)
    {
        return rcp(new TpetraVector(map, zeroOut));
    }

    if(map->lib() == UseEpetra)
    {
        return rcp(new EpetraIntVectorT<GlobalOrdinal, Node>(map, zeroOut));
    }

    XPETRA_FACTORY_END;
}

#endif      // #if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)



#if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)

RCP<Xpetra::Vector<int, int, long long, EpetraNode>>
VectorFactory<int, int, long long, EpetraNode>::
Build(const Teuchos::RCP<const Xpetra::Map<int, long long, EpetraNode>>& map, bool zeroOut)
{
    XPETRA_MONITOR("VectorFactory::Build");

    RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> bmap =
      Teuchos::rcp_dynamic_cast<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>(map);

    if(!bmap.is_null())
    {
        return rcp(new Xpetra::BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bmap, zeroOut));
    }

    if(map->lib() == UseTpetra)
    {
        return rcp(new TpetraVector(map, zeroOut));
    }

    if(map->lib() == UseEpetra)
    {
        return rcp(new EpetraIntVectorT<GlobalOrdinal, Node>(map, zeroOut));
    }

    XPETRA_FACTORY_END;
}

#endif      // #if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)



#endif // #if defined(HAVE_XPETRA_EPETRA)


}      // namespace Xpetra
