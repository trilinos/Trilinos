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
#include "Xpetra_MultiVectorFactory_decl.hpp"
#include "Xpetra_BlockedMultiVector.hpp"


namespace Xpetra {



// we need the Epetra specialization only if Epetra is enabled
#if defined(HAVE_XPETRA_EPETRA)

#if !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)


MultiVectorFactory<double, int, int, EpetraNode>::
MultiVectorFactory()
{
}


RCP<MultiVector<double, int, int, EpetraNode>>
MultiVectorFactory<double, int, int, EpetraNode>::
Build(const Teuchos::RCP<const Map<int, int, EpetraNode>>& map, size_t NumVectors, bool zeroOut)
{
    using BlockedMultiVector = Xpetra::BlockedMultiVector<double, int, int, EpetraNode>

      XPETRA_MONITOR("MultiVectorFactory::Build");

    RCP<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> bmap = Teuchos::rcp_dynamic_cast<const BlockedMap<int, int, EpetraNode>>(map);

    if(!bmap.is_null())
    {
        return rcp(new BlockedMultiVector(bmap, NumVectors, zeroOut));
    }

    if(map->lib() == UseTpetra)
    {
        return rcp(new TpetraMultiVector<double, int, int, EpetraNode>(map, NumVectors, zeroOut));
    }

    if(map->lib() == UseEpetra)
    {
        return rcp(new EpetraMultiVectorT<int, EpetraNode>(map, NumVectors, zeroOut));
    }

    XPETRA_FACTORY_END;
}


Teuchos::RCP<MultiVector<double, int, int, EpetraNode>>
MultiVectorFactory<double, int, int, EpetraNode>::
Build(const Teuchos::RCP<const Map<int, int, EpetraNode>>&              map,
      const Teuchos::ArrayView<const Teuchos::ArrayView<const double>>& ArrayOfPtrs,
      size_t                                                            NumVectors)
{
    XPETRA_MONITOR("MultiVectorFactory::Build");

    if(map->lib() == UseTpetra)
    {
        return rcp(new TpetraMultiVector<double, int, int, EpetraNode>(map, ArrayOfPtrs, NumVectors));
    }

    if(map->lib() == UseEpetra)
    {
        return rcp(new EpetraMultiVectorT<int, EpetraNode>(map, ArrayOfPtrs, NumVectors));
    }

    XPETRA_FACTORY_END;
}


Teuchos::RCP<MultiVector<double, int, int, EpetraNode>>
MultiVectorFactory<double, int, int, EpetraNode>::
Build(const Teuchos::RCP<const MultiVector<double, int, int, EpetraNode> > &source,
      Teuchos::DataAccess copyOrView)
{
    XPETRA_MONITOR("MultiVectorFactory::Build");

    if(source->getMap()->lib() == UseTpetra)
    {
      return rcp(new TpetraMultiVector<double, int, int, EpetraNode>(*source, copyOrView));
    }

    if(source->getMap()->lib() == UseEpetra)
    {
        return rcp(new EpetraMultiVectorT<int, EpetraNode>(*source, copyOrView));
    }

    XPETRA_FACTORY_END;
}


MultiVectorFactory<int, int, int, EpetraNode>::
MultiVectorFactory()
{
}


RCP<MultiVector<int, int, int, EpetraNode>>
MultiVectorFactory<int, int, int, EpetraNode>::
Build(const Teuchos::RCP<const Map<int, int, EpetraNode>>& map, size_t NumVectors, bool zeroOut)
{
    XPETRA_MONITOR("MultiVectorFactory::Build");

    RCP<const BlockedMap<int, int, EpetraNode>> bmap = Teuchos::rcp_dynamic_cast<const BlockedMap<int, int, EpetraNode>>(map);

    if(!bmap.is_null())
    {
        return rcp(new BlockedMultiVector<int, int, int, EpetraNode>(bmap, NumVectors, zeroOut));
    }

    if(map->lib() == UseTpetra)
    {
        return rcp(new TpetraMultiVector<int, int, int, EpetraNode>(map, NumVectors, zeroOut));
    }

    if(map->lib() == UseEpetra)
    {
        return rcp(new EpetraIntMultiVectorT<int, EpetraNode>(map, NumVectors, zeroOut));
    }

    XPETRA_FACTORY_END;
}


Teuchos::RCP<MultiVector<int, int, int, EpetraNode>>
MultiVectorFactory<int, int, int, EpetraNode>::
Build(const Teuchos::RCP<const Map<int, int, EpetraNode>>&           map,
      const Teuchos::ArrayView<const Teuchos::ArrayView<const int>>& ArrayOfPtrs,
      size_t                                                         NumVectors)
{
    XPETRA_MONITOR("MultiVectorFactory::Build");

    if(map->lib() == UseTpetra)
    {
        return rcp(new TpetraMultiVector<int, int, int, EpetraNode>(map, ArrayOfPtrs, NumVectors));
    }

    if(map->lib() == UseEpetra)
    {
        return rcp(new EpetraIntMultiVectorT<int, EpetraNode>(map, ArrayOfPtrs, NumVectors));
    }

    XPETRA_FACTORY_END;
}


Teuchos::RCP<MultiVector<int, int, int, EpetraNode>>
MultiVectorFactory<int, int, int, EpetraNode>::
Build(const Teuchos::RCP<const MultiVector<int, int, int, EpetraNode> > &source,
      Teuchos::DataAccess copyOrView)
{
    XPETRA_MONITOR("MultiVectorFactory::Build");

    if(source->getMap()->lib() == UseTpetra)
    {
      return rcp(new TpetraMultiVector<int, int, int, EpetraNode>(*source, copyOrView));
    }

    if(source->getMap()->lib() == UseEpetra)
    {
        return rcp(new EpetraIntMultiVectorT<int, EpetraNode>(*source, copyOrView));
    }

    XPETRA_FACTORY_END;
}


// we need the Epetra specialization only if Epetra is enabled
#if !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)



MultiVectorFactory<double, int, long long, EpetraNode>::MultiVectorFactory() {}


RCP<MultiVector<double, int, long long, EpetraNode>>
MultiVectorFactory<double, int, long long, EpetraNode>::
Build(const Teuchos::RCP<const Map<int, long long, EpetraNode>>& map,
      size_t                                                     NumVectors,
      bool                                                       zeroOut)
{
    XPETRA_MONITOR("MultiVectorFactory::Build");

    RCP<const BlockedMap<int, long long, EpetraNode>> bmap = Teuchos::rcp_dynamic_cast<const BlockedMap<int, long long, EpetraNode>>(map);

    if(!bmap.is_null())
    {
        return rcp(new BlockedMultiVector<double, int, long long, EpetraNode>(bmap, NumVectors, zeroOut));
    }

    if(map->lib() == UseTpetra)
    {
        return rcp(new TpetraMultiVector<double, int, long long, EpetraNode>(map, NumVectors, zeroOut));
    }

    if(map->lib() == UseEpetra)
    {
        return rcp(new EpetraMultiVectorT<long long, EpetraNode>(map, NumVectors, zeroOut));
    }

    XPETRA_FACTORY_END;
}


Teuchos::RCP<MultiVector<double, int, long long, EpetraNode>>
MultiVectorFactory<double, int, long long, EpetraNode>::
Build(const Teuchos::RCP<const Map<int, long long, EpetraNode>>&        map,
      const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar>>& ArrayOfPtrs,
      size_t                                                            NumVectors)
{
    XPETRA_MONITOR("MultiVectorFactory::Build");

    if(map->lib() == UseTpetra)
    {
        return rcp(new TpetraMultiVector<double, int, long long, EpetraNode>(map, ArrayOfPtrs, NumVectors));
    }

    if(map->lib() == UseEpetra)
    {
        return rcp(new EpetraMultiVectorT<long long, EpetraNode>(map, ArrayOfPtrs, NumVectors));
    }

    XPETRA_FACTORY_END;
}


Teuchos::RCP<MultiVector<double, int, long long, EpetraNode>>
MultiVectorFactory<double, int, long long, EpetraNode>::
Build(const Teuchos::RCP<const MultiVector<double, int, long long, EpetraNode> > &source,
      Teuchos::DataAccess copyOrView)
{
    XPETRA_MONITOR("MultiVectorFactory::Build");

    if(source->getMap()->lib() == UseTpetra)
    {
      return rcp(new TpetraMultiVector<double, int, long long, EpetraNode>(*source, copyOrView));
    }

    if(source->getMap()->lib() == UseEpetra)
    {
        return rcp(new EpetraMultiVectorT<long long, EpetraNode>(*source, copyOrView));
    }

    XPETRA_FACTORY_END;
}


MultiVectorFactory<int, int, long long, EpetraNode>::
MultiVectorFactory()
{
}


RCP<MultiVector<int, int, long long, EpetraNode>>
MultiVectorFactory<int, int, long long, EpetraNode>::
Build(const Teuchos::RCP<const Map<int, long long, EpetraNode>>& map,
      size_t                                                     NumVectors,
      bool                                                       zeroOut)
{
    XPETRA_MONITOR("MultiVectorFactory::Build");

    RCP<const BlockedMap<int, long long, EpetraNode>> bmap = Teuchos::rcp_dynamic_cast<const BlockedMap<int, long long, EpetraNode>>(map);

    if(!bmap.is_null())
    {
        return rcp(new BlockedMultiVector<int, int, long long, EpetraNode>(bmap, NumVectors, zeroOut));
    }

    if(map->lib() == UseTpetra)
    {
        return rcp(new TpetraMultiVector<int, int, long long, EpetraNode>(map, NumVectors, zeroOut));
    }

    if(map->lib() == UseEpetra)
    {
        return rcp(new EpetraIntMultiVectorT<long long, EpetraNode>(map, NumVectors, zeroOut));
    }

    XPETRA_FACTORY_END;
}


Teuchos::RCP<MultiVector<int, int, long long, EpetraNode>>
MultiVectorFactory<int, int, long long, EpetraNode>::
Build(const Teuchos::RCP<const Map<int, long long, Node>>&           map,
      const Teuchos::ArrayView<const Teuchos::ArrayView<const int>>& ArrayOfPtrs,
      size_t                                                         NumVectors)
{
    XPETRA_MONITOR("MultiVectorFactory::Build");

    if(map->lib() == UseTpetra)
    {
        return rcp(new TpetraMultiVector<int, int, long long, EpetraNode>(map, ArrayOfPtrs, NumVectors));
    }

    if(map->lib() == UseEpetra)
    {
        return rcp(new EpetraIntMultiVectorT<long long, EpetraNode>(map, ArrayOfPtrs, NumVectors));
    }

    XPETRA_FACTORY_END;
}


Teuchos::RCP<MultiVector<int, int, long long, EpetraNode>>
MultiVectorFactory<int, int, long long, EpetraNode>::
Build(const Teuchos::RCP<const MultiVector<int, int, long long, EpetraNode> > &source,
      Teuchos::DataAccess copyOrView)
{
    XPETRA_MONITOR("MultiVectorFactory::Build");

    if(source->getMap()->lib() == UseTpetra)
    {
      return rcp(new TpetraMultiVector<int, int, long long, EpetraNode>(*source, copyOrView));
    }

    if(source->getMap()->lib() == UseEpetra)
    {
        return rcp(new EpetraIntMultiVectorT<long long, EpetraNode>(*source, copyOrView));
    }

    XPETRA_FACTORY_END;
}


#endif      // END !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES)

#endif      // END !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES)

#endif      // END HAVE_XPETRA_EPETRA


}      // namespace Xpetra
