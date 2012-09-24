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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * Xpetra_MapExtractor.hpp
 *
 *  Created on: 08.08.2011
 *      Author: tobias
 */

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_MAPEXTRACTOR_HPP_
#define XPETRA_MAPEXTRACTOR_HPP_

#include <map>

#include <iostream>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Describable.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_Map.hpp>
#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraMap.hpp>
#include <Xpetra_EpetraImport.hpp>
#include <Xpetra_EpetraVector.hpp>
#include <Xpetra_EpetraMultiVector.hpp>
#endif

#ifdef HAVE_XPETRA_TPETRA
#include <Xpetra_TpetraMap.hpp>
#include <Xpetra_TpetraImport.hpp>
#include <Xpetra_TpetraVector.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#endif

#include <Xpetra_ImportFactory.hpp>


namespace Xpetra
{
  template <class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class MapExtractor : public Teuchos::Describable
  {
    typedef Xpetra::Map<LocalOrdinal,GlobalOrdinal, Node> MapClass;
    typedef Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> VectorClass;
    typedef Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MultiVectorClass;
    typedef Xpetra::Import<LocalOrdinal,GlobalOrdinal,Node> ImportClass;

  public:

    /** \name Extract subblocks from full map */
    //@{
    virtual void ExtractVector(Teuchos::RCP<const VectorClass>& full, size_t block, Teuchos::RCP<VectorClass>& partial) const = 0;
    virtual void ExtractVector(Teuchos::RCP<      VectorClass>& full, size_t block, Teuchos::RCP<VectorClass>& partial) const = 0;
    virtual Teuchos::RCP<VectorClass> ExtractVector(Teuchos::RCP<const VectorClass>& full, size_t block) const = 0;
    virtual Teuchos::RCP<VectorClass> ExtractVector(Teuchos::RCP<      VectorClass>& full, size_t block) const = 0;

    virtual void ExtractVector(Teuchos::RCP<const MultiVectorClass>& full, size_t block, Teuchos::RCP<MultiVectorClass>& partial) const = 0;
    virtual void ExtractVector(Teuchos::RCP<      MultiVectorClass>& full, size_t block, Teuchos::RCP<MultiVectorClass>& partial) const = 0;
    virtual Teuchos::RCP<MultiVectorClass> ExtractVector(Teuchos::RCP<const MultiVectorClass>& full, size_t block) const = 0;
    virtual Teuchos::RCP<MultiVectorClass> ExtractVector(Teuchos::RCP<      MultiVectorClass>& full, size_t block) const = 0;
    //@}

    //virtual Teuchos::RCP<Xpetra::Vector<LocalOrdinal,GlobalOrdinal,Node> >ExtractVector(Teuchos::RCP<const Xpetra::Vector<LocalOrdinal,GlobalOrdinal,Node> > full, int block) const {};

    /** \name Insert subblocks into full map */
    //@{
    virtual void InsertVector(Teuchos::RCP<const VectorClass>& partial, size_t block, Teuchos::RCP<VectorClass>& full) const = 0;
    virtual void InsertVector(Teuchos::RCP<      VectorClass>& partial, size_t block, Teuchos::RCP<VectorClass>& full) const = 0;

    virtual void InsertVector(Teuchos::RCP<const MultiVectorClass>& partial, size_t block, Teuchos::RCP<MultiVectorClass>& full) const = 0;
    virtual void InsertVector(Teuchos::RCP<      MultiVectorClass>& partial, size_t block, Teuchos::RCP<MultiVectorClass>& full) const = 0;

    //@}

    virtual Teuchos::RCP<VectorClass> getVector(size_t i) const = 0;
    virtual Teuchos::RCP<MultiVectorClass> getVector(size_t i, size_t numvec) const = 0;

    /** \name Maps */
    //@{

    /// number of partial maps
    virtual size_t NumMaps() const = 0;

    /// get the map
    virtual const Teuchos::RCP<const MapClass> getMap(size_t i) const = 0;

    /// the full map
    virtual const Teuchos::RCP<const MapClass> getFullMap() const = 0;

    /// returns map index in map extractor which contains GID or -1 otherwise
    virtual size_t getMapIndexForGID(GlobalOrdinal gid) const = 0;

    //@}

  private:
    virtual bool CheckConsistency() const = 0;
  };
}


#endif /* XPETRA_MAPEXTRACTOR_HPP_ */
