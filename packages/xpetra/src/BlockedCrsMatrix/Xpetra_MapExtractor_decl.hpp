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
#ifndef XPETRA_MAPEXTRACTOR_DECL_HPP_
#define XPETRA_MAPEXTRACTOR_DECL_HPP_

#include <map>

#include <iostream>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Describable.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_Map.hpp>

#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MapUtils.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Vector.hpp>

namespace Xpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of BlockedMultiVector, needed to prevent circular inclusions
  template<class S, class LO, class GO, class N> class BlockedMultiVector;

  // forward declaration of BlockedMap, needed because some functions take them as parameters 
  // (This should go away when BlockedMap is converted to ETI)
  template<class LO, class GO, class N> class BlockedMap;
#endif



  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class MapExtractor : public Teuchos::Describable
  {

  public:
    typedef Scalar scalar_type;
    typedef LocalOrdinal local_ordinal_type;
    typedef GlobalOrdinal global_ordinal_type;
    typedef Node node_type;

  private:
#undef XPETRA_MAPEXTRACTOR_SHORT
#include "Xpetra_UseShortNames.hpp"

  public:

    //! MapExtractor basic constructor
    //!
    //! @param[in] fullmap Full map containing all GIDs throughout the full vector. This parameter is only important if bThyraMode == false (see below)
    //! @param[in] maps    Vector containing submaps. The set of all GIDs stored in the submaps should be the same than stored in fullmap, if bThyraMode == false. In Thyra mode, the submaps should contain consecutive GIDs starting with 0 in each submap.
    //! @param[in] bThyraMode Flag which allows to switch between generating a MapExtractor in Thyra mode or Xpetra mode
    //!
    //! In Thyra mode, fullmap is not important as a fullmap with unique blocked GIDs is automatically generated which map the GIDs of the submaps
    //! to uniquely defined GIDs in the fullmap. The user has to provide a fullmap in Thyra mode to specify the underlying linear algebra library
    //! (Epetra or Tpetra) and some other map information (e.g. indexBase). This could be fixed.
    //!
    //! In Xpetra mode, the fullmap has to be the same as the union of the GIDs stored in the submaps in maps. The intersection of the GIDs of the sub-
    //! maps in maps must be empty.
    MapExtractor(const RCP<const Map>& fullmap, const std::vector<RCP<const Map> >& maps, bool bThyraMode = false);


    //! Expert constructor for Thyra maps
    MapExtractor(const std::vector<RCP<const Map> >& maps, const std::vector<RCP<const Map> >& thyramaps);


    /*!
     * Constructor which accepts a const version
     * of a blocked map
     *
     * \param map BlockedMap defining the block structure of the multi vector
     */
    MapExtractor(const Teuchos::RCP< const Xpetra::BlockedMap<LocalOrdinal,GlobalOrdinal,Node> > &map);


    //! copy constructor
    MapExtractor(const MapExtractor& input);


    //! Destructor.
    virtual ~MapExtractor();


    /** \name Extract subblocks from full map */
    //@{
    void             ExtractVector(const                  Vector& full, size_t block,           Vector& partial) const;
    void             ExtractVector(const             MultiVector& full, size_t block,      MultiVector& partial) const;
    void             ExtractVector(RCP<const             Vector>& full, size_t block, RCP<     Vector>& partial) const;
    void             ExtractVector(RCP<                  Vector>& full, size_t block, RCP<     Vector>& partial) const;
    void             ExtractVector(RCP<const        MultiVector>& full, size_t block, RCP<MultiVector>& partial) const;
    void             ExtractVector(RCP<             MultiVector>& full, size_t block, RCP<MultiVector>& partial) const;

    RCP<     Vector> ExtractVector(RCP<const             Vector>& full, size_t block,             bool  bThyraMode = false) const;
    RCP<     Vector> ExtractVector(RCP<                  Vector>& full, size_t block,             bool  bThyraMode = false) const;
    RCP<MultiVector> ExtractVector(RCP<const        MultiVector>& full, size_t block,             bool  bThyraMode = false) const;
    RCP<MultiVector> ExtractVector(RCP<             MultiVector>& full, size_t block,             bool  bThyraMode = false) const;

    RCP<MultiVector> ExtractVector(RCP<const Xpetra::BlockedMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& full, size_t block, bool bThyraMode = false) const;
    RCP<MultiVector> ExtractVector(RCP<      Xpetra::BlockedMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& full, size_t block, bool bThyraMode = false) const;

    //@}

    /** \name Insert subblocks into full map */
    //@{
    void InsertVector(const          Vector& partial, size_t block,          Vector& full, bool bThyraMode = false) const;
    void InsertVector(const     MultiVector& partial, size_t block,     MultiVector& full, bool bThyraMode = false) const;
    void InsertVector(RCP<const      Vector> partial, size_t block, RCP<     Vector> full, bool bThyraMode = false) const;
    void InsertVector(RCP<           Vector> partial, size_t block, RCP<     Vector> full, bool bThyraMode = false) const;
    void InsertVector(RCP<const MultiVector> partial, size_t block, RCP<MultiVector> full, bool bThyraMode = false) const;
    void InsertVector(RCP<      MultiVector> partial, size_t block, RCP<MultiVector> full, bool bThyraMode = false) const;
    void InsertVector(RCP<const MultiVector> partial, size_t block, RCP<Xpetra::BlockedMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > full, bool bThyraMode = false) const;
    void InsertVector(RCP<      MultiVector> partial, size_t block, RCP<Xpetra::BlockedMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > full, bool bThyraMode = false) const;

    //@}

    RCP<     Vector> getVector(size_t i, bool bThyraMode = false, bool bZero= true) const;
    RCP<MultiVector> getVector(size_t i, size_t numvec,   bool bThyraMode = false, bool bZero = true) const;

    /// returns true, if sub maps are stored in Thyra-style numbering
    bool getThyraMode() const;

    /** \name Maps */
    //@{

    /// number of partial maps
    size_t NumMaps() const;

    /// get the map
    /// returns the sub map i from list of sub maps
    /// depending on the parameter bThyraMode the sub map that is returned uses Thyra or Xpetra numbering
    /// Note: Thyra-numbering is only allowed if the MapExtractor is also constructed using Thyra numbering
    const RCP<const Map> getMap(size_t i, bool bThyraMode = false) const;

    /// get the underlying BlockedMap object (as Map)
    const RCP<const Map> getMap() const;

    /// get the underlying BlockedMap object (as BlockedMap)
    const RCP<const Xpetra::BlockedMap<LocalOrdinal,GlobalOrdinal,Node>> getBlockedMap() const;

    /// the full map
    const RCP<const Map> getFullMap() const;

    /// returns map index in map extractor which contains GID
    size_t getMapIndexForGID(GlobalOrdinal gid) const;

    //@}

  private:
    Teuchos::RCP<const Xpetra::BlockedMap<LocalOrdinal,GlobalOrdinal,Node>> map_;         ///< blocked map containing the sub block maps (either thyra or xpetra mode)
  };


} // namespace xpetra

#define XPETRA_MAPEXTRACTOR_SHORT
#endif /* XPETRA_MAPEXTRACTOR_HPP_ */
