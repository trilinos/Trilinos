// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
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
// ***********************************************************************
// @HEADER

#ifndef TPETRA_DISTOBJECT_HPP
#define TPETRA_DISTOBJECT_HPP

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CombineMode.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Distributor.hpp"

#include <Teuchos_Describable.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace Tpetra {

  //! A base class for distributed objects that support import and export operations.

  /*! The DistObject is a base class for all Tpetra distributed global objects.  It provides the basic
      mechanisms and interface specifications for importing and exporting operations using Tpetra::Import and
      Tpetra::Export objects.
    
    <b> Distributed Global vs. Replicated Local.</b>
    
    <ul>
    <li> Distributed Global objects - In most instances, a distributed object will be partitioned
    across multiple memory images associated with multiple processors.  In this case, there is 
    a unique copy of each element and elements are spread across all images specified by 
    the Teuchos::Comm object.
    <li> Replicated Local Objects - Some algorithms use objects that are too small to
    be distributed across all processors, the Hessenberg matrix in a GMRES
    computation.  In other cases, such as with block iterative methods,  block dot product 
    functions produce small dense matrices that are required by all images.  
    Replicated local objects handle these types of situation.
    </ul>
  */

  template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal>
  class DistObject : public Teuchos::Describable {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! constructor
    DistObject(const Map<LocalOrdinal,GlobalOrdinal> &map);

    //! copy constructor
    DistObject(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> &source);

    //! destructor
    virtual ~DistObject();

    //@}

    //! @name Import/Export Methods
    //@{ 

    //! Import
    void doImport(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> &source, 
                  const Import<LocalOrdinal,GlobalOrdinal> &importer, CombineMode CM);

    //! Export
    void doExport(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> &dest, 
                  const Export<LocalOrdinal,GlobalOrdinal> &exporter, CombineMode CM);

    //! Import (using an Exporter)
    void doImport(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> &source,
                  const Export<LocalOrdinal,GlobalOrdinal>& exporter, CombineMode CM);

    //! Export (using an Importer)
    void doExport(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> &dest,
                  const Import<LocalOrdinal,GlobalOrdinal>& importer, CombineMode CM);

    //@}

    //! @name Attribute Accessor Methods
    //@{ 

    //! Accessor for whether or not this is a global object
    inline bool isDistributed() const;

    //! Access function for the Tpetra::Map this DistObject was constructed with.
    inline const Map<LocalOrdinal,GlobalOrdinal> & getMap() const;

    //@}

    //! @name I/O methods
    //@{ 

    //! Print method.

    void print(std::ostream &os) const;

    //@} 

  protected:

    //! Perform transfer (redistribution) of data across memory images.
    virtual void doTransfer(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> &source,
                            CombineMode CM,
                            Teuchos_Ordinal numSameIDs,
                            const Teuchos::ArrayView<const LocalOrdinal> &permuteToLIDs,
                            const Teuchos::ArrayView<const LocalOrdinal> &permuteFromLIDs,
                            const Teuchos::ArrayView<const LocalOrdinal> &remoteLIDs,
                            const Teuchos::ArrayView<const LocalOrdinal> &exportLIDs,
                            Distributor &distor,
                            bool doReverse);

    // The following four methods must be implemented by the derived class

    //! Allows the source and target (\e this) objects to be compared for compatibility.
    /*! Return true if they are compatible, return false if they aren't. Also return the number of Scalar variables representing an entry. */ 
    virtual bool checkSizes(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> & source, Teuchos_Ordinal &packetSize) = 0;

    //! Perform copies and permutations that are local to this image.
    /*!
      \param source In
             On entry, the DistObject that we are importing from.
      \param numSameIDs In
             On entry, the number of elements that are the same on the source and dest objects.
         (i.e. The element is owned by the same image in both source and dest, 
         and no permutation occurs.)
      \param numPermuteIDs In
             On entry, the number of elements that are locally permuted between source and dest objects.
      \param permuteToLIDs In
             On entry, contains a list of the elements that are permuted. (Listed by their LID in the
         destination DistObject.)
      \param permuteFromLIDs In
             On entry, contains a list of the elements that are permuted. (Listed by their LID in the
         source DistObject.)
    */
    virtual void copyAndPermute(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> & source,
                                      Teuchos_Ordinal numSameIDs,
                                const Teuchos::ArrayView<const LocalOrdinal> &permuteToLIDs,
                                const Teuchos::ArrayView<const LocalOrdinal> &permuteFromLIDs) = 0;

    //! Perform any packing or preparation required for communication.
    /*!
      \param source In
             On entry, the DistObject that we are importing from.
      \param exportLIDs In
             On entry, a list of the entries we will be sending to other images.
             (Listed by their LID in the source DistObject.)
      \param exports Out
             On exit, buffer for data we will be sending out.
      \param distor In
             On entry, contains the Distributor object we are using.         
    */
    virtual void packAndPrepare(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> & source,
                                const Teuchos::ArrayView<const LocalOrdinal> &exportLIDs,
                                const Teuchos::ArrayView<Scalar> &exports,
                                Distributor &distor) = 0;

    //! Perform any unpacking and combining after communication.
    /*!
      \param importLIDs In
             On entry, a list of the entries we received from other images.
             (Listed by their LID in the target DistObject.)
      \param imports In
             Buffer containing data we received.
      \param distor In
             The Distributor object we are using.
      \param CM In
             The Tpetra::CombineMode to use when combining the imported entries with existing entries.
    */
    virtual void unpackAndCombine(const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
                                  const Teuchos::ArrayView<const Scalar> &imports,
                                  Distributor &distor,
                                  CombineMode CM) = 0;

  private:

    const Map<LocalOrdinal,GlobalOrdinal> map_;
    // buffers into which packed data is imported
    Teuchos::Array<Scalar> imports_;
    // buffers from which packed data is exported
    Teuchos::Array<Scalar> exports_;

  }; // class DistObject

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  DistObject<Scalar,LocalOrdinal,GlobalOrdinal>::DistObject(const Map<LocalOrdinal,GlobalOrdinal>& map, Teuchos::RCP<const Teuchos::Comm<int> > comm)
  : map_(map)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  DistObject<Scalar,LocalOrdinal,GlobalOrdinal>::DistObject(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal>& source)
  : map_(source.map_)
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  DistObject<Scalar,LocalOrdinal,GlobalOrdinal>::~DistObject() 
  {}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void DistObject<Scalar,LocalOrdinal,GlobalOrdinal>::doImport(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> & A, 
                                            const Import<LocalOrdinal,GlobalOrdinal> & importer, CombineMode CM) 
  {
    TEST_FOR_EXCEPTION( getMap() != importer.getTargetMap(), std::runtime_error, "Target Maps don't match.");
    TEST_FOR_EXCEPTION( A.getMap() != importer.getSourceMap(), std::runtime_error, "Source Maps don't match.");
    Teuchos_Ordinal numSameIDs = importer.getNumSameIDs();
    const Teuchos::ArrayView<const LocalOrdinal> exportLIDs      = importer.getExportLIDs();
    const Teuchos::ArrayView<const LocalOrdinal> remoteLIDs      = importer.getRemoteLIDs();
    const Teuchos::ArrayView<const LocalOrdinal> permuteToLIDs   = importer.getPermuteToLIDs();
    const Teuchos::ArrayView<const LocalOrdinal> permuteFromLIDs = importer.getPermuteFromLIDs();
    this->doTransfer(A, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
                     importer.getDistributor(), false);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void DistObject<Scalar,LocalOrdinal,GlobalOrdinal>::doExport(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> & A, 
                                            const Export<LocalOrdinal,GlobalOrdinal> & exporter, CombineMode CM) 
  {
    TEST_FOR_EXCEPTION( getMap() != exporter.getTargetMap(), std::runtime_error, "Target Maps don't match.");
    TEST_FOR_EXCEPTION( A.getMap() != exporter.getSourceMap(), std::runtime_error, "Source Maps don't match.");
    Teuchos_Ordinal numSameIDs = exporter.getNumSameIDs();
    Teuchos::ArrayView<const LocalOrdinal> exportLIDs      = exporter.getExportLIDs();
    Teuchos::ArrayView<const LocalOrdinal> remoteLIDs      = exporter.getRemoteLIDs();
    Teuchos::ArrayView<const LocalOrdinal> permuteToLIDs   = exporter.getPermuteToLIDs();
    Teuchos::ArrayView<const LocalOrdinal> permuteFromLIDs = exporter.getPermuteFromLIDs();
    doTransfer(A, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
               exporter.getDistributor(), false);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void DistObject<Scalar,LocalOrdinal,GlobalOrdinal>::doImport(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> & A,
                                            const Export<LocalOrdinal,GlobalOrdinal> & exporter, CombineMode CM)
  {
    TEST_FOR_EXCEPTION( getMap() != exporter.getSourceMap(), std::runtime_error, "Target Maps don't match.");
    TEST_FOR_EXCEPTION( A.getMap() != exporter.getTargetMap(), std::runtime_error, "Source Maps don't match.");
    Teuchos_Ordinal numSameIDs = exporter.getNumSameIDs();
    Teuchos::ArrayView<const LocalOrdinal> exportLIDs      = exporter.getRemoteLIDs();
    Teuchos::ArrayView<const LocalOrdinal> remoteLIDs      = exporter.getExportLIDs();
    Teuchos::ArrayView<const LocalOrdinal> permuteToLIDs   = exporter.getPermuteFromLIDs();
    Teuchos::ArrayView<const LocalOrdinal> permuteFromLIDs = exporter.getPermuteToLIDs();
    doTransfer(A, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
               exporter.getDistributor(), true);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void DistObject<Scalar,LocalOrdinal,GlobalOrdinal>::doExport(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> & A,
                                            const Import<LocalOrdinal,GlobalOrdinal> & importer, CombineMode CM)
  {
    TEST_FOR_EXCEPTION( getMap() != importer.getSourceMap(), std::runtime_error, "Target Maps don't match.");
    TEST_FOR_EXCEPTION( A.getMap() != importer.getTargetMap(), std::runtime_error, "Source Maps don't match.");
    Teuchos_Ordinal numSameIDs = importer.getNumSameIDs();
    Teuchos::ArrayView<const LocalOrdinal> exportLIDs      = importer.getRemoteLIDs();
    Teuchos::ArrayView<const LocalOrdinal> remoteLIDs      = importer.getExportLIDs();
    Teuchos::ArrayView<const LocalOrdinal> permuteToLIDs   = importer.getPermuteFromLIDs();
    Teuchos::ArrayView<const LocalOrdinal> permuteFromLIDs = importer.getPermuteToLIDs();
    doTransfer(A, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
               importer.getDistributor(), true);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  bool DistObject<Scalar,LocalOrdinal,GlobalOrdinal>::isDistributed() const 
  {
    return map_.isDistributed();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  const Map<LocalOrdinal,GlobalOrdinal> & DistObject<Scalar,LocalOrdinal,GlobalOrdinal>::getMap() const 
  {
    return map_;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void DistObject<Scalar,LocalOrdinal,GlobalOrdinal>::doTransfer(
      const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> & source,
      CombineMode CM,
      Teuchos_Ordinal numSameIDs, 
      const Teuchos::ArrayView<const LocalOrdinal> &permuteToLIDs, 
      const Teuchos::ArrayView<const LocalOrdinal> &permuteFromLIDs,
      const Teuchos::ArrayView<const LocalOrdinal> &remoteLIDs,    
      const Teuchos::ArrayView<const LocalOrdinal> &exportLIDs,
      Distributor &distor, bool doReverse) 
  {
    Teuchos_Ordinal packetSize;
    TEST_FOR_EXCEPTION( checkSizes(source,packetSize) == false, std::runtime_error, 
        "Tpetra::DistObject::doTransfer(): checkSizes() indicates that DistOjbects are not size-compatible.");
    Teuchos_Ordinal sbufLen = exportLIDs.size()*packetSize;
    Teuchos_Ordinal rbufLen = remoteLIDs.size()*packetSize;
    exports_.resize(sbufLen);
    imports_.resize(rbufLen);
    if (numSameIDs + permuteToLIDs.size()) {
      copyAndPermute(source,numSameIDs,permuteToLIDs,permuteFromLIDs);
    }
    packAndPrepare(source,exportLIDs,exports_(),distor);
    if ((isDistributed() && doReverse) || (source.isDistributed() && !doReverse)) 
    {
      // call one of the doPostsAndWaits functions
      if (doReverse) {
        distor.doReversePostsAndWaits(exports_().getConst(),packetSize,imports_());
      }
      else {
        distor.doPostsAndWaits(exports_().getConst(),packetSize,imports_());
      }
      unpackAndCombine(remoteLIDs,imports_(),distor,CM);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void DistObject<Scalar,LocalOrdinal,GlobalOrdinal>::print(std::ostream &os) const
  {
    using std::endl;
    os << "Tpetra::DistObject" << endl
       << " export buffer size: " << exports_.size() << endl
       << " import buffer size: " << imports_.size() << endl
       << "Map:" << endl
       << map_;
  }

} // namespace Tpetra

#endif /* TPETRA_DISTOBJECT_HPP */
