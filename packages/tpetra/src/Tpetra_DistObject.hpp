// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
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

#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>

namespace Tpetra {

  //! Tpetra::DistObject: A class for constructing and using dense multi-vectors, vectors and matrices in parallel.

  /*! The DistObject is a base class for all Tpetra distributed global objects.  It provides the basic
      mechanisms and interface specifications for importing and exporting operations using Tpetra::Import and
      Tpetra::Export objects.
    
    <b> Distributed Global vs. Replicated Local.</b>
    
    <ul>
    <li> Distributed Global objects - In most instances, a distributed object will be partitioned
    across multiple memory images associated with multiple processors.  In this case, there is 
    a unique copy of each element and elements are spread across all images specified by 
    the Tpetra::Platform object.
    <li> Replicated Local Objects - Some algorithms use objects that are too small to
    be distributed across all processors, the Hessenberg matrix in a GMRES
    computation.  In other cases, such as with block iterative methods,  block dot product 
    functions produce small dense matrices that are required by all images.  
    Replicated local objects handle these types of situation.
    </ul>
  */

  template<typename Ordinal, typename Scalar>
  class DistObject : public Teuchos::Object {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! constructor
    DistObject(const Map<Ordinal>& map, Teuchos::RCP<const Teuchos::Comm<Ordinal> > comm);

    //! constructor, taking label
    DistObject(const Map<Ordinal>& map, 
           Teuchos::RCP<const Teuchos::Comm<Ordinal> > comm,
           const std::string & Label);

    //! copy constructor
    DistObject(const DistObject<Ordinal, Scalar>& source);

    //! destructor
    virtual ~DistObject();

    //@}

    //! @name Import/Export Methods
    //@{ 

    //! Import
    void doImport(const DistObject<Ordinal, Scalar> & source, 
                  const Import<Ordinal> & importer, CombineMode CM);

    //! Export
    void doExport(const DistObject<Ordinal, Scalar> & source, 
                  const Export<Ordinal> & exporter, CombineMode CM);

    //! Import (using an Exporter)
    void doImport(const DistObject<Ordinal, Scalar> & source,
                  const Export<Ordinal> & exporter, CombineMode CM);

    //! Export (using an Importer)
    void doExport(const DistObject<Ordinal, Scalar> & source,
                  const Import<Ordinal> & importer, CombineMode CM);

    //@}

    //! @name Attribute Accessor Methods
    //@{ 

    //! Accessor for whether or not this is a global object
    inline bool isDistributed() const;

    //! Access function for the Tpetra::Map this DistObject was constructed with.
    inline const Map<Ordinal> & getMap() const;

    //@}

    //! @name I/O methods
    //@{ 

    //! Print method.

    void print(std::ostream &os) const;

    //@} 

  protected:

    //! Perform transfer (redistribution) of data across memory images.
    virtual void doTransfer(const DistObject<Ordinal,Scalar> &source,
                            CombineMode CM,
                            Ordinal numSameIDs,
                            const Teuchos::ArrayView<const Ordinal> &permuteToLIDs,
                            const Teuchos::ArrayView<const Ordinal> &permuteFromLIDs,
                            const Teuchos::ArrayView<const Ordinal> &remoteLIDs,
                            const Teuchos::ArrayView<const Ordinal> &exportLIDs,
                            Distributor<Ordinal> &distor,
                            bool doReverse);

    // The following four methods must be implemented by the derived class

    //! Allows the source and target (\e this) objects to be compared for compatibility.
    /*! Return true if they are compatible, return false if they aren't. Also return the number of Scalar variables representing an entry. */ 
    virtual bool checkSizes(const DistObject<Ordinal, Scalar> & source, Ordinal &packetSize) = 0;

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
    virtual void copyAndPermute(const DistObject<Ordinal, Scalar> & source,
                                      Ordinal numSameIDs,
                                const Teuchos::ArrayView<const Ordinal> &permuteToLIDs,
                                const Teuchos::ArrayView<const Ordinal> &permuteFromLIDs) = 0;

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
    virtual void packAndPrepare(const DistObject<Ordinal,Scalar> & source,
                                const Teuchos::ArrayView<const Ordinal> &exportLIDs,
                                const Teuchos::ArrayView<Scalar> &exports,
                                Distributor<Ordinal> &distor) = 0;

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
    virtual void unpackAndCombine(const Teuchos::ArrayView<const Ordinal> &importLIDs,
                                  const Teuchos::ArrayView<const Scalar> &imports,
                                  Distributor<Ordinal> &distor,
                                  CombineMode CM) = 0;

  private:

    const Map<Ordinal> map_;
    Teuchos::RCP<const Teuchos::Comm<Ordinal> > Comm_;
    // buffers into which packed data is imported
    Teuchos::Array<Scalar> imports_;
    // buffers from which packed data is exported
    Teuchos::Array<Scalar> exports_;

  }; // class DistObject

  template <typename Ordinal, typename Scalar>
  DistObject<Ordinal,Scalar>::DistObject(const Map<Ordinal>& map, Teuchos::RCP<const Teuchos::Comm<Ordinal> > comm)
  : Teuchos::Object("Tpetra::DistObject")
  , map_(map)
  , Comm_(comm)
  {}

  template <typename Ordinal, typename Scalar>
  DistObject<Ordinal,Scalar>::DistObject(const Map<Ordinal>& map, 
      Teuchos::RCP<const Teuchos::Comm<Ordinal> > comm, const std::string & Label)
  : Teuchos::Object(Label.c_str())
  , map_(map)
  , Comm_(comm)
  {}

  template <typename Ordinal, typename Scalar>
  DistObject<Ordinal,Scalar>::DistObject(const DistObject<Ordinal, Scalar>& source)
  : Teuchos::Object(source.label())
  , map_(source.map_)
  , Comm_(source.Comm_)
  {}

  template <typename Ordinal, typename Scalar>
  DistObject<Ordinal,Scalar>::~DistObject() 
  {}

  template <typename Ordinal, typename Scalar>
  void DistObject<Ordinal,Scalar>::doImport(const DistObject<Ordinal,Scalar> & source, 
                                            const Import<Ordinal> & importer, CombineMode CM) 
  {
    TEST_FOR_EXCEPTION( getMap() != importer.getTargetMap(), std::runtime_error, "Target Maps don't match.");
    TEST_FOR_EXCEPTION( source.getMap() != importer.getSourceMap(), std::runtime_error, "Source Maps don't match.");
    Ordinal numSameIDs = importer.getNumSameIDs();
    const Teuchos::ArrayView<const Ordinal> exportLIDs      = importer.getExportLIDs();
    const Teuchos::ArrayView<const Ordinal> remoteLIDs      = importer.getRemoteLIDs();
    const Teuchos::ArrayView<const Ordinal> permuteToLIDs   = importer.getPermuteToLIDs();
    const Teuchos::ArrayView<const Ordinal> permuteFromLIDs = importer.getPermuteFromLIDs();
    this->doTransfer(source, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
                     importer.getDistributor(), false);
  }

  template <typename Ordinal, typename Scalar>
  void DistObject<Ordinal,Scalar>::doExport(const DistObject<Ordinal,Scalar> & source, 
                                            const Export<Ordinal> & exporter, CombineMode CM) 
  {
    TEST_FOR_EXCEPTION( getMap() != exporter.getTargetMap(), std::runtime_error, "Target Maps don't match.");
    TEST_FOR_EXCEPTION( source.getMap() != exporter.getSourceMap(), std::runtime_error, "Source Maps don't match.");
    Ordinal numSameIDs = exporter.getNumSameIDs();
    Teuchos::ArrayView<const Ordinal> exportLIDs      = exporter.getExportLIDs();
    Teuchos::ArrayView<const Ordinal> remoteLIDs      = exporter.getRemoteLIDs();
    Teuchos::ArrayView<const Ordinal> permuteToLIDs   = exporter.getPermuteToLIDs();
    Teuchos::ArrayView<const Ordinal> permuteFromLIDs = exporter.getPermuteFromLIDs();
    doTransfer(source, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
               exporter.getDistributor(), false);
  }

  template <typename Ordinal, typename Scalar>
  void DistObject<Ordinal,Scalar>::doImport(const DistObject<Ordinal,Scalar> & source,
                                            const Export<Ordinal> & exporter, CombineMode CM)
  {
    TEST_FOR_EXCEPTION( getMap() != exporter.getTargetMap(), std::runtime_error, "Target Maps don't match.");
    TEST_FOR_EXCEPTION( source.getMap() != exporter.getSourceMap(), std::runtime_error, "Source Maps don't match.");
    Ordinal numSameIDs = exporter.getNumSameIDs();
    Teuchos::ArrayView<const Ordinal> exportLIDs      = exporter.getRemoteLIDs();
    Teuchos::ArrayView<const Ordinal> remoteLIDs      = exporter.getExportLIDs();
    Teuchos::ArrayView<const Ordinal> permuteToLIDs   = exporter.getPermuteFromLIDs();
    Teuchos::ArrayView<const Ordinal> permuteFromLIDs = exporter.getPermuteToLIDs();
    doTransfer(source, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
               exporter.getDistributor(), true);
  }

  template <typename Ordinal, typename Scalar>
  void DistObject<Ordinal,Scalar>::doExport(const DistObject<Ordinal, Scalar> & source,
                                            const Import<Ordinal> & importer, CombineMode CM)
  {
    TEST_FOR_EXCEPTION( getMap() != importer.getTargetMap(), std::runtime_error, "Target Maps don't match.");
    TEST_FOR_EXCEPTION( source.getMap() != importer.getSourceMap(), std::runtime_error, "Source Maps don't match.");
    Ordinal numSameIDs = importer.getNumSameIDs();
    Teuchos::ArrayView<const Ordinal> exportLIDs      = importer.getRemoteLIDs();
    Teuchos::ArrayView<const Ordinal> remoteLIDs      = importer.getExportLIDs();
    Teuchos::ArrayView<const Ordinal> permuteToLIDs   = importer.getPermuteFromLIDs();
    Teuchos::ArrayView<const Ordinal> permuteFromLIDs = importer.getPermuteToLIDs();
    doTransfer(source, CM, numSameIDs, permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
               importer.getDistributor(), true);
  }

  template <typename Ordinal, typename Scalar>
  bool DistObject<Ordinal,Scalar>::isDistributed() const 
  {
    return map_.isDistributed();
  }

  template <typename Ordinal, typename Scalar>
  const Map<Ordinal> & DistObject<Ordinal,Scalar>::getMap() const 
  {
    return map_;
  }

  template <typename Ordinal, typename Scalar>
  void DistObject<Ordinal,Scalar>::doTransfer(
      const DistObject<Ordinal, Scalar> & source,
      CombineMode CM,
      Ordinal numSameIDs, 
      const Teuchos::ArrayView<const Ordinal> &permuteToLIDs, 
      const Teuchos::ArrayView<const Ordinal> &permuteFromLIDs,
      const Teuchos::ArrayView<const Ordinal> &remoteLIDs,    
      const Teuchos::ArrayView<const Ordinal> &exportLIDs,
      Distributor<Ordinal> &distor, bool doReverse) 
  {
    Ordinal packetSize;
    TEST_FOR_EXCEPTION( checkSizes(source,packetSize) == false, std::runtime_error, 
        "Tpetra::DistObject::doTransfer(): checkSizes() indicates that DistOjbects are not size-compatible.");
    Ordinal sbufLen = exportLIDs.size()*packetSize;
    Ordinal rbufLen = remoteLIDs.size()*packetSize;
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

  template <typename Ordinal, typename Scalar>
  void DistObject<Ordinal,Scalar>::print(std::ostream &os) const
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
