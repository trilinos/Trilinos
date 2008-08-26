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
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_RCP.hpp>
#include "Tpetra_Map.hpp"
#include "Tpetra_CombineMode.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
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

    DistObject error codes (positive for non-fatal, negative for fatal):
    <ol>
    <li> -1  importer's target Map does not match my Map
    <li> -2  importer's source Map does not match source's Map
    <li> -99 Internal DistObject error. Contact developer.
    </ol>
    
  */

  template<typename Ordinal, typename Scalar>
  class DistObject: public Teuchos::Object {

  public:

    //@{ \name Constructor/Destructor Methods

    //! constructor
    DistObject(const Map<Ordinal>& map, Teuchos::RCP< Teuchos::Comm<Ordinal> > comm);

    //! constructor, taking label
    DistObject(const Map<Ordinal>& map, 
           Teuchos::RCP< Teuchos::Comm<Ordinal> > comm,
           const std::string & Label);

    //! copy constructor
    DistObject(const DistObject<Ordinal, Scalar>& source);

    //! destructor
    virtual ~DistObject();

    //@}

    //@{ \name Import/Export Methods

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

    //@{ \name I/O methods

    //! print method
    virtual void print(std::ostream& os) const;

    //@}

    //@{ \name Attribute Accessor Methods

    //! Accessor for whether or not this is a global object
    bool isDistributed() const;

    //! Access function for the Tpetra::Map this DistObject was constructed with.
    const Map<Ordinal> & getMap() const;

    //@}

  protected:
 
    //! Perform actual transfer (redistribution) of data across memory images.
    virtual void doTransfer(const DistObject<Ordinal, Scalar> & source,
                 CombineMode CM,
                 Ordinal numSameIDs,
                 Ordinal numPermuteIDs,
                 Ordinal numRemoteIDs,
                 Ordinal numExportIDs,
                 const std::vector<Ordinal> & permuteToLIDs,
                 const std::vector<Ordinal> & permuteFromLIDs,
                 const std::vector<Ordinal> & remoteLIDs,
                 const std::vector<Ordinal> & exportLIDs,
                 std::vector<Scalar>& exports,
                 std::vector<Scalar>& imports,
                 Distributor<Ordinal> &distor,
                 bool doReverse);

    // The following four methods must be implemented by the derived class

    //! Allows the source and target (\e this) objects to be compared for compatibility.
    /*! Return true if they are compatible, return false if they aren't. */ 
    virtual bool checkSizes(const DistObject<Ordinal, Scalar> & source) = 0;

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
                   Ordinal numPermuteIDs,
                   const std::vector<Ordinal> & permuteToLIDs,
                   const std::vector<Ordinal> & permuteFromLIDs) = 0;

    //! Perform any packing or preparation required for communication.
    /*!
      \param source In
             On entry, the DistObject that we are importing from.
      \param numExportIDs In
             On entry, the number of elements we will be sending to other images.
      \param exportLIDs In
             On entry, a list of the elements we will be sending to other images.
         (Listed by their LID in the source DistObject.)
      \param exports Out
             On exit, buffer for data we will be sending out.
      \param packetSize Out
             On exit, will contain the number of Scalar variables used to pack
         a single element.
      \param distor In
             On entry, contains the Distributor object we are using.         
    */
    virtual void packAndPrepare(const DistObject<Ordinal, Scalar> & source,
                   Ordinal numExportIDs,
                   const std::vector<Ordinal> & exportLIDs,
                   std::vector<Scalar>& exports,
                   Ordinal & packetSize,
                   Distributor<Ordinal> &distor) = 0;
  
    //! Perform any unpacking and combining after communication.
    /*!
      \param numImportIDs In
             The number of elements we received from other images.
      \param importLIDs In
             On entry, a list of the elements we received from other images.
         (Listed by their LID in the target DistObject.)
      \param imports In
             Buffer for data we received.
      \param distor In
             The Distributor object we are using.
      \param CM In
             The Tpetra::CombineMode to use when combining the imported
         entries with existing entries.
    */
    virtual void unpackAndCombine(Ordinal numImportIDs,
                   const std::vector<Ordinal> & importLIDs,
                   const std::vector<Scalar> & imports,
                   Distributor<Ordinal> &distor,
                   CombineMode CM) = 0;

  private:
    
    const Map<Ordinal> map_;
    Teuchos::RCP< Teuchos::Comm<Ordinal> > Comm_;
    std::vector<Scalar> imports_;
    std::vector<Scalar> exports_;
    std::vector<Ordinal> sizes_;

  }; // class DistObject

  template <typename Ordinal, typename Scalar>
  DistObject<Ordinal,Scalar>::DistObject(const Map<Ordinal>& map, Teuchos::RCP< Teuchos::Comm<Ordinal> > comm)
  : Teuchos::Object("Tpetra::DistObject")
  , map_(map)
  , Comm_(comm)
  , imports_()
  , exports_()
  , sizes_()
  {}

  template <typename Ordinal, typename Scalar>
  DistObject<Ordinal,Scalar>::DistObject(const Map<Ordinal>& map, 
      Teuchos::RCP< Teuchos::Comm<Ordinal> > comm, const std::string & Label)
  : Teuchos::Object(Label.c_str())
  , map_(map)
  , Comm_(comm)
  , imports_()
  , exports_()
  , sizes_()
  {}

  template <typename Ordinal, typename Scalar>
  DistObject<Ordinal,Scalar>::DistObject(const DistObject<Ordinal, Scalar>& source)
  : Teuchos::Object(source.label())
  , map_(source.map_)
  , Comm_(source.Comm_)
  , imports_(source.imports_)
  , exports_(source.exports_)
  , sizes_(source.sizes_)
  {}

  template <typename Ordinal, typename Scalar>
  DistObject<Ordinal,Scalar>::~DistObject() 
  {}

  template <typename Ordinal, typename Scalar>
  void DistObject<Ordinal,Scalar>::doImport(const DistObject<Ordinal, Scalar> & source, 
                                           const Import<Ordinal> & importer, CombineMode CM) 
  {
    TEST_FOR_EXCEPTION( getMap() != importer.getTargetMap(), std::runtime_error,
        "Target Maps don't match.");
    TEST_FOR_EXCEPTION( source.getMap() != importer.getSourceMap(), std::runtime_error,
        "Source Maps don't match.");

    // copy variables from importer
    Ordinal numSameIDs = importer.getNumSameIDs();
    Ordinal numPermuteIDs = importer.getNumPermuteIDs();
    Ordinal numRemoteIDs = importer.getNumRemoteIDs();
    Ordinal numExportIDs = importer.getNumExportIDs();
    const std::vector<Ordinal> & exportLIDs = importer.getExportLIDs();
    const std::vector<Ordinal> & remoteLIDs = importer.getRemoteLIDs();
    const std::vector<Ordinal> & permuteToLIDs = importer.getPermuteToLIDs();
    const std::vector<Ordinal> & permuteFromLIDs = importer.getPermuteFromLIDs();

    // call doTransfer
    doTransfer(source, CM, numSameIDs, numPermuteIDs, numRemoteIDs, numExportIDs,
        permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
        exports_, imports_, importer.getDistributor(), false);
  }

  template <typename Ordinal, typename Scalar>
  void DistObject<Ordinal,Scalar>::doExport(const DistObject<Ordinal,Scalar> & source, 
                                           const Export<Ordinal> & exporter, CombineMode CM) 
  {
    TEST_FOR_EXCEPTION( getMap() != exporter.getTargetMap(), std::runtime_error,
        "Target Maps don't match.");
    TEST_FOR_EXCEPTION( source.getMap() != exporter.getSourceMap(), std::runtime_error,
        "Source Maps don't match.");

    // copy variables from exporter
    Ordinal numSameIDs = exporter.getNumSameIDs();
    Ordinal numPermuteIDs = exporter.getNumPermuteIDs();
    Ordinal numRemoteIDs = exporter.getNumRemoteIDs();
    Ordinal numExportIDs = exporter.getNumExportIDs();
    const std::vector<Ordinal> & exportLIDs = exporter.getExportLIDs();
    const std::vector<Ordinal> & remoteLIDs = exporter.getRemoteLIDs();
    const std::vector<Ordinal> & permuteToLIDs = exporter.getPermuteToLIDs();
    const std::vector<Ordinal> & permuteFromLIDs = exporter.getPermuteFromLIDs();

    // call doTransfer
    doTransfer(source, CM, numSameIDs, numPermuteIDs, numRemoteIDs, numExportIDs,
        permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
        exports_, imports_, exporter.getDistributor(), false);
  }

  template <typename Ordinal, typename Scalar>
  void DistObject<Ordinal,Scalar>::doImport(const DistObject<Ordinal,Scalar> & source,
                                           const Export<Ordinal> & exporter, CombineMode CM)
  {
    TEST_FOR_EXCEPTION( getMap() != exporter.getTargetMap(), std::runtime_error,
        "Target Maps don't match.");
    TEST_FOR_EXCEPTION( source.getMap() != exporter.getSourceMap(), std::runtime_error,
        "Source Maps don't match.");

    // copy variables from exporter
    // note that some of them are swapped
    Ordinal numSameIDs = exporter.getNumSameIDs();
    Ordinal numPermuteIDs = exporter.getNumPermuteIDs();
    Ordinal numRemoteIDs = exporter.getNumExportIDs();
    Ordinal numExportIDs = exporter.getNumRemoteIDs();
    const std::vector<Ordinal> & exportLIDs = exporter.getRemoteLIDs();
    const std::vector<Ordinal> & remoteLIDs = exporter.getExportLIDs();
    const std::vector<Ordinal> & permuteToLIDs = exporter.getPermuteFromLIDs();
    const std::vector<Ordinal> & permuteFromLIDs = exporter.getPermuteToLIDs();

    // call doTransfer
    doTransfer(source, CM, numSameIDs, numPermuteIDs, numRemoteIDs, numExportIDs,
        permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
        imports_, exports_, exporter.getDistributor(), true);
  }

  template <typename Ordinal, typename Scalar>
  void DistObject<Ordinal,Scalar>::doExport(const DistObject<Ordinal, Scalar> & source,
                                           const Import<Ordinal> & importer, CombineMode CM)
  {
    TEST_FOR_EXCEPTION( getMap() != importer.getTargetMap(), std::runtime_error,
        "Target Maps don't match.");
    TEST_FOR_EXCEPTION( source.getMap() != importer.getSourceMap(), std::runtime_error,
        "Source Maps don't match.");

    // copy variables from importer
    // note that some of them are swapped
    Ordinal numSameIDs = importer.getNumSameIDs();
    Ordinal numPermuteIDs = importer.getNumPermuteIDs();
    Ordinal numRemoteIDs = importer.getNumExportIDs();
    Ordinal numExportIDs = importer.getNumRemoteIDs();
    const std::vector<Ordinal> & exportLIDs = importer.getRemoteLIDs();
    const std::vector<Ordinal> & remoteLIDs = importer.getExportLIDs();
    const std::vector<Ordinal> & permuteToLIDs = importer.getPermuteFromLIDs();
    const std::vector<Ordinal> & permuteFromLIDs = importer.getPermuteToLIDs();

    // call doTransfer
    doTransfer(source, CM, numSameIDs, numPermuteIDs, numRemoteIDs, numExportIDs,
        permuteToLIDs, permuteFromLIDs, remoteLIDs, exportLIDs,
        imports_, exports_, importer.getDistributor(), true);
  }

  template <typename Ordinal, typename Scalar>
  void DistObject<Ordinal,Scalar>::print(std::ostream& /*os*/) const 
  {
    TEST_FOR_EXCEPT(true); // FINISH
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
      Ordinal numSameIDs, Ordinal numPermuteIDs, Ordinal numRemoteIDs, Ordinal numExportIDs,
      const std::vector<Ordinal> & permuteToLIDs, const std::vector<Ordinal> & permuteFromLIDs,
      const std::vector<Ordinal> & remoteLIDs,    const std::vector<Ordinal> & exportLIDs,
      std::vector<Scalar>& exports, std::vector<Scalar>& imports,
      Distributor<Ordinal> &distor, bool doReverse) 
  {
    const Ordinal zero = Teuchos::OrdinalTraits<Ordinal>::zero();

    checkSizes(source);

    if(numSameIDs + numPermuteIDs > zero)
      copyAndPermute(source, numSameIDs, numPermuteIDs, permuteToLIDs, permuteFromLIDs);

    // we don't have a "Zero" CombineMode like Epetra does, so we don't have to check for that

    Ordinal packetSize = zero; // dummy value
    bool varSizes = false;
    if((!sizes_.empty()) && (numExportIDs > zero))
      sizes_.resize(numExportIDs);
    packAndPrepare(source, numExportIDs, exportLIDs, exports, packetSize, distor);

    if((isDistributed() && doReverse) || (source.getMap().isDistributed() && !doReverse)) {
      // call one of the doPostsAndWaits functions
      if(doReverse) {
        if(varSizes)
          throw reportError("var-sized doReversePostsAndWaits not implemented yet", -99);
        else
          distor.doReversePostsAndWaits(exports, packetSize, imports);
      }
      else {
        if(varSizes)
          throw reportError("var-sized doPostsAndWaits not implemented yet", -99);
        else
          distor.doPostsAndWaits(exports, packetSize, imports);
      }
      unpackAndCombine(numRemoteIDs, remoteLIDs, imports, distor, CM);
    }
  }

} // namespace Tpetra

#endif /* TPETRA_DISTOBJECT_HPP */
