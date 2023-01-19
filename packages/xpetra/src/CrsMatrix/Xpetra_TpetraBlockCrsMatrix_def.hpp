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
#ifndef XPETRA_TPETRABLOCKCRSMATRIX_DEF_HPP
#define XPETRA_TPETRABLOCKCRSMATRIX_DEF_HPP


#include "Xpetra_TpetraBlockCrsMatrix_decl.hpp"
#include "Xpetra_TpetraCrsGraph.hpp"

namespace Xpetra {


    //! Constructor specifying fixed number of entries for each row (not implemented)
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, 
                         size_t maxNumEntriesPerRow, 
                         const Teuchos::RCP< Teuchos::ParameterList > &params)
    { 
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in"+std::string(__FILE__)+":"+std::to_string(__LINE__)); 
    }


    //! Constructor specifying (possibly different) number of entries in each row (not implemented)
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, 
                         const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, 
                         const Teuchos::RCP< Teuchos::ParameterList > &params)
    {
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in"+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }


    //! Constructor specifying column Map and fixed number of entries for each row (not implemented)
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, 
                         const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, 
                         size_t maxNumEntriesPerRow, 
                         const Teuchos::RCP< Teuchos::ParameterList > &params)
    {
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in"+std::string(__FILE__)+":"+std::to_string(__LINE__)); 
    }


    //! Constructor specifying column Map and number of entries in each row (not implemented)
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, 
                         const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, 
                         const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, 
                         const Teuchos::RCP< Teuchos::ParameterList > &params)
    { 
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in"+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }


    //! Constructor specifying a previously constructed graph ( not implemented )
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node> > &graph, 
                         const Teuchos::RCP< Teuchos::ParameterList > &params)
        // : mtx_(Teuchos::rcp(new Tpetra::BlockCrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, Node >(toTpetra(graph), params)))
        // * there is no Tpetra::BlockCrsMatrix(graph, params) c'tor.  We throw anyways here so no need to set mtx_.
    {
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in"+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }


    //! Constructor specifying a previously constructed graph & blocksize
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node> > &graph, 
                         const LocalOrdinal blockSize)
      : mtx_(Teuchos::rcp(new Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(*toTpetra(graph), blockSize))) 
    { }


    //! Constructor specifying a previously constructed graph, point maps & blocksize
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node> > &graph, 
                         const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& pointDomainMap,
                         const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& pointRangeMap,
                         const LocalOrdinal blockSize)
      : mtx_(Teuchos::rcp(new Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(*toTpetra(graph), *toTpetra(pointDomainMap), *toTpetra(pointRangeMap),blockSize)))
    { }


    //! Constructor for a fused import ( not implemented )
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP<const Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
                         const Import<LocalOrdinal,GlobalOrdinal,Node> & importer,
                         const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
                         const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
                         const Teuchos::RCP<Teuchos::ParameterList>& params)
    {
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in"+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }


    //! Constructor for a fused export (not implemented(
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP<const Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
                         const Export<LocalOrdinal,GlobalOrdinal,Node> & exporter,
                         const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
                         const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
                         const Teuchos::RCP<Teuchos::ParameterList>& params)
    {
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in"+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }


    //! Constructor for a fused import ( not implemented )
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP<const Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
                         const Import<LocalOrdinal,GlobalOrdinal,Node> & RowImporter,
                         const Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > DomainImporter,
                         const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
                         const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
                         const Teuchos::RCP<Teuchos::ParameterList>& params)
    {
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in"+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }
    

    //! Constructor for a fused export (not implemented(
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraBlockCrsMatrix(const Teuchos::RCP<const Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
                         const Export<LocalOrdinal,GlobalOrdinal,Node> & RowExporter,
                         const Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > DomainExporter,
                         const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
                         const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
                         const Teuchos::RCP<Teuchos::ParameterList>& params)
    {
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in"+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }


    //! Destructor.
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ~TpetraBlockCrsMatrix() {  }


    //@}


    //! Insert matrix entries, using global TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::IDs (not implemented)    
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    insertGlobalValues(GlobalOrdinal globalRow, 
                       const ArrayView< const GlobalOrdinal > &cols, 
                       const ArrayView< const Scalar > &vals)
    {
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in"+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }


    //! Insert matrix entries, using local IDs (not implemented)
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    insertLocalValues(LocalOrdinal localRow, 
                      const ArrayView< const LocalOrdinal > &cols, 
                      const ArrayView< const Scalar > &vals)
    {
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in"+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }


    //! Replace matrix entries, using global IDs (not implemented)
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceGlobalValues(GlobalOrdinal globalRow, 
                        const ArrayView< const GlobalOrdinal > &cols, 
                        const ArrayView< const Scalar > &vals)
    {
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in"+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }


    //! Replace matrix entries, using local IDs.
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceLocalValues (LocalOrdinal localRow,const ArrayView<const LocalOrdinal> &cols,const ArrayView<const Scalar> &vals)
    {
      XPETRA_MONITOR("TpetraBlockCrsMatrix::replaceLocalValues");
      mtx_->replaceLocalValues(localRow,cols.getRawPtr(),vals.getRawPtr(),cols.size());
    }


    //! Set all matrix entries equal to scalarThis.
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    setAllToScalar(const Scalar &alpha) 
    { 
      XPETRA_MONITOR("TpetraBlockCrsMatrix::setAllToScalar"); mtx_->setAllToScalar(alpha); 
    }


    //! Scale the current values of a matrix, this = alpha*this (not implemented)
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    scale(const Scalar &alpha)
    {
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in "+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }


    //! Allocates and returns ArrayRCPs of the Crs arrays --- This is an Xpetra-only routine.
    //** \warning This is an expert-only routine and should not be called from user code. (not implemented)
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    allocateAllValues(size_t numNonZeros,ArrayRCP<size_t> & rowptr, ArrayRCP<LocalOrdinal> & colind, ArrayRCP<Scalar> & values)
    {
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in "+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }


    //! Sets the 1D pointer arrays of the graph (not impelmented)
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    setAllValues(const ArrayRCP<size_t> & rowptr, const ArrayRCP<LocalOrdinal> & colind, const ArrayRCP<Scalar> & values)
    {
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in "+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }


    //! Gets the 1D pointer arrays of the graph (not implemented)
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getAllValues(ArrayRCP<const size_t>& rowptr, 
                 ArrayRCP<const LocalOrdinal>& colind, 
                 ArrayRCP<const Scalar>& values) const
    { 
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in "+std::string(__FILE__)+":"+std::to_string(__LINE__)); 
    }  

    //! Gets the 1D pointer arrays of the graph (not implemented)
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getAllValues(ArrayRCP<Scalar>& values) 
    { 
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in "+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }



    //@}
   
    // Transformational Methods
    //@{


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    resumeFill(const RCP< ParameterList > &params)
    { 
      /*noop*/ 
    }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    fillComplete(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &domainMap, 
                 const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rangeMap, 
                 const RCP< ParameterList > &params)
    { 
      /*noop*/ 
    }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    fillComplete(const RCP< ParameterList > &params)
    { 
      /*noop*/ 
    }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    replaceDomainMapAndImporter(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >& newDomainMap, 
                                Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > & newImporter)
    { 
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in "+std::string(__FILE__)+":"+std::to_string(__LINE__)); 
    }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    expertStaticFillComplete(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & domainMap,
                             const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & rangeMap,
                             const RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > &importer,
                             const RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > &exporter,
                             const RCP<ParameterList> &params)
    { 
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in "+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }  

    //@}


    //! @name Methods implementing RowMatrix


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getRowMap() const
    { 
      XPETRA_MONITOR("TpetraBlockCrsMatrix::getRowMap"); return toXpetra(mtx_->getRowMap()); 
    }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getColMap() const
    { 
      XPETRA_MONITOR("TpetraBlockCrsMatrix::getColMap"); return toXpetra(mtx_->getColMap()); 
    }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node> > 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getCrsGraph() const
    {
      XPETRA_MONITOR("TpetraBlockCrsMatrix::getCrsGraph"); 
      using G_t = Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>;
      using G_x = TpetraCrsGraph<LocalOrdinal,GlobalOrdinal,Node>;
      RCP<G_t> t_graph = Teuchos::rcp_const_cast<G_t>(Teuchos::rcpFromRef(mtx_->getCrsGraph()));
      RCP<const G_x> x_graph = rcp(new G_x(t_graph));
      return x_graph;
    }
    

    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    global_size_t 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getGlobalNumRows() const
    { XPETRA_MONITOR("TpetraBlockCrsMatrix::getGlobalNumRows"); return mtx_->getGlobalNumRows(); }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    global_size_t 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getGlobalNumCols() const
    { XPETRA_MONITOR("TpetraBlockCrsMatrix::getGlobalNumCols"); return mtx_->getGlobalNumCols(); }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    size_t 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getLocalNumRows() const
    { XPETRA_MONITOR("TpetraBlockCrsMatrix::getLocalNumRows"); return mtx_->getLocalNumRows(); }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    size_t 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getLocalNumCols() const
    { XPETRA_MONITOR("TpetraBlockCrsMatrix::getLocalNumCols"); return mtx_->getLocalNumCols(); }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    global_size_t 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getGlobalNumEntries() const
    { XPETRA_MONITOR("TpetraBlockCrsMatrix::getGlobalNumEntries"); return mtx_->getGlobalNumEntries(); }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    size_t 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getLocalNumEntries() const
    { XPETRA_MONITOR("TpetraBlockCrsMatrix::getLocalNumEntries"); return mtx_->getLocalNumEntries(); }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    size_t 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getNumEntriesInLocalRow(LocalOrdinal localRow) const
    { XPETRA_MONITOR("TpetraBlockCrsMatrix::getNumEntriesInLocalRow"); return mtx_->getNumEntriesInLocalRow(localRow); }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    size_t 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const
    { XPETRA_MONITOR("TpetraBlockCrsMatrix::getNumEntriesInGlobalRow"); return mtx_->getNumEntriesInGlobalRow(globalRow); }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    size_t TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getGlobalMaxNumRowEntries() const
    { XPETRA_MONITOR("TpetraBlockCrsMatrix::getGlobalMaxNumRowEntries"); return mtx_->getGlobalMaxNumRowEntries(); }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    size_t TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getLocalMaxNumRowEntries() const
    { XPETRA_MONITOR("TpetraBlockCrsMatrix::getLocalMaxNumRowEntries"); return mtx_->getLocalMaxNumRowEntries(); }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    bool TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isLocallyIndexed() const
    { XPETRA_MONITOR("TpetraBlockCrsMatrix::isLocallyIndexed"); return mtx_->isLocallyIndexed(); }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    bool TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isGloballyIndexed() const
    { XPETRA_MONITOR("TpetraBlockCrsMatrix::isGloballyIndexed"); return mtx_->isGloballyIndexed(); }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    bool TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isFillComplete() const
    { XPETRA_MONITOR("TpetraBlockCrsMatrix::isFillComplete"); return mtx_->isFillComplete(); }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    bool TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::isFillActive() const
    { XPETRA_MONITOR("TpetraBlockCrsMatrix::isFillActive"); return false; }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    typename ScalarTraits< Scalar >::magnitudeType TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getFrobeniusNorm() const
    { XPETRA_MONITOR("TpetraBlockCrsMatrix::getFrobeniusNorm"); return mtx_->getFrobeniusNorm(); }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    bool TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::supportsRowViews() const
    { XPETRA_MONITOR("TpetraBlockCrsMatrix::supportsRowViews"); return mtx_->supportsRowViews(); }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getLocalRowCopy(LocalOrdinal LocalRow, 
                    const ArrayView< LocalOrdinal > &Indices, 
                    const ArrayView< Scalar > &Values, 
                    size_t &NumEntries) const
    { 
        XPETRA_MONITOR("TpetraBlockCrsMatrix::getLocalRowCopy"); 
        typename Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::nonconst_local_inds_host_view_type indices("indices",Indices.size());
        typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::nonconst_values_host_view_type values("values", Values.size());

        mtx_->getLocalRowCopy(LocalRow, indices, values, NumEntries);
        for (size_t i=0;i<NumEntries; ++i) {
          Indices[i]=indices(i);
          Values[i]=values(i);
        }
    }
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getLocalRowView(LocalOrdinal LocalRow, ArrayView< const LocalOrdinal > &Indices,
                    ArrayView< const Scalar > &Values) const
    {
      XPETRA_MONITOR("TpetraBlockCrsMatrix::getLocalRowView");
      typename Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::local_inds_host_view_type indices;
      typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::values_host_view_type values;

      mtx_->getLocalRowView(LocalRow, indices, values);
      Indices = ArrayView<const LocalOrdinal> (indices.data(), indices.extent(0));
      Values = ArrayView<const Scalar> (reinterpret_cast<const Scalar*>(values.data()), values.extent(0));
    }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getGlobalRowView(GlobalOrdinal GlobalRow, 
                     ArrayView< const GlobalOrdinal > &Indices,
                     ArrayView< const Scalar > &Values) const
    { 
        XPETRA_MONITOR("TpetraBlockCrsMatrix::getGlobalRowView"); 
        typename Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::global_inds_host_view_type indices;
        typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::values_host_view_type values;

        mtx_->getGlobalRowView(GlobalRow, indices, values);
        Indices = ArrayView<const GlobalOrdinal> (indices.data(), indices.extent(0));
        Values = ArrayView<const Scalar> (reinterpret_cast<const Scalar*>(values.data()), values.extent(0));
    }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getGlobalRowCopy(GlobalOrdinal GlobalRow, 
                     const ArrayView< GlobalOrdinal > &Indices,
                     const ArrayView< Scalar > &Values,
                     size_t &NumEntries) const
    { 
        XPETRA_MONITOR("TpetraBlockCrsMatrix::getGlobalRowCopy"); 
        typename Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::nonconst_global_inds_host_view_type indices("indices",Indices.size());
        typename Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::nonconst_values_host_view_type values("values", Values.size());

        mtx_->getGlobalRowCopy(GlobalRow, indices, values, NumEntries);
        for (size_t i=0;i<NumEntries; ++i) {
          Indices[i]=indices(i);
          Values[i]=values(i);
        }
    }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    bool 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    haveGlobalConstants() const
    { return true; }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    apply(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &X, 
          MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &Y, 
          Teuchos::ETransp mode, 
          Scalar alpha, 
          Scalar beta) const
    { XPETRA_MONITOR("TpetraBlockCrsMatrix::apply"); mtx_->apply(toTpetra(X), toTpetra(Y), mode, alpha, beta); }

    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    apply(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &X,
          MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &Y,
          Teuchos::ETransp mode,
          Scalar alpha,
          Scalar beta,
          bool sumInterfaceValues,
          const RCP<Import<LocalOrdinal, GlobalOrdinal, Node> >& regionInterfaceImporter,
          const Teuchos::ArrayRCP<LocalOrdinal>& regionInterfaceLIDs) const
    { }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getDomainMap() const
    { XPETRA_MONITOR("TpetraBlockCrsMatrix::getDomainMap"); return toXpetra(mtx_->getDomainMap()); }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getRangeMap() const
    { XPETRA_MONITOR("TpetraBlockCrsMatrix::getRangeMap"); return toXpetra(mtx_->getRangeMap()); }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    std::string 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    description() const
    { XPETRA_MONITOR("TpetraBlockCrsMatrix::description"); return mtx_->description(); }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    describe(Teuchos::FancyOStream &out, 
             const Teuchos::EVerbosityLevel verbLevel) const
    { 
        XPETRA_MONITOR("TpetraBlockCrsMatrix::describe"); 
        mtx_->describe(out, verbLevel); 
    }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    setObjectLabel( const std::string &objectLabel )
    {
        XPETRA_MONITOR("TpetraCrsMatrix::setObjectLabel");
        Teuchos::LabeledObject::setObjectLabel(objectLabel);
        mtx_->setObjectLabel(objectLabel);
    }



    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getLocalDiagCopy(Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &diag) const
    {
        XPETRA_MONITOR("TpetraBlockCrsMatrix::getLocalDiagCopy");
        XPETRA_DYNAMIC_CAST(TpetraVectorClass, 
                            diag, 
                            tDiag, 
                            "Xpetra::TpetraBlockCrsMatrix.getLocalDiagCopy() only accept Xpetra::TpetraVector as input arguments.");
        mtx_->getLocalDiagCopy(*tDiag.getTpetra_Vector());
    }


    //! Get a copy of the diagonal entries owned by this node, with local row indices.
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getLocalDiagCopy(Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &diag, 
                     const Teuchos::ArrayView<const size_t> &offsets) const
    {
        throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in "+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }



    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getLocalDiagOffsets(Teuchos::ArrayRCP<size_t> &offsets) const
    {
        XPETRA_MONITOR("TpetraBlockCrsMatrix::getLocalDiagOffsets");

        const size_t lclNumRows = mtx_->getGraph()->getLocalNumRows();
        if (static_cast<size_t>(offsets.size()) < lclNumRows) 
        {
            offsets.resize(lclNumRows);
        }

        // The input ArrayRCP must always be a host pointer.  Thus, if
        // device_type::memory_space is Kokkos::HostSpace, it's OK for us
        // to write to that allocation directly as a Kokkos::View.
        typedef typename Node::device_type device_type;
        typedef typename device_type::memory_space memory_space;
        if (std::is_same<memory_space, Kokkos::HostSpace>::value) 
        {
            // It is always syntactically correct to assign a raw host
            // pointer to a device View, so this code will compile correctly
            // even if this branch never runs.
            typedef Kokkos::View<size_t*, device_type, Kokkos::MemoryUnmanaged> output_type;
            output_type offsetsOut (offsets.getRawPtr(), offsets.size());
            mtx_->getLocalDiagOffsets(offsetsOut);
        }
        else 
        {
            Kokkos::View<size_t*, device_type> offsetsTmp ("diagOffsets", offsets.size());
            mtx_->getLocalDiagOffsets(offsetsTmp);
            typedef Kokkos::View<size_t*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> output_type;
            output_type offsetsOut(offsets.getRawPtr(), offsets.size());
            Kokkos::deep_copy(offsetsOut, offsetsTmp);
        }
    }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    replaceDiag(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &diag)
    {
        throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix::replaceDiag: function not implemented in "+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    leftScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x)
    {
        throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in "+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    rightScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x)
    {
        throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in "+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    getMap() const
    { 
        XPETRA_MONITOR("TpetraBlockCrsMatrix::getMap"); return rcp( new TpetraMap< LocalOrdinal, GlobalOrdinal, Node >(mtx_->getMap()) ); 
    }


    //! Import.
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    doImport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &source,
             const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM)
    {
        throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in "+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }


    //! Export.
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    doExport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &dest,
                  const Import< LocalOrdinal, GlobalOrdinal, Node >& importer, CombineMode CM)
    {
        throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in "+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }


    //! Import (using an Exporter).
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    doImport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &source,
                  const Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM)
    {   
        throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in "+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }


    //! Export (using an Importer).
    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    doExport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &dest,
                  const Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM)
    {
        throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in "+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }


    template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void 
    TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
    removeEmptyProcessesInPlace (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& newMap)
    {
        throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix function not implemented in "+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool 
TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
hasMatrix() const
{ 
    return !mtx_.is_null();
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
TpetraBlockCrsMatrix(const Teuchos::RCP<Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &mtx) 
: mtx_(mtx)
{  }


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > 
TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getTpetra_BlockCrsMatrix() const
{ 
    return mtx_; 
}


// TODO: remove
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > 
TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getTpetra_BlockCrsMatrixNonConst() const
{ 
    return mtx_; 
} 

// was:     typedef typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type local_matrix_type;
//using local_matrix_type = typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type;

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_type
TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getLocalMatrixDevice () const
{
    throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix does not support getLocalMatrix due to missing Kokkos::CrsMatrix in Tpetra's experimental implementation in "+std::string(__FILE__)+":"+std::to_string(__LINE__));

#ifndef __NVCC__
    local_matrix_type ret;
#endif  // __NVCC__

    TEUCHOS_UNREACHABLE_RETURN(ret);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::local_matrix_type::HostMirror
TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getLocalMatrixHost () const
{
    throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix does not support getLocalMatrix due to missing Kokkos::CrsMatrix in Tpetra's experimental implementation in "+std::string(__FILE__)+":"+std::to_string(__LINE__));

#ifndef __NVCC__
    typename local_matrix_type::HostMirror ret;
#endif  // __NVCC__

    TEUCHOS_UNREACHABLE_RETURN(ret);
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void 
TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
setAllValues (const typename local_matrix_type::row_map_type& ptr,
              const typename local_matrix_type::StaticCrsGraphType::entries_type::non_const_type& ind,
              const typename local_matrix_type::values_type& val)
{
    throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix does not support setAllValues due to missing Kokkos::CrsMatrix in Tpetra's experimental implementation in "+std::string(__FILE__)+":"+std::to_string(__LINE__));
}


#ifdef HAVE_XPETRA_EPETRA

#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_INT))))


  // specialization of TpetraBlockCrsMatrix for GO=LO=int and Node=EpetraNode
  template <class Scalar>
  class TpetraBlockCrsMatrix<Scalar,int,int,EpetraNode>
    : public CrsMatrix<Scalar,int,int,EpetraNode>//, public TpetraRowMatrix<Scalar,int,int,Node>
  {

    // The following typedef are used by the XPETRA_DYNAMIC_CAST() macro.
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef EpetraNode Node;
    typedef TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> TpetraBlockCrsMatrixClass;
    typedef TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TpetraVectorClass;
    typedef TpetraImport<LocalOrdinal,GlobalOrdinal,Node> TpetraImportClass;
    typedef TpetraExport<LocalOrdinal,GlobalOrdinal,Node> TpetraExportClass;

  public:

    //! @name Constructor/Destructor Methods

    //! Constructor specifying fixed number of entries for each row (not implemented)
    TpetraBlockCrsMatrix(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, size_t maxNumEntriesPerRow, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "int", typeid(EpetraNode).name() );
    }

    //! Constructor specifying (possibly different) number of entries in each row (not implemented)
    TpetraBlockCrsMatrix(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "int", typeid(EpetraNode).name() );
    }

    //! Constructor specifying column Map and fixed number of entries for each row (not implemented)
    TpetraBlockCrsMatrix(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, size_t maxNumEntriesPerRow, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "int", typeid(EpetraNode).name() );
    }

    //! Constructor specifying column Map and number of entries in each row (not implemented)
    TpetraBlockCrsMatrix(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "int", typeid(EpetraNode).name() );
    }

    //! Constructor specifying a previously constructed graph ( not implemented )
    TpetraBlockCrsMatrix(const Teuchos::RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node> > &graph, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "int", typeid(EpetraNode).name() );
    }

    //! Constructor specifying a previously constructed graph & blocksize
    TpetraBlockCrsMatrix(const Teuchos::RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node> > &graph, const LocalOrdinal blockSize) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "int", typeid(EpetraNode).name() );
    }

    //! Constructor for a fused import ( not implemented )
    TpetraBlockCrsMatrix(const Teuchos::RCP<const Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
                    const Import<LocalOrdinal,GlobalOrdinal,Node> & importer,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap = Teuchos::null,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap = Teuchos::null,
       const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
    { XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "int", typeid(EpetraNode).name() );}

    //! Constructor for a fused export (not implemented(
    TpetraBlockCrsMatrix(const Teuchos::RCP<const Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
                    const Export<LocalOrdinal,GlobalOrdinal,Node> & exporter,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap = Teuchos::null,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap = Teuchos::null,
                    const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
    { XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "int", typeid(EpetraNode).name() );}

    //! Constructor for a fused import ( not implemented )
    TpetraBlockCrsMatrix(const Teuchos::RCP<const Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
                    const Import<LocalOrdinal,GlobalOrdinal,Node> & RowImporter,
                    const Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > DomainImporter,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
       const Teuchos::RCP<Teuchos::ParameterList>& params)
    { XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "int", typeid(EpetraNode).name() );}

    //! Constructor for a fused export (not implemented(
    TpetraBlockCrsMatrix(const Teuchos::RCP<const Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
                    const Export<LocalOrdinal,GlobalOrdinal,Node> & RowExporter,
                    const Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > DomainExporter,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
                    const Teuchos::RCP<Teuchos::ParameterList>& params)
    { 
        XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "int", typeid(EpetraNode).name() );
    }

    //! Destructor.
    ~TpetraBlockCrsMatrix() {  }


    //! @name Insertion/Removal Methods

    //! Insert matrix entries, using global IDs (not implemented)
    void insertGlobalValues(GlobalOrdinal globalRow, const ArrayView< const GlobalOrdinal > &cols, const ArrayView< const Scalar > &vals)
    {}

    //! Insert matrix entries, using local IDs (not implemented)
    void insertLocalValues(LocalOrdinal localRow, const ArrayView< const LocalOrdinal > &cols, const ArrayView< const Scalar > &vals)
    {}

    //! Replace matrix entries, using global IDs (not implemented)
    void replaceGlobalValues(GlobalOrdinal globalRow, const ArrayView< const GlobalOrdinal > &cols, const ArrayView< const Scalar > &vals)
    {}

    //! Replace matrix entries, using local IDs.
    void replaceLocalValues (LocalOrdinal localRow,const ArrayView<const LocalOrdinal> &cols,const ArrayView<const Scalar> &vals)
    {}

    //! Set all matrix entries equal to scalarThis.
    void setAllToScalar(const Scalar &alpha) {}

    //! Scale the current values of a matrix, this = alpha*this (not implemented)
    void scale(const Scalar &alpha)
    {}

    //! Allocates and returns ArrayRCPs of the Crs arrays --- This is an Xpetra-only routine.
    //** \warning This is an expert-only routine and should not be called from user code. (not implemented)
    void allocateAllValues(size_t numNonZeros,ArrayRCP<size_t> & rowptr, ArrayRCP<LocalOrdinal> & colind, ArrayRCP<Scalar> & values)
    {}

    //! Sets the 1D pointer arrays of the graph (not impelmented)
    void setAllValues(const ArrayRCP<size_t> & rowptr, const ArrayRCP<LocalOrdinal> & colind, const ArrayRCP<Scalar> & values)
    {}

    //! Gets the 1D pointer arrays of the graph (not implemented)
    void getAllValues(ArrayRCP<const size_t>& rowptr, ArrayRCP<const LocalOrdinal>& colind, ArrayRCP<const Scalar>& values) const
    {}
    
    
    //! Gets the 1D pointer arrays of the graph (not implemented)
    void getAllValues(ArrayRCP<Scalar>& values) 
    {}

    //! @name Transformational Methods

    //!
    void resumeFill(const RCP< ParameterList > &params=null) { /*noop*/ }

    //! Signal that data entry is complete, specifying domain and range maps.
    void fillComplete(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &domainMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rangeMap, const RCP< ParameterList > &params=null) { /*noop*/ }

    //! Signal that data entry is complete.
    void fillComplete(const RCP< ParameterList > &params=null) { /*noop*/ }


    //!  Replaces the current domainMap and importer with the user-specified objects.
    void replaceDomainMapAndImporter(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >& newDomainMap, Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> >  & newImporter)
    {}

    //! Expert static fill complete
    void expertStaticFillComplete(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & domainMap,
                                  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & rangeMap,
                                  const RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > &importer=Teuchos::null,
                                  const RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > &exporter=Teuchos::null,
                                  const RCP<ParameterList> &params=Teuchos::null)
    {}


    //! @name Methods implementing RowMatrix

    //! Returns the Map that describes the row distribution in this matrix.
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getRowMap() const { return Teuchos::null; }

    //! Returns the Map that describes the column distribution in this matrix.
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getColMap() const { return Teuchos::null; }

    //! Returns the CrsGraph associated with this matrix.
    RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node> > getCrsGraph() const
    {return Teuchos::null;}

    //! Number of global elements in the row map of this matrix.
    global_size_t getGlobalNumRows() const { return 0; }

    //! Number of global columns in the matrix.
    global_size_t getGlobalNumCols() const { return 0; }

    //! Returns the number of matrix rows owned on the calling node.
    size_t getLocalNumRows() const { return 0; }

    //! Returns the number of columns connected to the locally owned rows of this matrix.
    size_t getLocalNumCols() const { return 0; }

    //! Returns the global number of entries in this matrix.
    global_size_t getGlobalNumEntries() const { return 0; }

    //! Returns the local number of entries in this matrix.
    size_t getLocalNumEntries() const { return 0; }

    //! Returns the current number of entries on this node in the specified local row.
    size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const { return 0; }

    //! Returns the current number of entries in the (locally owned) global row.
    size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const { return 0; }

    //! Returns the maximum number of entries across all rows/columns on all nodes.
    size_t getGlobalMaxNumRowEntries() const { return 0; }

    //! Returns the maximum number of entries across all rows/columns on this node.
    size_t getLocalMaxNumRowEntries() const { return 0; }

    //! If matrix indices are in the local range, this function returns true. Otherwise, this function returns false.
    bool isLocallyIndexed() const { return false; }

    //! If matrix indices are in the global range, this function returns true. Otherwise, this function returns false.
    bool isGloballyIndexed() const { return false; }

    //! Returns true if the matrix is in compute mode, i.e. if fillComplete() has been called.
    bool isFillComplete() const { return false; }

    //! Returns true if the matrix is in edit mode.
    bool isFillActive() const { return false; }

    //! Returns the Frobenius norm of the matrix.
    typename ScalarTraits< Scalar >::magnitudeType getFrobeniusNorm() const { return Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero()); }

    //! Returns true if getLocalRowView() and getGlobalRowView() are valid for this class.
    bool supportsRowViews() const { return false; }

    //! Extract a list of entries in a specified local row of the matrix. Put into storage allocated by calling routine.
    void getLocalRowCopy(LocalOrdinal LocalRow, const ArrayView< LocalOrdinal > &Indices, const ArrayView< Scalar > &Values, size_t &NumEntries) const {  }

    //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
    void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView< const GlobalOrdinal > &indices, ArrayView< const Scalar > &values) const {  }

    //! Extract a list of entries in a specified global row of this matrix. Put into pre-allocated storage.
    void getGlobalRowCopy(GlobalOrdinal GlobalRow, const ArrayView< GlobalOrdinal > &indices, const ArrayView< Scalar > &values, size_t &numEntries) const {  }

    //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
    void getLocalRowView(LocalOrdinal LocalRow, ArrayView< const LocalOrdinal > &indices, ArrayView< const Scalar > &values) const {  }

    //! Returns true if globalConstants have been computed; false otherwise
    bool haveGlobalConstants() const {return false;}


    //! @name Methods implementing Operator

    //! Computes the sparse matrix-multivector multiplication.
    void apply(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &X, MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &Y, Teuchos::ETransp mode=Teuchos::NO_TRANS, Scalar alpha=ScalarTraits< Scalar >::one(), Scalar beta=ScalarTraits< Scalar >::zero()) const {  }

    //! Returns the Map associated with the domain of this operator. This will be null until fillComplete() is called.
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getDomainMap() const { return Teuchos::null; }

    //!
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getRangeMap() const { return Teuchos::null; }


    //! @name Overridden from Teuchos::Describable

    //! A simple one-line description of this object.
    std::string description() const { return std::string(""); }

    //! Print the object with some verbosity level to an FancyOStream object.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const {  }


    //! Deep copy constructor
    TpetraBlockCrsMatrix(const TpetraBlockCrsMatrix& matrix) {}

    //! Get a copy of the diagonal entries owned by this node, with local row idices 
    void getLocalDiagCopy(Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &diag) const {    }

    //! Get offsets of the diagonal entries in the matrix.
    void getLocalDiagOffsets(Teuchos::ArrayRCP<size_t> &offsets) const {    }

    //! Get a copy of the diagonal entries owned by this node, with local row indices.
    void getLocalDiagCopy(Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &diag, const Teuchos::ArrayView<const size_t> &offsets) const
    {}

    void replaceDiag(const Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &diag) {    }

    void leftScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) { }
    void rightScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) { }


    //! Implements DistObject interface

    //! Access function for the Tpetra::Map this DistObject was constructed with.
    Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > getMap() const { return Teuchos::null; }

    //! Import.
    void doImport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &source,
                  const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM)
    {}

    //! Export.
    void doExport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &dest,
                  const Import< LocalOrdinal, GlobalOrdinal, Node >& importer, CombineMode CM)
    {}

    //! Import (using an Exporter).
    void doImport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &source,
                  const Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM)
    {}

    //! Export (using an Importer).
    void doExport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &dest,
                  const Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM)
    {}

    void removeEmptyProcessesInPlace (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& newMap)
    {}



    //! @name Xpetra specific

    //! Does this have an underlying matrix
    bool hasMatrix() const { return false; }

    //! TpetraBlockCrsMatrix constructor to wrap a Tpetra::BlockCrsMatrix object
    TpetraBlockCrsMatrix(const Teuchos::RCP<Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &mtx) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "int", typeid(EpetraNode).name() );
    }

    //! Get the underlying Tpetra matrix
    RCP<const Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getTpetra_BlockCrsMatrix() const { return Teuchos::null; }

    //! Get the underlying Tpetra matrix
    RCP<Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getTpetra_BlockCrsMatrixNonConst() const { return Teuchos::null; }

    typedef typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type local_matrix_type;

    local_matrix_type getLocalMatrix () const {
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix does not support getLocalMatrix due to missing Kokkos::CrsMatrix in Tpetra's experimental implementation in "+std::string(__FILE__)+":"+std::to_string(__LINE__));
      local_matrix_type ret;
      return ret; // make compiler happy
    }

    void setAllValues (const typename local_matrix_type::row_map_type& ptr,
                       const typename local_matrix_type::StaticCrsGraphType::entries_type::non_const_type& ind,
                       const typename local_matrix_type::values_type& val)
    {
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix does not support setAllValues due to missing Kokkos::CrsMatrix in Tpetra's experimental implementation in "+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }

    }; // TpetraBlockCrsMatrix class


#endif  // #if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_INT))) 




#if ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))) || \
    (!defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_SERIAL) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))))

  // specialization of TpetraBlockCrsMatrix for GO=long long and Node=EpetraNode
  template <class Scalar>
  class TpetraBlockCrsMatrix<Scalar,int,long long,EpetraNode>
    : public CrsMatrix<Scalar,int,long long,EpetraNode>//, public TpetraRowMatrix<Scalar,int,int,Node>
  {

    // The following typedef are used by the XPETRA_DYNAMIC_CAST() macro.
    typedef int LocalOrdinal;
    typedef long long GlobalOrdinal;
    typedef EpetraNode Node;
    typedef TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> TpetraBlockCrsMatrixClass;
    typedef TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> TpetraVectorClass;
    typedef TpetraImport<LocalOrdinal,GlobalOrdinal,Node> TpetraImportClass;
    typedef TpetraExport<LocalOrdinal,GlobalOrdinal,Node> TpetraExportClass;

  public:

    //! @name Constructor/Destructor Methods

    //! Constructor specifying fixed number of entries for each row (not implemented)
    TpetraBlockCrsMatrix(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, size_t maxNumEntriesPerRow, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)
    {XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "long long", typeid(EpetraNode).name() );}

    //! Constructor specifying (possibly different) number of entries in each row (not implemented)
    TpetraBlockCrsMatrix(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)
    {XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "long long", typeid(EpetraNode).name() );}

    //! Constructor specifying column Map and fixed number of entries for each row (not implemented)
    TpetraBlockCrsMatrix(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, size_t maxNumEntriesPerRow, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)
    {XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "long long", typeid(EpetraNode).name() );}

    //! Constructor specifying column Map and number of entries in each row (not implemented)
    TpetraBlockCrsMatrix(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rowMap, const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &colMap, const ArrayRCP< const size_t > &NumEntriesPerRowToAlloc, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)
    {XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "long long", typeid(EpetraNode).name() );}

    //! Constructor specifying a previously constructed graph ( not implemented )
    TpetraBlockCrsMatrix(const Teuchos::RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node> > &graph, const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null)
    {XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "long long", typeid(EpetraNode).name() );}

    //! Constructor specifying a previously constructed graph & blocksize
    TpetraBlockCrsMatrix(const Teuchos::RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node> > &graph, const LocalOrdinal blockSize)
    {XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "long long", typeid(EpetraNode).name() );}




    //! Constructor for a fused import (not implemented)
    TpetraBlockCrsMatrix(const Teuchos::RCP<const Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
                    const Import<LocalOrdinal,GlobalOrdinal,Node> & importer,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap = Teuchos::null,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap = Teuchos::null,
       const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
    {XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "long long", typeid(EpetraNode).name() );}

    //! Constructor for a fused export (not implemented)
    TpetraBlockCrsMatrix(const Teuchos::RCP<const Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
                    const Export<LocalOrdinal,GlobalOrdinal,Node> & exporter,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap = Teuchos::null,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap = Teuchos::null,
                    const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null)
    {XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "long long", typeid(EpetraNode).name() );}

    //! Constructor for a fused import (not implemented)
    TpetraBlockCrsMatrix(const Teuchos::RCP<const Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
                    const Import<LocalOrdinal,GlobalOrdinal,Node> & RowImporter,
                    const Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > DomainImporter,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
       const Teuchos::RCP<Teuchos::ParameterList>& params)
    {XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "long long", typeid(EpetraNode).name() );}

    //! Constructor for a fused export (not implemented)
    TpetraBlockCrsMatrix(const Teuchos::RCP<const Tpetra::BlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& sourceMatrix,
                    const Export<LocalOrdinal,GlobalOrdinal,Node> & RowExporter,
                    const Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > DomainExporter,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& domainMap,
                    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rangeMap,
                    const Teuchos::RCP<Teuchos::ParameterList>& params)
    {XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "long long", typeid(EpetraNode).name() );}

    //! Destructor.
    ~TpetraBlockCrsMatrix() {  }


    //! @name Insertion/Removal Methods

    //! Insert matrix entries, using global IDs (not implemented)
    void insertGlobalValues(GlobalOrdinal globalRow, const ArrayView< const GlobalOrdinal > &cols, const ArrayView< const Scalar > &vals)
    {}

    //! Insert matrix entries, using local IDs (not implemented)
    void insertLocalValues(LocalOrdinal localRow, const ArrayView< const LocalOrdinal > &cols, const ArrayView< const Scalar > &vals)
    {}

    //! Replace matrix entries, using global IDs (not implemented)
    void replaceGlobalValues(GlobalOrdinal globalRow, const ArrayView< const GlobalOrdinal > &cols, const ArrayView< const Scalar > &vals)
    {}

    //! Replace matrix entries, using local IDs.
    void replaceLocalValues (LocalOrdinal localRow,const ArrayView<const LocalOrdinal> &cols,const ArrayView<const Scalar> &vals)
    {}

    //! Set all matrix entries equal to scalarThis.
    void setAllToScalar(const Scalar &alpha) {}

    //! Scale the current values of a matrix, this = alpha*this (not implemented)
    void scale(const Scalar &alpha)
    {}

    //! Allocates and returns ArrayRCPs of the Crs arrays --- This is an Xpetra-only routine.
    //** \warning This is an expert-only routine and should not be called from user code. (not implemented)
    void allocateAllValues(size_t numNonZeros,ArrayRCP<size_t> & rowptr, ArrayRCP<LocalOrdinal> & colind, ArrayRCP<Scalar> & values)
    {}

    //! Sets the 1D pointer arrays of the graph (not impelmented)
    void setAllValues(const ArrayRCP<size_t> & rowptr, const ArrayRCP<LocalOrdinal> & colind, const ArrayRCP<Scalar> & values)
    {}

    //! Gets the 1D pointer arrays of the graph (not implemented)
    void getAllValues(ArrayRCP<const size_t>& rowptr, ArrayRCP<const LocalOrdinal>& colind, ArrayRCP<const Scalar>& values) const
    {}

    
    //! Gets the 1D pointer arrays of the graph (not implemented)
    void getAllValues(ArrayRCP<Scalar>& values) 
    {}


    //! @name Transformational Methods

    //!
    void resumeFill(const RCP< ParameterList > &params=null) { /*noop*/ }

    //! Signal that data entry is complete, specifying domain and range maps.
    void fillComplete(const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &domainMap, const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &rangeMap, const RCP< ParameterList > &params=null) { /*noop*/ }

    //! Signal that data entry is complete.
    void fillComplete(const RCP< ParameterList > &params=null) { /*noop*/ }


    //!  Replaces the current domainMap and importer with the user-specified objects.
    void replaceDomainMapAndImporter(const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >& newDomainMap, Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> >  & newImporter)
    {}

    //! Expert static fill complete
    void expertStaticFillComplete(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & domainMap,
                                  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & rangeMap,
                                  const RCP<const Import<LocalOrdinal,GlobalOrdinal,Node> > &importer=Teuchos::null,
                                  const RCP<const Export<LocalOrdinal,GlobalOrdinal,Node> > &exporter=Teuchos::null,
                                  const RCP<ParameterList> &params=Teuchos::null)
    {}


    //! @name Methods implementing RowMatrix

    //! Returns the Map that describes the row distribution in this matrix.
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getRowMap() const { return Teuchos::null; }

    //! Returns the Map that describes the column distribution in this matrix.
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getColMap() const { return Teuchos::null; }

    //! Returns the CrsGraph associated with this matrix.
    RCP< const CrsGraph< LocalOrdinal, GlobalOrdinal, Node> > getCrsGraph() const
    {return Teuchos::null;}

    //! Number of global elements in the row map of this matrix.
    global_size_t getGlobalNumRows() const { return 0; }

    //! Number of global columns in the matrix.
    global_size_t getGlobalNumCols() const { return 0; }

    //! Returns the number of matrix rows owned on the calling node.
    size_t getLocalNumRows() const { return 0; }

    //! Returns the number of columns connected to the locally owned rows of this matrix.
    size_t getLocalNumCols() const { return 0; }

    //! Returns the global number of entries in this matrix.
    global_size_t getGlobalNumEntries() const { return 0; }

    //! Returns the local number of entries in this matrix.
    size_t getLocalNumEntries() const { return 0; }

    //! Returns the current number of entries on this node in the specified local row.
    size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const { return 0; }

    //! Returns the current number of entries in the (locally owned) global row.
    size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const { return 0; }

    //! Returns the maximum number of entries across all rows/columns on all nodes.
    size_t getGlobalMaxNumRowEntries() const { return 0; }

    //! Returns the maximum number of entries across all rows/columns on this node.
    size_t getLocalMaxNumRowEntries() const { return 0; }

    //! If matrix indices are in the local range, this function returns true. Otherwise, this function returns false.
    bool isLocallyIndexed() const { return false; }

    //! If matrix indices are in the global range, this function returns true. Otherwise, this function returns false.
    bool isGloballyIndexed() const { return false; }

    //! Returns true if the matrix is in compute mode, i.e. if fillComplete() has been called.
    bool isFillComplete() const { return false; }

    //! Returns true if the matrix is in edit mode.
    bool isFillActive() const { return false; }

    //! Returns the Frobenius norm of the matrix.
    typename ScalarTraits< Scalar >::magnitudeType getFrobeniusNorm() const { return Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero()); }

    //! Returns true if getLocalRowView() and getGlobalRowView() are valid for this class.
    bool supportsRowViews() const { return false; }

    //! Extract a list of entries in a specified local row of the matrix. Put into storage allocated by calling routine.
    void getLocalRowCopy(LocalOrdinal LocalRow, const ArrayView< LocalOrdinal > &Indices, const ArrayView< Scalar > &Values, size_t &NumEntries) const {  }

    //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
    void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView< const GlobalOrdinal > &indices, ArrayView< const Scalar > &values) const {  }

    //! Extract a list of entries in a specified global row of this matrix. Put into pre-allocated storage.
    void getGlobalRowCopy(GlobalOrdinal GlobalRow, const ArrayView< GlobalOrdinal > &indices, const ArrayView< Scalar > &values, size_t &numEntries) const {  }

    //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
    void getLocalRowView(LocalOrdinal LocalRow, ArrayView< const LocalOrdinal > &indices, ArrayView< const Scalar > &values) const {  }

    //! Returns true if globalConstants have been computed; false otherwise
    bool haveGlobalConstants() const {return true;}


    //! @name Methods implementing Operator

    //! Computes the sparse matrix-multivector multiplication.
    void apply(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &X, MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &Y, Teuchos::ETransp mode=Teuchos::NO_TRANS, Scalar alpha=ScalarTraits< Scalar >::one(), Scalar beta=ScalarTraits< Scalar >::zero()) const {  }

    //! Returns the Map associated with the domain of this operator. This will be null until fillComplete() is called.
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getDomainMap() const { return Teuchos::null; }

    //!
    const RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > >  getRangeMap() const { return Teuchos::null; }


    //! @name Overridden from Teuchos::Describable

    //! A simple one-line description of this object.
    std::string description() const { return std::string(""); }

    //! Print the object with some verbosity level to an FancyOStream object.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const {  }

    //! Deep copy constructor
    TpetraBlockCrsMatrix(const TpetraBlockCrsMatrix& matrix) {}

    //! Get a copy of the diagonal entries owned by this node, with local row idices 
    void getLocalDiagCopy(Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &diag) const {    }

    //! Get offsets of the diagonal entries in the matrix.
    void getLocalDiagOffsets(Teuchos::ArrayRCP<size_t> &offsets) const {    }

    //! Get a copy of the diagonal entries owned by this node, with local row indices.
    void getLocalDiagCopy(Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &diag, const Teuchos::ArrayView<const size_t> &offsets) const
    {}

    void replaceDiag(Vector< Scalar, LocalOrdinal, GlobalOrdinal, Node > &diag) const {    }

    void leftScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) { }
    void rightScale (const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) { }

    //! Implements DistObject interface

    //! Access function for the Tpetra::Map this DistObject was constructed with.
    Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > getMap() const { return Teuchos::null; }

    //! Import.
    void doImport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &source,
                  const Import< LocalOrdinal, GlobalOrdinal, Node > &importer, CombineMode CM)
    {}

    //! Export.
    void doExport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &dest,
                  const Import< LocalOrdinal, GlobalOrdinal, Node >& importer, CombineMode CM)
    {}

    //! Import (using an Exporter).
    void doImport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &source,
                  const Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM)
    {}

    //! Export (using an Importer).
    void doExport(const DistObject<char, LocalOrdinal, GlobalOrdinal, Node> &dest,
                  const Export< LocalOrdinal, GlobalOrdinal, Node >& exporter, CombineMode CM)
    {}

    void removeEmptyProcessesInPlace (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& newMap)
    {}



    //! @name Xpetra specific

    //! Does this have an underlying matrix
    bool hasMatrix() const { return false; }

    //! TpetraBlockCrsMatrix constructor to wrap a Tpetra::BlockCrsMatrix object
    TpetraBlockCrsMatrix(const Teuchos::RCP<Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &mtx) {
      XPETRA_TPETRA_ETI_EXCEPTION( typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name() , typeid(TpetraBlockCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,EpetraNode>).name(), "long long", typeid(EpetraNode).name() );
    }

    //! Get the underlying Tpetra matrix
    RCP<const Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getTpetra_BlockCrsMatrix() const { return Teuchos::null; }

    //! Get the underlying Tpetra matrix
    RCP<Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getTpetra_BlockCrsMatrixNonConst() const { return Teuchos::null; }

    typedef typename Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type local_matrix_type;

    local_matrix_type getLocalMatrix () const {
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix does not support getLocalMatrix due to missing Kokkos::CrsMatrix in Tpetra's experimental implementation in "+std::string(__FILE__)+":"+std::to_string(__LINE__));
      local_matrix_type ret;
      TEUCHOS_UNREACHABLE_RETURN(ret);
    }

    void setAllValues (const typename local_matrix_type::row_map_type& ptr,
                       const typename local_matrix_type::StaticCrsGraphType::entries_type::non_const_type& ind,
                       const typename local_matrix_type::values_type& val)
    {
      throw std::runtime_error("Xpetra::TpetraBlockCrsMatrix does not support setAllValues due to missing Kokkos::CrsMatrix in Tpetra's experimental implementation in "+std::string(__FILE__)+":"+std::to_string(__LINE__));
    }

    }; // TpetraBlockCrsMatrix class


#endif  // IF ((defined(EPETRA_HAVE_OMP) && (!defined(HAVE_TPETRA_INST_OPENMP) || !defined(HAVE_TPETRA_INST_INT_LONG_LONG))) 




#endif // HAVE_XPETRA_EPETRA


} // Xpetra namespace

#endif // XPETRA_TPETRABLOCKCRSMATRIX_DEF_HPP


