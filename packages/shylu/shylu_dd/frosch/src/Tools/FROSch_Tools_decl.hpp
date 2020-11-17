//@HEADER
// ************************************************************************
//
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
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
// Questions? Contact Alexander Heinlein (alexander.heinlein@uni-koeln.de)
//
// ************************************************************************
//@HEADER

#ifndef _FROSCH_TOOLS_DECL_HPP
#define _FROSCH_TOOLS_DECL_HPP

#include <ShyLU_DDFROSch_config.h>

#include <FROSch_Output.h>
#include <FROSch_Timers.h>

#include <Tpetra_Distributor.hpp>

#include <MatrixMarket_Tpetra.hpp>

#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_ExportFactory.hpp>

#ifdef HAVE_SHYLU_DDFROSCH_ZOLTAN2
#include <Zoltan2_MatrixAdapter.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#endif


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    #if defined HAVE_XPETRA_EPETRA || defined HAVE_TPETRA_INT_INT
    typedef int DefaultGlobalOrdinal;
    #elif !defined HAVE_TPETRA_INT_LONG_LONG
    typedef long DefaultGlobalOrdinal;
    #else
    typedef long long DefaultGlobalOrdinal;
    #endif

    enum DofOrdering {NodeWise=0,DimensionWise=1,Custom=2};

    enum NullSpace {LaplaceNullSpace=0,LinearElasticityNullSpace=1};

        enum Verbosity {None=0,All=1};

    template <typename LO,
              typename GO>
    class OverlappingData {

    protected:

        using IntVec        = Array<int>;

        using LOVec         = Array<LO>;

    public:

        OverlappingData(GO gid,
                        int pid,
                        LO lid);

        int Merge(const RCP<OverlappingData<LO,GO> > od) const;


        GO GID_;

        mutable IntVec PIDs_;

        mutable LOVec LIDs_;

    };

    template <typename LO,typename GO>
    int MergeList(Array<RCP<OverlappingData<LO,GO> > > &odList);

    template <typename LO,
              typename GO,
              typename NO>
    class LowerPIDTieBreak : public Tpetra::Details::TieBreak<LO,GO> {

    protected:

        using CommPtr                   = RCP<const Comm<int> >;

        using XMap                      = Map<LO,GO,NO>;
        using XMapPtr                   = RCP<XMap>;
        using ConstXMapPtr              = RCP<const XMap>;
        using XMapPtrVecPtr             = ArrayRCP<XMapPtr>;
        using ConstXMapPtrVecPtr        = ArrayRCP<ConstXMapPtr>;
        using XMapPtrVecPtr2D           = ArrayRCP<XMapPtrVecPtr>;
        using ConstXMapPtrVecPtr2D      = ArrayRCP<ConstXMapPtrVecPtr>;

        using OverlappingDataPtr        = RCP<OverlappingData<LO,GO> >;
        using OverlappingDataPtrVec     = Array<OverlappingDataPtr>;

        using UN                        = unsigned;

        using IntVec                    = Array<int>;
        using IntVecVecPtr              = ArrayRCP<IntVec>;

        using LOVec                     = Array<LO>;

        using GOVec                     = Array<GO>;
        using GOVecPtr                  = ArrayRCP<GO>;
        using GOVecVec                  = Array<GOVec>;
        using GOVecVecPtr               = ArrayRCP<GOVec>;

    public:
        LowerPIDTieBreak(CommPtr comm,
                         ConstXMapPtr originalMap,
                         UN dimension,
                         UN levelID = 1); // This is in order to estimate the length of SendImageIDs_ and ExportEntries_ in advance

        virtual bool mayHaveSideEffects() const {
            return false;
        }

        IntVecVecPtr& getComponents()
        {
            return ComponentsSubdomains_;
        }

        int sendDataToOriginalMap();

        virtual size_t selectedIndex(GO GID,
                                          const vector<pair<int,LO> > & pid_and_lid) const;

    protected:

        CommPtr MpiComm_;

        ConstXMapPtr OriginalMap_;

        mutable LO ElementCounter_; // This is mutable such that it can be modified in selectedIndex()

        mutable OverlappingDataPtrVec OverlappingDataList_; // This is mutable such that it can be modified in selectedIndex()

        IntVecVecPtr ComponentsSubdomains_; // This is mutable such that it can be modified in selectedIndex()

        UN LevelID_ = 1;
    };

    template <class SC, class LO, class GO, class NO>
    void writeMM(std::string fileName, Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > &matrix_);

    template<class SC, class LO, class GO, class NO>
    void writeMM(Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> >& matrix_,std::string fileName);

    template <class SC, class LO, class GO, class NO>
    void readMM(std::string fileName, Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > &matrix_,RCP<const Comm<int> > &comm);

    template <class LO,class GO,class NO>
    RCP<const Map<LO,GO,NO> > BuildUniqueMap(const RCP<const Map<LO,GO,NO> > map,
                                             bool useCreateOneToOneMap = true,
                                             RCP<Tpetra::Details::TieBreak<LO,GO> > tieBreak = null);

    template <class SC,class LO,class GO,class NO>
    ArrayRCP<RCP<const Map<LO,GO,NO> > > BuildRepeatedSubMaps(RCP<const Matrix<SC,LO,GO,NO> > matrix,
                                                              ArrayRCP<const RCP<Map<LO,GO,NO> > > subMaps);

    template <class SC,class LO,class GO,class NO>
    ArrayRCP<RCP<const Map<LO,GO,NO> > > BuildRepeatedSubMaps(RCP<const CrsGraph<LO,GO,NO> > graph,
                                                              ArrayRCP<const RCP<Map<LO,GO,NO> > > subMaps);

    template <class SC,class LO,class GO,class NO>
    RCP<Map<LO,GO,NO> > BuildRepeatedMapNonConstOld(RCP<const Matrix<SC,LO,GO,NO> > matrix);

    template <class SC,class LO,class GO,class NO>
    RCP<const Map<LO,GO,NO> > BuildRepeatedMapOld(RCP<const Matrix<SC,LO,GO,NO> > matrix);

    template <class LO,class GO,class NO>
    RCP<Map<LO,GO,NO> > BuildRepeatedMapNonConstOld(RCP<const CrsGraph<LO,GO,NO> > graph);

    template <class LO,class GO,class NO>
    RCP<const Map<LO,GO,NO> > BuildRepeatedMapOld(RCP<const CrsGraph<LO,GO,NO> > graph);

    template <class SC,class LO,class GO,class NO>
    RCP<Map<LO,GO,NO> > BuildRepeatedMapNonConst(RCP<const Matrix<SC,LO,GO,NO> > matrix);

    template <class SC,class LO,class GO,class NO>
    RCP<const Map<LO,GO,NO> > BuildRepeatedMap(RCP<const Matrix<SC,LO,GO,NO> > matrix);

    template <class LO,class GO,class NO>
    RCP<Map<LO,GO,NO> > BuildRepeatedMapNonConst(RCP<const CrsGraph<LO,GO,NO> > graph);

    template <class LO,class GO,class NO>
    RCP<const Map<LO,GO,NO> > BuildRepeatedMap(RCP<const CrsGraph<LO,GO,NO> > graph);

    template <class LO,class GO,class NO>
    Teuchos::RCP<Xpetra::Map<LO,GO,NO> > BuildMapFromNodeMapRepeated(Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > &nodesMap,
                                                                     unsigned dofsPerNode,
                                                                     unsigned dofOrdering);

    template <class SC,class LO,class GO,class NO>
    int ExtendOverlapByOneLayer_Old(RCP<const Matrix<SC,LO,GO,NO> > inputMatrix,
                                    RCP<const Map<LO,GO,NO> > inputMap,
                                    RCP<const Matrix<SC,LO,GO,NO> > &outputMatrix,
                                    RCP<const Map<LO,GO,NO> > &outputMap);

    template <class SC,class LO,class GO,class NO>
    int ExtendOverlapByOneLayer(RCP<const Matrix<SC,LO,GO,NO> > inputMatrix,
                                RCP<const Map<LO,GO,NO> > inputMap,
                                RCP<const Matrix<SC,LO,GO,NO> > &outputMatrix,
                                RCP<const Map<LO,GO,NO> > &outputMap);

    template <class LO,class GO,class NO>
    int ExtendOverlapByOneLayer(RCP<const CrsGraph<LO,GO,NO> > inputGraph,
                                RCP<const Map<LO,GO,NO> > inputMap,
                                RCP<const CrsGraph<LO,GO,NO> > &outputGraph,
                                RCP<const Map<LO,GO,NO> > &outputMap);

    /*! \brief Sort the Xpetra::Map by the global IDs \c x
     * \param[in] inputMap Unsorted input map
     */
    template <class LO,class GO,class NO>
    RCP<const Map<LO,GO,NO> > SortMapByGlobalIndex(RCP<const Map<LO,GO,NO> > inputMap);

    template <class LO,class GO,class NO>
    RCP<Map<LO,GO,NO> > AssembleMaps(ArrayView<RCP<const Map<LO,GO,NO> > > mapVector,
                                     ArrayRCP<ArrayRCP<LO> > &partMappings);

    template <class LO,class GO,class NO>
    RCP<Map<LO,GO,NO> > AssembleMapsNonConst(ArrayView<RCP<Map<LO,GO,NO> > > mapVector,
                                             ArrayRCP<ArrayRCP<LO> > &partMappings);

    template <class LO,class GO,class NO>
    RCP<Map<LO,GO,NO> > AssembleSubdomainMap(unsigned numberOfBlocks,
                                             ArrayRCP<ArrayRCP<RCP<const Map<LO,GO,NO> > > > dofsMaps,
                                             ArrayRCP<unsigned> dofsPerNode);

    template <class LO,class GO,class NO>
    RCP<Map<LO,GO,NO> > MergeMapsNonConst(ArrayRCP<RCP<const Map<LO,GO,NO> > > mapVector);

    template <class LO,class GO,class NO>
    RCP<const Map<LO,GO,NO> > MergeMaps(ArrayRCP<RCP<const Map<LO,GO,NO> > > mapVector);

    template <class LO,class GO,class NO>
    int BuildDofMaps(const RCP<const Map<LO,GO,NO> > map,
                     unsigned dofsPerNode,
                     unsigned dofOrdering,
                     RCP<const Map<LO,GO,NO> > &nodesMap,
                     ArrayRCP<RCP<const Map<LO,GO,NO> > > &dofMaps,
                     GO offset = 0);

    template <class LO,class GO,class NO>
    int BuildDofMapsVec(const ArrayRCP<RCP<const Map<LO,GO,NO> > > mapVec,
                        ArrayRCP<unsigned> dofsPerNodeVec,
                        ArrayRCP<FROSch::DofOrdering> dofOrderingVec,
                        ArrayRCP<RCP<const Map<LO,GO,NO> > > &nodesMapVec,
                        ArrayRCP<ArrayRCP<RCP<const Map<LO,GO,NO> > > >&dofMapsVec);


    template <class LO,class GO,class NO>
    RCP<Map<LO,GO,NO> > BuildMapFromDofMaps(const ArrayRCP<RCP<Map<LO,GO,NO> > > &dofMaps,
                                            unsigned dofsPerNode,
                                            unsigned dofOrdering);

    template <class LO,class GO,class NO>
    RCP<Map<LO,GO,NO> > BuildMapFromNodeMap(RCP<const Map<LO,GO,NO> > &nodesMap,
                                            unsigned dofsPerNode,
                                            unsigned dofOrdering);

    template <class LO,class GO,class NO>
    ArrayRCP<RCP<const Map<LO,GO,NO> > > BuildNodeMapsFromDofMaps(ArrayRCP<ArrayRCP<RCP<const Map<LO,GO,NO> > > >dofsMapsVecVec,
                                                                  ArrayRCP<unsigned> dofsPerNodeVec,
                                                                  ArrayRCP<DofOrdering> dofOrderingVec);

    template <class LO,class GO,class NO>
    ArrayRCP<RCP<Map<LO,GO,NO> > > BuildSubMaps(RCP<const Map<LO,GO,NO> > &fullMap,
                                                ArrayRCP<GO> maxSubGIDVec);

    template <class SC,class LO,class GO,class NO>
    ArrayRCP<GO> FindOneEntryOnlyRowsGlobal(RCP<const Matrix<SC,LO,GO,NO> > matrix,
                                            RCP<const Map<LO,GO,NO> > repeatedMap);

    template <class LO,class GO,class NO>
    ArrayRCP<GO> FindOneEntryOnlyRowsGlobal(RCP<const CrsGraph<LO,GO,NO> > graph,
                                            RCP<const Map<LO,GO,NO> > repeatedMap);

    template <class SC,class LO>
    bool ismultiple(ArrayView<SC> A,
                    ArrayView<SC> B);

    template<class T>
    inline void sort(T &v);

    template<class T>
    inline void sortunique(T &v);

    template <class SC, class LO,class GO,class NO>
    RCP<MultiVector<SC,LO,GO,NO> > ModifiedGramSchmidt(RCP<const MultiVector<SC,LO,GO,NO> > multiVector,
                                                       ArrayView<unsigned> zero = ArrayView<unsigned>());

    template <class SC, class LO,class GO,class NO>
    RCP<const MultiVector<SC,LO,GO,NO> > BuildNullSpace(unsigned dimension,
                                                        unsigned nullSpaceType,
                                                        RCP<const Map<LO,GO,NO> > repeatedMap,
                                                        unsigned dofsPerNode,
                                                        ArrayRCP<RCP<const Map<LO,GO,NO> > > dofsMaps,
                                                        RCP<const MultiVector<SC,LO,GO,NO> > nodeList = null);

#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
    template <class SC,class LO,class GO,class NO>
    struct ConvertToXpetra {

    public:

        static RCP<Map<LO,GO,NO> > ConvertMap(UnderlyingLib lib,
                                              const Epetra_BlockMap &map,
                                              RCP<const Comm<int> > comm);

        static RCP<Matrix<SC,LO,GO,NO> > ConvertMatrix(UnderlyingLib lib,
                                                       Epetra_CrsMatrix &matrix,
                                                       RCP<const Comm<int> > comm);

        static RCP<MultiVector<SC,LO,GO,NO> > ConvertMultiVector(UnderlyingLib lib,
                                                                 Epetra_MultiVector &vector,
                                                                 RCP<const Comm<int> > comm);
    };

    template <class SC,class LO,class NO>
    struct ConvertToXpetra<SC,LO,int,NO> {

    public:

        static RCP<Map<LO,int,NO> > ConvertMap(UnderlyingLib lib,
                                               const Epetra_BlockMap &map,
                                               RCP<const Comm<int> > comm);

        static RCP<Matrix<SC,LO,int,NO> > ConvertMatrix(UnderlyingLib lib,
                                                        Epetra_CrsMatrix &matrix,
                                                        RCP<const Comm<int> > comm);

        static RCP<MultiVector<SC,LO,int,NO> > ConvertMultiVector(UnderlyingLib lib,
                                                                  Epetra_MultiVector &vector,
                                                                  RCP<const Comm<int> > comm);
    };

    template <class SC,class LO,class NO>
    struct ConvertToXpetra<SC,LO,long long,NO> {

    public:

        static RCP<Map<LO,long long,NO> > ConvertMap(UnderlyingLib lib,
                                                     const Epetra_BlockMap &map,
                                                     RCP<const Comm<int> > comm);

        static RCP<Matrix<SC,LO,long long,NO> > ConvertMatrix(UnderlyingLib lib,
                                                              Epetra_CrsMatrix &matrix,
                                                              RCP<const Comm<int> > comm);

        static RCP<MultiVector<SC,LO,long long,NO> > ConvertMultiVector(UnderlyingLib lib,
                                                                        Epetra_MultiVector &vector,
                                                                        RCP<const Comm<int> > comm);
    };
#endif

    template <class Type>
    RCP<Type> ExtractPtrFromParameterList(ParameterList& paramList,
                                          string namePtr="Ptr");

    template <class Type>
    ArrayRCP<Type> ExtractVectorFromParameterList(ParameterList& paramList,
                                                  string nameVector="Vector");

#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
    template <class LO,class GO,class NO>
    RCP<Epetra_Map> ConvertToEpetra(const Map<LO,GO,NO> &map,
                                    RCP<Epetra_Comm> epetraComm);

    template <class SC,class LO,class GO,class NO>
    RCP<Epetra_MultiVector> ConvertToEpetra(const MultiVector<SC,LO,GO,NO> &vector,
                                            RCP<Epetra_Comm> epetraComm);

    template <class SC,class LO,class GO,class NO>
    RCP<Epetra_CrsMatrix> ConvertToEpetra(const Matrix<SC,LO,GO,NO> &matrix,
                                          RCP<Epetra_Comm> epetraComm);
#endif

    template <class LO>
    Array<LO> GetIndicesFromString(string string);

#ifdef HAVE_SHYLU_DDFROSCH_ZOLTAN2
    template <class SC,class LO,class GO,class NO>
    int RepartionMatrixZoltan2(RCP<Matrix<SC,LO,GO,NO> > &crsMatrix,
                               RCP<ParameterList> parameterList);

    template <class LO,class GO, class NO>
    int BuildRepMapZoltan(RCP<CrsGraph<LO,GO,NO> > Xgraph,
                          RCP<CrsGraph<LO,GO,NO> > B,
                          RCP<ParameterList> parameterList,
                          Teuchos::RCP<const Teuchos::Comm<int> > TeuchosComm,
                          RCP<Map<LO,GO,NO> > &RepeatedMap);
#endif

    /*!
    \brief Throw runtime error due to missing package in build configuration

    As many packages are optional, we might detect only at runtime that are certain package
    is not included into the build configuration, but still is used by FROSch.
    Use this routine to throw a generic error message with some information for the user
    and provide details how to fix it.

    \param[in] forschObj FROSch object that is asking for the missing package
    \param[in] packageName Name of the missing package
    */
    inline void ThrowErrorMissingPackage(const string& froschObj,
                                         const string& packageName)
    {
        // Create the error message
        stringstream errMsg;
        errMsg << froschObj << " is asking for the Trilinos packate '"<< packageName << "', "
        "but this package is not included in your build configuration. "
        "Please enable '" << packageName << "' in your build configuration to be used with ShyLU_DDFROSch.";

        // Throw the error
        FROSCH_ASSERT(false, errMsg.str());

        return;
    }
}

#endif
