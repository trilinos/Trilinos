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

#define FROSCH_ASSERT(A,S) if(!(A)) { std::cerr<<"Assertion failed. "<<S<<std::endl; std::cout.flush(); throw std::out_of_range("Assertion.");};

#include <ShyLU_DDFROSch_config.h>

#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_ExportFactory.hpp>

#ifdef HAVE_SHYLU_DDFROSCH_ZOLTAN2
#include <Zoltan2_MatrixAdapter.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#endif

namespace FROSch {

    enum DofOrdering {NodeWise=0,DimensionWise=1,Custom=2};

    enum NullSpace {LaplaceNullSpace=0,LinearElasticityNullSpace=1};

        enum Verbosity {None=0,All=1};

    template <typename LO,
              typename GO>
    class OverlappingData {

    protected:

        using IntVec        = Teuchos::Array<int>;

        using LOVec         = Teuchos::Array<LO>;

    public:

        OverlappingData(GO gid,
                        int pid,
                        LO lid);

        int Merge(const Teuchos::RCP<OverlappingData<LO,GO> > od) const;

        GO GID_;

        mutable IntVec PIDs_;

        mutable LOVec LIDs_;

    };

    template <typename LO,typename GO>
    int MergeList(Teuchos::Array<Teuchos::RCP<OverlappingData<LO,GO> > > &odList);

    template <typename LO,
              typename GO,
              typename NO>
    class LowerPIDTieBreak : public Tpetra::Details::TieBreak<LO,GO> {

    protected:

        using CommPtr                   = Teuchos::RCP<const Teuchos::Comm<int> >;

        using ConstMapPtr               = Teuchos::RCP<const Xpetra::Map<LO,GO,NO> >;

        using OverlappingDataPtr        = Teuchos::RCP<OverlappingData<LO,GO> >;
        using OverlappingDataPtrVec     = Teuchos::Array<OverlappingDataPtr>;

        using IntVec                    = Teuchos::Array<int>;
        using IntVecVecPtr              = Teuchos::ArrayRCP<IntVec>;

        using LOVec                     = Teuchos::Array<LO>;

        using GOVec                     = Teuchos::Array<GO>;
        using GOVecPtr                  = Teuchos::ArrayRCP<GO>;
        using GOVecVec                  = Teuchos::Array<GOVec>;
        using GOVecVecPtr               = Teuchos::ArrayRCP<GOVec>;

    public:
        LowerPIDTieBreak(CommPtr comm,
                         ConstMapPtr originalMap); // This is in order to estimate the length of SendImageIDs_ and ExportEntries_ in advance

        virtual bool mayHaveSideEffects() const {
            return false;
        }

        IntVecVecPtr& getComponents()
        {
            return ComponentsSubdomains_;
        }

        int sendDataToOriginalMap();

        virtual std::size_t selectedIndex(GO GID,
                                          const std::vector<std::pair<int,LO> > & pid_and_lid) const;

    protected:

        CommPtr MpiComm_;

        ConstMapPtr OriginalMap_;

        mutable LO ElementCounter_; // This is mutable such that it can be modified in selectedIndex()

        mutable OverlappingDataPtrVec OverlappingDataList_; // This is mutable such that it can be modified in selectedIndex()

        IntVecVecPtr ComponentsSubdomains_; // This is mutable such that it can be modified in selectedIndex()
    };

    template <class LO,class GO,class NO>
    Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > BuildUniqueMap(const Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > map,
                                                              bool useCreateOneToOneMap = true,
                                                              Teuchos::RCP<Tpetra::Details::TieBreak<LO,GO> > tieBreak = Teuchos::null);

    template <class SC,class LO,class GO,class NO>
    Teuchos::ArrayRCP<Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > > BuildRepeatedSubMaps(Teuchos::RCP<const Xpetra::Matrix<SC,LO,GO,NO> > matrix,
                                                                                  Teuchos::ArrayRCP<const Teuchos::RCP<Xpetra::Map<LO,GO,NO> > > subMaps);

    template <class SC,class LO,class GO,class NO>
    Teuchos::ArrayRCP<Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > > BuildRepeatedSubMaps(Teuchos::RCP<const Xpetra::CrsGraph<LO,GO,NO> > graph,
                                                                                  Teuchos::ArrayRCP<const Teuchos::RCP<Xpetra::Map<LO,GO,NO> > > subMaps);

    template <class SC,class LO,class GO,class NO>
    Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > BuildRepeatedMap(Teuchos::RCP<const Xpetra::Matrix<SC,LO,GO,NO> > matrix);

    template <class LO,class GO,class NO>
    Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > BuildRepeatedMap(Teuchos::RCP<const Xpetra::CrsGraph<LO,GO,NO> > graph);

    template <class SC,class LO,class GO,class NO>
    int ExtendOverlapByOneLayer_Old(Teuchos::RCP<const Xpetra::Matrix<SC,LO,GO,NO> > inputMatrix,
                                    Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > inputMap,
                                    Teuchos::RCP<const Xpetra::Matrix<SC,LO,GO,NO> > &outputMatrix,
                                    Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > &outputMap);

    template <class SC,class LO,class GO,class NO>
    int ExtendOverlapByOneLayer(Teuchos::RCP<const Xpetra::Matrix<SC,LO,GO,NO> > inputMatrix,
                                Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > inputMap,
                                Teuchos::RCP<const Xpetra::Matrix<SC,LO,GO,NO> > &outputMatrix,
                                Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > &outputMap);

    template <class LO,class GO,class NO>
    int ExtendOverlapByOneLayer(Teuchos::RCP<const Xpetra::CrsGraph<LO,GO,NO> > inputGraph,
                                Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > inputMap,
                                Teuchos::RCP<const Xpetra::CrsGraph<LO,GO,NO> > &outputGraph,
                                Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > &outputMap);

    template <class LO,class GO,class NO>
    Teuchos::RCP<Xpetra::Map<LO,GO,NO> > AssembleMaps(Teuchos::ArrayView<Teuchos::RCP<Xpetra::Map<LO,GO,NO> > > mapVector,
                                                      Teuchos::ArrayRCP<Teuchos::ArrayRCP<LO> > &partMappings);

    template <class LO,class GO,class NO>
    Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > MergeMaps(Teuchos::ArrayRCP<Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > > mapVector);

    template <class LO,class GO,class NO>
    int BuildDofMaps(const Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > map,
                     unsigned dofsPerNode,
                     unsigned dofOrdering,
                     Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > &nodesMap,
                     Teuchos::ArrayRCP<Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > > &dofMaps,
                     GO offset = 0);

    template <class LO,class GO,class NO>
    int BuildDofMapsVec(const Teuchos::ArrayRCP<Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > > mapVec,
                        Teuchos::ArrayRCP<unsigned> dofsPerNodeVec,
                        Teuchos::ArrayRCP<FROSch::DofOrdering> dofOrderingVec,
                        Teuchos::ArrayRCP<Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > > &nodesMapVec,
                        Teuchos::ArrayRCP<Teuchos::ArrayRCP<Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > > >&dofMapsVec);


    template <class LO,class GO,class NO>
    Teuchos::RCP<Xpetra::Map<LO,GO,NO> > BuildMapFromDofMaps(const Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::Map<LO,GO,NO> > > &dofMaps,
                                                             unsigned dofsPerNode,
                                                             unsigned dofOrdering);

    template <class LO,class GO,class NO>
    Teuchos::RCP<Xpetra::Map<LO,GO,NO> > BuildMapFromNodeMap(Teuchos::RCP<Xpetra::Map<LO,GO,NO> > &nodesMap,
                                                             unsigned dofsPerNode,
                                                             unsigned dofOrdering);

    template <class LO,class GO,class NO>
    Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::Map<LO,GO,NO> > > BuildSubMaps(Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > &fullMap,
                                                                          Teuchos::ArrayRCP<GO> maxSubGIDVec);

    template <class SC,class LO,class GO,class NO>
    Teuchos::ArrayRCP<GO> FindOneEntryOnlyRowsGlobal(Teuchos::RCP<const Xpetra::Matrix<SC,LO,GO,NO> > matrix,
                                                     Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > repeatedMap);

    template <class SC,class LO>
    bool ismultiple(Teuchos::ArrayView<SC> A,
                    Teuchos::ArrayView<SC> B);

    template<class T>
    inline void sortunique(T &v);

    template <class SC, class LO,class GO,class NO>
    Teuchos::RCP<Xpetra::MultiVector<SC,LO,GO,NO> > ModifiedGramSchmidt(Teuchos::RCP<const Xpetra::MultiVector<SC,LO,GO,NO> > multiVector,
                                                                        Teuchos::ArrayView<unsigned> zero = Teuchos::ArrayView<unsigned>());

    template <class SC, class LO,class GO,class NO>
    Teuchos::RCP<const Xpetra::MultiVector<SC,LO,GO,NO> > BuildNullSpace(unsigned dimension,
                                                                         unsigned nullSpaceType,
                                                                         Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > repeatedMap,
                                                                         unsigned dofsPerNode,
                                                                         Teuchos::ArrayRCP<Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > > dofsMaps,
                                                                         Teuchos::RCP<const Xpetra::MultiVector<SC,LO,GO,NO> > nodeList = Teuchos::null);

#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
    template <class SC,class LO,class GO,class NO>
    struct ConvertToXpetra {
    
    public:
        
        static Teuchos::RCP<Xpetra::Map<LO,GO,NO> > ConvertMap(Xpetra::UnderlyingLib lib,
                                                               const Epetra_BlockMap &map,
                                                               Teuchos::RCP<const Teuchos::Comm<int> > comm);
        
        static Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > ConvertMatrix(Xpetra::UnderlyingLib lib,
                                                                        Epetra_CrsMatrix &matrix,
                                                                        Teuchos::RCP<const Teuchos::Comm<int> > comm);
        
        static Teuchos::RCP<Xpetra::MultiVector<SC,LO,GO,NO> > ConvertMultiVector(Xpetra::UnderlyingLib lib,
                                                                                  Epetra_MultiVector &vector,
                                                                                  Teuchos::RCP<const Teuchos::Comm<int> > comm);
    };
    
    template <class SC,class LO,class NO>
    struct ConvertToXpetra<SC,LO,int,NO> {
        
    public:
        
        static Teuchos::RCP<Xpetra::Map<LO,int,NO> > ConvertMap(Xpetra::UnderlyingLib lib,
                                                                const Epetra_BlockMap &map,
                                                                Teuchos::RCP<const Teuchos::Comm<int> > comm);
        
        static Teuchos::RCP<Xpetra::Matrix<SC,LO,int,NO> > ConvertMatrix(Xpetra::UnderlyingLib lib,
                                                                         Epetra_CrsMatrix &matrix,
                                                                         Teuchos::RCP<const Teuchos::Comm<int> > comm);
        
        static Teuchos::RCP<Xpetra::MultiVector<SC,LO,int,NO> > ConvertMultiVector(Xpetra::UnderlyingLib lib,
                                                                                   Epetra_MultiVector &vector,
                                                                                   Teuchos::RCP<const Teuchos::Comm<int> > comm);
    };
    
    template <class SC,class LO,class NO>
    struct ConvertToXpetra<SC,LO,long long,NO> {
        
    public:
        
        static Teuchos::RCP<Xpetra::Map<LO,long long,NO> > ConvertMap(Xpetra::UnderlyingLib lib,
                                                                           const Epetra_BlockMap &map,
                                                                           Teuchos::RCP<const Teuchos::Comm<int> > comm);
        
        static Teuchos::RCP<Xpetra::Matrix<SC,LO,long long,NO> > ConvertMatrix(Xpetra::UnderlyingLib lib,
                                                                               Epetra_CrsMatrix &matrix,
                                                                               Teuchos::RCP<const Teuchos::Comm<int> > comm);
        
        static Teuchos::RCP<Xpetra::MultiVector<SC,LO,long long,NO> > ConvertMultiVector(Xpetra::UnderlyingLib lib,
                                                                                         Epetra_MultiVector &vector,
                                                                                         Teuchos::RCP<const Teuchos::Comm<int> > comm);
    };
#endif

    template <class Type>
    Teuchos::RCP<Type> ExtractPtrFromParameterList(Teuchos::ParameterList& paramList,
                                                   std::string namePtr="Ptr");

    template <class Type>
    Teuchos::ArrayRCP<Type> ExtractVectorFromParameterList(Teuchos::ParameterList& paramList,
                                                           std::string nameVector="Vector");

#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
    template <class LO,class GO,class NO>
    Teuchos::RCP<Epetra_Map> ConvertToEpetra(const Xpetra::Map<LO,GO,NO> &map,
                                             Teuchos::RCP<Epetra_Comm> epetraComm);

    template <class SC, class LO,class GO, class NO>
    Teuchos::RCP<Epetra_MultiVector> ConvertToEpetra(const Xpetra::MultiVector<SC,LO,GO,NO> &vector,
                                                     Teuchos::RCP<Epetra_Comm> epetraComm);

    template <class SC, class LO,class GO, class NO>
    Teuchos::RCP<Epetra_CrsMatrix> ConvertToEpetra(const Xpetra::Matrix<SC,LO,GO,NO> &matrix,
                                                   Teuchos::RCP<Epetra_Comm> epetraComm);
#endif

    template <class LO>
    Teuchos::Array<LO> GetIndicesFromString(std::string string);

#ifdef HAVE_SHYLU_DDFROSCH_ZOLTAN2
    template <class SC, class LO,class GO,class NO>
    int RepartionMatrixZoltan2(Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > &crsMatrix, Teuchos::RCP<Teuchos::ParameterList> parameterList);
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
    void ThrowErrorMissingPackage(const std::string& froschObj, const std::string& packageName);
}

#endif
