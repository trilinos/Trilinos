// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Test for FROSch::ExtendOverlapByLayers_Threaded.
//
// For every case the result is compared, per rank and GID for GID, with
//   (A) ExtendOverlapByOneLayer applied overlap times (original graph as input
//       graph every time), followed by SortMapByGlobalIndex, and
//   (B) a plain breadth-first search done in this test on the gathered global
//       graph: layer 0 = start GIDs, layer i+1 = union of the column GIDs of the
//       rows in layer i (layer i is not kept implicitly).

#include <ShyLU_DDFROSch_config.h>

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <map>
#include <set>
#include <vector>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_StackedTimer.hpp>

// Galeri::Xpetra
#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"
#include "Galeri_XpetraParameters.hpp"
#include "Galeri_XpetraUtils.hpp"
#include "Galeri_XpetraMaps.hpp"

// Xpetra include
#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MatrixFactory.hpp>

#include <FROSch_Tools_def.hpp>


using UN    = unsigned;
using SC    = double;
using LO    = int;
using GO    = FROSch::DefaultGlobalOrdinal;
using NO    = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType;

using namespace std;
using namespace Teuchos;
using namespace Xpetra;

using XMap          = Map<LO,GO,NO>;
using XGraph        = CrsGraph<LO,GO,NO>;
using ConstXMapPtr  = RCP<const XMap>;
using ConstXGraphPtr= RCP<const XGraph>;
using Adjacency     = std::map<GO,std::vector<GO> >;

namespace {

    const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();

    std::vector<GO> toVector(ArrayView<const GO> view)
    {
        std::vector<GO> vec(view.size());
        for (size_t i=0; i<(size_t) view.size(); i++) vec[i] = view[i];
        return vec;
    }

    ConstXMapPtr buildMap(const std::vector<GO> &gids,
                          RCP<const Comm<int> > comm)
    {
        Array<GO> gidsArray(gids.begin(),gids.end());
        return MapFactory<LO,GO,NO>::Build(UseTpetra,INVALID,gidsArray(),0,comm).getConst();
    }

    // Graph with the given owned rows; the columns of row g are adj.at(g) (possibly empty, possibly without g)
    ConstXGraphPtr buildGraph(const std::vector<GO> &ownedRows,
                              const Adjacency &adj,
                              RCP<const Comm<int> > comm)
    {
        ConstXMapPtr rowMap = buildMap(ownedRows,comm);
        size_t maxNnz = 1;
        for (size_t i=0; i<ownedRows.size(); i++) {
            maxNnz = std::max(maxNnz,adj.at(ownedRows[i]).size());
        }
        RCP<XGraph> graph = CrsGraphFactory<LO,GO,NO>::Build(rowMap,maxNnz);
        for (size_t i=0; i<ownedRows.size(); i++) {
            const std::vector<GO> &cols = adj.at(ownedRows[i]);
            if (cols.size()>0) {
                Array<GO> colsArray(cols.begin(),cols.end());
                graph->insertGlobalIndices(ownedRows[i],colsArray());
            }
        }
        graph->fillComplete();
        return graph.getConst();
    }

    // Every rank receives the complete (small) global graph: row GID -> column GIDs
    Adjacency gatherGlobalGraph(ConstXGraphPtr graph)
    {
        RCP<const Comm<int> > comm = graph->getComm();
        ConstXMapPtr rowMap = graph->getRowMap();
        ConstXMapPtr colMap = graph->getColMap();

        Array<GO> sendBuffer;
        for (LO i=0; i<(LO) graph->getLocalNumRows(); i++) {
            ArrayView<const LO> indices;
            graph->getLocalRowView(i,indices);
            sendBuffer.push_back(rowMap->getGlobalElement(i));
            sendBuffer.push_back(GO(indices.size()));
            for (LO k=0; k<(LO) indices.size(); k++) {
                sendBuffer.push_back(colMap->getGlobalElement(indices[k]));
            }
        }

        int numProcs = comm->getSize();
        int length = (int) sendBuffer.size();
        Array<int> lengths(numProcs,0);
        gatherAll<int,int>(*comm,1,&length,numProcs,lengths.getRawPtr());
        int maxLength = 1;
        for (int p=0; p<numProcs; p++) maxLength = std::max(maxLength,lengths[p]);
        sendBuffer.resize(maxLength,GO(0));
        Array<GO> recvBuffer(maxLength*numProcs);
        gatherAll<int,GO>(*comm,maxLength,sendBuffer.getRawPtr(),maxLength*numProcs,recvBuffer.getRawPtr());

        Adjacency adj;
        for (int p=0; p<numProcs; p++) {
            int pos = p*maxLength;
            int end = pos+lengths[p];
            while (pos<end) {
                GO row = recvBuffer[pos++];
                GO numCols = recvBuffer[pos++];
                std::vector<GO> cols((size_t) numCols);
                for (GO k=0; k<numCols; k++) cols[k] = recvBuffer[pos++];
                adj[row] = cols;
            }
        }
        return adj;
    }

    // Plain BFS; returns layers 0..overlap, each sorted and unique
    std::vector<std::vector<GO> > referenceLayers(const Adjacency &adj,
                                                  ArrayView<const GO> start,
                                                  int overlap,
                                                  bool &referenceOk)
    {
        std::vector<std::vector<GO> > layers(overlap+1);
        std::vector<GO> startVec = toVector(start);
        std::set<GO> current(startVec.begin(),startVec.end());
        layers[0].assign(current.begin(),current.end());
        for (int l=0; l<overlap; l++) {
            std::set<GO> next;
            for (typename std::set<GO>::const_iterator it=current.begin(); it!=current.end(); ++it) {
                Adjacency::const_iterator row = adj.find(*it);
                if (row==adj.end()) {
                    referenceOk = false; // GID is not a row of the graph: test setup error
                    continue;
                }
                next.insert(row->second.begin(),row->second.end());
            }
            current.swap(next);
            layers[l+1].assign(current.begin(),current.end());
        }
        return layers;
    }

    template <class T>
    std::string listToString(const T &list)
    {
        std::ostringstream os;
        os << "[";
        size_t n = 0;
        for (typename T::const_iterator it=list.begin(); it!=list.end(); ++it, ++n) {
            if (n>=20) {os << " ..."; break;}
            os << (n>0 ? " " : "") << *it;
        }
        os << "] (size " << list.size() << ")";
        return os.str();
    }

    // One collective check; returns the local pass flag
    bool checkOverlap(const std::string &label,
                      ConstXGraphPtr graph,
                      ConstXMapPtr startMap,
                      int overlap,
                      const Adjacency &adj)
    {
        RCP<const Comm<int> > comm = graph->getComm();
        const int rank = comm->getRank();
        bool ok = true;
        std::ostringstream reason;

        // Function under test; layerSizes is pre-filled with junk on purpose
        Array<LO> layerSizes(3,LO(-7));
        ConstXMapPtr result = FROSch::ExtendOverlapByLayers_Threaded<LO,GO,NO>(graph,startMap,overlap,layerSizes);

        // Reference A: ExtendOverlapByOneLayer overlap times, then SortMapByGlobalIndex
        ConstXMapPtr mapA = startMap;
        for (int i=0; i<overlap; i++) {
            ConstXGraphPtr outputGraph;
            ConstXMapPtr outputMap;
            FROSch::ExtendOverlapByOneLayer<LO,GO,NO>(graph,mapA,outputGraph,outputMap);
            mapA = outputMap;
        }
        mapA = FROSch::SortMapByGlobalIndex<LO,GO,NO>(mapA);

        // Reference B: plain BFS
        bool referenceOk = true;
        std::vector<std::vector<GO> > layersB = referenceLayers(adj,startMap->getLocalElementList(),overlap,referenceOk);
        if (!referenceOk) {
            ok = false; reason << " reference B: a GID in a layer is not a row of the graph (test setup);";
        }

        if (result.is_null()) {
            ok = false; reason << " returned map is null;";
        } else {
            std::vector<GO> res = toVector(result->getLocalElementList());
            std::vector<GO> refA = toVector(mapA->getLocalElementList());
            const std::vector<GO> &refB = layersB[overlap];

            if (refA!=refB) {
                ok = false; reason << " references disagree: A=" << listToString(refA) << " B=" << listToString(refB) << ";";
            }
            if (res!=refA) {
                ok = false; reason << " result != reference A: result=" << listToString(res) << " A=" << listToString(refA) << ";";
            }
            if (res!=refB) {
                ok = false; reason << " result != reference B: result=" << listToString(res) << " B=" << listToString(refB) << ";";
            }
            for (size_t i=1; i<res.size(); i++) {
                if (!(res[i-1]<res[i])) {
                    ok = false; reason << " result not strictly ascending at local index " << i << ";";
                    break;
                }
            }
            if (result->getComm().is_null() || result->getComm()->getSize()!=comm->getSize() || result->getComm()->getRank()!=rank) {
                ok = false; reason << " result map comm differs from graph comm;";
            }
        }

        if ((int) layerSizes.size()!=overlap) {
            ok = false; reason << " layerSizes.size()=" << layerSizes.size() << " expected " << overlap << ";";
        } else {
            for (int i=0; i<overlap; i++) {
                if (layerSizes[i]!=(LO) layersB[i+1].size()) {
                    ok = false; reason << " layerSizes[" << i << "]=" << layerSizes[i] << " expected " << layersB[i+1].size() << ";";
                }
            }
        }

        if (!ok) {
            cout << "FAILED [" << label << ", overlap " << overlap << "] on rank " << rank << ":" << reason.str() << endl;
        }

        int localPass = ok ? 1 : 0, globalPass = 0;
        reduceAll(*comm,REDUCE_MIN,localPass,ptr(&globalPass));
        if (rank==0) {
            cout << (globalPass ? "passed " : "FAILED ") << "[" << label << ", overlap " << overlap << "]" << endl;
        }
        return ok;
    }

    // Contiguous partition of 0..n-1 with weights; returns [begin,end) of this rank
    void weightedPartition(GO n,
                           const std::vector<int> &weights,
                           int rank,
                           GO &begin,
                           GO &end)
    {
        long long total = 0, before = 0;
        for (size_t p=0; p<weights.size(); p++) {
            total += weights[p];
            if ((int) p<rank) before += weights[p];
        }
        begin = GO((n*before)/total);
        end = GO((n*(before+weights[rank]))/total);
    }

    std::vector<GO> contiguousRows(GO n,
                                   RCP<const Comm<int> > comm,
                                   const std::vector<int> &weights = std::vector<int>())
    {
        std::vector<int> w = weights;
        if (w.empty()) w.assign(comm->getSize(),1);
        GO begin, end;
        weightedPartition(n,w,comm->getRank(),begin,end);
        std::vector<GO> rows;
        for (GO g=begin; g<end; g++) rows.push_back(g);
        return rows;
    }

    Adjacency chainAdjacency(GO n,
                             GO offset = 0)
    {
        Adjacency adj;
        for (GO i=0; i<n; i++) {
            std::vector<GO> cols;
            if (i>0) cols.push_back(offset+i-1);
            cols.push_back(offset+i);
            if (i<n-1) cols.push_back(offset+i+1);
            adj[offset+i] = cols;
        }
        return adj;
    }

    Adjacency grid2DAdjacency(GO nx,
                              GO ny)
    {
        Adjacency adj;
        for (GO y=0; y<ny; y++) {
            for (GO x=0; x<nx; x++) {
                std::vector<GO> cols;
                if (y>0) cols.push_back(x+nx*(y-1));
                if (x>0) cols.push_back(x-1+nx*y);
                cols.push_back(x+nx*y);
                if (x<nx-1) cols.push_back(x+1+nx*y);
                if (y<ny-1) cols.push_back(x+nx*(y+1));
                adj[x+nx*y] = cols;
            }
        }
        return adj;
    }

    bool runGraphCase(const std::string &label,
                      ConstXGraphPtr graph,
                      ConstXMapPtr startMap,
                      const std::vector<int> &overlaps)
    {
        Adjacency adj = gatherGlobalGraph(graph);
        bool ok = true;
        for (size_t i=0; i<overlaps.size(); i++) {
            ok = checkOverlap(label,graph,startMap,overlaps[i],adj) && ok;
        }
        return ok;
    }

    // 1: 1D chain, diagonal stored, contiguous partition
    bool caseChain(RCP<const Comm<int> > comm)
    {
        const GO n = 12; // diameter 11
        std::vector<GO> rows = contiguousRows(n,comm);
        ConstXGraphPtr graph = buildGraph(rows,chainAdjacency(n),comm);
        int ks[] = {0,1,2,5,11,12,20};
        bool ok = runGraphCase("chain",graph,graph->getRowMap(),std::vector<int>(ks,ks+7));

        // Same with the start GIDs given in descending order: the result must still be sorted
        std::vector<GO> reversed(rows.rbegin(),rows.rend());
        int ksRev[] = {0,1,2};
        ok = runGraphCase("chain, start map in descending order",graph,buildMap(reversed,comm),std::vector<int>(ksRev,ksRev+3)) && ok;
        return ok;
    }

    // 2: 2D 5-point grid, uneven contiguous partition (weights 1,2,3,...)
    bool caseGrid(RCP<const Comm<int> > comm)
    {
        const GO nx = 7, ny = 5;
        std::vector<int> weights(comm->getSize());
        for (int p=0; p<comm->getSize(); p++) weights[p] = p+1;
        std::vector<GO> rows = contiguousRows(nx*ny,comm,weights);
        ConstXGraphPtr graph = buildGraph(rows,grid2DAdjacency(nx,ny),comm);
        int ks[] = {1,3};
        return runGraphCase("grid2D uneven partition",graph,graph->getRowMap(),std::vector<int>(ks,ks+2));
    }

    // 3: rank 0 owns no rows
    bool caseEmptyRank(RCP<const Comm<int> > comm)
    {
        const GO nx = 6, ny = 6;
        std::vector<int> weights(comm->getSize(),1);
        if (comm->getSize()>1) weights[0] = 0;
        std::vector<GO> rows = contiguousRows(nx*ny,comm,weights);
        ConstXGraphPtr graph = buildGraph(rows,grid2DAdjacency(nx,ny),comm);
        int ks[] = {2};
        bool ok = runGraphCase("grid2D, rank 0 owns no rows",graph,graph->getRowMap(),std::vector<int>(ks,ks+1));

        // Rank 0 owns no rows but starts from a GID owned elsewhere (all its rows must be fetched)
        std::vector<GO> start = rows;
        if (comm->getRank()==0 && comm->getSize()>1) start.assign(1,GO(20));
        int ks2[] = {0,1,2};
        ok = runGraphCase("grid2D, rank 0 owns no rows, start GID 20 on rank 0",graph,buildMap(start,comm),std::vector<int>(ks2,ks2+3)) && ok;
        return ok;
    }

    // 4: two uncoupled chains (GIDs 0..9 and 10..19)
    bool caseDisconnected(RCP<const Comm<int> > comm)
    {
        Adjacency adj = chainAdjacency(10,0);
        Adjacency adj2 = chainAdjacency(10,10);
        adj.insert(adj2.begin(),adj2.end());
        std::vector<GO> rows = contiguousRows(20,comm);
        ConstXGraphPtr graph = buildGraph(rows,adj,comm);
        int ks[] = {1,3,9,12,25};
        return runGraphCase("two disconnected chains",graph,graph->getRowMap(),std::vector<int>(ks,ks+5));
    }

    // 5: chain where row 4 has no diagonal entry and row 7 is empty
    bool caseNoDiagonal(RCP<const Comm<int> > comm)
    {
        const GO n = 10;
        Adjacency adj = chainAdjacency(n);
        adj[4] = std::vector<GO>(); adj[4].push_back(3); adj[4].push_back(5);
        adj[7] = std::vector<GO>();
        std::vector<GO> rows = contiguousRows(n,comm);
        ConstXGraphPtr graph = buildGraph(rows,adj,comm);
        int ks[] = {1,2,3};
        bool ok = runGraphCase("chain without diagonal in row 4, empty row 7",graph,graph->getRowMap(),std::vector<int>(ks,ks+3));

        // Start from single GIDs: {4} on rank 0 (4 must not stay in layer 1), {7} on the last rank (layers become empty)
        std::vector<GO> start = rows;
        if (comm->getRank()==0) start.assign(1,GO(4));
        if (comm->getRank()==comm->getSize()-1 && comm->getSize()>1) start.assign(1,GO(7));
        ok = runGraphCase("chain without diagonal, start {4} on rank 0 and {7} on last rank",graph,buildMap(start,comm),std::vector<int>(ks,ks+3)) && ok;
        return ok;
    }

    // 6: large, gapped GIDs (base 2^33, stride 3), round-robin ownership
    bool caseGapped(RCP<const Comm<int> > comm)
    {
        if (sizeof(GO)<8) {
            if (comm->getRank()==0) cout << "skipped [gapped GIDs]: GlobalOrdinal has " << sizeof(GO) << " bytes" << endl;
            return true;
        }
        const long long base = 8589934592LL; // 2^33
        const long long n = 40;
        Adjacency adj;
        std::vector<GO> rows;
        for (long long i=0; i<n; i++) {
            std::vector<GO> cols;
            if (i>0) cols.push_back(GO(base+3*(i-1)));
            cols.push_back(GO(base+3*i));
            if (i<n-1) cols.push_back(GO(base+3*(i+1)));
            adj[GO(base+3*i)] = cols;
            if (i%comm->getSize()==comm->getRank()) rows.push_back(GO(base+3*i));
        }
        ConstXGraphPtr graph = buildGraph(rows,adj,comm);
        int ks[] = {1,2,3,40};
        return runGraphCase("gapped GIDs 2^33+3i, round-robin",graph,graph->getRowMap(),std::vector<int>(ks,ks+4));
    }

    // 7: long chain, many layers
    bool caseRehash(RCP<const Comm<int> > comm)
    {
        const GO n = 2000;
        std::vector<GO> rows = contiguousRows(n,comm);
        ConstXGraphPtr graph = buildGraph(rows,chainAdjacency(n),comm);
        int ks[] = {200};
        bool ok = runGraphCase("chain 2000 rows",graph,graph->getRowMap(),std::vector<int>(ks,ks+1));

        // Start from one GID per rank so that the set grows from 1 to ~1400 local GIDs
        std::vector<GO> start;
        if (rows.size()>0) start.push_back(rows[rows.size()/2]);
        int ksSingle[] = {200,700};
        ok = runGraphCase("chain 2000 rows, one start GID per rank",graph,buildMap(start,comm),std::vector<int>(ksSingle,ksSingle+2)) && ok;
        return ok;
    }

    // 8: Galeri Laplace2D/3D with DPN dofs per node (NodeWise), start map = repeated map
    bool caseGaleri(RCP<const Comm<int> > comm,
                    int dimension,
                    int dofsPerNode,
                    int M)
    {
        int N = 0;
        if (dimension==2) {
            N = (int) (pow(comm->getSize(),1/2.) + 100*numeric_limits<double>::epsilon());
            if (N*N!=comm->getSize()) {
                if (comm->getRank()==0) cout << "FAILED: galeri case DIM=2 needs a square number of ranks" << endl;
                return false;
            }
        } else if (dimension==3) {
            N = (int) (pow(comm->getSize(),1/3.) + 100*numeric_limits<double>::epsilon());
            if (N*N*N!=comm->getSize()) {
                if (comm->getRank()==0) cout << "FAILED: galeri case DIM=3 needs a cubic number of ranks" << endl;
                return false;
            }
        } else {
            if (comm->getRank()==0) cout << "FAILED: DIM must be 2 or 3" << endl;
            return false;
        }

        UnderlyingLib xpetraLib = UseTpetra;
        ParameterList GaleriList;
        GaleriList.set("nx", GO(N*M));
        GaleriList.set("ny", GO(N*M));
        GaleriList.set("nz", GO(N*M));
        GaleriList.set("mx", GO(N));
        GaleriList.set("my", GO(N));
        GaleriList.set("mz", GO(N));

        RCP<const Map<LO,GO,NO> > UniqueMapTmp;
        RCP<Matrix<SC,LO,GO,NO> > KTmp;
        if (dimension==2) {
            UniqueMapTmp = Galeri::Xpetra::CreateMap<LO,GO,NO>(xpetraLib,"Cartesian2D",comm,GaleriList);
            RCP<Galeri::Xpetra::Problem<Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> > > Problem = Galeri::Xpetra::BuildProblem<SC,LO,GO,Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("Laplace2D",UniqueMapTmp,GaleriList);
            KTmp = Problem->BuildMatrix();
        } else {
            UniqueMapTmp = Galeri::Xpetra::CreateMap<LO,GO,NO>(xpetraLib,"Cartesian3D",comm,GaleriList);
            RCP<Galeri::Xpetra::Problem<Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> > > Problem = Galeri::Xpetra::BuildProblem<SC,LO,GO,Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("Laplace3D",UniqueMapTmp,GaleriList);
            KTmp = Problem->BuildMatrix();
        }

        // NodeWise dof numbering, as in test/Overlap
        const UN dpn = (UN) dofsPerNode;
        Array<GO> uniqueMapArray(dpn*UniqueMapTmp->getLocalNumElements());
        for (LO i=0; i<(LO) UniqueMapTmp->getLocalNumElements(); i++) {
            for (UN j=0; j<dpn; j++) {
                uniqueMapArray[dpn*i+j] = dpn*UniqueMapTmp->getGlobalElement(i)+j;
            }
        }
        RCP<Map<LO,GO,NO> > UniqueMap = MapFactory<LO,GO,NO>::Build(xpetraLib,INVALID,uniqueMapArray(),0,comm);
        RCP<Matrix<SC,LO,GO,NO> > K = MatrixFactory<SC,LO,GO,NO>::Build(UniqueMap,KTmp->getGlobalMaxNumRowEntries());
        for (LO i=0; i<(LO) UniqueMapTmp->getLocalNumElements(); i++) {
            ArrayView<const LO> indices;
            ArrayView<const SC> values;
            KTmp->getLocalRowView(i,indices,values);
            for (UN j=0; j<dpn; j++) {
                Array<GO> indicesArray(indices.size());
                for (LO k=0; k<indices.size(); k++) {
                    indicesArray[k] = dpn*KTmp->getColMap()->getGlobalElement(indices[k])+j;
                }
                K->insertGlobalValues(dpn*KTmp->getRowMap()->getGlobalElement(i)+j,indicesArray(),values);
            }
        }
        K->fillComplete();

        ConstXGraphPtr graph = K->getCrsGraph();
        ConstXMapPtr repeatedMap = FROSch::BuildRepeatedMap<LO,GO,NO>(graph);

        std::ostringstream label;
        label << "galeri Laplace" << dimension << "D DPN " << dofsPerNode;
        int ks[] = {1,2};
        bool ok = runGraphCase(label.str()+", repeated start map",graph,repeatedMap,std::vector<int>(ks,ks+2));
        ok = runGraphCase(label.str()+", unique start map",graph,graph->getRowMap(),std::vector<int>(ks,ks+2)) && ok;
        return ok;
    }

}

int main(int argc, char *argv[])
{
    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > CommWorld = DefaultPlatform::getDefaultPlatform().getComm();

    CommandLineProcessor My_CLP;

    RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();

    string testCase = "chain";
    My_CLP.setOption("case",&testCase,"Test case: chain, grid, empty, disconnected, nodiag, gapped, rehash, galeri.");
    int Dimension = 2;
    My_CLP.setOption("DIM",&Dimension,"Dimension (galeri case).");
    int DofsPerNode = 1;
    My_CLP.setOption("DPN",&DofsPerNode,"Dofs per node (galeri case).");
    int M = 5;
    My_CLP.setOption("M",&M,"H / h (galeri case).");

    My_CLP.recogniseAllOptions(true);
    My_CLP.throwExceptions(false);
    CommandLineProcessor::EParseCommandLineReturn parseReturn = My_CLP.parse(argc,argv);
    if (parseReturn == CommandLineProcessor::PARSE_HELP_PRINTED) {
        return(EXIT_SUCCESS);
    }

    CommWorld->barrier();
    RCP<StackedTimer> stackedTimer = rcp(new StackedTimer("Overlap Layers Test"));
    TimeMonitor::setStackedTimer(stackedTimer);

    if (CommWorld->getRank()==0) cout << "###############################\n# Overlap Layers Test: " << testCase << "\n###############################" << endl;

    bool localOk = true;
    if (!testCase.compare("chain")) {
        localOk = caseChain(CommWorld);
    } else if (!testCase.compare("grid")) {
        localOk = caseGrid(CommWorld);
    } else if (!testCase.compare("empty")) {
        localOk = caseEmptyRank(CommWorld);
    } else if (!testCase.compare("disconnected")) {
        localOk = caseDisconnected(CommWorld);
    } else if (!testCase.compare("nodiag")) {
        localOk = caseNoDiagonal(CommWorld);
    } else if (!testCase.compare("gapped")) {
        localOk = caseGapped(CommWorld);
    } else if (!testCase.compare("rehash")) {
        localOk = caseRehash(CommWorld);
    } else if (!testCase.compare("galeri")) {
        localOk = caseGaleri(CommWorld,Dimension,DofsPerNode,M);
    } else {
        if (CommWorld->getRank()==0) cout << "FAILED: unknown case " << testCase << endl;
        localOk = false;
    }

    int localPass = localOk ? 1 : 0, globalPass = 0;
    reduceAll(*CommWorld,REDUCE_MIN,localPass,ptr(&globalPass));

    CommWorld->barrier();
    stackedTimer->stop("Overlap Layers Test");
    StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    stackedTimer->report(*out,CommWorld,options);

    if (CommWorld->getRank()==0) cout << (globalPass ? "End Result: TEST PASSED" : "End Result: TEST FAILED") << endl;
    return globalPass ? EXIT_SUCCESS : EXIT_FAILURE;
}
