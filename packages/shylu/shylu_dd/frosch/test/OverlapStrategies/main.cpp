// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Test for the "Adding Layers Strategy" = "Threaded" option of
// FROSch::AlgebraicOverlappingOperator. The operator is built with the
// strategies "CrsGraph" (reference), "CrsMatrix" and "Threaded" on the same
// Galeri Laplace matrix, for "Sort Overlapping Map" = true and false, and the
// overlapping maps and the results of apply() are compared.

#include <ShyLU_DDFROSch_config.h>

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_StackedTimer.hpp>

// Galeri::Xpetra
#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"
#include "Galeri_XpetraParameters.hpp"
#include "Galeri_XpetraUtils.hpp"
#include "Galeri_XpetraMaps.hpp"

// Xpetra include
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MatrixFactory.hpp>

#include <FROSch_AlgebraicOverlappingOperator_def.hpp>
#include <FROSch_Tools_def.hpp>


using UN    = unsigned;
using SC    = double;
using LO    = int;
using GO    = FROSch::DefaultGlobalOrdinal;
using NO    = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType;

using namespace std;
using namespace Teuchos;
using namespace Xpetra;

namespace {

    const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();

    // Test-only subclass that exposes the protected overlapping map
    class ExposedAlgebraicOverlappingOperator : public FROSch::AlgebraicOverlappingOperator<SC,LO,GO,NO> {
    public:
        ExposedAlgebraicOverlappingOperator(RCP<const Matrix<SC,LO,GO,NO> > k,
                                            RCP<ParameterList> parameterList) :
        FROSch::AlgebraicOverlappingOperator<SC,LO,GO,NO> (k,parameterList)
        {
        }

        RCP<const Map<LO,GO,NO> > getOverlappingMapForTest() const
        {
            return this->OverlappingMap_;
        }
    };

    struct StrategyResult {
        std::vector<GO> gids;
        RCP<MultiVector<SC,LO,GO,NO> > y;
    };

    std::vector<GO> toVector(ArrayView<const GO> view)
    {
        std::vector<GO> vec(view.size());
        for (size_t i=0; i<(size_t) view.size(); i++) vec[i] = view[i];
        return vec;
    }

    bool strictlyAscending(const std::vector<GO> &v)
    {
        for (size_t i=1; i<v.size(); i++) {
            if (!(v[i-1]<v[i])) return false;
        }
        return true;
    }

    StrategyResult runStrategy(RCP<const Matrix<SC,LO,GO,NO> > K,
                               RCP<const Map<LO,GO,NO> > repeatedMap,
                               int overlap,
                               const string &strategy,
                               bool sortOverlappingMap,
                               RCP<const MultiVector<SC,LO,GO,NO> > x)
    {
        RCP<ParameterList> parameterList = rcp(new ParameterList("AlgebraicOverlappingOperator"));
        parameterList->set("Adding Layers Strategy",strategy);
        parameterList->set("Sort Overlapping Map",sortOverlappingMap);
        parameterList->sublist("Solver").set("SolverType",string("Amesos2"));
        parameterList->sublist("Solver").set("Solver",string("Klu"));

        RCP<ExposedAlgebraicOverlappingOperator> op = rcp(new ExposedAlgebraicOverlappingOperator(K,parameterList));
        op->initialize(overlap,repeatedMap);
        op->compute();

        StrategyResult result;
        result.gids = toVector(op->getOverlappingMapForTest()->getLocalElementList());
        result.y = MultiVectorFactory<SC,LO,GO,NO>::Build(K->getRangeMap(),1);
        op->apply(*x,*result.y,true);
        return result;
    }

}

int main(int argc, char *argv[])
{
    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > CommWorld = DefaultPlatform::getDefaultPlatform().getComm();

    CommandLineProcessor My_CLP;

    RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();

    int M = 5;
    My_CLP.setOption("M",&M,"H / h.");
    int Dimension = 2;
    My_CLP.setOption("DIM",&Dimension,"Dimension.");
    int Overlap = 1;
    My_CLP.setOption("O",&Overlap,"Overlap.");
    int DofsPerNode = 1;
    My_CLP.setOption("DPN",&DofsPerNode,"Dofs per node.");
    double Tolerance = 1.e-12;
    My_CLP.setOption("TOL",&Tolerance,"Tolerance for the relative difference of apply().");

    My_CLP.recogniseAllOptions(true);
    My_CLP.throwExceptions(false);
    CommandLineProcessor::EParseCommandLineReturn parseReturn = My_CLP.parse(argc,argv);
    if (parseReturn == CommandLineProcessor::PARSE_HELP_PRINTED) {
        return(EXIT_SUCCESS);
    }

    CommWorld->barrier();
    RCP<StackedTimer> stackedTimer = rcp(new StackedTimer("Overlap Strategies Test"));
    TimeMonitor::setStackedTimer(stackedTimer);

    const int rank = CommWorld->getRank();
    bool ok = true;

    int N = 0;
    if (Dimension == 2) {
        N = (int) (pow(CommWorld->getSize(),1/2.) + 100*numeric_limits<double>::epsilon()); // 1/H
        if (N*N!=CommWorld->getSize()) ok = false;
    } else if (Dimension == 3) {
        N = (int) (pow(CommWorld->getSize(),1/3.) + 100*numeric_limits<double>::epsilon()); // 1/H
        if (N*N*N!=CommWorld->getSize()) ok = false;
    } else {
        ok = false;
    }
    if (!ok) {
        if (rank==0) cout << "FAILED: DIM must be 2 or 3 and the number of ranks a square (2D) or cube (3D)" << endl;
        return(EXIT_FAILURE);
    }

    UnderlyingLib xpetraLib = UseTpetra;

    CommWorld->barrier(); if (rank==0) cout << "############\n# Assembly #\n############\n" << endl;

    ParameterList GaleriList;
    GaleriList.set("nx", GO(N*M));
    GaleriList.set("ny", GO(N*M));
    GaleriList.set("nz", GO(N*M));
    GaleriList.set("mx", GO(N));
    GaleriList.set("my", GO(N));
    GaleriList.set("mz", GO(N));

    RCP<const Map<LO,GO,NO> > UniqueMapTmp;
    RCP<Matrix<SC,LO,GO,NO> > KTmp;
    if (Dimension==2) {
        UniqueMapTmp = Galeri::Xpetra::CreateMap<LO,GO,NO>(xpetraLib,"Cartesian2D",CommWorld,GaleriList);
        RCP<Galeri::Xpetra::Problem<Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> > > Problem = Galeri::Xpetra::BuildProblem<SC,LO,GO,Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("Laplace2D",UniqueMapTmp,GaleriList);
        KTmp = Problem->BuildMatrix();
    } else {
        UniqueMapTmp = Galeri::Xpetra::CreateMap<LO,GO,NO>(xpetraLib,"Cartesian3D",CommWorld,GaleriList);
        RCP<Galeri::Xpetra::Problem<Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> > > Problem = Galeri::Xpetra::BuildProblem<SC,LO,GO,Map<LO,GO,NO>,CrsMatrixWrap<SC,LO,GO,NO>,MultiVector<SC,LO,GO,NO> >("Laplace3D",UniqueMapTmp,GaleriList);
        KTmp = Problem->BuildMatrix();
    }

    // NodeWise dof numbering, as in test/Overlap
    const UN dpn = (UN) DofsPerNode;
    Array<GO> uniqueMapArray(dpn*UniqueMapTmp->getLocalNumElements());
    for (LO i=0; i<(LO) UniqueMapTmp->getLocalNumElements(); i++) {
        for (UN j=0; j<dpn; j++) {
            uniqueMapArray[dpn*i+j] = dpn*UniqueMapTmp->getGlobalElement(i)+j;
        }
    }
    RCP<Map<LO,GO,NO> > UniqueMap = MapFactory<LO,GO,NO>::Build(xpetraLib,INVALID,uniqueMapArray(),0,CommWorld);
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

    RCP<const Map<LO,GO,NO> > RepeatedMap = FROSch::BuildRepeatedMap<LO,GO,NO>(K->getCrsGraph());

    // Input vector: deterministic pseudo-random values from a hash of the GID (independent of the partition)
    RCP<MultiVector<SC,LO,GO,NO> > x = MultiVectorFactory<SC,LO,GO,NO>::Build(K->getRangeMap(),1);
    {
        ArrayRCP<SC> data = x->getDataNonConst(0);
        for (LO i=0; i<(LO) K->getRangeMap()->getLocalNumElements(); i++) {
            unsigned long long h = (unsigned long long) K->getRangeMap()->getGlobalElement(i);
            h = h*6364136223846793005ULL + 1442695040888963407ULL;
            h ^= (h >> 33);
            data[i] = SC((h % 1000003ULL)/1000003.0 - 0.5);
        }
    }

    const string reference = "CrsGraph";
    const string others[] = {"CrsMatrix","Threaded"};
    const bool sortModes[] = {true,false};

    for (int s=0; s<2; s++) {
        const bool sortMap = sortModes[s];
        CommWorld->barrier(); if (rank==0) cout << "\n####################################\n# Sort Overlapping Map = " << (sortMap ? "true " : "false") << "     #\n####################################" << endl;

        StrategyResult ref = runStrategy(K,RepeatedMap,Overlap,reference,sortMap,x.getConst());
        Array<SC> refNorm(1);
        ref.y->norm2(refNorm());
        if (!(refNorm[0]>0.0)) {
            ok = false;
            if (rank==0) cout << "FAILED [sort " << sortMap << "]: reference apply() returned a zero vector" << endl;
        }
        if (sortMap && !strictlyAscending(ref.gids)) {
            ok = false;
            cout << "FAILED [sort " << sortMap << ", CrsGraph] on rank " << rank << ": reference map not sorted" << endl;
        }
        std::vector<GO> refSorted = ref.gids;
        std::sort(refSorted.begin(),refSorted.end());

        for (int o=0; o<2; o++) {
            const string &strategy = others[o];
            StrategyResult res = runStrategy(K,RepeatedMap,Overlap,strategy,sortMap,x.getConst());

            bool localOk = true;
            std::ostringstream reason;
            if (sortMap) {
                if (res.gids!=ref.gids) {
                    localOk = false; reason << " overlapping map differs from CrsGraph (local sizes " << res.gids.size() << " vs " << ref.gids.size() << ");";
                }
            } else {
                std::vector<GO> resSorted = res.gids;
                std::sort(resSorted.begin(),resSorted.end());
                if (resSorted!=refSorted) {
                    localOk = false; reason << " overlapping GID set differs from CrsGraph (local sizes " << res.gids.size() << " vs " << ref.gids.size() << ");";
                }
                if (!strategy.compare("Threaded") && !strictlyAscending(res.gids)) {
                    localOk = false; reason << " Threaded map is not sorted although it always should be;";
                }
            }
            if (!localOk) {
                cout << "FAILED [sort " << sortMap << ", " << strategy << "] on rank " << rank << ":" << reason.str() << endl;
            }

            // apply(): relative difference to CrsGraph (collective)
            res.y->update(-1.0,*ref.y,1.0);
            Array<SC> diffNorm(1);
            res.y->norm2(diffNorm());
            const double relDiff = (refNorm[0]>0.0) ? diffNorm[0]/refNorm[0] : diffNorm[0];
            if (!(relDiff<=Tolerance)) {
                localOk = false;
                if (rank==0) cout << "FAILED [sort " << sortMap << ", " << strategy << "]: apply() relative difference " << relDiff << " > " << Tolerance << endl;
            }

            int localPass = localOk ? 1 : 0, globalPass = 0;
            reduceAll(*CommWorld,REDUCE_MIN,localPass,ptr(&globalPass));
            if (rank==0) cout << (globalPass ? "passed " : "FAILED ") << "[sort " << sortMap << ", " << strategy << " vs CrsGraph, overlap " << Overlap << "], apply() relative difference " << relDiff << endl;
            ok = ok && localOk;
        }
    }

    int localPass = ok ? 1 : 0, globalPass = 0;
    reduceAll(*CommWorld,REDUCE_MIN,localPass,ptr(&globalPass));

    CommWorld->barrier();
    stackedTimer->stop("Overlap Strategies Test");
    StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    stackedTimer->report(*out,CommWorld,options);

    if (rank==0) cout << (globalPass ? "End Result: TEST PASSED" : "End Result: TEST FAILED") << endl;
    return globalPass ? EXIT_SUCCESS : EXIT_FAILURE;
}
