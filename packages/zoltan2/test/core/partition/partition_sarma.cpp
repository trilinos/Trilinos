// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 *   Created by mbenlioglu on Nov 10, 2020.
 */

#include <Zoltan2_config.h>

#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

typedef Zoltan2::BasicUserTypes <zscalar_t, zlno_t, zgno_t> myTypes_t;

#ifdef HAVE_ZOLTAN2_SARMA

struct Parameters {
    std::string alg = "pal";
    sarma::Order order_type = sarma::Order::NAT;
    std::string order_str = "nat";
    Ordinal row_parts = 8, col_parts = 0;
    Value max_load = 0;
    int seed = 2147483647;
    double sparsify = 1.0;
    bool triangular = false, use_data = false;
};

const static std::vector<std::pair<Parameters, std::vector<Ordinal> > > expected = {
        {{.alg="opal"}, {0, 48, 94, 143, 175, 218, 257, 295, 332, 0, 48, 94, 143, 175, 218, 257, 295, 332}}
};

auto testFromFile(const RCP<const Teuchos::Comm<int> > &comm, int nparts, std::string &filename, bool doRemap,
                  const Parameters sarma) {
    int me = comm->getRank();

    UserInputForTests userInputForTests(testDataFilePath, filename, comm, false, true);
    RCP<tcrsMatrix_t> matrix = userInputForTests.getUITpetraCrsMatrix();
    RCP<const tcrsMatrix_t> matrixConst = rcp_const_cast<const tcrsMatrix_t>(matrix);


    xCM_tCM_t matrixAdapter(matrixConst, 0);

    #ifdef HAVE_ZOLTAN2_MPI
    double begin = MPI_Wtime();
    #endif

    Teuchos::ParameterList params("test params");
    params.set("debug_level", "basic_status");
    params.set("num_global_parts", nparts);
    params.set("algorithm", "sarma");
    if (doRemap) params.set("remap_parts", true);
    Teuchos::ParameterList &zparams = params.sublist("zoltan_parameters", false);
    zparams.set("DEBUG_LEVEL", "0");

    // Sarma params
    Teuchos::ParameterList &sparams = params.sublist("sarma_parameters", false);
    sparams.set("alg", sarma.alg);
    sparams.set("order", sarma.order_str);
    sparams.set("row_parts", sarma.row_parts);
    sparams.set("col_parts", sarma.col_parts);
    sparams.set("max_load", sarma.max_load);
    sparams.set("sparsify", sarma.sparsify);
    sparams.set("triangular", sarma.triangular);
    sparams.set("use_data", sarma.use_data);
    sparams.set("seed", sarma.seed);

    #ifdef HAVE_ZOLTAN2_MPI
    Zoltan2::PartitioningProblem<xCM_tCM_t> problem(&matrixAdapter, &params, comm);
    #else
    Zoltan2::PartitioningProblem<xCM_tCM_t> problem(&matrixAdapter, &params);
    #endif
    if (me == 0) std::cout << "Problem constructed" << std::endl;

    problem.solve();

    #ifdef HAVE_ZOLTAN2_MPI
    if (me == 0)
        std::cout << "Run-time:\t" << MPI_Wtime() - begin << std::endl;
    #endif
    return problem.getSolution();
}
#endif

int main(int argc, char *argv[]) {
    #ifdef HAVE_ZOLTAN2_SARMA
    Tpetra::ScopeGuard tscope(&argc, &argv);
    RCP<const Teuchos::Comm<int> > tcomm = Tpetra::getDefaultComm();

    Parameters params;

    int nParts = 8;
    bool doRemap = false;
    std::string filename = "USAir97";

    // Run-time options

    for (auto test : expected){

        auto zsoln = testFromFile(tcomm, nParts, filename, doRemap, test.first);

        // compare zoltan vs raw

        const int *zparts = zsoln.getPartListView();
        for (unsigned i = 0; i < test.second.size(); ++i) {
            if (zparts[i] != (int) test.second[i]) {
                std::cout << "FAIL" << std::endl;
                return EXIT_FAILURE;
            }
        }
    }
    #endif
    std::cout << "PASS" << std::endl;

    return EXIT_SUCCESS;
}
