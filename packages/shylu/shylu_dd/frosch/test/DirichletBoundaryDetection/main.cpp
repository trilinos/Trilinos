// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Solve a Q1 Laplace problem with homogeneous Dirichlet boundary conditions
// on four square subdomains. The entire domain is [0,8] x [0,8].
// The bottom-left subdomain is [0,4] x [0,4].
// After one GMRES iteration from zero, check every DOF at the shared
// boundary node (4, 0).
// Boundary detection removes the node from the coarse basis support, so the
// value stays zero when detection is enabled. With detection disabled, the
// coarse correction makes it nonzero.
// As a second test, reset the solution to zero and verify convergence in a
// normal solve.
// Two fundamentally different tests are implemented:
// * scalar, as described above
// * two blocks, and each block contains discretized Laplacians with two dofs per node
// The block case has essentially four independent block matrices.

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Tpetra_Core.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ThyraUtils.hpp>

#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>
#include <Stratimikos_FROSch_def.hpp>
#include <FROSch_Tools_def.hpp>
#include <FROSch_Types.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

using Scalar = double;
using LO = int;
using GO = FROSch::DefaultGlobalOrdinal;
using Node = Tpetra::Map<>::node_type;
using Map = Xpetra::Map<LO, GO, Node>;
using MapFactory = Xpetra::MapFactory<LO, GO, Node>;
using MatrixFactory = Xpetra::MatrixFactory<Scalar, LO, GO, Node>;
using MultiVectorFactory = Xpetra::MultiVectorFactory<Scalar, LO, GO, Node>;
using ThyraUtils = Xpetra::ThyraUtils<Scalar, LO, GO, Node>;

constexpr int subdomainWidth = 4; // 4 x 4 elements per subdomain
constexpr int nodesPerSide = 2 * subdomainWidth + 1;
constexpr int nodesPerBlock = nodesPerSide * nodesPerSide;

bool isBoundary(const int x, const int y)
{
    // The entire boundary of the square is used as a Dirichlet boundary.
    return x == 0 || y == 0 || x == nodesPerSide - 1 || y == nodesPerSide - 1;
}

} // namespace

int main(int argc, char* argv[])
{
    Tpetra::ScopeGuard scope(&argc, &argv);
    const auto comm = Tpetra::getDefaultComm();
    const int rank = comm->getRank();

    // Read method/solver to be used and if boundary detection should be en- or disabled.
    std::string method;
    std::string detection;
    Teuchos::CommandLineProcessor options(false);
    options.setOption("method", &method, "FROSch preconditioner to test", true);
    options.setOption("detection", &detection, "Boundary detection setting (enabled or disabled)", true);
    const auto parseResult = options.parse(argc, argv);
    if (parseResult != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
        return parseResult == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    // Make sure a solver was selected.
    const bool use_block_prec = (method == "TwoLevelBlockPreconditioner_GDSW");
    const bool use_nonblock_prec = (method == "TwoLevelPreconditioner_GDSW");
    const bool use_rgdsw_prec = (method == "RGDSWPreconditioner");
    const bool use_gdsw_prec = (method == "GDSWPreconditioner");
    const bool detectBoundary = detection == "enabled";
    if ((!use_block_prec && !use_nonblock_prec && !use_rgdsw_prec && !use_gdsw_prec) ||
            (!detectBoundary && detection != "disabled") ||
            (comm->getSize() != 4)) {
        if (rank == 0) {
            std::cerr << "Expected a supported method, enabled or disabled detection, and four MPI ranks.\n";
        }
        return EXIT_FAILURE;
    }

    // Assign IDs to mesh nodes.
    const int numBlocks = use_block_prec ? 2 : 1;
    const int dofsPerNode = use_block_prec ? 2 : 1;
    const int dofsPerBlock = nodesPerBlock * dofsPerNode;
    const int xStart = (rank % 2) * subdomainWidth; // start of x coordinate for current MPI rank
    const int yStart = (rank / 2) * subdomainWidth; // start of y coordinate for current MPI rank
    Teuchos::Array<GO> repeatedIds;
    Teuchos::Array<GO> ownedIds;
    for (int y = yStart; y <= yStart + subdomainWidth; ++y) {
        for (int x = xStart; x <= xStart + subdomainWidth; ++x) {
            const GO id = y * nodesPerSide + x;
            repeatedIds.push_back(id);
            // Interface nodes belong to the lowest-numbered adjacent rank.
            if ((xStart == 0 || x > xStart) && (yStart == 0 || y > yStart)) {
                ownedIds.push_back(id);
            }
        }
    }

    // Build repeated node map
    const auto invalid = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
    Teuchos::RCP<const Map> repeatedNodeMap = MapFactory::Build(Xpetra::UseTpetra, invalid, repeatedIds(), 0, comm);

    // Build repeated dof map
    // Block case: Both blocks use the same node-wise DOF numbering within the block.
    //             So the repeated dof map can be used for both blocks.
    const auto repeatedDofMap = FROSch::BuildMapFromNodeMap(repeatedNodeMap, dofsPerNode, FROSch::NodeWise);

    // Build a repeated map for the entire system.
    Teuchos::Array<GO> systemIds;
    for (int block = 0; block < numBlocks; ++block) {
        for (const GO id : ownedIds) {
            for (int component = 0; component < dofsPerNode; ++component) {
                systemIds.push_back(block * dofsPerBlock + dofsPerNode * id + component);
            }
        }
    }
    const auto map = MapFactory::Build(Xpetra::UseTpetra, numBlocks * dofsPerBlock, systemIds(), 0, comm);

    // Assemble matrix, right-hand side and initialize solution.
    const auto matrix = MatrixFactory::Build(map, 9);
    const auto rhs = MultiVectorFactory::Build(map, 1);
    const auto solution = MultiVectorFactory::Build(map, 1);
    solution->putScalar(0.0);
    auto rhsView = rhs->getDataNonConst(0);
    for (LO localRow = 0; localRow < systemIds.size(); ++localRow) {
        const GO row = systemIds[localRow];
        const GO node = (row % dofsPerBlock) / dofsPerNode;
        const GO offset = row - dofsPerNode * node; // block offset
        const int x = node % nodesPerSide;          // x-coordinate within the domain
        const int y = node / nodesPerSide;          // y-coordinate within the domain
        const bool boundary = isBoundary(x, y);
        rhsView[localRow] = boundary ? 0.0 : 1.0; // Set Dirichlet boundary rows to unit vectors

        Teuchos::Array<GO> columns;
        Teuchos::Array<Scalar> values;
        // The assembled Q1 stencil on a uniform square mesh has eight equal off-diagonals.
        for (int ny = std::max(0, y - 1); ny <= std::min(nodesPerSide - 1, y + 1); ++ny) {
            for (int nx = std::max(0, x - 1); nx <= std::min(nodesPerSide - 1, x + 1); ++nx) {
                const GO column = offset + dofsPerNode * (ny * nodesPerSide + nx);
                columns.push_back(column);
                if (row == column) {
                    // Entries on diagonal
                    values.push_back(boundary ? 1.0 : 8.0 / 3.0);
                } else {
                    // Entries on off-diagonals
                    // Check whether row or column belongs to the Dirichlet boundary.
                    // Store zero entries explicitly to test inspection of boundary detection via matrix values
                    values.push_back(boundary || isBoundary(nx, ny) ? 0.0 : -1.0 / 3.0);
                }
            }
        }
        matrix->insertGlobalValues(row, columns(), values());
    }
    rhsView = Teuchos::null;
    matrix->fillComplete();

    // Get and set parameters
    const auto parameters = Teuchos::getParametersFromXmlFile("Parameters.xml");
    const int fullIterationLimit = parameters->sublist("Linear Solver Types")
                                           .sublist("Belos")
                                           .sublist("Solver Types")
                                           .sublist("Block GMRES")
                                           .get<int>("Maximum Iterations");
    auto &paramListFrosch = parameters->sublist("Preconditioner Types").sublist("FROSch");
    const std::string preconditioner =
            use_nonblock_prec ? "TwoLevelPreconditioner" : (use_block_prec ? "TwoLevelBlockPreconditioner" : method);
    paramListFrosch.set("FROSch Preconditioner Type", preconditioner);
    if (use_nonblock_prec || use_block_prec) paramListFrosch.set("CoarseOperator Type", "GDSWCoarseOperator");
    paramListFrosch.sublist("GDSWCoarseOperator")
            .sublist("Dirichlet boundary rows in global matrix")
            .set("Find Dirichlet boundary rows", detectBoundary);
    if (use_rgdsw_prec) {
        paramListFrosch.sublist("RGDSWCoarseOperator").setParameters(paramListFrosch.sublist("GDSWCoarseOperator"));
    }

    // Pass parameters and maps to FROSch.
    if (use_block_prec) {
        // Pass the repeated dof map "num blocks" times to FROSch.
        paramListFrosch.set("Repeated Map Vector", Teuchos::ArrayRCP<Teuchos::RCP<Map>>(numBlocks, repeatedDofMap));

        paramListFrosch.set("DofsPerNode Vector", Teuchos::ArrayRCP<unsigned>(numBlocks, dofsPerNode));
        paramListFrosch.set("DofOrdering Vector", Teuchos::ArrayRCP<FROSch::DofOrdering>(numBlocks, FROSch::NodeWise));
    } else {
        paramListFrosch.set("Repeated Map", repeatedDofMap);
    }

    // Thyra and Stratimikos boilerplate
    using CrsMatrixWrap = Xpetra::CrsMatrixWrap<Scalar, LO, GO, Node>;
    const auto operatorThyra =
            ThyraUtils::toThyra(Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(matrix, true)->getCrsMatrix());
    const auto rhsThyra = ThyraUtils::toThyraMultiVector(rhs);
    const auto solutionThyra =
            Teuchos::rcp_const_cast<Thyra::MultiVectorBase<Scalar>>(ThyraUtils::toThyraMultiVector(solution));
    Stratimikos::LinearSolverBuilder<Scalar> builder;
    Stratimikos::enableFROSch<LO, GO, Node>(builder);
    builder.setParameterList(parameters);
    const auto solverFactory = builder.createLinearSolveStrategy("");
    const auto solver = Thyra::linearOpWithSolve<Scalar>(*solverFactory, operatorThyra);

    // First test: Do only one iteration with Belos and then check the result.
    Thyra::SolveCriteria<Scalar> criteria;
    criteria.extraParameters = Teuchos::parameterList();
    criteria.extraParameters->set("Maximum Iterations", 1);
    const auto probe = solver->solve(Thyra::NOTRANS, *rhsThyra, solutionThyra.ptr(), Teuchos::ptr(&criteria));
    constexpr const char* iterationKey = "Belos/Iteration Count";
    const int probeIterations = !probe.extraParameters.is_null() && probe.extraParameters->isType<int>(iterationKey)
                                ? probe.extraParameters->get<int>(iterationKey) : -1;
    bool probePassed = (probeIterations == 1);
    if (!probePassed) {
        std::cerr << method << ", detection " << detection << ", rank " << rank
                  << ", expected one probe iteration, observed " << probeIterations << '\n';
    }

    // Check the values at the node at the very right side of the bottom-left subdomain.
    // They should be zero if boundary detection is enabled and not zero otherwise.
    {
        constexpr Scalar zeroTolerance = 1e-12;
        constexpr Scalar nonzeroMinimum = 1e-3;
        const auto values = solution->getData(0);
        for (int block = 0; block < numBlocks; ++block) {
            for (int component = 0; component < dofsPerNode; ++component) {
                // dof at position (x,y) = (4,0)
                const GO dof = block * dofsPerBlock + dofsPerNode * subdomainWidth + component;
                const LO localDof = map->getLocalElement(dof);

                // Check if the dof belongs to the current rank.
                if (localDof == Teuchos::OrdinalTraits<LO>::invalid()) {
                    continue;
                }

                // Verify that the value is as expected.
                const Scalar value = values[localDof];
                const Scalar magnitude = std::abs(value);
                const bool boundaryPassed = std::isfinite(value) &&
                    (detectBoundary ? magnitude <= zeroTolerance : magnitude >= nonzeroMinimum);
                std::cout << method << ", detection " << detection << ", block " << block + 1
                          << ", component " << component + 1
                          << ", node (" << subdomainWidth << ", 0) after one iteration = " << value
                          << ", expected magnitude " << (detectBoundary ? "<= " : ">= ")
                          << (detectBoundary ? zeroTolerance : nonzeroMinimum) << '\n';
                probePassed = probePassed && boundaryPassed;
            }
        }
    }

    // Test 2: Restart the iteration and let it run to verify it converges.
    solution->putScalar(0.0);
    criteria.extraParameters->set("Maximum Iterations", fullIterationLimit);
    const auto status = solver->solve(Thyra::NOTRANS, *rhsThyra, solutionThyra.ptr(), Teuchos::ptr(&criteria));
    const bool converged = status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED;
    const int passed = probePassed && converged;
    if (rank == 0 || !converged) {
        std::cout << method << ", detection " << detection << ", rank " << rank
                  << ", full solve " << (converged ? "converged" : "did not converge") << '\n';
    }
    int allPassed = 0;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &passed, &allPassed);
    return allPassed ? EXIT_SUCCESS : EXIT_FAILURE;
}
