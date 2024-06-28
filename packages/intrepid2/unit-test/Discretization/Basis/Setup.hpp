// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef INTREPID2_UNIT_TEST_DISCRETIZATION_SETUP_HPP
#define INTREPID2_UNIT_TEST_DISCRETIZATION_SETUP_HPP

#include <iostream>

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid2::Test
{

/**
 * @brief Setup test output stream.
 *
 * It will also print device and host configurations, as well as a test header
 * based on @p test_name.
 */
template <typename DeviceType>
auto setup_output_stream(const bool verbose, const std::string& test_name, const std::vector<std::string>& test_content)
{
    Teuchos::RCP<std::ostream> out_stream;

    if (verbose)
        out_stream = Teuchos::rcp(&std::cout, false);
    else
        out_stream = Teuchos::make_rcp<Teuchos::oblackholestream>();

    using      execution_space = typename DeviceType::execution_space;
    using host_execution_space = Kokkos::DefaultHostExecutionSpace;

    *out_stream << "DeviceSpace::  ";      execution_space().print_configuration(*out_stream, false);
    *out_stream << "HostSpace::    "; host_execution_space().print_configuration(*out_stream, false);

    const std::string sep   = "================================================================================";
    const std::string empty = "|                                                                              |";
    const std::string unitt = "|                 Unit Test (";

    *out_stream
        << sep   << std::endl
        << empty << std::endl
        << unitt << test_name << ")" << std::string(sep.size() - unitt.size() - test_name.size() - 2, ' ') << "|" << std::endl
        << empty << std::endl;

    for(const auto& step : test_content)
    {
        *out_stream << "|     " << step << std::string(sep.size() - step.size() - 7, ' ') << "|" << std::endl;
    }

    *out_stream
        << empty << std::endl
        << "|  Questions? Contact  Pavel Bochev  (pbboche@sandia.gov),                     |\n"
        << "|                      Robert Kirby  (robert.c.kirby@ttu.edu),                 |\n"
        << "|                      Denis Ridzal  (dridzal@sandia.gov),                     |\n"
        << "|                      Kara Peterson (dridzal@sandia.gov),                     |\n"
        << "|                      Mauro Perego  (mperego@sandia.gov),                     |\n"
        << "|                      Kyungjoo Kim  (kyukim@sandia.gov).                      |\n"
        << empty << std::endl
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid            |\n"
        << "|  Trilinos website:   http://trilinos.sandia.gov                              |\n"
        << empty << std::endl
        << sep   << std::endl;

    return out_stream;
}

} // namespace Intrepid2::Test

#endif // INTREPID2_UNIT_TEST_DISCRETIZATION_SETUP_HPP
