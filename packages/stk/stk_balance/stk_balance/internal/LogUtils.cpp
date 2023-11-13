// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "LogUtils.hpp"
#include <Kokkos_Core.hpp>
#include "stk_util/environment/Env.hpp"
#include <stk_util/diag/PrintTable.hpp>
#include <stk_util/registry/ProductRegistry.hpp>
#include <stk_util/diag/Platform.hpp>
#include <stk_util/registry/RegisterProduct.hpp>
#include "stk_util/environment/EnvData.hpp"
#include "stk_util/parallel/Parallel.hpp"
#include "Trilinos_version.h"

namespace stk {
namespace balance {

void initialize_environment(MPI_Comm communicator, const char** argv)
{
  sierra::Env::set_mpi_communicator(communicator);

  stk::EnvData& env_data = stk::EnvData::instance();
  env_data.m_executablePath = argv[0];

  size_t pathSeparator = env_data.m_executablePath.find_last_of("\\/");
  if (pathSeparator != std::string::npos) {
    env_data.m_productName = env_data.m_executablePath.substr(pathSeparator + 1);
  }
  else {
    env_data.m_productName = env_data.m_executablePath;
  }

  ProductRegistry::instance().setProductName(env_data.m_productName);

  ProductRegistry::AttributeMap& product_attributes = ProductRegistry::instance().getProductAttributeMap(env_data.m_productName);
  product_attributes[ProductRegistry::BUILD_TIME]   = __DATE__ " " __TIME__;
  product_attributes[ProductRegistry::VERSION]      = stk::ProductRegistry::version();
  product_attributes[ProductRegistry::TITLE]        = env_data.m_productName;

  stk::register_product();

  ProductRegistry::AttributeMap& attr_map = ProductRegistry::instance().addProduct(sierra::Env::osname().c_str());
  attr_map[ProductRegistry::VERSION] = sierra::Env::osversion().c_str();

#if TRILINOS_MAJOR_MINOR_VERSION == 0
  std::string version = "Dev";
#else
  std::ostringstream os;
  os << TRILINOS_MAJOR_MINOR_VERSION / 10000 << "." << TRILINOS_MAJOR_MINOR_VERSION / 100 % 100 << "."
     << TRILINOS_MAJOR_MINOR_VERSION % 100;
  std::string version(os.str());
#endif

  ProductRegistry::instance().addTPL("Trilinos", version);
}

void print_banner(std::ostream & outputStream)
{
  stk::PrintTable banner_environment;
  stk::PrintTable banner_support;

  banner_environment.setAutoEndCol(false);
  banner_support.setAutoEndCol(false);

  banner_environment << justify(stk::PrintTable::Cell::RIGHT) << "Directory:" << end_col
                     << justify(stk::PrintTable::Cell::LEFT) << sierra::Env::working_directory() << end_col << end_row
                     << justify(stk::PrintTable::Cell::RIGHT) << "Executable:" << end_col
                     << justify(stk::PrintTable::Cell::LEFT) << sierra::Env::executable_file() << end_col << end_row
                     << justify(stk::PrintTable::Cell::RIGHT) << "Built:" << end_col
                     << justify(stk::PrintTable::Cell::LEFT) << ProductRegistry::executable_date() << end_col << end_row
                     << justify(stk::PrintTable::Cell::RIGHT) << "Run Started:" << end_col
                     << justify(stk::PrintTable::Cell::LEFT) << sierra::format_time(sierra::Env::start_time()) << end_col << end_row
                     << justify(stk::PrintTable::Cell::RIGHT) << "User:" << end_col
                     << justify(stk::PrintTable::Cell::LEFT) << sierra::Env::username() << end_col << end_row
                     << justify(stk::PrintTable::Cell::RIGHT) << "Architecture:" << end_col
                     << justify(stk::PrintTable::Cell::LEFT) << (sierra::Env::architecture().empty() ? sierra::Env::hostname()
                                                                                                     : sierra::Env::architecture().c_str())
                                                             << end_col << end_row
                     << justify(stk::PrintTable::Cell::RIGHT) << "Host:" << end_col
                     << justify(stk::PrintTable::Cell::LEFT) << sierra::Env::hostname() << end_col << end_row
                     << justify(stk::PrintTable::Cell::RIGHT) << "Hardware:" << end_col
                     << justify(stk::PrintTable::Cell::LEFT) << sierra::Env::hardware() << end_col << end_row
                     << justify(stk::PrintTable::Cell::RIGHT) << "Operating System:" << end_col
                     << justify(stk::PrintTable::Cell::LEFT) << sierra::Env::osname() << end_col << end_row
                     << justify(stk::PrintTable::Cell::RIGHT) << "Processors:" << end_col
                     << justify(stk::PrintTable::Cell::LEFT) << sierra::Env::parallel_size() << end_col << end_row;

#ifdef KOKKOS_ENABLE_CUDA
  banner_environment << justify(stk::PrintTable::Cell::RIGHT) << "Accelerator:" << end_col
                     << justify(stk::PrintTable::Cell::LEFT) << "May be running on Nvidia GPUs" << end_col << end_row;
#endif

#ifdef KOKKOS_ENABLE_HIP
  banner_environment << justify(stk::PrintTable::Cell::RIGHT) << "Accelerator:" << end_col
                     << justify(stk::PrintTable::Cell::LEFT) << "May be running on AMD GPUs" << end_col << end_row;
#endif

#ifdef KOKKOS_ENABLE_OPENMP
  banner_environment << justify(stk::PrintTable::Cell::RIGHT) << "Accelerator:" << end_col
                     << justify(stk::PrintTable::Cell::LEFT) << "Using " << Kokkos::DefaultExecutionSpace().concurrency()
                        << " OpenMP Threads" << end_col << end_row;
//                     << justify(stk::PrintTable::Cell::LEFT) << "Using " << omp_get_max_threads() << " OpenMP Threads" << end_col << end_row;
#endif

  banner_support << justify(stk::PrintTable::Cell::CENTER) << "Product" << end_col
                 << justify(stk::PrintTable::Cell::CENTER) << "Version" << end_col << end_header;

  for (auto& product : ProductRegistry::instance().productMap()) {
    const std::string& name = product.first;
    ProductRegistry::AttributeMap& attribute_map = product.second;
    banner_support << justify(stk::PrintTable::Cell::RIGHT) << name << end_col;

    ProductRegistry::AttributeMap::iterator ver_it = attribute_map.find(ProductRegistry::VERSION);
    if (ver_it != attribute_map.end() && !(*ver_it).second.empty()) {
      banner_support << justify(stk::PrintTable::Cell::LEFT) << (*ver_it).second;
    }
    banner_support << end_col << end_row;
  }

  outputStream << std::endl << banner_environment << banner_support
               << "------------------------------------------------------------------------------------------"
               << std::endl;

}

}
}
