//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
//
// Created by Poliakoff, David Zoeller on 4/27/21.
//

#ifndef KOKKOSKERNELS_PERFTESTUTILITIES_HPP
#define KOKKOSKERNELS_PERFTESTUTILITIES_HPP
#include "KokkosKernels_default_types.hpp"
#include "KokkosKernels_config.h"
#include "KokkosKernels_IOUtils.hpp"
#include <common/RunParams.hpp>
#include <common/QuickKernelBase.hpp>
#include <common/KernelBase.hpp>
#include <dirent.h>

namespace test {
void set_input_data_path(const std::string &path_to_data);

std::string get_input_data_path();

}  // namespace test

namespace KokkosSparse {

template <class Scalar, class Ordinal, class ExecutionSpace, class, class Offset>
class CrsMatrix;
}

// helper function for get_directories
inline bool isDirectory(std::string path) {
  DIR *dirp;  // Pointer to a directory
  dirp = opendir(path.c_str());
  // bool var indicating that dirp is not NULL, i.e., a true statement
  bool isDir = dirp != NULL;
  if (dirp != NULL) closedir(dirp);
  return isDir;
}

/** TODO: fix in C++17 */
inline std::vector<std::string> get_directories(std::string path) {
  DIR *d;
  std::vector<std::string> paths;
  struct dirent *dir;
  d = opendir(path.c_str());
  if (d) {
    while ((dir = readdir(d)) != NULL) {
      std::string nname = std::string(dir->d_name);
      // Check to see if item is a directory
      // if (isDirectory(path + '/' + nname))
      if (nname != "." && nname != ".." && isDirectory(path + '/' + dir->d_name))
        // std::vector::emplace_back: insert a new element to the end of vector
        paths.emplace_back(dir->d_name);
    }
    closedir(d);
  }
  return paths;
}

namespace readers {

template <class Scalar, class Ordinal, class ExecutionSpace, class Offset>
using matrix_type = KokkosSparse::CrsMatrix<Scalar, Ordinal, ExecutionSpace, void, Offset>;

template <class>
struct test_reader;

template <class Scalar, class Ordinal, class ExecutionSpace, class Offset>
struct test_reader<matrix_type<Scalar, Ordinal, ExecutionSpace, Offset>> {
  static matrix_type<Scalar, Ordinal, ExecutionSpace, Offset> read(const std::string &filename) {
    return KokkosKernels::Impl::read_kokkos_crst_matrix<matrix_type<Scalar, Ordinal, ExecutionSpace, Offset>>(
        filename.c_str());
  }
};

}  // namespace readers
template <class... SubComponents>
struct data_retriever {
  std::string root_path;
  std::string sub_path;
  struct test_case {
    std::string filename;
    std::tuple<SubComponents...> test_data;
  };
  std::vector<test_case> test_cases;
  std::string make_full_path_to_data_file(std::string repo, std::string path_to_data, std::string dataset,
                                          std::string filename) {
    return root_path + "/" + repo + "/" + path_to_data + dataset + "/" + filename;
  }
  template <class... Locations>
  data_retriever(std::string path_to_data, Locations... locations) : sub_path(path_to_data) {
    root_path = test::get_input_data_path();

    // TODO: way to list the directories in the root path
    std::vector<std::string> data_repos = get_directories(root_path + "/");
    // TODO: list directories in subpaths
    for (auto repo : data_repos) {
      std::vector<std::string> datasets = get_directories(root_path + "/" + repo + "/" + path_to_data + "/");
      for (auto dataset : datasets) {
        test_cases.push_back(test_case{repo + "/" + dataset,
                                       std::make_tuple(readers::test_reader<SubComponents>::read(
                                           make_full_path_to_data_file(repo, path_to_data, dataset, locations))...)});
      }
    }
  }
};

using test_list = std::vector<rajaperf::KernelBase *>;

#endif  // KOKKOSKERNELS_PERFTESTUTILITIES_HPP
