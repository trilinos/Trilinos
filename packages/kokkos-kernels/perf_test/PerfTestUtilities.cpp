// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
//
// Created by Poliakoff, David Zoeller on 4/27/21.
//

#include <string>

namespace test {

std::string inputDataPath;

void set_input_data_path(const std::string& path_to_data) { inputDataPath = path_to_data; }
std::string get_input_data_path() { return inputDataPath; }
}  // namespace test
