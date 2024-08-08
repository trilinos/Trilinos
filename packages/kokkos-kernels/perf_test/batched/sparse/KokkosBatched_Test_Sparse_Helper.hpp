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

template <class XType>
void writeArrayToMM(std::string name, const XType x) {
  std::ofstream myfile;
  myfile.open(name);

  typename XType::HostMirror x_h = Kokkos::create_mirror_view(x);

  Kokkos::deep_copy(x_h, x);

  myfile << "%% MatrixMarket 2D Array\n%" << std::endl;
  myfile << x_h.extent(0) << " " << x_h.extent(1) << std::endl;

  for (size_t i = 0; i < x_h.extent(0); ++i) {
    for (size_t j = 0; j < x_h.extent(1); ++j) {
      myfile << std::setprecision(15) << x_h(i, j) << " ";
    }
    myfile << std::endl;
  }

  myfile.close();
}

void readSizesFromMM(std::string name, int &nrows, int &ncols, int &nnz, int &N) {
  std::ifstream input(name);
  while (input.peek() == '%') input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  std::string line_sizes;

  getline(input, line_sizes);

  std::stringstream iss(line_sizes);

  int number;
  std::vector<int> sizes;
  while (iss >> number) sizes.push_back(number);

  nrows = sizes[0];
  ncols = sizes[1];

  nnz = 0;
  N   = 0;

  if (sizes.size() >= 3) nnz = sizes[2];

  if (sizes.size() == 4) N = sizes[3];
}

template <class XType>
void readArrayFromMM(std::string name, const XType &x) {
  std::ifstream input(name);

  while (input.peek() == '%') input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  typename XType::HostMirror x_h = Kokkos::create_mirror_view(x);

  for (size_t i = 0; i < x_h.extent(0); ++i)
    for (size_t j = 0; j < x_h.extent(1); ++j) input >> x_h(i, j);

  input.close();

  Kokkos::deep_copy(x, x_h);
}

template <class AType>
void readDenseFromMM(std::string name, const AType &A) {
  std::ifstream input(name);

  while (input.peek() == '%') input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  typename AType::HostMirror A_h = Kokkos::create_mirror_view(A);

  Kokkos::deep_copy(A_h, 0.);

  int read_row;
  int read_col;

  int N, Blk, nnz, nrows;
  readSizesFromMM(name, Blk, nrows, nnz, N);

  for (int i = 0; i < nnz; ++i) {
    input >> read_row >> read_col;
    --read_row;
    --read_col;

    for (int j = 0; j < N; ++j) input >> A_h(j, read_row, read_col);
  }

  input.close();

  Kokkos::deep_copy(A, A_h);
}

template <class VType, class IntType>
void readCRSFromMM(std::string name, const VType &V, const IntType &r, const IntType &c) {
  std::ifstream input(name);

  while (input.peek() == '%') input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  typename VType::HostMirror V_h   = Kokkos::create_mirror_view(V);
  typename IntType::HostMirror r_h = Kokkos::create_mirror_view(r);
  typename IntType::HostMirror c_h = Kokkos::create_mirror_view(c);

  int current_row = 0;
  int read_row;

  size_t nnz = c_h.extent(0);
  int nrows  = r_h.extent(0) - 1;

  r_h(0) = 0;

  for (size_t i = 0; i < nnz; ++i) {
    input >> read_row >> c_h(i);
    --read_row;
    --c_h(i);
    for (int tmp_row = current_row + 1; tmp_row <= read_row; ++tmp_row) r_h(tmp_row) = i;
    current_row = read_row;

    // if (VType::rank == 1)
    //  input >> V_h(i);
    if (VType::rank == 2)
      for (size_t j = 0; j < V_h.extent(0); ++j) input >> V_h(j, i);
  }

  r_h(nrows) = nnz;

  input.close();

  Kokkos::deep_copy(V, V_h);
  Kokkos::deep_copy(r, r_h);
  Kokkos::deep_copy(c, c_h);
}

template <class VType, class IntType>
void getInvDiagFromCRS(const VType &V, const IntType &r, const IntType &c, const VType &diag) {
  auto diag_values_host = Kokkos::create_mirror_view(diag);
  auto values_host      = Kokkos::create_mirror_view(V);
  auto row_ptr_host     = Kokkos::create_mirror_view(r);
  auto colIndices_host  = Kokkos::create_mirror_view(c);

  Kokkos::deep_copy(values_host, V);
  Kokkos::deep_copy(row_ptr_host, r);
  Kokkos::deep_copy(colIndices_host, c);

  int current_index;
  int N       = diag.extent(0);
  int BlkSize = diag.extent(1);

  for (int i = 0; i < BlkSize; ++i) {
    for (current_index = row_ptr_host(i); current_index < row_ptr_host(i + 1); ++current_index) {
      if (colIndices_host(current_index) == i) break;
    }
    for (int j = 0; j < N; ++j) diag_values_host(j, i) = 1. / values_host(j, current_index);
  }

  Kokkos::deep_copy(diag, diag_values_host);
}
