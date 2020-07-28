// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Ioss_DatabaseIO.h"   // for DatabaseIO
#include "Ioss_ElementBlock.h" // for ElementBlock
#include "Ioss_Property.h"     // for Property
#include <algorithm>           // for max_element, min_element
#include <cstddef>             // for size_t, nullptr
#include <fmt/ostream.h>
#include <iomanip>     // for operator<<, setw
#include <iostream>    // for operator<<, basic_ostream, etc
#include <string>      // for char_traits, operator<<
#include <sys/types.h> // for int64_t
#include <vector>      // for vector

namespace {
  void comp_grad12x(double *const grad_ptr, const double *const x, const double *const y,
                    const double *const z)
  {
    double R42 = (z[3] - z[1]);
    double R52 = (z[4] - z[1]);
    double R54 = (z[4] - z[3]);

    double R63 = (z[5] - z[2]);
    double R83 = (z[7] - z[2]);
    double R86 = (z[7] - z[5]);

    double R31 = (z[2] - z[0]);
    double R61 = (z[5] - z[0]);
    double R74 = (z[6] - z[3]);

    double R72 = (z[6] - z[1]);
    double R75 = (z[6] - z[4]);
    double R81 = (z[7] - z[0]);

    double t1 = (R63 + R54);
    double t2 = (R61 + R74);
    double t3 = (R72 + R81);

    double t4 = (R86 + R42);
    double t5 = (R83 + R52);
    double t6 = (R75 + R31);

    grad_ptr[0] =
        (y[1] * t1) - (y[2] * R42) - (y[3] * t5) + (y[4] * t4) + (y[5] * R52) - (y[7] * R54);
    grad_ptr[1] =
        (y[2] * t2) + (y[3] * R31) - (y[0] * t1) - (y[5] * t6) + (y[6] * R63) - (y[4] * R61);
    grad_ptr[2] =
        (y[3] * t3) + (y[0] * R42) - (y[1] * t2) - (y[6] * t4) + (y[7] * R74) - (y[5] * R72);
    grad_ptr[3] =
        (y[0] * t5) - (y[1] * R31) - (y[2] * t3) + (y[7] * t6) + (y[4] * R81) - (y[6] * R83);
    grad_ptr[4] =
        (y[5] * t3) + (y[6] * R86) - (y[7] * t2) - (y[0] * t4) - (y[3] * R81) + (y[1] * R61);
    grad_ptr[5] =
        (y[6] * t5) - (y[4] * t3) - (y[7] * R75) + (y[1] * t6) - (y[0] * R52) + (y[2] * R72);
    grad_ptr[6] =
        (y[7] * t1) - (y[5] * t5) - (y[4] * R86) + (y[2] * t4) - (y[1] * R63) + (y[3] * R83);
    grad_ptr[7] =
        (y[4] * t2) - (y[6] * t1) + (y[5] * R75) - (y[3] * t6) - (y[2] * R74) + (y[0] * R54);

    R42 = (x[3] - x[1]);
    R52 = (x[4] - x[1]);
    R54 = (x[4] - x[3]);

    R63 = (x[5] - x[2]);
    R83 = (x[7] - x[2]);
    R86 = (x[7] - x[5]);

    R31 = (x[2] - x[0]);
    R61 = (x[5] - x[0]);
    R74 = (x[6] - x[3]);

    R72 = (x[6] - x[1]);
    R75 = (x[6] - x[4]);
    R81 = (x[7] - x[0]);

    t1 = (R63 + R54);
    t2 = (R61 + R74);
    t3 = (R72 + R81);

    t4 = (R86 + R42);
    t5 = (R83 + R52);
    t6 = (R75 + R31);

    grad_ptr[8] =
        (z[1] * t1) - (z[2] * R42) - (z[3] * t5) + (z[4] * t4) + (z[5] * R52) - (z[7] * R54);
    grad_ptr[9] =
        (z[2] * t2) + (z[3] * R31) - (z[0] * t1) - (z[5] * t6) + (z[6] * R63) - (z[4] * R61);
    grad_ptr[10] =
        (z[3] * t3) + (z[0] * R42) - (z[1] * t2) - (z[6] * t4) + (z[7] * R74) - (z[5] * R72);
    grad_ptr[11] =
        (z[0] * t5) - (z[1] * R31) - (z[2] * t3) + (z[7] * t6) + (z[4] * R81) - (z[6] * R83);
    grad_ptr[12] =
        (z[5] * t3) + (z[6] * R86) - (z[7] * t2) - (z[0] * t4) - (z[3] * R81) + (z[1] * R61);
    grad_ptr[13] =
        (z[6] * t5) - (z[4] * t3) - (z[7] * R75) + (z[1] * t6) - (z[0] * R52) + (z[2] * R72);
    grad_ptr[14] =
        (z[7] * t1) - (z[5] * t5) - (z[4] * R86) + (z[2] * t4) - (z[1] * R63) + (z[3] * R83);
    grad_ptr[15] =
        (z[4] * t2) - (z[6] * t1) + (z[5] * R75) - (z[3] * t6) - (z[2] * R74) + (z[0] * R54);

    R42 = (y[3] - y[1]);
    R52 = (y[4] - y[1]);
    R54 = (y[4] - y[3]);

    R63 = (y[5] - y[2]);
    R83 = (y[7] - y[2]);
    R86 = (y[7] - y[5]);

    R31 = (y[2] - y[0]);
    R61 = (y[5] - y[0]);
    R74 = (y[6] - y[3]);

    R72 = (y[6] - y[1]);
    R75 = (y[6] - y[4]);
    R81 = (y[7] - y[0]);

    t1 = (R63 + R54);
    t2 = (R61 + R74);
    t3 = (R72 + R81);

    t4 = (R86 + R42);
    t5 = (R83 + R52);
    t6 = (R75 + R31);

    grad_ptr[16] =
        (x[1] * t1) - (x[2] * R42) - (x[3] * t5) + (x[4] * t4) + (x[5] * R52) - (x[7] * R54);
    grad_ptr[17] =
        (x[2] * t2) + (x[3] * R31) - (x[0] * t1) - (x[5] * t6) + (x[6] * R63) - (x[4] * R61);
    grad_ptr[18] =
        (x[3] * t3) + (x[0] * R42) - (x[1] * t2) - (x[6] * t4) + (x[7] * R74) - (x[5] * R72);
    grad_ptr[19] =
        (x[0] * t5) - (x[1] * R31) - (x[2] * t3) + (x[7] * t6) + (x[4] * R81) - (x[6] * R83);
    grad_ptr[20] =
        (x[5] * t3) + (x[6] * R86) - (x[7] * t2) - (x[0] * t4) - (x[3] * R81) + (x[1] * R61);
    grad_ptr[21] =
        (x[6] * t5) - (x[4] * t3) - (x[7] * R75) + (x[1] * t6) - (x[0] * R52) + (x[2] * R72);
    grad_ptr[22] =
        (x[7] * t1) - (x[5] * t5) - (x[4] * R86) + (x[2] * t4) - (x[1] * R63) + (x[3] * R83);
    grad_ptr[23] =
        (x[4] * t2) - (x[6] * t1) + (x[5] * R75) - (x[3] * t6) - (x[2] * R74) + (x[0] * R54);
  }

  inline double dot8(const double *const x1, const double *const x2)
  {
    double d1 = x1[0] * x2[0] + x1[1] * x2[1];
    double d2 = x1[2] * x2[2] + x1[3] * x2[3];
    double d4 = x1[4] * x2[4] + x1[5] * x2[5];
    double d5 = x1[6] * x2[6] + x1[7] * x2[7];
    return d1 + d2 + d4 + d5;
  }
} // namespace

template <typename T>
void hex_volume_internal(Ioss::ElementBlock *block, const std::vector<double> &coordinates,
                         const std::vector<T> &connectivity, size_t nelem)
{
  const double one12th = 1.0 / 12.0;

  double              gradop12x[24];
  double              x[8];
  double              y[8];
  double              z[8];
  std::vector<double> volume(nelem);

  size_t t1 = Ioss::Utils::timer();

  size_t count = 0;
  for (size_t ielem = 0; ielem < nelem; ++ielem) {
    if (count++ >= nelem / 100) {
      fmt::print(Ioss::OUTPUT(), ".");
      count = 0;
    }
    for (size_t j = 0; j < 8; j++) {
      size_t node = connectivity[ielem * 8 + j] - 1;
      x[j]        = coordinates[node * 3 + 0];
      y[j]        = coordinates[node * 3 + 1];
      z[j]        = coordinates[node * 3 + 2];
    }

    comp_grad12x(&gradop12x[0], x, y, z);
    const double volume12x = dot8(x, &gradop12x[0]);
    volume[ielem]          = volume12x * one12th;
  }
  size_t t2 = Ioss::Utils::timer();

  if (nelem > 0) {
    fmt::print(Ioss::OUTPUT(),
               "\n{:12}\tMin volume = {:12}  Max volume = {:12}  Elements = {:12n}  Time/Elem = "
               "{:5.3f} micro-sec.\n",
               block->name(), *std::min_element(volume.begin(), volume.end()),
               *std::max_element(volume.begin(), volume.end()), nelem,
               1000000 * double(t2 - t1) / nelem);
  }
}

void hex_volume(Ioss::ElementBlock *block, const std::vector<double> &coordinates)
{

  size_t nelem = block->entity_count();
  if (block->get_database()->int_byte_size_api() == 8) {
    std::vector<int64_t> connectivity(nelem);
    block->get_field_data("connectivity_raw", connectivity);
    hex_volume_internal(block, coordinates, connectivity, nelem);
  }
  else {
    std::vector<int> connectivity(nelem);
    block->get_field_data("connectivity_raw", connectivity);
    hex_volume_internal(block, coordinates, connectivity, nelem);
  }
}
