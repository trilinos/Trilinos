// Copyright(C) 2021, 2022, 2023, 2025, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#include <exception>
#include <fstream>
#include <iterator>
#include <string>

#include "add_to_log.h"
#include "fmt/chrono.h"
#include "fmt/color.h"
#include "fmt/format.h"
#include "fmt/ostream.h"
#include "time_stamp.h"
#include "tokenize.h"

#include <exodusII.h>

#include <Ionit_Initializer.h>
#include <Ioss_MemoryUtils.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_SmartAssert.h>
#include <Ioss_Utils.h>

#include "Grid.h"
#include "ZE_SystemInterface.h"

#ifdef SEACAS_HAVE_MPI
#include <mpi.h>
#endif

//! \file

namespace {
  Grid define_lattice(SystemInterface &interFace, Ioss::ParallelUtils &pu);
} // namespace

std::string  tsFormat    = "[{:%H:%M:%S}]";
unsigned int debug_level = 0;

template <typename INT> double zellij(SystemInterface &interFace, INT dummy);

int main(int argc, char *argv[])
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);

#endif

  Ioss::ParallelUtils pu{};
  int                 my_rank = pu.parallel_rank();

  try {
    if (my_rank == 0) {
      SystemInterface::show_version();
#ifdef SEACAS_HAVE_MPI
      fmt::print("\tParallel Capability Enabled.\n");
#else
      fmt::print("\tParallel Capability Not Enabled.\n");
#endif
    }
    Ioss::Init::Initializer io;

    SystemInterface interFace(my_rank);
    bool            ok = interFace.parse_options(argc, argv);

    if (!ok) {
      fmt::print(stderr, fmt::fg(fmt::color::red),
                 "\nERROR: Problems parsing command line arguments.\n\n");
      exit(EXIT_FAILURE);
    }

    debug_level = interFace.debug();

    if (debug_level & 1) {
      ex_opts(EX_VERBOSE | EX_DEBUG);
    }
    else {
      ex_opts(0);
    }

    double time = 0.0;
    if (interFace.ints32bit()) {
      time = zellij(interFace, 0);
    }
    else {
      time = zellij(interFace, static_cast<int64_t>(0));
    }

    if (my_rank == 0) {
      fmt::print("\n Zellij execution successful.\n");
      add_to_log(argv[0], time);
    }

#ifdef SEACAS_HAVE_MPI
    MPI_Finalize();
#endif
  }
  catch (std::exception &e) {
    fmt::print(stderr, fmt::fg(fmt::color::red), "ERROR: Standard exception: {}\n", e.what());
  }
}

template <typename INT> double zellij(SystemInterface &interFace, INT /*dummy*/)
{
  double              begin = Ioss::Utils::timer();
  Ioss::ParallelUtils pu{};

  if (debug_level & 2) {
    fmt::print(stderr, "{} Begin Execution\n", time_stamp(tsFormat));
  }

  auto grid = define_lattice(interFace, pu);

  grid.generate_sidesets();
  grid.set_coordinate_offsets();
  grid.decompose(interFace.decomp_method());

  if (debug_level & 2) {
    fmt::print(stderr, "{} Lattice Decomposed\n", time_stamp(tsFormat));
  }

  // All unit cells have been mapped into the IxJ grid, now calculate all node / element offsets
  // and the global node and element counts...
  //
  // Iterate through the grid starting with (0,0) and accumulate node and element counts...
  grid.process(interFace, (INT)0);

  /*************************************************************************/
  // EXIT program
  if (debug_level & 2) {
    fmt::print(stderr, "{} Execution Complete\n", time_stamp(tsFormat));
  }

  double end = Ioss::Utils::timer();
  double hwm = (double)Ioss::MemoryUtils::get_hwm_memory_info() / 1024.0 / 1024.0;
  if (pu.parallel_rank() == 0) {
    fmt::print("\n Total Execution Time     = {:.5} seconds.\n", end - begin);
    fmt::print(" High-Water Memory Use    = {:.3} MiBytes.\n", hwm);
  }
  return (end - begin);
}

namespace {
  Grid define_lattice(SystemInterface &interFace, Ioss::ParallelUtils &pu)
  {
    int my_rank = pu.parallel_rank();

    Grid grid(interFace);

    if (debug_level & 2) {
      pu.progress("Defining Unit Cells...");
    }
    std::string filename = interFace.lattice();

    std::ifstream input(filename, std::ios::in);
    if (!input.is_open()) {
      fmt::print(stderr, fmt::fg(fmt::color::red),
                 "\nERROR: Could not open/access lattice file '{}'\n", filename);
      exit(EXIT_FAILURE);
    }

    bool        in_dictionary{false};
    bool        in_lattice{false};
    std::string line;
    while (getline(input, line)) {
      if (line.empty()) {
        continue;
      }

      auto tokens = Ioss::tokenize(line, " ");
      if (tokens[0] == "BEGIN_DICTIONARY") {
        SMART_ASSERT(!in_lattice && !in_dictionary);
        in_dictionary = true;
      }
      else if (tokens[0] == "END_DICTIONARY") {
        SMART_ASSERT(!in_lattice && in_dictionary);
        in_dictionary = false;
        if (debug_level & 2) {
          pu.progress("Unit Cells Defined...");
        }
      }
      else if (in_dictionary) {
        if (tokens.size() != 2) {
          fmt::print(stderr, fmt::fg(fmt::color::red),
                     "\nERROR: There are {} entries on a lattice dictionary line; there should be "
                     "only 2:\n\t'{}'.\n\n",
                     tokens.size(), line);
          exit(EXIT_FAILURE);
        }

        grid.add_unit_cell(tokens[0], tokens[1], interFace.ints32bit());
      }
      else if (tokens[0] == "BEGIN_LATTICE") {
        SMART_ASSERT(!in_lattice && !in_dictionary);
        in_lattice = true;
        break;
      }
    }

    if (!in_lattice) {
      fmt::print(
          stderr, fmt::fg(fmt::color::red),
          "\nERROR: Reached end of input file without finding a 'BEGIN_LATTICE' command\n\n");
      exit(EXIT_FAILURE);
    }

    // Tokenize line to get I J K size of lattice
    auto tokens = Ioss::tokenize(line, " ");
    SMART_ASSERT(tokens[0] == "BEGIN_LATTICE")(tokens[0])(line);

    if (tokens.size() != 4) {
      fmt::print(stderr, fmt::fg(fmt::color::red),
                 "\nERROR: The 'BEGIN_LATTICE' line has incorrect syntax.  It should be "
                 "'BEGIN_LATTICE I J K'\n"
                 "\tThe line was '{}'\n\n",
                 line);
      exit(EXIT_FAILURE);
    }

    SMART_ASSERT(tokens.size() == 4)(tokens.size())(line);
    int II = std::stoi(tokens[1]);
    int JJ = std::stoi(tokens[2]);
    int KK = std::stoi(tokens[3]);
    SMART_ASSERT(KK == 1);

    if (interFace.repeat() > 1) {
      II *= interFace.repeat();
      JJ *= interFace.repeat();
    }

    grid.set_extent(II, JJ, KK);
    grid.handle_file_count();

    if (my_rank == 0) {
      fmt::print("\n Lattice:\tUnit Cells: {},\tGrid Size:  {} x {} x {}\n",
                 fmt::group_digits(grid.unit_cells().size()), fmt::group_digits(II),
                 fmt::group_digits(JJ), fmt::group_digits(KK));
    }
    if (interFace.ranks() > 1) {
      fmt::print("         \t[{}] Ranks: {}, Outputting {} ranks starting at rank {}.\n", my_rank,
                 fmt::group_digits(interFace.ranks()), fmt::group_digits(interFace.rank_count()),
                 fmt::group_digits(interFace.start_rank()));
    }

    // Now process the lattice portion of the lattice file...
    size_t row{0};
    while (getline(input, line)) {
      if (line.empty()) {
        continue;
      }

      tokens = Ioss::tokenize(line, " ");
      if (tokens[0] == "END_LATTICE") {
        SMART_ASSERT(in_lattice && !in_dictionary);
        in_lattice = false;

        // Check row count to make sure matches 'J' size of lattice
        if (row != grid.JJ()) {
          fmt::print(stderr, fmt::fg(fmt::color::red),
                     "\nERROR: Only {} rows of the {} x {} lattice were defined.\n\n", row, II, JJ);
          exit(EXIT_FAILURE);
        }
        break;
      }
      if (in_lattice) {
        // TODO: Currently assumes that each row in the lattice is defined on a single row;
        //       This will need to be relaxed since a lattice of 5000x5000 would result in
        //       lines that are too long and would be easier to split a row over multiple lines...
        if (tokens.size() * interFace.repeat() != grid.II()) {
          fmt::print(
              stderr, fmt::fg(fmt::color::red),
              "\nERROR: Line {} of the lattice definition has {} entries.  It should have {}.\n\n",
              row + 1, tokens.size() * interFace.repeat(), grid.II());
          exit(EXIT_FAILURE);
        }

        if (row >= grid.JJ()) {
          fmt::print(stderr, fmt::fg(fmt::color::red),
                     "\nERROR: There are too many rows in the lattice definition. The lattice is "
                     "{} x {}.\n\n",
                     grid.II(), grid.JJ());
          exit(EXIT_FAILURE);
        }

        if (tokens.size() * interFace.repeat() != grid.II()) {
          fmt::print(stderr, fmt::fg(fmt::color::red),
                     "\nERROR: In row {}, there is an incorrect number of entries.  There should "
                     "be {}, but found {}.\n",
                     row + 1, grid.II(), tokens.size() * interFace.repeat());
          exit(EXIT_FAILURE);
        }

        auto repeat = interFace.repeat();
        for (int j = 0; j < repeat; j++) {
          size_t col = 0;
          for (auto &key : tokens) {

            for (int i = 0; i < repeat; i++) {
              if (!grid.initialize(col++, row, key)) {
                fmt::print(
                    stderr, fmt::fg(fmt::color::red),
                    "\nERROR: In row {}, column {}, the lattice specifies a unit cell ({}) that "
                    "has not been defined.\n\n",
                    row + 1, col + 1, key);
                exit(EXIT_FAILURE);
              }
            }
          }
          row++;
        }
      }
    }
    if (debug_level & 2) {
      fmt::print(stderr, "{} Lattice Defined\n", time_stamp(tsFormat));
    }
    return grid;
  }
} // namespace
