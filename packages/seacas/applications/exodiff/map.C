// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <algorithm>
#include <cfloat>
#include <cstdlib>
#include <iomanip>
#include <numeric>

#include "ED_SystemInterface.h"
#include "Tolerance.h"
#include "exoII_read.h"
#include "exo_block.h"
#include "fmt/ostream.h"
#include "iqsort.h"
#include "smart_assert.h"
#include "util.h"

namespace {
  double find_range(const double *x, size_t num_nodes);

  template <typename INT>
  int64_t Find(double x0, double y0, double z0, const std::vector<double> &x,
               const std::vector<double> &y, const std::vector<double> &z,
               const std::vector<INT> &id, int dim, bool ignore_dups);

  template <typename INT>
  void Compute_Node_Map(std::vector<INT> &node_map, ExoII_Read<INT> &file1, ExoII_Read<INT> &file2);
} // namespace

template <typename INT>
void Compute_Maps(std::vector<INT> &node_map, std::vector<INT> &elmt_map, ExoII_Read<INT> &file1,
                  ExoII_Read<INT> &file2)
{
  SMART_ASSERT(file1.Open());
  SMART_ASSERT(file2.Open());

  size_t num_nodes = file1.Num_Nodes();
  size_t num_elmts = file1.Num_Elements();
  int    dim       = file1.Dimension();

  //  ********************  elements  ********************  //

  // Load global ids (0-offset) into id array.
  std::vector<INT> id(num_elmts);
  std::iota(id.begin(), id.end(), 0);

  // Get map storage.
  node_map.resize(num_nodes);
  std::fill(node_map.begin(), node_map.end(), -1);

  elmt_map.resize(num_elmts);
  std::fill(elmt_map.begin(), elmt_map.end(), -1);

  // Create storage for midpoints.
  std::vector<double> x2;
  x2.reserve(num_elmts);
  std::vector<double> y2;
  y2.reserve(dim > 1 ? num_elmts : 0);
  std::vector<double> z2;
  z2.reserve(dim > 2 ? num_elmts : 0);

  // Load coordinates for file 2 and get pointers to them.
  file2.Load_Nodal_Coordinates();
  const auto *x2_f = file2.X_Coords();
  const auto *y2_f = file2.Y_Coords();
  const auto *z2_f = file2.Z_Coords();

  // Load connectivities for all blocks in second file.
  file2.Load_Element_Block_Descriptions();

  {
    // Compute midpoints of each element and place into x,y,z arrays.
    size_t num_blocks = file2.Num_Element_Blocks();
    for (size_t b = 0; b < num_blocks; ++b) {
      const Exo_Block<INT> *block              = file2.Get_Element_Block_by_Index(b);
      size_t                num_elmts_in_block = block->Size();
      size_t                num_nodes_per_elmt = block->Num_Nodes_per_Element();
      for (size_t i = 0; i < num_elmts_in_block; ++i) {
        const INT *conn  = block->Connectivity(i); // Connectivity for element i.
        double     sum_x = 0.0;
        double     sum_y = 0.0;
        double     sum_z = 0.0;
        for (size_t j = 0; j < num_nodes_per_elmt; ++j) {
          sum_x += x2_f[conn[j] - 1];
          if (dim > 1) {
            sum_y += y2_f[conn[j] - 1];
          }
          if (dim > 2) {
            sum_z += z2_f[conn[j] - 1];
          }
        }
        x2.push_back(sum_x / static_cast<double>(num_nodes_per_elmt));
        if (dim > 1) {
          y2.push_back(sum_y / static_cast<double>(num_nodes_per_elmt));
        }
        if (dim > 2) {
          z2.push_back(sum_z / static_cast<double>(num_nodes_per_elmt));
        }
      }
    }
  }

  // Sort by x value.
  index_qsort(Data(x2), Data(id), num_elmts);

#if 0
  fmt::print("******************  elmts  ******************** \n");
  {for (size_t i = 0; i < num_elmts; ++i)
      fmt::print("{})\t{}\t{}\t{}\t{}\n"
                 i, x2[id[i]], y2[id[i]], z2[id[i]], id[i]);}
  fmt::print("******************  elmts  ******************** \n");
#endif
  //  Load and get nodal coordinates for first file.
  file1.Load_Nodal_Coordinates();
  const auto *x1_f = file1.X_Coords();
  const auto *y1_f = file1.Y_Coords();
  const auto *z1_f = file1.Z_Coords();

  // Cannot ignore the comparisons, so make sure the coord_tol_type
  // is not -1 which is "ignore"
  ToleranceMode save_tolerance_type = interFace.coord_tol.type;
  if (save_tolerance_type == ToleranceMode::IGNORE_) {
    interFace.coord_tol.type = ToleranceMode::ABSOLUTE_;
  }

  // Match elmts in first file to their corresponding elmts in second.
  size_t num_blocks = file1.Num_Element_Blocks();
  size_t e1         = 0;

  for (size_t b = 0; b < num_blocks; ++b) {
    const Exo_Block<INT> *block1 = file1.Get_Element_Block_by_Index(b);
    file1.Load_Element_Block_Description(b);
    size_t num_elmts_in_block = block1->Size();
    size_t num_nodes_per_elmt = block1->Num_Nodes_per_Element();
    for (size_t i = 0; i < num_elmts_in_block; ++i) {
      // Connectivity for element i.
      const INT *conn1 = block1->Connectivity(i);

      // Compute midpoint.
      double mid_x = 0.0;
      double mid_y = 0.0;
      double mid_z = 0.0;

      for (size_t j = 0; j < num_nodes_per_elmt; ++j) {
        SMART_ASSERT(conn1[j] >= 1 && conn1[j] <= (INT)num_nodes);
        mid_x += x1_f[conn1[j] - 1];
        if (dim > 1) {
          mid_y += y1_f[conn1[j] - 1];
        }
        if (dim > 2) {
          mid_z += z1_f[conn1[j] - 1];
        }
      }
      mid_x /= static_cast<double>(num_nodes_per_elmt);
      if (dim > 1) {
        mid_y /= static_cast<double>(num_nodes_per_elmt);
      }
      if (dim > 2) {
        mid_z /= static_cast<double>(num_nodes_per_elmt);
      }

      // Locate midpoint in sorted array.
      int64_t sort_idx = Find(mid_x, mid_y, mid_z, x2, y2, z2, id, dim, interFace.ignore_dups);

      if (sort_idx < 0) {
        Error(fmt::format("Files are different (couldn't match element {} from block {} from first "
                          "file to second)\n",
                          i + 1, file1.Block_Id(b)));
      }
      size_t e2 = id[sort_idx];

      // Assign element map for this element.
      elmt_map[e1] = e2;

      {
        // Determine the block and elmt index of matched element.
        auto bl_idx = file2.Global_to_Block_Local(e2 + 1);

        const Exo_Block<INT> *block2 = file2.Get_Element_Block_by_Index(bl_idx.first);
        SMART_ASSERT(block2 != nullptr);

        // Check that the element types are the same.
        if (num_nodes_per_elmt != block2->Num_Nodes_per_Element()) {
          Error(fmt::format("Files are different.\n"
                            " In File 1: Element {} in Block {} has {}  and\n"
                            " In File 2: Element {} in Block {} has {}\n",
                            fmt::group_digits(i + 1), file1.Block_Id(b), num_nodes_per_elmt,
                            fmt::group_digits(bl_idx.second + 1), file2.Block_Id(bl_idx.first),
                            block2->Num_Nodes_per_Element()));
        }

        // Get connectivity for file2 element.
        const INT *conn2 = block2->Connectivity(bl_idx.second);

        // Match each node in the first elmt with a node in the second
        // and assign node_map.
        for (size_t ln1 = 0; ln1 < num_nodes_per_elmt; ++ln1) {
          // Grab coordinate of node in first file.
          double x1_val = x1_f[conn1[ln1] - 1];
          double y1_val = dim > 1 ? y1_f[conn1[ln1] - 1] : 0.0;
          double z1_val = dim > 2 ? z1_f[conn1[ln1] - 1] : 0.0;

          size_t found = 0;
          for (size_t ln2 = 0; ln2 < num_nodes_per_elmt; ++ln2) {
            // Grab coordinate of node in second file.
            double x2_val = x2_f[conn2[ln2] - 1];
            double y2_val = dim > 1 ? y2_f[conn2[ln2] - 1] : 0.0;
            double z2_val = dim > 2 ? z2_f[conn2[ln2] - 1] : 0.0;

            if (!interFace.coord_tol.Diff(x1_val, x2_val) &&
                !interFace.coord_tol.Diff(y1_val, y2_val) &&
                !interFace.coord_tol.Diff(z1_val, z2_val)) {
              // assert that if this node has been given a map
              // previously, that it agrees with the latest
              // assignment.
              if (node_map[conn1[ln1] - 1] >= 0 && node_map[conn1[ln1] - 1] != conn2[ln2] - 1) {

                if (!interFace.ignore_dups) {
                  // Node in file 1.
                  INT    node1 = conn1[ln1];
                  double x1a   = x1_f[node1 - 1];
                  double y1a   = dim >= 2 ? y1_f[node1 - 1] : 0.0;
                  double z1a   = dim >= 3 ? z1_f[node1 - 1] : 0.0;

                  // Node in file 2 that was already mapped to node 1 in file 1
                  INT    n1  = node_map[conn1[ln1] - 1] + 1;
                  double x2a = x2_f[n1 - 1];
                  double y2a = dim >= 2 ? y2_f[n1 - 1] : 0.0;
                  double z2a = dim >= 3 ? z2_f[n1 - 1] : 0.0;

                  // Node in file 2 that is now being mapped to node 1 in file 1
                  INT    n2  = conn2[ln2];
                  double x2b = x2_f[n2 - 1];
                  double y2b = dim >= 2 ? y2_f[n2 - 1] : 0.0;
                  double z2b = dim >= 3 ? z2_f[n2 - 1] : 0.0;

                  SMART_ASSERT(!interFace.coord_tol.Diff(x2a, x2b) &&
                               !interFace.coord_tol.Diff(y2a, y2b) &&
                               !interFace.coord_tol.Diff(z2a, z2b));
                  Error(fmt::format("No unique node mapping possible.\n"
                                    "\tFile 1, Node {} at ({}, {}, {}) maps to both:\n"
                                    "\tFile 2, Node {} at ({}, {}, {}) and\n"
                                    "\tFile 2, Node {} at ({}, {}, {})\n\n",
                                    fmt::group_digits(node1), x1a, y1a, z1a, fmt::group_digits(n1),
                                    x2a, y2a, z2a, fmt::group_digits(n2), x2b, y2b, z2b));
                }
                found = 1;
                break;
              }
              node_map[conn1[ln1] - 1] = conn2[ln2] - 1;
              found                    = 1;
              break;
            }
          }
          if (!found) {
            std::ostringstream out;
            fmt::print(out,
                       "\nCannot find a match for node at position {} in first element.\n"
                       "\tFile 1: Element {} in Block {} nodes:\n",
                       ln1 + 1, fmt::group_digits(i + 1), file1.Block_Id(b));
            for (size_t l1 = 0; l1 < num_nodes_per_elmt; ++l1) {
              double x_val = x1_f[conn1[l1] - 1];
              double y_val = dim > 1 ? y1_f[conn1[l1] - 1] : 0.0;
              double z_val = dim > 2 ? z1_f[conn1[l1] - 1] : 0.0;
              fmt::print(out, "\t({})\t{}\t{:.9e}\t{:.9e}\t{:.9e}\n", l1 + 1,
                         fmt::group_digits(conn1[l1]), x_val, y_val, z_val);
            }
            fmt::print(out, "\tFile 2: Element {} in Block {} nodes:\n",
                       fmt::group_digits(bl_idx.second + 1), file1.Block_Id(b));
            for (size_t l3 = 0; l3 < num_nodes_per_elmt; ++l3) {
              double x_val = x2_f[conn2[l3] - 1];
              double y_val = dim > 1 ? y2_f[conn2[l3] - 1] : 0.0;
              double z_val = dim > 2 ? z2_f[conn2[l3] - 1] : 0.0;
              fmt::print(out, "\t({})\t{}\t{:.9e}\t{:.9e}\t{:.9e}\n", l3 + 1,
                         fmt::group_digits(conn2[l3]), x_val, y_val, z_val);
            }
            fmt::print(out, "Coordinates compared using tolerance: {} ({}), floor: {}\n",
                       interFace.coord_tol.value, interFace.coord_tol.typestr(),
                       interFace.coord_tol.floor);
            Error(out);
          }
        } // End of local node loop on file1's element.
      } // End of local node search block.

      ++e1;

    } // End of loop on elements in file1 element block.

    file1.Free_Element_Block(b);

  } // End of loop on file1 blocks.

  // Check that all nodes in the file have been matched...  If any
  // unmatched nodes are found, then perform a node-based matching
  // algorithm...
  for (size_t i = 0; i < num_nodes; i++) {
    if (node_map[i] < 0) {
      Compute_Node_Map(node_map, file1, file2);
      break;
    }
  }

  file1.Free_Nodal_Coordinates();
  file2.Free_Nodal_Coordinates();
  file2.Free_Element_Blocks();

  interFace.coord_tol.type = save_tolerance_type;
}

template <typename INT>
void Compute_Partial_Maps(std::vector<INT> &node_map, std::vector<INT> &elmt_map,
                          ExoII_Read<INT> &file1, ExoII_Read<INT> &file2)
{
  SMART_ASSERT(file1.Open());
  SMART_ASSERT(file2.Open());

  size_t num_nodes1 = file1.Num_Nodes();
  size_t num_elmts1 = file1.Num_Elements();

  size_t num_nodes2 = file2.Num_Nodes();
  size_t num_elmts2 = file2.Num_Elements();
  int    dim        = file1.Dimension();
  SMART_ASSERT(dim == file2.Dimension());

  //  ********************  elements  ********************  //

  // Load global ids (0-offset) into id array.
  std::vector<INT> id2(num_elmts2);
  std::iota(id2.begin(), id2.end(), 0);

  // Get map storage.
  node_map.resize(num_nodes1);
  std::fill(node_map.begin(), node_map.end(), -1);

  elmt_map.resize(num_elmts1);
  std::fill(elmt_map.begin(), elmt_map.end(), -1);

  // Create storage for midpoints.
  std::vector<double> x2;
  std::vector<double> y2;
  std::vector<double> z2;
  x2.reserve(num_elmts2);
  if (dim > 1) {
    y2.reserve(num_elmts2);
  }
  if (dim > 2) {
    z2.reserve(num_elmts2);
  }

  // Load coordinates for file 2 and get pointers to them.
  file2.Load_Nodal_Coordinates();
  const auto *x2_f = file2.X_Coords();
  const auto *y2_f = file2.Y_Coords();
  const auto *z2_f = file2.Z_Coords();

  // Load connectivities for all blocks in second file.
  file2.Load_Element_Block_Descriptions();

  {
    // Compute midpoints of each element and place into x,y,z arrays.
    size_t num_blocks2 = file2.Num_Element_Blocks();
    double sum_x;
    double sum_y;
    double sum_z;
    for (size_t b = 0; b < num_blocks2; ++b) {
      const Exo_Block<INT> *block              = file2.Get_Element_Block_by_Index(b);
      size_t                num_elmts_in_block = block->Size();
      size_t                num_nodes_per_elmt = block->Num_Nodes_per_Element();
      for (size_t i = 0; i < num_elmts_in_block; ++i) {
        const INT *conn = block->Connectivity(i); // Connectivity for element i.
        sum_x           = 0.0;
        sum_y           = 0.0;
        sum_z           = 0.0;
        for (size_t j = 0; j < num_nodes_per_elmt; ++j) {
          sum_x += x2_f[conn[j] - 1];
          if (dim > 1) {
            sum_y += y2_f[conn[j] - 1];
          }
          if (dim > 2) {
            sum_z += z2_f[conn[j] - 1];
          }
        }
        x2.push_back(sum_x / static_cast<double>(num_nodes_per_elmt));
        if (dim > 1) {
          y2.push_back(sum_y / static_cast<double>(num_nodes_per_elmt));
        }
        if (dim > 2) {
          z2.push_back(sum_z / static_cast<double>(num_nodes_per_elmt));
        }
      }
    }
  }

  // Sort by x value.
  index_qsort(Data(x2), Data(id2), num_elmts2);

#if 0
  fmt::print("******************  elmts  ******************** \n");
  {for (size_t i = 0; i < num_elmts; ++i)
      fmt::print("{})\t{}\t{}\t{}\t{}\n"
                 i, x2[id[i]], y2[id[i]], z2[id[i]], id[i]);}
  fmt::print("******************  elmts  ******************** \n");
#endif
  //  Load and get nodal coordinates for first file.
  file1.Load_Nodal_Coordinates();
  const auto *x1_f = file1.X_Coords();
  const auto *y1_f = file1.Y_Coords();
  const auto *z1_f = file1.Z_Coords();

  // Cannot ignore the comparisons, so make sure the coord_tol_type
  // is not -1 which is "ignore"
  ToleranceMode save_tolerance_type = interFace.coord_tol.type;
  if (save_tolerance_type == ToleranceMode::IGNORE_) {
    interFace.coord_tol.type = ToleranceMode::ABSOLUTE_;
  }

  // Match elmts in first file to their corresponding elmts in second.
  size_t num_blocks1 = file1.Num_Element_Blocks();
  size_t e1          = 0;

  bool   first     = true;
  size_t unmatched = 0;
  for (size_t b = 0; b < num_blocks1; ++b) {
    const Exo_Block<INT> *block1 = file1.Get_Element_Block_by_Index(b);
    file1.Load_Element_Block_Description(b);
    size_t num_elmts_in_block = block1->Size();
    size_t num_nodes_per_elmt = block1->Num_Nodes_per_Element();
    for (size_t i = 0; i < num_elmts_in_block; ++i) {
      // Connectivity for element i.
      const INT *conn1 = block1->Connectivity(i);

      // Compute midpoint.
      double mid_x = 0.0;
      double mid_y = 0.0;
      double mid_z = 0.0;

      for (size_t j = 0; j < num_nodes_per_elmt; ++j) {
        SMART_ASSERT(conn1[j] >= 1 && conn1[j] <= (INT)num_nodes1);
        mid_x += x1_f[conn1[j] - 1];
        if (dim > 1) {
          mid_y += y1_f[conn1[j] - 1];
        }
        if (dim > 2) {
          mid_z += z1_f[conn1[j] - 1];
        }
      }
      mid_x /= static_cast<double>(num_nodes_per_elmt);
      if (dim > 1) {
        mid_y /= static_cast<double>(num_nodes_per_elmt);
      }
      if (dim > 2) {
        mid_z /= static_cast<double>(num_nodes_per_elmt);
      }

      // Locate midpoint in sorted array.
      int64_t sort_idx = Find(mid_x, mid_y, mid_z, x2, y2, z2, id2, dim, interFace.ignore_dups);
      if (sort_idx < 0) {
        unmatched++;
        if (first && interFace.show_unmatched) {
          fmt::print("exodiff: Doing Partial Comparison: No Match for (b.e):\n");
        }
        first = false;
        if (interFace.show_unmatched) {
          fmt::print("{}.{}, ", file1.Block_Id(b), (i + 1));
        }
      }
      else {
        size_t e2    = id2[sort_idx];
        elmt_map[e1] = e2;

        // Assign element map for this element.

        // Determine the block and elmt index of matched element.
        auto bl_idx = file2.Global_to_Block_Local(e2 + 1);

        const Exo_Block<INT> *block2 = file2.Get_Element_Block_by_Index(bl_idx.first);
        SMART_ASSERT(block2 != nullptr);

        // Check that the element types are the same.
        if (num_nodes_per_elmt != block2->Num_Nodes_per_Element()) {
          Error(fmt::format("Files are different.\n"
                            " In File 1: Element {} in Block {} has {}  and\n"
                            " In File 2: Element {} in Block {} has {}\n",
                            fmt::group_digits(i + 1), file1.Block_Id(b), num_nodes_per_elmt,
                            fmt::group_digits(bl_idx.second + 1), file2.Block_Id(bl_idx.first),
                            block2->Num_Nodes_per_Element()));
        }

        // Get connectivity for file2 element.
        const INT *conn2 = block2->Connectivity(bl_idx.second);

        // Match each node in the first elmt with a node in the second
        // and assign node_map.
        for (size_t ln1 = 0; ln1 < num_nodes_per_elmt; ++ln1) {
          // Grab coordinate of node in first file.
          double x1_val = x1_f[conn1[ln1] - 1];
          double y1_val = dim > 1 ? y1_f[conn1[ln1] - 1] : 0.0;
          double z1_val = dim > 2 ? z1_f[conn1[ln1] - 1] : 0.0;
          size_t found  = 0;
          for (size_t ln2 = 0; ln2 < num_nodes_per_elmt; ++ln2) {
            // Grab coordinate of node in second file.
            double x2_val = x2_f[conn2[ln2] - 1];
            double y2_val = dim > 1 ? y2_f[conn2[ln2] - 1] : 0.0;
            double z2_val = dim > 2 ? z2_f[conn2[ln2] - 1] : 0.0;

            if (!interFace.coord_tol.Diff(x1_val, x2_val) &&
                !interFace.coord_tol.Diff(y1_val, y2_val) &&
                !interFace.coord_tol.Diff(z1_val, z2_val)) {
              node_map[conn1[ln1] - 1] = conn2[ln2] - 1;
              found                    = 1;
              break;
            }
          }
          if (!found) {
            std::ostringstream out;
            fmt::print(out,
                       "\nCannot find a match for node at position {} in first element.\n"
                       "\tFile 1: Element {} in Block {} nodes:\n",
                       ln1 + 1, fmt::group_digits(i + 1), file1.Block_Id(b));
            for (size_t l1 = 0; l1 < num_nodes_per_elmt; ++l1) {
              double x_val = x1_f[conn1[l1] - 1];
              double y_val = dim > 1 ? y1_f[conn1[l1] - 1] : 0.0;
              double z_val = dim > 2 ? z1_f[conn1[l1] - 1] : 0.0;
              fmt::print(out, "\t({})\t{}\t{:.9e}\t{:.9e}\t{:.9e}\n", l1 + 1,
                         fmt::group_digits(conn1[l1]), x_val, y_val, z_val);
            }
            fmt::print(out, "\tFile 2: Element {} in Block {} nodes:\n",
                       fmt::group_digits(bl_idx.second + 1), file1.Block_Id(b));
            for (size_t l3 = 0; l3 < num_nodes_per_elmt; ++l3) {
              double x_val = x2_f[conn2[l3] - 1];
              double y_val = dim > 1 ? y2_f[conn2[l3] - 1] : 0.0;
              double z_val = dim > 2 ? z2_f[conn2[l3] - 1] : 0.0;
              fmt::print(out, "\t({})\t{}\t{:.9e}\t{:.9e}\t{:.9e}\n", l3 + 1,
                         fmt::group_digits(conn2[l3]), x_val, y_val, z_val);
            }
            fmt::print(out, "Coordinates compared using tolerance: {} ({}), floor: {}\n",
                       interFace.coord_tol.value, interFace.coord_tol.typestr(),
                       interFace.coord_tol.floor);
            Error(out);
          }
        } // End of local node loop on file1's element.
      } // End of local node search block.

      ++e1;

    } // End of loop on elements in file1 element block.
    file1.Free_Element_Block(b);

  } // End of loop on file1 blocks.
  if (!first) {
    fmt::print("\nPartial Map selected -- {} elements unmatched\n", fmt::group_digits(unmatched));
  }
  else {
    if (num_elmts1 == num_elmts2 && num_nodes1 == num_nodes2) {
      fmt::print(
          "exodiff: INFO .. Partial Map was specified, but not needed.  All elements matched.\n");
    }
  }

  // Check that all nodes in the file have been matched...  If any
  // unmatched nodes are found, then perform a node-based matching
  // algorithm...
  //   for (size_t i=0; i < num_nodes; i++) {
  //     if (node_map[i] < 0) {
  //       Compute_Node_Map(node_map, file1, file2);
  //       break;
  //     }
  //   }

  file1.Free_Nodal_Coordinates();
  file2.Free_Nodal_Coordinates();
  file2.Free_Element_Blocks();

  interFace.coord_tol.type = save_tolerance_type;
}

namespace {
  template <typename INT> bool check_sort(const INT *map, size_t count)
  {
    for (size_t i = 1; i < count; i++) {
      if (map[i - 1] > map[i]) {
        return true;
      }
    }
    return false;
  }

  template <typename INT>
  bool internal_compute_maps(std::vector<INT> &map, const INT *file1_id_map,
                             const INT *file2_id_map, size_t count, const char *type)
  {
    std::vector<INT> id1;
    id1.reserve(count);
    std::vector<INT> id2;
    id2.reserve(count);
    for (size_t i = 0; i < count; i++) {
      id1.push_back(i);
      id2.push_back(i);
    }

    // Check whether sorting needed...
    bool sort1_needed = check_sort(file1_id_map, count);
    if (sort1_needed) {
      index_qsort(file1_id_map, Data(id1), count);
    }

    bool sort2_needed = check_sort(file2_id_map, count);
    if (sort2_needed) {
      index_qsort(file2_id_map, Data(id2), count);
    }

    for (size_t i = 0; i < count; i++) {
      if (file1_id_map[id1[i]] == file2_id_map[id2[i]]) {
        map[id1[i]] = id2[i];
      }
      else {
        Error(fmt::format("Unable to match {0} {1} in first file with {0} in second file.\n", type,
                          fmt::group_digits(file1_id_map[id1[i]])));
      }
    }

    // See if there is any mapping happening...
    bool mapped = false;
    for (INT i = 0; i < (INT)count; i++) {
      if (i != map[i]) {
        mapped = true;
        break;
      }
    }
    return mapped;
  }
} // namespace

template <typename INT>
void Compute_FileId_Maps(std::vector<INT> &node_map, std::vector<INT> &elmt_map,
                         ExoII_Read<INT> &file1, ExoII_Read<INT> &file2)
{
  // Compute map of nodes and elements in file1 to nodes and elements in file2
  // Use the internal exodus node and element number maps in file1 and file2 to
  // do the matching.  Currently assume (and verify) that number of nodes and
  // elements match in the two files.

  SMART_ASSERT(file1.Open());
  SMART_ASSERT(file2.Open());

  {
    size_t num_nodes = file1.Num_Nodes();
    SMART_ASSERT(num_nodes == file2.Num_Nodes());

    node_map.resize(num_nodes);
    file1.Load_Node_Map();
    file2.Load_Node_Map();
    const INT *node_id_map1 = file1.Get_Node_Map();
    const INT *node_id_map2 = file2.Get_Node_Map();

    if (!internal_compute_maps(node_map, node_id_map1, node_id_map2, num_nodes, "node")) {
      node_map.clear();
    }
  }

  {
    size_t num_elmts = file1.Num_Elements();
    SMART_ASSERT(num_elmts == file2.Num_Elements());
    elmt_map.resize(num_elmts);
    file1.Load_Element_Map();
    file2.Load_Element_Map();
    const INT *elem_id_map1 = file1.Get_Element_Map();
    const INT *elem_id_map2 = file2.Get_Element_Map();

    if (!internal_compute_maps(elmt_map, elem_id_map1, elem_id_map2, num_elmts, "element")) {
      elmt_map.clear();
    }
  }
}

template <typename INT>
void Dump_Maps(const std::vector<INT> &node_map, const std::vector<INT> &elmt_map,
               ExoII_Read<INT> &file1)
{
  size_t ijk;
  fmt::print("\n=== node number map (file1 -> file2) local ids\n");
  bool one_to_one = true;
  if (!node_map.empty()) {
    for (ijk = 0; ijk < file1.Num_Nodes(); ++ijk) {
      if ((INT)ijk != node_map[ijk]) {
        one_to_one = false;
        break;
      }
    }
  }
  if (!one_to_one) {
    for (ijk = 0; ijk < file1.Num_Nodes(); ++ijk) {
      fmt::print("{} -> {}\n", (ijk + 1), (node_map[ijk] + 1));
    }
  }
  else {
    fmt::print(" *** Node map is one-to-one\n");
  }

  fmt::print("\n=== element number map (file1 -> file2) local ids\n");
  one_to_one = true;
  if (!elmt_map.empty()) {
    for (ijk = 0; ijk < file1.Num_Elements(); ++ijk) {
      if ((INT)ijk != elmt_map[ijk]) {
        one_to_one = false;
        break;
      }
    }
  }
  if (!one_to_one) {
    for (ijk = 0; ijk < file1.Num_Elements(); ++ijk) {
      fmt::print("{} -> {}\n", (ijk + 1), (elmt_map[ijk] + 1));
    }
  }
  else {
    fmt::print(" *** Element map is one-to-one\n");
  }
  fmt::print("===\n");
}

namespace {
  template <typename INT>
  void Compute_Node_Map(std::vector<INT> &node_map, ExoII_Read<INT> &file1, ExoII_Read<INT> &file2)
  {
    // This function is called if and only if there are nodes that were
    // not matched in the Compute_Map function.  This is typically the
    // case if there are 'free' nodes which are not connected to any
    // elements.

    size_t           num_nodes = file1.Num_Nodes();
    std::vector<INT> mapped_2(num_nodes, -1);

    // Cannot ignore the comparisons, so make sure the coord_tol_type
    // is not -1 which is "ignore"
    ToleranceMode save_tolerance_type = interFace.coord_tol.type;
    if (save_tolerance_type == ToleranceMode::IGNORE_) {
      interFace.coord_tol.type = ToleranceMode::ABSOLUTE_;
    }

    // Find unmapped nodes in file2; count the unmapped nodes in file_1.
    // The code below sets the 'mapped_2' entry to 1 for each node in
    // file2 which has been mapped to a node in file1
    size_t count_1 = 0;
    for (size_t i = 0; i < num_nodes; i++) {
      if (node_map[i] != -1) {
        mapped_2[node_map[i]] = 1;
      }
      else {
        count_1++;
      }
    }

    // Get list of all unmapped nodes in file1 and file2. A file1
    // unmapped node will have a '-1' entry in 'node_map' and a file2
    // unmapped node will have a '-1' entry in 'mapped_2'.  Reuse the
    // 'mapped_2' array to hold the list.
    std::vector<INT> mapped_1(count_1);
    size_t           count_2 = 0;
    count_1                  = 0;
    for (size_t i = 0; i < num_nodes; i++) {
      if (node_map[i] == -1) {
        mapped_1[count_1++] = i;
      }
      if (mapped_2[i] == -1) {
        mapped_2[count_2++] = i;
      }
    }

    // check that umnapped node counts are equal.  If not, output
    // message and exit.
    if (count_1 != count_2) {
      Error(fmt::format("Files are different (free node count in file1 is "
                        "{} but file2 free node count is {})\n",
                        fmt::group_digits(count_1), fmt::group_digits(count_2)));
    }

    // Now, need to match all nodes in 'mapped_1' with nodes in
    // 'mapped_2'
    // Get pointer to coordinates...
    const auto *x1_f = file1.X_Coords();
    const auto *y1_f = file1.Y_Coords();
    const auto *z1_f = file1.Z_Coords();

    const auto *x2_f = file2.X_Coords();
    const auto *y2_f = file2.Y_Coords();
    const auto *z2_f = file2.Z_Coords();

    // For now, we will try a brute force matching with the hopes that
    // the number of unmatched nodes is 'small'.  If this proves to be a
    // bottleneck, replace with a better algorithm; perhaps the sorted
    // matching process used in gjoin...
    size_t matched = 0;
    int    dim     = file1.Dimension();
    size_t j;
    for (size_t i = 0; i < count_1; i++) {
      size_t id_1 = mapped_1[i];
      for (j = 0; j < count_2; j++) {
        if (mapped_2[j] >= 0) {
          size_t id_2 = mapped_2[j];
          if ((dim == 1 && !interFace.coord_tol.Diff(x1_f[id_1], x2_f[id_2])) ||
              (dim == 2 && !interFace.coord_tol.Diff(x1_f[id_1], x2_f[id_2]) &&
               !interFace.coord_tol.Diff(y1_f[id_1], y2_f[id_2])) ||
              (dim == 3 && !interFace.coord_tol.Diff(x1_f[id_1], x2_f[id_2]) &&
               !interFace.coord_tol.Diff(y1_f[id_1], y2_f[id_2]) &&
               !interFace.coord_tol.Diff(z1_f[id_1], z2_f[id_2]))) {
            node_map[id_1] = id_2;
            mapped_2[j]    = -1;
            matched++;
            break;
          }
        }
      }
    }

    // Check that all nodes were matched.
    if (matched != count_1) {
      Error(fmt::format("Unable to match all free nodes in the model.  There are {}"
                        " unmatched nodes remaining.\n",
                        fmt::group_digits(count_1 - matched)));
    }
    interFace.coord_tol.type = save_tolerance_type;
  }

  template <typename INT>
  int64_t Find(double x0, double y0, double z0, const std::vector<double> &x,
               const std::vector<double> &y, const std::vector<double> &z,
               const std::vector<INT> &id, int dim, bool ignore_dups)
  {
    if (x.empty()) {
      return -1;
    }

    // Cannot ignore the comparisons, so make sure the coord_tol_type
    // is not -1 which is "ignore"
    ToleranceMode save_tolerance_type = interFace.coord_tol.type;
    if (save_tolerance_type == ToleranceMode::IGNORE_) {
      interFace.coord_tol.type = ToleranceMode::ABSOLUTE_;
    }

    // Find the index such that x0 > x[0,1,...,low-1] and x0 >= x[low]
    // where x[N] is infinity.
    auto   N    = x.size();
    size_t low  = 0;
    size_t high = N;
    while (low < high) {
      size_t mid = (low + high) / 2;
      if (x[id[mid]] < x0) {
        low = mid + 1;
      }
      else {
        high = mid;
      }
    }

    int64_t i = low == N ? N - 1 : low; // Make sure index falls within array bounds.

    if (i == 0 && interFace.coord_tol.Diff(x[id[i]], x0)) {
      // Could not find an index within tolerance on x coordinate.
      return -1;
    }

    // Drop to first index before which the tolerance fails.
    while (i > 0 && !interFace.coord_tol.Diff(x[id[i - 1]], x0)) {
      --i;
    }

    // Search until tolerance between the x coordinate fails or a match is found.
    // If a match is found, the loop continues in order to check for dups.

    int64_t index = -1;
    do {
      if (dim == 1 || (dim == 2 && !interFace.coord_tol.Diff(y[id[i]], y0)) ||
          (dim == 3 && !interFace.coord_tol.Diff(y[id[i]], y0) &&
           !interFace.coord_tol.Diff(z[id[i]], z0))) {
        if (index >= 0) {
          if (ignore_dups) {
            return index;
          }

          double x1 = x[id[i]];
          double y1 = dim > 1 ? y[id[i]] : 0.0;
          double z1 = dim > 2 ? z[id[i]] : 0.0;

          double x2 = x[id[index]];
          double y2 = dim > 1 ? y[id[index]] : 0.0;
          double z2 = dim > 2 ? z[id[index]] : 0.0;

          Warning(fmt::format("Two elements in file 2 have the same midpoint (within tolerance).\n"
                              "\tLocal element {} at ({}, {}, {}) and\n"
                              "\tLocal element {} at ({}, {}, {})\n"
                              "\tNo unique element mapping possible.\n",
                              fmt::group_digits(id[i] + 1), x1, y1, z1,
                              fmt::group_digits(id[index] + 1), x2, y2, z2));
          return -1;
        }

        index = i;
      }
    } while (++i < (int64_t)N && !interFace.coord_tol.Diff(x[id[i]], x0));

    interFace.coord_tol.type = save_tolerance_type;
    return index;
  }

  inline double dist_sqrd(double x1, double x2) { return (x2 - x1) * (x2 - x1); }

  inline double dist_sqrd(double x1, double y1, double x2, double y2)
  {
    double d1 = x2 - x1;
    d1 *= d1;
    double d2 = y2 - y1;
    d2 *= d2;
    return (d1 + d2);
  }

  inline double dist_sqrd(double x1, double y1, double z1, double x2, double y2, double z2)
  {
    double d1 = x2 - x1;
    d1 *= d1;
    double d2 = y2 - y1;
    d2 *= d2;
    double d3 = z2 - z1;
    d3 *= d3;
    return (d1 + d2 + d3);
  }

  double find_range(const double *x, size_t num_nodes)
  {
    auto range = std::minmax_element(x, x + num_nodes);
    return *range.second - *range.first;
  }
} // namespace

template <typename INT> double Find_Min_Coord_Sep(ExoII_Read<INT> &file)
{
  size_t num_nodes = file.Num_Nodes();
  if (num_nodes < 2) {
    return 0.0;
  }

  file.Load_Nodal_Coordinates();
  const auto *x = file.X_Coords();
  const auto *y = file.Y_Coords();
  const auto *z = file.Z_Coords();

  std::vector<INT> indx(num_nodes);
  std::iota(indx.begin(), indx.end(), 0);

  // Find coordinate with largest range...
  const double *r     = x;
  double        range = find_range(x, num_nodes);
  if (file.Dimension() > 1) {
    double yrange = find_range(y, num_nodes);
    if (yrange > range) {
      range = yrange;
      r     = y;
    }
  }

  if (file.Dimension() > 2) {
    double zrange = find_range(z, num_nodes);
    if (zrange > range) {
      range = zrange;
      r     = z;
    }
  }

  // Sort based on coordinate with largest range...
  index_qsort(r, Data(indx), num_nodes);

  double min = DBL_MAX;
  switch (file.Dimension()) {
  case 1: {
    for (size_t i = 0; i < num_nodes; i++) {
      for (size_t j = i + 1; j < num_nodes; j++) {
        double tmp = dist_sqrd(x[indx[i]], x[indx[j]]);
        if (tmp >= min) {
          break;
        }
        min = tmp;
      }
    }
    break;
  }

  case 2: {
    for (size_t i = 0; i < num_nodes; i++) {
      for (size_t j = i + 1; j < num_nodes; j++) {
        double delr = dist_sqrd(r[indx[i]], r[indx[j]]);
        if (delr > min) {
          break;
        }
        double tmp = dist_sqrd(x[indx[i]], y[indx[i]], x[indx[j]], y[indx[j]]);
        min        = min < tmp ? min : tmp;
      }
    }
    break;
  }

  case 3: {
    for (size_t i = 0; i < num_nodes; i++) {
      for (size_t j = i + 1; j < num_nodes; j++) {
        double delr = dist_sqrd(r[indx[i]], r[indx[j]]);
        if (delr > min) {
          break;
        }
        double tmp =
            dist_sqrd(x[indx[i]], y[indx[i]], z[indx[i]], x[indx[j]], y[indx[j]], z[indx[j]]);
        min = min < tmp ? min : tmp;
      }
    }
    break;
  }
  }
  return sqrt(min);
}

template <typename INT>
bool Compare_Maps_Internal(const std::vector<INT> &entity_map, bool partial_flag,
                           const INT *entity_id_map1, const INT *entity_id_map2,
                           size_t num_entities1, size_t num_entities2, const char *type)
{
  bool diff       = false;
  int  warn_count = 0;

  if (!entity_map.empty()) {
    if (!interFace.dump_mapping) {
      // There is a map between file1 and file2, but all entities are
      // used in both files.
      for (size_t i = 0; i < num_entities1; i++) {
        size_t idx = entity_map[i];
        if (idx >= num_entities2) {
          continue;
        }
        if (entity_id_map1[i] != entity_id_map2[idx]) {
          if (!(entity_id_map2[idx] == 0 &&
                partial_flag)) { // Don't output diff if non-matched and partial
            fmt::print(stderr,
                       "exodiff: WARNING .. The local {} {} with global id {} in file1 has the "
                       "global id "
                       "{} in file2.\n",
                       type, i + 1, entity_id_map1[i], entity_id_map2[idx]);
            diff = true;
            warn_count++;
            if (warn_count >= interFace.max_warnings) {
              fmt::print(stderr, "exodiff: WARNING .. Too many warnings, skipping remainder...\n");
              break;
            }
          }
        }
      }
    }
  }
  else {
    // No entity mapping between file1 and file2 -- do a straight compare.
    for (size_t i = 0; i < num_entities1; i++) {
      if (i >= num_entities2) {
        break;
      }
      if (entity_id_map1[i] != entity_id_map2[i]) {
        if (!(entity_id_map2[i] == 0 &&
              partial_flag)) { // Don't output diff if non-matched and partial
          fmt::print(
              stderr,
              "exodiff: WARNING .. The local {} {} with global id {} in file1 has the global id "
              "{} in file2.\n",
              type, i + 1, entity_id_map1[i], entity_id_map2[i]);
          diff = true;
          warn_count++;
          if (warn_count >= interFace.max_warnings) {
            fmt::print(stderr, "exodiff: WARNING .. Too many warnings, skipping remainder...\n");
            break;
          }
        }
      }
    }
  }
  return diff;
}

template <typename INT>
bool Compare_Maps(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, const std::vector<INT> &node_map,
                  const std::vector<INT> &elmt_map, bool partial_flag)
{
  // Check whether the node and element number maps from both file1
  // and file2 match which indicates that we are comparing the same
  // element and node in each file.

  size_t num_nodes1 = file1.Num_Nodes();
  size_t num_nodes2 = file2.Num_Nodes();

  // NOTE: file1 maps are already loaded...
  file2.Load_Node_Map();

  const INT *node_id_map1 = file1.Get_Node_Map();
  const INT *node_id_map2 = file2.Get_Node_Map();

  bool diff_nodes = Compare_Maps_Internal(node_map, partial_flag, node_id_map1, node_id_map2,
                                          num_nodes1, num_nodes2, "node");
  file2.Free_Node_Map();

  size_t num_elmts1 = file1.Num_Elements();
  size_t num_elmts2 = file2.Num_Elements();

  // NOTE: file1 maps are already loaded...
  file2.Load_Element_Map();

  const INT *elem_id_map1 = file1.Get_Element_Map();
  const INT *elem_id_map2 = file2.Get_Element_Map();

  bool diff_elems = Compare_Maps_Internal(elmt_map, partial_flag, elem_id_map1, elem_id_map2,
                                          num_elmts1, num_elmts2, "element");

  file2.Free_Element_Map();

  if (diff_nodes || diff_elems) {
    fmt::print("\n");
  }
  return diff_nodes || diff_elems;
}

template void Compute_Maps(std::vector<int> &node_map, std::vector<int> &elmt_map,
                           ExoII_Read<int> &file1, ExoII_Read<int> &file2);
template bool Compare_Maps(ExoII_Read<int> &file1, ExoII_Read<int> &file2,
                           const std::vector<int> &node_map, const std::vector<int> &elmt_map,
                           bool partial_flag);

template void   Compute_Partial_Maps(std::vector<int> &node_map, std::vector<int> &elmt_map,
                                     ExoII_Read<int> &file1, ExoII_Read<int> &file2);
template void   Compute_FileId_Maps(std::vector<int> &node_map, std::vector<int> &elmt_map,
                                    ExoII_Read<int> &file1, ExoII_Read<int> &file2);
template void   Dump_Maps(const std::vector<int> &node_map, const std::vector<int> &elmt_map,
                          ExoII_Read<int> &file1);
template double Find_Min_Coord_Sep(ExoII_Read<int> &file);

template void Compute_Maps(std::vector<int64_t> &node_map, std::vector<int64_t> &elmt_map,
                           ExoII_Read<int64_t> &file1, ExoII_Read<int64_t> &file2);
template bool Compare_Maps(ExoII_Read<int64_t> &file1, ExoII_Read<int64_t> &file2,
                           const std::vector<int64_t> &node_map,
                           const std::vector<int64_t> &elmt_map, bool partial_flag);

template void Compute_Partial_Maps(std::vector<int64_t> &node_map, std::vector<int64_t> &elmt_map,
                                   ExoII_Read<int64_t> &file1, ExoII_Read<int64_t> &file2);
template void Compute_FileId_Maps(std::vector<int64_t> &node_map, std::vector<int64_t> &elmt_map,
                                  ExoII_Read<int64_t> &file1, ExoII_Read<int64_t> &file2);
template void Dump_Maps(const std::vector<int64_t> &node_map, const std::vector<int64_t> &elmt_map,
                        ExoII_Read<int64_t> &file1);
template double Find_Min_Coord_Sep(ExoII_Read<int64_t> &file);
