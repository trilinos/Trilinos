// Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once
#include "exo_read.h"

enum class MapType { FILE_ORDER = 0, PARTIAL, USE_FILE_IDS, DISTANCE };

template <typename INT>
void Compute_Maps(std::vector<INT> &node_map, std::vector<INT> &elmt_map, Exo_Read<INT> &file1,
                  Exo_Read<INT> &file2);

template <typename INT>
void Compute_Partial_Maps(std::vector<INT> &node_map, std::vector<INT> &elmt_map,
                          Exo_Read<INT> &file1, Exo_Read<INT> &file2);

template <typename INT>
void Compute_FileId_Maps(std::vector<INT> &node_map, std::vector<INT> &elmt_map,
                         Exo_Read<INT> &file1, Exo_Read<INT> &file2);

template <typename INT>
void Dump_Maps(const std::vector<INT> &node_map, const std::vector<INT> &elmt_map,
               Exo_Read<INT> &file1);

template <typename INT>
bool Check_Maps(const std::vector<INT> &node_map, const std::vector<INT> &elmt_map,
                const Exo_Read<INT> &file1, const Exo_Read<INT> &file2);

template <typename INT>
bool Compare_Maps(Exo_Read<INT> &file1, Exo_Read<INT> &file2, const std::vector<INT> &node_map,
                  const std::vector<INT> &elmt_map, bool partial_flag);
