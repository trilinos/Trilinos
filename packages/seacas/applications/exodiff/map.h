// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#ifndef EXODIFF_MAP_H
#define EXODIFF_MAP_H
#include "exoII_read.h"

enum MAP_TYPE_enum { FILE_ORDER = 0, PARTIAL, USE_FILE_IDS, DISTANCE };

template <typename INT>
void Compute_Maps(INT *&node_map, INT *&elmt_map, ExoII_Read<INT> &file1, ExoII_Read<INT> &file2);

template <typename INT>
void Compute_Partial_Maps(INT *&node_map, INT *&elmt_map, ExoII_Read<INT> &file1,
                          ExoII_Read<INT> &file2);

template <typename INT>
void Compute_FileId_Maps(INT *&node_map, INT *&elmt_map, ExoII_Read<INT> &file1,
                         ExoII_Read<INT> &file2);

template <typename INT>
void Dump_Maps(const INT *node_map, const INT *elmt_map, ExoII_Read<INT> &file1);

template <typename INT>
bool Check_Maps(const INT *node_map, const INT *elmt_map, const ExoII_Read<INT> &file1,
                const ExoII_Read<INT> &file2);

template <typename INT>
bool Compare_Maps(ExoII_Read<INT> &file1, ExoII_Read<INT> &file2, const INT *node_map,
                  const INT *elmt_map, bool partial_flag);

#endif
