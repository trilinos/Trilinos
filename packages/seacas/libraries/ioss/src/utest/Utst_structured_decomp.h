/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef Utst_structured_decomp_h
#define Utst_structured_decomp_h

#include <cgns/Iocgns_StructuredZoneData.h>
#include <vector>

void cleanup(std::vector<Iocgns::StructuredZoneData *> &zones);
void check_split_assign(std::vector<Iocgns::StructuredZoneData *> &zones,
                        double load_balance_tolerance, size_t proc_count, double min_toler = 0.9,
                        double max_toler = 1.0);

#endif
