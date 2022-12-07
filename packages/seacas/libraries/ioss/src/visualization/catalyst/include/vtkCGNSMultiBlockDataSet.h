// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef __VTK_CGNS_MULTI_BLOCK_DATA_SET_H
#define __VTK_CGNS_MULTI_BLOCK_DATA_SET_H

#include "vtkMultiBlockDataSet.h"
#include <map>

class vtkCGNSMultiBlockDataSet : public vtkMultiBlockDataSet
{
public:
  vtkTypeMacro(vtkCGNSMultiBlockDataSet, vtkMultiBlockDataSet);
  void PrintSelf(ostream &os, vtkIndent indent);

  static vtkCGNSMultiBlockDataSet *New();

  void CreateBase(int base_id, const std::string &base_name);

  void AddStructuredZoneData(int base_id, int zone_id, const std::string &zone_name,
                             const std::string &data_name, int ni, int nj, int nk, int comp_count,
                             bool is_cell_field, char field_suffix_separator, double *data,
                             int size);

protected:
  vtkCGNSMultiBlockDataSet();
  ~vtkCGNSMultiBlockDataSet();

private:
  vtkCGNSMultiBlockDataSet(const vtkCGNSMultiBlockDataSet &); // Not implemented.
  void operator=(const vtkCGNSMultiBlockDataSet &);           // Not implemented.

  struct base
  {
    int                base_location;
    std::map<int, int> zone_id_to_zone_location_map;
  };

  std::map<int, base> base_id_to_base_map;
};

#endif /* __VTK_CGNS_MULTI_BLOCK_DATA_SET_H */
