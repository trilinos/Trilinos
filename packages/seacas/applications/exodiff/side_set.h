// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef SIDE_SET_H
#define SIDE_SET_H

#include "exo_entity.h"
#include <iostream>

template <typename INT> class ExoII_Read;

template <typename INT> class Side_Set : public Exo_Entity
{
public:
  Side_Set();
  Side_Set(int file_id, size_t id);
  Side_Set(int file_id, size_t id, size_t ns, size_t ndf = 0);
  ~Side_Set() override;

  void                apply_map(const INT *elmt_map);
  const INT *         Elements() const;
  const INT *         Sides() const;
  std::pair<INT, INT> Side_Id(size_t position) const;
  size_t              Side_Index(size_t position) const;

  std::pair<INT, INT> Distribution_Factor_Range(size_t side) const;
  const double *      Distribution_Factors() const;
  void                Free_Distribution_Factors() const;

  int    Check_State() const;
  size_t Distribution_Factor_Count() const { return num_dist_factors; }

private:
  Side_Set(const Side_Set &);                  // Not written.
  const Side_Set &operator=(const Side_Set &); // Not written.

  void load_sides(const INT *elmt_map = nullptr) const;
  void load_df() const;
  void entity_load_params() override;

  EXOTYPE     exodus_type() const override;
  const char *label() const override { return "Sideset"; }
  const char *short_label() const override { return "sideset"; }

  size_t num_dist_factors{0};

  mutable INT *   elmts{nullptr};
  mutable INT *   sides{nullptr};
  mutable INT *   sideIndex{nullptr};
  mutable INT *   dfIndex{nullptr};
  mutable double *dist_factors{nullptr};

  friend class ExoII_Read<INT>;
};

#endif
