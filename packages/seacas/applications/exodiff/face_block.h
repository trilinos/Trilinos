// Copyright(C) 1999-2022 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef FACE_BLOCK_H
#define FACE_BLOCK_H

#include "exo_entity.h"
#include <iostream>

template <typename INT> class ExoII_Read;

template <typename INT> class Face_Block : public Exo_Entity
{
public:
  Face_Block();
  Face_Block(int file_id, size_t id);
  Face_Block(int file_id, size_t id, size_t ne);
  ~Face_Block() override;

  size_t Face_Index(size_t position) const;

  int Check_State() const;

private:
  Face_Block(const Face_Block &)                  = delete; // Not written.
  const Face_Block &operator=(const Face_Block &) = delete; // Not written.

  void entity_load_params() override;

  EXOTYPE     exodus_type() const override;
  const char *label() const override { return "Faceblock"; }
  const char *short_label() const override { return "faceblock"; }

  std::string elmt_type{};
  int         num_faces_per_elmt{-1};

  friend class ExoII_Read<INT>;
};

#endif
