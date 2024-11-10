// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details
#pragma once

#include <exodusII.h>
#include <iostream>
#include <string>
#include <vector>

#if defined(EX_API_VERS_NODOT)
#if EX_API_VERS_NODOT > 467
using EXOTYPE = ex_entity_type;
#else
using EXOTYPE = int;
#endif
#else
using EXOTYPE = int;
#endif

template <typename INT> class ExoII_Read;

class Exo_Entity
{
public:
  Exo_Entity() = default;
  Exo_Entity(int file_id, size_t id);
  Exo_Entity(int file_id, size_t id, size_t nnodes);
  virtual ~Exo_Entity();
  Exo_Entity(const Exo_Entity &)                  = delete;
  const Exo_Entity &operator=(const Exo_Entity &) = delete;

  virtual size_t Size() const { return numEntity; }

  size_t Id() const { return id_; }
  size_t Index() const { return index_; }

  void initialize(int file_id, size_t id);

  bool        is_valid_var(size_t var_index) const;
  size_t      var_count() const { return numVars; }
  std::string Load_Results(int time_step, int var_index);
  std::string Load_Results(int t1, int t2, double proportion, int var_index); // Interpolation

  const double *Get_Results(int var_index) const;
  void          Free_Results();

  int           attr_count() const { return numAttr; }
  std::string   Load_Attributes(int attr_index);
  const double *Get_Attributes(int attr_index) const;
  void          Free_Attributes();

  const std::string              &Get_Attribute_Name(int attr_index) const;
  const std::string              &Name() const { return name_; }
  const std::vector<std::string> &Attribute_Names() const { return attributeNames; }
  int                             Find_Attribute_Index(const std::string &name) const;

  // Return "Element Block", "Nodeset", "Sideset, depending on underlying type.
  virtual const char *label() const = 0;

  // Return "block", "nodelist", "surface", depending on underlying type.
  virtual const char *short_label() const = 0;

  bool generatedName_{true};

protected:
  std::string  name_{};
  int          fileId{-1};
  ex_entity_id id_{EX_INVALID_ID};
  size_t       index_{0};    // 0-offset index into Exodus nodeset list.
  size_t       numEntity{0}; // Number of items (nodes, sides, elements)

private:
  // Return EX_ELEM_BLOCK, EX_NODE_SET, ... of underlying type
  virtual EXOTYPE exodus_type() const = 0;

  virtual int  Check_State() const  = 0;
  virtual void entity_load_params() = 0;
  void         internal_load_params();

  void get_truth_table() const;

  mutable int *truth_{nullptr};        // Array; holds local truth table for this entity
  int          currentStep{0};         // Time step number of the current results.
  int          numVars{0};             // Total number of variables in the file.
  double     **results_{nullptr};      // Array of pointers (length numVars)
                                       // to arrays of results (length num_entity).
  int                   numAttr{0};    // Total number of attributes in the file.
  std::vector<double *> attributes_{}; // Array of pointers (length numAttr)
                                       // to arrays of attributes (length num_entity).
  std::vector<std::string> attributeNames{};

  template <typename INT> friend class ExoII_Read;
};
