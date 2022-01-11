// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef AKRI_INTERFACEID_H_
#define AKRI_INTERFACEID_H_

#include <set>
#include <map>
#include <vector>
#include <ostream>

#include <Akri_TypeDefs.hpp>

namespace krino {

class InterfaceID;
typedef std::map<InterfaceID, double> CrossingMap;
typedef std::map<InterfaceID, int> CrossingSignMap;

class InterfaceID
{
public:
  InterfaceID() : ls_1(-1), ls_2(-1) {}
  InterfaceID(const int first, const int second) : ls_1(std::min(first, second)),
                                                   ls_2(std::max(first, second)) {}

  bool is_single_ls() const { return ls_1 == ls_2; }
  int first_ls() const { return ls_1; }
  int second_ls() const {return ls_2; }
  void fill_sorted_domains(std::vector<int> & sortedDomains) const;

  static std::vector<InterfaceID> all_phase_pairs(const std::set<int> & phases);
private:
  int ls_1;
  int ls_2;
};

inline bool operator<(const InterfaceID lhs, const InterfaceID rhs)
{
  if(lhs.first_ls() < rhs.first_ls()) return true;
  if(rhs.first_ls() < lhs.first_ls()) return false;
  if(lhs.second_ls() < rhs.second_ls()) return true;
  return false;
}

inline bool operator==(const InterfaceID lhs, const InterfaceID rhs)
{
  return (lhs.first_ls() == rhs.first_ls()) && (lhs.second_ls() == rhs.second_ls());
}

inline bool operator!=(const InterfaceID lhs, const InterfaceID rhs)
{
  return !(lhs == rhs);
}

inline std::ostream & operator<<(std::ostream &os, const InterfaceID lhs)
{
  os << "InterfaceID(" << lhs.first_ls() << ", " << lhs.second_ls() << ")";
  return os;
}

inline std::vector<InterfaceID> InterfaceID::all_phase_pairs(const std::set<int> & phases)
{
  std::vector<InterfaceID> result;
  std::set<int>::const_iterator end_phase = phases.end();
  for(std::set<int>::const_iterator it = phases.begin(); it != end_phase; ++it)
  {
    std::set<int>::const_iterator it2 = it;
    for(++it2; it2 != end_phase; ++it2)
    {
      result.push_back(InterfaceID(*it, *it2));
    }
  }
  return result;
}

inline void InterfaceID::fill_sorted_domains(std::vector<int> & sortedDomains) const
{
  sortedDomains.clear();
  if (first_ls() == second_ls())
  {
    sortedDomains.push_back(first_ls());
  }
  else
  {
    sortedDomains.push_back(first_ls());
    sortedDomains.push_back(second_ls());
  }
}

} // namespace krino

#endif /* AKRI_INTERFACEID_H_ */
