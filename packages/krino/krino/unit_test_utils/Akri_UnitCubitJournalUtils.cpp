// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_UnitCubitJournalUtils.hpp>
#include <regex>
#include <sstream>

namespace krino {

template<size_t N>
static std::array<int,N> journal_node_indices(const std::string & line)
{
  static const std::regex nodeRef(R"(\{krinoNode0\+(\d+)\})");
  std::array<int,N> indices;
  size_t n = 0;
  for (auto it = std::sregex_iterator(line.begin(), line.end(), nodeRef); it != std::sregex_iterator() && n < N; ++it)
    indices[n++] = std::stoi((*it)[1].str()) - 1;
  return indices;
}

static stk::math::Vector3d journal_coordinates(const std::string & values)
{
  std::istringstream in(values);
  stk::math::Vector3d x;
  in >> x[0] >> x[1] >> x[2];
  return x;
}

ParsedCubitJournal parse_cubit_journal(const std::string & journal)
{
  ParsedCubitJournal parsed;
  std::istringstream lines(journal);
  std::string line;
  auto starts_with = [&](const std::string & start) { return line.rfind(start, 0) == 0; };
  while (std::getline(lines, line))
  {
    if (starts_with("create node location "))
      parsed.nodes.push_back(journal_coordinates(line.substr(21)));
    else if (starts_with("create vertex "))
      parsed.vertices.push_back(journal_coordinates(line.substr(14)));
    else if (starts_with("create edge node "))
      parsed.edges.push_back(journal_node_indices<2>(line));
    else if (starts_with("create tri node "))
      parsed.tris.push_back(journal_node_indices<3>(line));
    else if (starts_with("create tet node "))
      parsed.tets.push_back(journal_node_indices<4>(line));
  }
  return parsed;
}

double total_area_of_tris(const ParsedCubitJournal & journal)
{
  double area = 0.;
  for (auto && tri : journal.tris)
    area += 0.5*Cross(journal.nodes[tri[1]]-journal.nodes[tri[0]], journal.nodes[tri[2]]-journal.nodes[tri[0]]).length();
  return area;
}

double total_length_of_edges(const ParsedCubitJournal & journal)
{
  double length = 0.;
  for (auto && edge : journal.edges)
    length += (journal.nodes[edge[1]]-journal.nodes[edge[0]]).length();
  return length;
}

}
