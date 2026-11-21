// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_STLFileParsing.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_Facet.hpp>

#include <stk_util/util/ReportHandler.hpp>

#include <cctype>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

namespace krino
{

static void flip_facet_with_negative_sign(std::array<stk::math::Vector3d, 3> & pts, const int sign)
{
  if (sign == -1)
  {
    std::swap(pts[1][0], pts[2][0]);
    std::swap(pts[1][1], pts[2][1]);
    std::swap(pts[1][2], pts[2][2]);
  }
}

static bool is_finite_normal(const stk::math::Vector3d normal)
{
  bool finiteFlag = true;
  if (!std::isfinite(normal[0]) || !std::isfinite(normal[1]) || !std::isfinite(normal[2]))
  {
    finiteFlag = false;
    krinolog << "STL file has an invalid facet, which is skipped." << std::endl;
    krinolog << "  facet normal " << normal[0] << " " << normal[1] << " " << normal[2] << std::endl;
  }
  return finiteFlag;
}

static std::string trim_whitespace(const std::string & s)
{
  size_t start = 0;
  while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start])))
  {
    ++start;
  }

  size_t end = s.size();
  while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1])))
  {
    --end;
  }

  return s.substr(start, end - start);
}

static bool starts_with(const std::string & s, const std::string & prefix)
{
  return s.size() >= prefix.size() && s.compare(0, prefix.size(), prefix) == 0;
}

static void parse_ascii_vector3d(const std::string & s, stk::math::Vector3d & v)
{
  std::istringstream iss(s);
  if (!(iss >> v[0] >> v[1] >> v[2])) throw std::runtime_error("Failed to parse 3D vector");
}

static void parse_ascii_prefix_and_vector3d(
    const std::string & line, const std::string & prefix, stk::math::Vector3d & v)
{
  if (!starts_with(line, prefix)) throw std::runtime_error("Expected '" + prefix + " x y z'");

  parse_ascii_vector3d(line.substr(prefix.length()), v);
}

static void parse_ascii_prefix_and_vector3d(
    std::istream & in, const std::string & prefix, stk::math::Vector3d & v)
{
  std::string line;
  if (!std::getline(in, line))
    throw std::runtime_error("Unexpected EOF before '" + prefix + " x y z'");
  line = trim_whitespace(line);
  parse_ascii_prefix_and_vector3d(line, prefix, v);
}

static void parse_and_expect_token(std::istream & in, const std::string & expectedToken)
{
  std::string line;
  if (!std::getline(in, line))
    throw std::runtime_error("Unexpected EOF before " + expectedToken + "'");
  line = trim_whitespace(line);
  if (line != expectedToken) throw std::runtime_error("Expected '" + expectedToken + "'");
}

static void parse_ascii_facet(std::istream & in,
    const std::string & firstLine,
    stk::math::Vector3d & facetNormal,
    std::array<stk::math::Vector3d, 3> & facetPts)
{
  parse_ascii_prefix_and_vector3d(firstLine, "facet normal", facetNormal);

  parse_and_expect_token(in, "outer loop");

  for (int i = 0; i < 3; ++i)
    parse_ascii_prefix_and_vector3d(in, "vertex", facetPts[i]);

  parse_and_expect_token(in, "endloop");
  parse_and_expect_token(in, "endfacet");
}

static void read_ascii_facet(std::istream & in,
    const std::string & firstLine,
    const stk::math::Vector3d & scale,
    const int sign,
    std::array<stk::math::Vector3d, 3> & facetPts,
    bool & isFacetValid)
{
  stk::math::Vector3d facetNormal;
  parse_ascii_facet(in, firstLine, facetNormal, facetPts);

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      facetPts[i][j] *= scale[j];
  isFacetValid = is_finite_normal(facetNormal) && !Facet3d::is_degenerate(facetPts);
  flip_facet_with_negative_sign(facetPts, sign);
}

static float read_binary_float(std::istream & input)
{
  float val = 0;
  input.read((char *)&val, sizeof(val));
  return val;
}

static stk::math::Vector3d read_binary_vector(std::istream & input)
{
  float x = read_binary_float(input);
  float y = read_binary_float(input);
  float z = read_binary_float(input);
  return stk::math::Vector3d(x, y, z);
}

static void read_binary_facet(std::istream & input,
    const stk::math::Vector3d & scale,
    const int sign,
    std::array<stk::math::Vector3d, 3> & facetPts,
    bool & isFacetValid)
{
  const stk::math::Vector3d normal = read_binary_vector(input);

  for (auto && pt : facetPts)
  {
    const stk::math::Vector3d vec = read_binary_vector(input);
    pt = stk::math::Vector3d(vec[0] * scale[0], vec[1] * scale[1], vec[2] * scale[2]);
  }

  char dummy[2];
  input.read(dummy, 2);

  isFacetValid = is_finite_normal(normal) && !Facet3d::is_degenerate(facetPts);
  flip_facet_with_negative_sign(facetPts, sign);
}

STLFacetReader::STLFacetReader(
    std::istream & input, const int sign, const stk::math::Vector3d & scale)
    : myInput(input), mySign(sign), myScale(scale)
{
}

void STLFacetReader::read_header()
{
  std::string line;
  std::getline(myInput, line);
  myIsAscii = starts_with(line, "solid");

  if (!myIsAscii)
  {
    char header_info[80] = "";
    myInput.clear();
    myInput.seekg(0);
    myInput.read(header_info, 80);
    STK_ThrowErrorMsgIf(
        !myInput.good(), "Problem reading STL file, cannot read binary file header.");
    myInput.read((char *)&myNumBinaryFacets, sizeof(myNumBinaryFacets));
    STK_ThrowErrorMsgIf(!myInput.good(), "Problem reading STL file, cannot read number of facets.");
    krinolog << "Reading binary STL file with " << myNumBinaryFacets << " facets."
             << stk::diag::dendl;
  }
}

bool STLFacetReader::read_next_valid_facet(std::array<stk::math::Vector3d, 3> & facetPts)
{
  bool isFacetValid = false;
  if (!myIsAscii)
  {
    while (myNumBinaryFacetsRead < myNumBinaryFacets)
    {
      ++myNumBinaryFacetsRead;
      read_binary_facet(myInput, myScale, mySign, facetPts, isFacetValid);
      if (isFacetValid) return true;
    }
    return false;
  }

  std::string line;
  while (std::getline(myInput, line))
  {
    line = trim_whitespace(line);
    if (starts_with(line, "facet"))
    {
      read_ascii_facet(myInput, line, myScale, mySign, facetPts, isFacetValid);
      if (isFacetValid) return true;
    }
    else
    {
      if (starts_with(line, "solid"))
        krinolog << "STL file contains multiple solids which are being read together. This will "
                    "not behave well if the surfaces intersect."
                 << stk::diag::dendl;
      else if (!line.empty() && !starts_with(line, "endsolid"))
        throw std::runtime_error("Unexpected line '" + line + "'");
    }
  }
  return false;
}

} // namespace krino
