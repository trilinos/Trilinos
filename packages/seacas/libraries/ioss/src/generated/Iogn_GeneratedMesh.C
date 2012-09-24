// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
//         
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <generated/Iogn_GeneratedMesh.h>
#include <tokenize.h>

#include <cmath>
#include <cstring>
#include <cstdlib>

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <assert.h>

namespace Iogn {
  GeneratedMesh::GeneratedMesh(int64_t num_x, int64_t num_y, int64_t num_z,
			       int proc_count, int my_proc) :
    numX(num_x), numY(num_y), numZ(num_z), myNumZ(num_z), myStartZ(0),
    processorCount(proc_count), myProcessor(my_proc),
    offX(0), offY(0), offZ(0),
    sclX(1), sclY(1), sclZ(1),
    doRotation(false)
  {
    initialize();
  }

  GeneratedMesh::GeneratedMesh(const std::string &parameters,
			       int proc_count, int my_proc) :
    numX(0), numY(0), numZ(0), myNumZ(0), myStartZ(0),
    processorCount(proc_count), myProcessor(my_proc),
    offX(0), offY(0), offZ(0),
    sclX(1), sclY(1), sclZ(1),
    doRotation(false)
  {
    // Possible that the 'parameters' has the working directory path
    // prepended to the parameter list.  Strip off everything in front
    // of the last '/' (if any)...
    std::vector<std::string> params;
    Ioss::tokenize(parameters, "/", params);
    
    std::vector<std::string> groups;
    Ioss::tokenize(params[params.size()-1], "|+", groups);

    // First 'group' is the interval specification -- IxJxK
    std::vector<std::string> tokens;
    Ioss::tokenize(groups[0], "x", tokens);
    assert(tokens.size() == 3);
    numX = std::strtol(tokens[0].c_str(), NULL, 10);
    numY = std::strtol(tokens[1].c_str(), NULL, 10);
    numZ = std::strtol(tokens[2].c_str(), NULL, 10);

    initialize();
    parse_options(groups);

  }

  GeneratedMesh::~GeneratedMesh() {}

  void GeneratedMesh::initialize()
  {
    if (processorCount > numZ) {
      if (myProcessor == 0) {
	std::cerr << "ERROR: (Iogn::GeneratedMesh::initialize)\n"
		  << "       The number of mesh intervals in the Z direction (" << numZ << ")\n"
		  << "       must be at least as large as the number of processors (" << processorCount << ").\n"
		  << "       The current parameters do not meet that requirement. Execution will terminate.\n";
      }
      std::exit(EXIT_FAILURE);
    }

    if (processorCount > 1) {
      myNumZ = numZ / processorCount;
      if (myProcessor < (numZ % processorCount)) myNumZ++;

      // Determine myStartZ for this processor...
      size_t extra = numZ % processorCount;
      if (extra > myProcessor)
	extra = myProcessor;
      size_t per_proc  = numZ / processorCount;
      myStartZ = myProcessor * per_proc + extra;
    } else {
      myNumZ = numZ;
    }

    for (int i=0; i < 3; i++) {
      for (int j=0; j < 3; j++) {
	rotmat[i][j] = 0.0;
      }
      rotmat[i][i] = 1.0;
    }
  }

  int64_t GeneratedMesh::add_shell_block(ShellLocation loc)
  {
    shellBlocks.push_back(loc);
    return shellBlocks.size();
  }

  int64_t GeneratedMesh::add_nodeset(ShellLocation loc)
  {
    nodesets.push_back(loc);
    return nodesets.size();
  }

  int64_t GeneratedMesh::add_sideset(ShellLocation loc)
  {
    sidesets.push_back(loc);
    return sidesets.size();
  }

  void GeneratedMesh::set_bbox(double xmin, double ymin, double zmin,
			       double xmax, double ymax, double zmax)
  {
    // NOTE: All calculations are based on the currently
    // active interval settings. If scale or offset or zdecomp
    // specified later in the option list, you may not get the
    // desired bounding box.
    double x_range = xmax - xmin;
    double y_range = ymax - ymin;
    double z_range = zmax - zmin;

    sclX = x_range / static_cast<double>(numX);
    sclY = y_range / static_cast<double>(numY);
    sclZ = z_range / static_cast<double>(numZ);

    offX = xmin;
    offY = ymin;
    offZ = zmin;
  }

  void GeneratedMesh::set_scale(double scl_x, double scl_y, double scl_z)
  {
    sclX = scl_x;
    sclY = scl_y;
    sclZ = scl_z;
  }

  void GeneratedMesh::set_offset(double off_x, double off_y, double off_z)
  {
    offX = off_x;
    offY = off_y;
    offZ = off_z;
  }

  void GeneratedMesh::parse_options(const std::vector<std::string> &groups)
  {
    for (size_t i=1; i < groups.size(); i++) {
      std::vector<std::string> option;
      Ioss::tokenize(groups[i], ":", option);
      // option[0] is the type of the option and option[1] is the argument to the option.

      if (option[0] == "shell") {
	// Option of the form  "shell:xXyYzZ"
	// The argument specifies whether there is a shell block
	// at the location. 'x' is minX, 'X' is maxX, etc.
	size_t length = option[1].size();
	for (size_t j=0; j < length; j++) {
	  switch (option[1][j]) {
	  case 'x':
	    add_shell_block(MX);
	    break;
	  case 'X':
	    add_shell_block(PX);
	    break;
	  case 'y':
	    add_shell_block(MY);
	    break;
	  case 'Y':
	    add_shell_block(PY);
	    break;
	  case 'z':
	    add_shell_block(MZ);
	    break;
	  case 'Z':
	    add_shell_block(PZ);
	    break;
	  default:
	    std::cerr << "ERROR: Unrecognized shell location option '"
		      << option[1][j]
		      << "'.";
	  }
	}
      }
      else if (option[0] == "nodeset") {
	// Option of the form  "nodeset:xXyYzZ"
	// The argument specifies whether there is a nodeset
	// at the location. 'x' is minX, 'X' is maxX, etc.
	size_t length = option[1].size();
	for (size_t j=0; j < length; j++) {
	  switch (option[1][j]) {
	  case 'x':
	    add_nodeset(MX);
	    break;
	  case 'X':
	    add_nodeset(PX);
	    break;
	  case 'y':
	    add_nodeset(MY);
	    break;
	  case 'Y':
	    add_nodeset(PY);
	    break;
	  case 'z':
	    add_nodeset(MZ);
	    break;
	  case 'Z':
	    add_nodeset(PZ);
	    break;
	  default:
	    std::cerr << "ERROR: Unrecognized nodeset location option '"
		      << option[1][j]
		      << "'.";
	  }
	}
      }
      else if (option[0] == "sideset") {
	// Option of the form  "sideset:xXyYzZ"
	// The argument specifies whether there is a sideset
	// at the location. 'x' is minX, 'X' is maxX, etc.
	size_t length = option[1].size();
	for (size_t j=0; j < length; j++) {
	  switch (option[1][j]) {
	  case 'x':
	    add_sideset(MX);
	    break;
	  case 'X':
	    add_sideset(PX);
	    break;
	  case 'y':
	    add_sideset(MY);
	    break;
	  case 'Y':
	    add_sideset(PY);
	    break;
	  case 'z':
	    add_sideset(MZ);
	    break;
	  case 'Z':
	    add_sideset(PZ);
	    break;
	  default:
	    std::cerr << "ERROR: Unrecognized sideset location option '"
		      << option[1][j]
		      << "'.";
	  }
	}
      }
      else if (option[0] == "scale") {
	// Option of the form  "scale:xs,ys,zs
	std::vector<std::string> tokens;
	Ioss::tokenize(option[1], ",", tokens);
	assert(tokens.size() == 3);
	sclX = std::strtod(tokens[0].c_str(), NULL);
	sclY = std::strtod(tokens[1].c_str(), NULL);
	sclZ = std::strtod(tokens[2].c_str(), NULL);
      }

      else if (option[0] == "offset") {
	// Option of the form  "offset:xo,yo,zo
	std::vector<std::string> tokens;
	Ioss::tokenize(option[1], ",", tokens);
	assert(tokens.size() == 3);
	offX = std::strtod(tokens[0].c_str(), NULL);
	offY = std::strtod(tokens[1].c_str(), NULL);
	offZ = std::strtod(tokens[2].c_str(), NULL);
      }

      else if (option[0] == "zdecomp") {
	// Option of the form  "zdecomp:1,1,2,2,1,2,...
	// Specifies the number of intervals in the z direction
	// for each processor.  The number of tokens must match
	// the number of processors.  Note that the new numZ will
	// be the sum of the intervals specified in this command.
	std::vector<std::string> tokens;
	Ioss::tokenize(option[1], ",", tokens);
	assert(tokens.size() == processorCount);
	Int64Vector Zs;
	numZ = 0;
	for (size_t j = 0; j < processorCount; j++) {
	  Zs.push_back(std::strtol(tokens[j].c_str(), NULL, 10));
	  numZ += Zs[j];
	}
	myNumZ = Zs[myProcessor];
	myStartZ = 0;
	for (size_t j=0; j < myProcessor; j++) {
	  myStartZ += Zs[j];
	}
      }

      else if (option[0] == "bbox") {
	// Bounding-Box Option of the form  "bbox:xmin,ymin,zmin,xmax,ymax,zmaxo
	std::vector<std::string> tokens;
	Ioss::tokenize(option[1], ",", tokens);
	assert(tokens.size() == 6);
	double xmin = std::strtod(tokens[0].c_str(), NULL);
	double ymin = std::strtod(tokens[1].c_str(), NULL);
	double zmin = std::strtod(tokens[2].c_str(), NULL);
	double xmax = std::strtod(tokens[3].c_str(), NULL);
	double ymax = std::strtod(tokens[4].c_str(), NULL);
	double zmax = std::strtod(tokens[5].c_str(), NULL);

	set_bbox(xmin, ymin, zmin,  xmax, ymax, zmax);
      }

      else if (option[0] == "rotate") {
	// Rotate Option of the form  "rotate:axis,angle,axis,angle,...
	std::vector<std::string> tokens;
	Ioss::tokenize(option[1], ",", tokens);
	assert(tokens.size() %2 == 0);
	for (size_t ir=0; ir < tokens.size();) {
	  std::string axis = tokens[ir++];
	  double angle_degree = std::strtod(tokens[ir++].c_str(), NULL);
	  set_rotation(axis, angle_degree);
	}
      }

      else if (option[0] == "help") {
	std::cerr << "\nValid Options for GeneratedMesh parameter string:\n"
		  << "\tIxJxK -- specifies intervals; must be first option. Ex: 4x10x12\n"
		  << "\toffset:xoff, yoff, zoff\n"
		  << "\tscale: xscl, yscl, zscl\n"
		  << "\tzdecomp:n1,n2,n3,...,n#proc\n"
		  << "\tbbox: xmin, ymin, zmin, xmax, ymax, zmax\n"
		  << "\trotate: axis,angle,axis,angle,...\n"
		  << "\tshell:xXyYzZ (specifies which plane to apply shell)\n"
	          << "\tnodeset:xXyXzZ (specifies which plane to apply nodeset)\n"
	          << "\tsideset:xXyXzZ (specifies which plane to apply sideset)\n"
		  << "\tshow -- show mesh parameters\n"
		  << "\thelp -- show this list\n\n";
      }

      else if (option[0] == "show") {
	show_parameters();
      }

      else {
	std::cerr << "ERROR: Unrecognized option '" << option[0]
		  << "'.  It will be ignored.\n";
      }
    }
  }

  void GeneratedMesh::show_parameters() const
  {
    if (myProcessor == 0) {
      std::cerr << "\nMesh Parameters:\n"
		<< "\tIntervals: " << numX << " by " << numY << " by " << numZ << "\n"
		<< "\tX = " << sclX << " * (0.." << numX << ") + " << offX
		<< "\tRange: " << offX << " <= X <= "  << offX + numX * sclX << "\n"
		<< "\tY = " << sclY << " * (0.." << numY << ") + " << offY
		<< "\tRange: " << offY << " <= Y <= "  << offY + numY * sclY << "\n"
		<< "\tZ = " << sclZ << " * (0.." << numZ << ") + " << offZ
		<< "\tRange: " << offZ << " <= Z <= "  << offZ + numZ * sclZ << "\n\n"
		<< "\tNode Count (total)    = " << std::setw(9) << node_count() << "\n"
		<< "\tElement Count (total) = " << std::setw(9) << element_count() << "\n"
		<< "\tBlock Count           = " << std::setw(9) << block_count() << "\n"
		<< "\tNodeset Count         = " << std::setw(9) << nodeset_count() << "\n"
		<< "\tSideset Count         = " << std::setw(9) << sideset_count() << "\n\n";
      if (doRotation) {
	std::cerr << "\tRotation Matrix: \n\t" << std::scientific ;
	for (int ii=0; ii < 3; ii++) {
	  for (int jj=0; jj < 3; jj++) {
	    std::cerr << std::setw(14) << rotmat[ii][jj] << "\t";
	  }
	  std::cerr << "\n\t";
	}
	std::cerr << std::fixed << "\n";
      }
    }
  }

  int64_t GeneratedMesh::node_count() const
  {
    return (numX+1) * (numY+1) * (numZ+1);
  }

  int64_t GeneratedMesh::node_count_proc() const
  {
    return (numX+1) * (numY+1) * (myNumZ+1);
  }

  int64_t GeneratedMesh::block_count() const
  {
    return shellBlocks.size() + 1;
  }

  int64_t GeneratedMesh::nodeset_count() const
  {
    return nodesets.size();
  }

  int64_t GeneratedMesh::sideset_count() const
  {
    return sidesets.size();
  }

  int64_t GeneratedMesh::element_count() const
  {
    int64_t count = element_count(1);
    for (size_t i=0; i < shellBlocks.size(); i++) {
      count += element_count(i+2);
    }
    return count;
  }

  int64_t GeneratedMesh::element_count_proc() const
  {
    int64_t count = 0;
    for (int64_t i=0; i < block_count(); i++) {
      count += element_count_proc(i+1);
    }
    return count;
  }

  int64_t GeneratedMesh::element_count(int64_t block_number) const
  {
    assert(block_number <= block_count());

    if (block_number == 1) {
      return numX * numY * numZ;
    } else {
      ShellLocation loc = shellBlocks[block_number-2];
      return shell_element_count(loc);
    }
  }

  int64_t GeneratedMesh::shell_element_count(ShellLocation loc) const
  {
    switch (loc) {
    case MX:
    case PX:
      return numY * numZ;
    case MY:
    case PY:
      return numX * numZ;
    case MZ:
    case PZ:
      return numX * numY;
    }
    return 0;
  }

  int64_t GeneratedMesh::element_count_proc(int64_t block_number) const
  {
    assert(block_number <= block_count());

    if (block_number == 1) {
      return numX * numY * myNumZ;
    } else {
      ShellLocation loc = shellBlocks[block_number-2];
      return shell_element_count_proc(loc);
    }
  }

  int64_t GeneratedMesh::shell_element_count_proc(ShellLocation loc) const
  {
    switch (loc) {
    case MX:
    case PX:
      return numY * myNumZ;
    case MY:
    case PY:
      return numX * myNumZ;
    case MZ:
      if (myProcessor == 0)
	return numX * numY;
      else
	return 0;
    case PZ:
      if (myProcessor == processorCount -1)
	return numX * numY;
      else
	return 0;
    }
    return 0;
  }

  int64_t GeneratedMesh::nodeset_node_count(int64_t id) const
  {
    // id is position in nodeset list + 1
    assert(id > 0 && (size_t)id <= nodesets.size());
    ShellLocation loc = nodesets[id-1];
    switch (loc) {
    case MX:
    case PX:
      return (numY+1) * (numZ+1);
    case MY:
    case PY:
      return (numX+1) * (numZ+1);
    case MZ:
    case PZ:
      return (numX+1) * (numY+1);
    }
    return 0;
  }

  int64_t GeneratedMesh::nodeset_node_count_proc(int64_t id) const
  {
    // id is position in nodeset list + 1
    assert(id > 0 && (size_t)id <= nodesets.size());
    ShellLocation loc = nodesets[id-1];
    switch (loc) {
    case MX:
    case PX:
      return (numY+1) * (myNumZ+1);
    case MY:
    case PY:
      return (numX+1) * (myNumZ+1);
    case MZ:
      if (myProcessor == 0)
	return (numX+1) * (numY+1);
      else
	return 0;
    case PZ:
      if (myProcessor == processorCount -1)
	return (numX+1) * (numY+1);
      else
	return 0;
    }
    return 0;
  }

  int64_t GeneratedMesh::sideset_side_count(int64_t id) const
  {
    // id is position in sideset list + 1
    assert(id > 0 && (size_t)id <= sidesets.size());
    ShellLocation loc = sidesets[id-1];
    switch (loc) {
    case MX:
    case PX:
      return numY * numZ;
    case MY:
    case PY:
      return numX * numZ;
    case MZ:
    case PZ:
      return numX * numY;
    }
    return 0;
  }

  int64_t GeneratedMesh::sideset_side_count_proc(int64_t id) const
  {
    // id is position in sideset list + 1
    assert(id > 0 && (size_t)id <= sidesets.size());
    ShellLocation loc = sidesets[id-1];
    switch (loc) {
    case MX:
    case PX:
      return numY * myNumZ;
    case MY:
    case PY:
      return numX * myNumZ;
    case MZ:
      if (myProcessor == 0)
	return numX * numY;
      else
	return 0;
    case PZ:
      if (myProcessor == processorCount -1)
	return numX * numY;
      else
	return 0;
    }
    return 0;
  }

  std::pair<std::string,int> GeneratedMesh::topology_type(int64_t block_number) const
  {
    assert(block_number <= block_count() && block_number > 0);

    if (block_number == 1) {
      return std::make_pair(std::string("hex8"), 8);
    } else {
      return std::make_pair(std::string("shell4"), 4);
    }
  }

  void GeneratedMesh::node_map(Int64Vector &map)
  {
    int64_t count = node_count_proc();
    map.reserve(count);
    int64_t offset = myStartZ * (numX+1) * (numY+1);
    for (int64_t i=0; i < count; i++) {
      map.push_back(offset + i + 1);
    }
  }

  void GeneratedMesh::node_map(IntVector &map)
  {
    int count = node_count_proc();
    map.resize(count);
    int offset = myStartZ * (numX+1) * (numY+1);
    for (int i=0; i < count; i++) {
      map[i] = offset + i + 1;
    }
  }

  int64_t GeneratedMesh::communication_node_count_proc() const
  {
    int64_t count = (numX+1) * (numY+1);
    if (myProcessor != 0 && myProcessor != processorCount-1)
      count *= 2;
    
    return count;
  }

  void GeneratedMesh::node_communication_map(Int64Vector &map, std::vector<int> &proc)
  {
    int64_t count = (numX+1) * (numY+1);
    int64_t slab = count;
    if (myProcessor != 0 && myProcessor != processorCount-1)
      count *= 2;

    map.resize(count);
    proc.resize(count);
    int64_t j = 0;
    if (myProcessor != 0) {
      int64_t offset = myStartZ * (numX+1) * (numY+1);
      for (int64_t i=0; i < slab; i++) {
	map[j] = offset + i + 1;
	proc[j++] = static_cast<int>(myProcessor-1);
      }
    }
    if (myProcessor != processorCount-1) {
      int64_t offset = (myStartZ + myNumZ) * (numX+1) * (numY+1);
      for (int64_t i=0; i < slab; i++) {
	map[j] = offset + i + 1;
	proc[j++] = static_cast<int>(myProcessor+1);
      }
    }
  }

  void GeneratedMesh::element_map(int64_t block_number, Int64Vector &map) const
  {
    assert(block_number <= block_count() && block_number > 0);

    int64_t count = element_count_proc(block_number);
    map.reserve(count);

    if (block_number == 1) {
      // Hex block...
      count = element_count_proc(1);
      int64_t offset = myStartZ * numX * numY;
      for (int64_t i=0; i < count; i++) {
	map.push_back(offset + i + 1);
      }
    } else {
      int64_t start = element_count(1);

      // Shell blocks...
      for (int64_t ib=0; (size_t)ib < shellBlocks.size(); ib++) {
	count = element_count_proc(ib+2);
	if (block_number == ib + 2) {
	  int64_t offset = 0;
	  ShellLocation loc = shellBlocks[ib];

	  switch (loc) {
	  case MX:
	  case PX:
	    offset = myStartZ * numY;
	    break;

	  case MY:
	  case PY:
	    offset = myStartZ * numX;
	    break;

	  case MZ:
	  case PZ:
	    offset = 0;
	    break;
	  }
	  for (int64_t i=0; i < count; i++) {
	    map.push_back(start + offset + i + 1);
	  }
	} else {
	  start += element_count(ib+2);
	}
      }
    }
  }

  void GeneratedMesh::element_map(int block_number, IntVector &map) const
  {
    assert(block_number <= block_count() && block_number > 0);

    int count = element_count_proc(block_number);
    map.resize(count);

    if (block_number == 1) {
      // Hex block...
      count = element_count_proc(1);
      int offset = myStartZ * numX * numY;
      for (int i=0; i < count; i++) {
	map[i] = offset + i + 1;
      }
    } else {
      int start = element_count(1);

      // Shell blocks...
      for (int ib=0; (size_t)ib < shellBlocks.size(); ib++) {
	count = element_count_proc(ib+2);
	if (block_number == ib + 2) {
	  int offset = 0;
	  ShellLocation loc = shellBlocks[ib];

	  switch (loc) {
	  case MX:
	  case PX:
	    offset = myStartZ * numY;
	    break;

	  case MY:
	  case PY:
	    offset = myStartZ * numX;
	    break;

	  case MZ:
	  case PZ:
	    offset = 0;
	    break;
	  }
	  for (int i=0; i < count; i++) {
	    map[i] = start + offset + i + 1;
	  }
	} else {
	  start += element_count(ib+2);
	}
      }
    }
  }

  void GeneratedMesh::element_map(Int64Vector &map) const
  {
    int64_t count = element_count_proc();
    map.reserve(count);

    // Hex block...
    count = element_count_proc(1);
    int64_t offset = myStartZ * numX * numY;
    for (int64_t i=0; i < count; i++) {
      map.push_back(offset + i + 1);
    }

    int64_t start = element_count(1);

    // Shell blocks...
    for (int64_t ib=0; (size_t)ib < shellBlocks.size(); ib++) {
      count = element_count_proc(ib+2);
      offset = 0;
      ShellLocation loc = shellBlocks[ib];

      switch (loc) {
      case MX:
      case PX:
	offset = myStartZ * numY;
	break;

      case MY:
      case PY:
	offset = myStartZ * numX;
	break;

      case MZ:
      case PZ:
	offset = 0;
	break;
      }
      for (int64_t i=0; i < count; i++) {
	map.push_back(start + offset + i + 1);
      }
      start += element_count(ib+2);
    }
  }

  void GeneratedMesh::element_map(IntVector &map) const
  {
    int count = element_count_proc();
    map.resize(count);

    int k = 0;
    // Hex block...
    count = element_count_proc(1);
    int offset = myStartZ * numX * numY;
    for (int i=0; i < count; i++) {
      map[k++] = offset + i + 1;
    }

    int start = element_count(1);

    // Shell blocks...
    for (int ib=0; (size_t)ib < shellBlocks.size(); ib++) {
      count = element_count_proc(ib+2);
      offset = 0;
      ShellLocation loc = shellBlocks[ib];

      switch (loc) {
      case MX:
      case PX:
	offset = myStartZ * numY;
	break;

      case MY:
      case PY:
	offset = myStartZ * numX;
	break;

      case MZ:
      case PZ:
	offset = 0;
	break;
      }
      for (int i=0; i < count; i++) {
	map[k++] = start + offset + i + 1;
      }
      start += element_count(ib+2);
    }
  }

  void GeneratedMesh::element_surface_map(ShellLocation loc, Int64Vector &map) const
  {
    int64_t count = shell_element_count_proc(loc);
    map.resize(2*count);

    int64_t index  = 0;
    int64_t offset = 0;

    switch (loc) {
    case MX:
      offset = myStartZ * numX * numY + 1;  // 1-based elem id
      for (size_t k = 0; k < myNumZ; ++k) {
	for (size_t j = 0; j < numY; ++j) {
	  map[index++] = offset; 
	  map[index++] = 3; // 0-based local face id
	  offset += numX;
	}
      }
      break;

    case PX:
      offset = myStartZ * numX * numY + numX;
      for (size_t k = 0; k < myNumZ; ++k) {
	for (size_t j = 0; j < numY; ++j) {
	  map[index++] = offset; // 1-based elem id
	  map[index++] = 1; // 0-based local face id
	  offset += numX;
	}
      }
      break;
      
    case MY:
      offset = myStartZ * numX * numY + 1;
      for (size_t k = 0; k < myNumZ; ++k) {
	for (size_t i = 0; i < numX; ++i) {
	  map[index++] = offset++;
	  map[index++] = 0; // 0-based local face id
	}
	offset+= numX * (numY-1);
      }
      break;

    case PY:
      offset = myStartZ * numX * numY + numX * (numY-1) +1;
      for (size_t k = 0; k < myNumZ; ++k) {
	for (size_t i = 0; i < numX; ++i) {
	  map[index++] = offset++;
	  map[index++] = 2; // 0-based local face id
	}
	offset+= numX * (numY-1);
      }
      break;

    case MZ:
      if (myProcessor == 0) {
	offset = 1;
	for (size_t i=0; i < numY; i++) {
	  for (size_t j=0; j < numX; j++) {
	    map[index++] = offset++;
	    map[index++] = 4;
	  }
	}
      }
      break;

    case PZ:
      if (myProcessor == processorCount-1) {
	offset = (numZ-1)*numX*numY + 1;
	for (size_t i=0, k=0; i < numY; i++) {
	  for (size_t j=0; j < numX; j++, k++) {
	    map[index++] = offset++;
	    map[index++] = 5;
	  }
	}
      }
      break;
    }
  }

  void GeneratedMesh::coordinates(std::vector<double> &coord) const
  {
    /* create global coordinates */
    int64_t count = node_count_proc();
    coord.resize(count * 3);
    coordinates(&coord[0]);
  }

  void GeneratedMesh::coordinates(double *coord) const
  {
    /* create global coordinates */
    int64_t count = node_count_proc();

    int64_t k = 0;
    for (size_t m=myStartZ; m < myStartZ+myNumZ+1; m++) {
      for (size_t i=0; i < numY+1; i++) {
	for (size_t j=0; j < numX+1; j++) {
	  coord[k++] = sclX * static_cast<double>(j) + offX;
	  coord[k++] = sclY * static_cast<double>(i) + offY;
	  coord[k++] = sclZ * static_cast<double>(m) + offZ;
	}
      }
    }

    if (doRotation) {
      for (int64_t i=0; i < count*3; i+=3) {
	double xn = coord[i+0];
	double yn = coord[i+1];
	double zn = coord[i+2];
	coord[i+0] = xn * rotmat[0][0] + yn * rotmat[1][0] + zn * rotmat[2][0];
	coord[i+1] = xn * rotmat[0][1] + yn * rotmat[1][1] + zn * rotmat[2][1];
	coord[i+2] = xn * rotmat[0][2] + yn * rotmat[1][2] + zn * rotmat[2][2];
      }
    }
  }

  void GeneratedMesh::coordinates(std::vector<double> &x,
				  std::vector<double> &y,
				  std::vector<double> &z) const
  {
    /* create global coordinates */
    int64_t count = node_count_proc();
    x.reserve(count);
    y.reserve(count);
    z.reserve(count);

    int64_t k = 0;
    for (size_t m=myStartZ; m < myStartZ+myNumZ+1; m++) {
      for (size_t i=0; i < numY+1; i++) {
	for (size_t j=0; j < numX+1; j++) {
	  x.push_back(sclX * static_cast<double>(j) + offX);
	  y.push_back(sclY * static_cast<double>(i) + offY);
	  z.push_back(sclZ * static_cast<double>(m) + offZ);
	  ++k;
	}
      }
    }
    if (doRotation) {
      for (int64_t i=0; i < count; i++) {
	double xn = x[i];
	double yn = y[i];
	double zn = z[i];
	x.push_back(xn * rotmat[0][0] + yn * rotmat[1][0] + zn * rotmat[2][0]);
	y.push_back(xn * rotmat[0][1] + yn * rotmat[1][1] + zn * rotmat[2][1]);
	z.push_back(xn * rotmat[0][2] + yn * rotmat[1][2] + zn * rotmat[2][2]);
      }
    }
  }

  void GeneratedMesh::coordinates(int component,
				  std::vector<double> &xyz) const
  {
    assert(!doRotation);
    /* create global coordinates */
    size_t count = node_count_proc();
    xyz.reserve(count);

    if (component == 1) {
      for (size_t m=myStartZ; m < myStartZ+myNumZ+1; m++) {
	for (size_t i=0; i < numY+1; i++) {
	  for (size_t j=0; j < numX+1; j++) {
	    xyz.push_back(sclX * static_cast<double>(j) + offX);
	  }
	}
      }
    }
    else if (component == 2) {
      for (size_t m=myStartZ; m < myStartZ+myNumZ+1; m++) {
	for (size_t i=0; i < numY+1; i++) {
	  for (size_t j=0; j < numX+1; j++) {
	    xyz.push_back(sclY * static_cast<double>(i) + offY);
	  }
	}
      }
    }
    else if (component == 3) {
      for (size_t m=myStartZ; m < myStartZ+myNumZ+1; m++) {
	for (size_t i=0; i < numY+1; i++) {
	  for (size_t j=0; j < numX+1; j++) {
	    xyz.push_back(sclZ * static_cast<double>(m) + offZ);
	  }
	}
      }
    }    
  }

  void GeneratedMesh::connectivity(int64_t block_number, Int64Vector &connect) const
  {
    if (block_number == 1) {  // HEX Element Block
      connect.resize(element_count_proc(block_number)*8);
    } else {
      connect.resize(element_count_proc(block_number)*4);
    }
    connectivity(block_number, &connect[0]);
  }

  void GeneratedMesh::connectivity(int64_t block_number, int64_t *connect) const
  {
    assert(block_number <= block_count());

    int64_t xp1yp1 = (numX+1) * (numY+1);

    /* build connectivity array (node list) for mesh */
    if (block_number == 1) {  // HEX Element Block

      size_t cnt = 0;
      for (size_t m=myStartZ; m < myNumZ+myStartZ; m++) {
	for (size_t i=0, k=0; i < numY; i++) {
	  for (size_t j=0; j < numX; j++, k++) {
	    size_t base = (m*xp1yp1) + k + i + 1;
	    ;
	    connect[cnt++] = base;
	    connect[cnt++] = base+1;
	    connect[cnt++] = base+numX+2;
	    connect[cnt++] = base+numX+1;

	    connect[cnt++] = xp1yp1 + base;
	    connect[cnt++] = xp1yp1 + base+1;
	    connect[cnt++] = xp1yp1 + base+numX+2;
	    connect[cnt++] = xp1yp1 + base+numX+1;
	  }
	}
      }
    } else { // Shell blocks....
      ShellLocation loc = shellBlocks[block_number-2];

      size_t cnt = 0;
      switch (loc) {
      case MX:  // Minumum X Face
	for (size_t i=0; i < myNumZ; i++) {
	  size_t layer_off = i * xp1yp1;
	  for (size_t j=0; j < numY; j++) {
	    size_t base = layer_off + j * (numX+1) + 1 + myStartZ * xp1yp1;
	    connect[cnt++] = base;
	    connect[cnt++] = base + xp1yp1;
	    connect[cnt++] = base + xp1yp1 + (numX+1);
	    connect[cnt++] = base + (numX+1);
	  }
	}
	break;
      case PX: // Maximum X Face
	for (size_t i=0; i < myNumZ; i++) {
	  size_t layer_off = i * xp1yp1;
	  for (size_t j=0; j < numY; j++) {
	    size_t base = layer_off + j * (numX+1) + numX + 1 + myStartZ * xp1yp1;
	    connect[cnt++] = base;
	    connect[cnt++] = base + (numX+1);
	    connect[cnt++] = base + xp1yp1 +(numX+1);
	    connect[cnt++] = base + xp1yp1;
	  }
	}
	break;
      case MY: // Minumum Y Face
	for (size_t i=0; i < myNumZ; i++) {
	  size_t layer_off = i * xp1yp1;
	  for (size_t j=0; j < numX; j++) {
	    size_t base = layer_off + j + 1 + myStartZ * xp1yp1;
	    connect[cnt++] = base;
	    connect[cnt++] = base + 1;
	    connect[cnt++] = base + xp1yp1 + 1;
	    connect[cnt++] = base + xp1yp1;
	  }
	}
	break;
      case PY: // Maximum Y Face
	for (size_t i=0; i < myNumZ; i++) {
	  size_t layer_off = i * xp1yp1;
	  for (size_t j=0; j < numX; j++) {
	    size_t base = layer_off + (numX+1)*(numY) + j + 1 + myStartZ * xp1yp1;
	    connect[cnt++] = base;
	    connect[cnt++] = base + xp1yp1;
	    connect[cnt++] = base + xp1yp1 + 1;
	    connect[cnt++] = base + 1;
	  }
	}
	break;
      case MZ: // Minumum Z Face
	if (myProcessor == 0) {
	  for (size_t i=0, k=0; i < numY; i++) {
	    for (size_t j=0; j < numX; j++, k++) {
	      size_t base = i + k + 1 + myStartZ * xp1yp1;
	      connect[cnt++] = base;
	      connect[cnt++] = base+numX+1;
	      connect[cnt++] = base+numX+2;
	      connect[cnt++] = base+1;
	    }
	  }
	}
	break;
      case PZ: // Maximum Z Face
	if (myProcessor == processorCount-1) {
	  for (size_t i=0, k=0; i < numY; i++) {
	    for (size_t j=0; j < numX; j++, k++) {
	      size_t base = xp1yp1 * (numZ - myStartZ) + k + i + 1 + myStartZ * xp1yp1;
	      connect[cnt++] = base;
	      connect[cnt++] = base+1;
	      connect[cnt++] = base+numX+2;
	      connect[cnt++] = base+numX+1;
	    }
	  }
	}
	break;
      }
      assert(cnt == size_t(4 * element_count_proc(block_number)));
    }
    return;
  }

  void GeneratedMesh::connectivity(int64_t block_number, IntVector &connect) const
  {
    if (block_number == 1) {  // HEX Element Block
      connect.resize(element_count_proc(block_number)*8);
    } else {
      connect.resize(element_count_proc(block_number)*4);
    }
    connectivity(block_number, &connect[0]);
  }

  void GeneratedMesh::connectivity(int64_t block_number, int *connect) const
  {
    assert(block_number <= block_count());

    int xp1yp1 = (numX+1) * (numY+1);

    /* build connectivity array (node list) for mesh */
    if (block_number == 1) {  // HEX Element Block

      int cnt = 0;
      for (size_t m=myStartZ; m < myNumZ+myStartZ; m++) {
	for (size_t i=0, k=0; i < numY; i++) {
	  for (size_t j=0; j < numX; j++, k++) {
	    size_t base = (m*xp1yp1) + k + i + 1;
	    ;
	    connect[cnt++] = base;
	    connect[cnt++] = base+1;
	    connect[cnt++] = base+numX+2;
	    connect[cnt++] = base+numX+1;

	    connect[cnt++] = xp1yp1 + base;
	    connect[cnt++] = xp1yp1 + base+1;
	    connect[cnt++] = xp1yp1 + base+numX+2;
	    connect[cnt++] = xp1yp1 + base+numX+1;
	  }
	}
      }
    } else { // Shell blocks....
      ShellLocation loc = shellBlocks[block_number-2];

      size_t cnt = 0;
      switch (loc) {
      case MX:  // Minumum X Face
	for (size_t i=0; i < myNumZ; i++) {
	  size_t layer_off = i * xp1yp1;
	  for (size_t j=0; j < numY; j++) {
	    size_t base = layer_off + j * (numX+1) + 1 + myStartZ * xp1yp1;
	    connect[cnt++] = base;
	    connect[cnt++] = base + xp1yp1;
	    connect[cnt++] = base + xp1yp1 + (numX+1);
	    connect[cnt++] = base + (numX+1);
	  }
	}
	break;
      case PX: // Maximum X Face
	for (size_t i=0; i < myNumZ; i++) {
	  size_t layer_off = i * xp1yp1;
	  for (size_t j=0; j < numY; j++) {
	    size_t base = layer_off + j * (numX+1) + numX + 1 + myStartZ * xp1yp1;
	    connect[cnt++] = base;
	    connect[cnt++] = base + (numX+1);
	    connect[cnt++] = base + xp1yp1 +(numX+1);
	    connect[cnt++] = base + xp1yp1;
	  }
	}
	break;
      case MY: // Minumum Y Face
	for (size_t i=0; i < myNumZ; i++) {
	  size_t layer_off = i * xp1yp1;
	  for (size_t j=0; j < numX; j++) {
	    size_t base = layer_off + j + 1 + myStartZ * xp1yp1;
	    connect[cnt++] = base;
	    connect[cnt++] = base + 1;
	    connect[cnt++] = base + xp1yp1 + 1;
	    connect[cnt++] = base + xp1yp1;
	  }
	}
	break;
      case PY: // Maximum Y Face
	for (size_t i=0; i < myNumZ; i++) {
	  size_t layer_off = i * xp1yp1;
	  for (size_t j=0; j < numX; j++) {
	    size_t base = layer_off + (numX+1)*(numY) + j + 1 + myStartZ * xp1yp1;
	    connect[cnt++] = base;
	    connect[cnt++] = base + xp1yp1;
	    connect[cnt++] = base + xp1yp1 + 1;
	    connect[cnt++] = base + 1;
	  }
	}
	break;
      case MZ: // Minumum Z Face
	if (myProcessor == 0) {
	  for (size_t i=0, k=0; i < numY; i++) {
	    for (size_t j=0; j < numX; j++, k++) {
	      size_t base = i + k + 1 + myStartZ * xp1yp1;
	      connect[cnt++] = base;
	      connect[cnt++] = base+numX+1;
	      connect[cnt++] = base+numX+2;
	      connect[cnt++] = base+1;
	    }
	  }
	}
	break;
      case PZ: // Maximum Z Face
	if (myProcessor == processorCount-1) {
	  for (size_t i=0, k=0; i < numY; i++) {
	    for (size_t j=0; j < numX; j++, k++) {
	      size_t base = xp1yp1 * (numZ - myStartZ) + k + i + 1 + myStartZ * xp1yp1;
	      connect[cnt++] = base;
	      connect[cnt++] = base+1;
	      connect[cnt++] = base+numX+2;
	      connect[cnt++] = base+numX+1;
	    }
	  }
	}
	break;
      }
      assert(cnt == size_t(4 * element_count_proc(block_number)));
    }
    return;
  }

  void GeneratedMesh::nodeset_nodes(int64_t id, Int64Vector &nodes) const
  {
    // id is position in nodeset list + 1
    assert(id > 0 && (size_t)id <= nodesets.size());
    ShellLocation loc = nodesets[id-1];
    nodes.resize(nodeset_node_count_proc(id));

    size_t xp1yp1 = (numX+1) * (numY+1);
    size_t k = 0;

    switch (loc) {
    case MX:  // Minumum X Face
      for (size_t i=0; i < myNumZ+1; i++) {
	size_t layer_off = myStartZ * xp1yp1 + i * xp1yp1;
	for (size_t j=0; j < numY+1; j++) {
	  nodes[k++] = layer_off + j * (numX+1) + 1;
	}
      }
      break;
    case PX: // Maximum X Face
      for (size_t i=0; i < myNumZ+1; i++) {
	size_t layer_off = myStartZ * xp1yp1 + i * xp1yp1;
	for (size_t j=0; j < numY+1; j++) {
	  nodes[k++] = layer_off + j * (numX+1) + numX + 1;
	}
      }
      break;
    case MY: // Minumum Y Face
      for (size_t i=0; i < myNumZ+1; i++) {
	size_t layer_off = myStartZ * xp1yp1 + i * xp1yp1;
	for (size_t j=0; j < numX+1; j++) {
	  nodes[k++] = layer_off + j + 1;
	}
      }
      break;
    case PY: // Maximum Y Face
      for (size_t i=0; i < myNumZ+1; i++) {
	size_t layer_off = myStartZ * xp1yp1 + i * xp1yp1;
	for (size_t j=0; j < numX+1; j++) {
	  nodes[k++] = layer_off + (numX+1)*(numY) + j + 1;
	}
      }
      break;
    case MZ: // Minumum Z Face
      if (myProcessor == 0) {
	for (size_t i=0; i < (numY+1) * (numX+1); i++) {
	  nodes[i] = i+1;
	}
      }
      break;
    case PZ: // Maximum Z Face
      if (myProcessor == processorCount-1) {
	size_t offset = (numY+1) * (numX+1) * numZ;
	for (size_t i=0; i < (numY+1) * (numX+1); i++) {
	  nodes[i] = offset + i+1;
	}
      }
      break;
    }
  }

  void GeneratedMesh::sideset_elem_sides(int64_t id, Int64Vector &elem_sides) const
  {
    // id is position in sideset list + 1
    assert(id > 0 && (size_t)id <= sidesets.size());
    ShellLocation loc = sidesets[id-1];
    
    // If there is a shell block on this face, then the sideset is
    // applied to the shell block; if not, it is applied to the
    // underlying hex elements.
    bool underlying_shell = false;
    int64_t shell_block = 0;
    for (size_t i = 0; i < shellBlocks.size(); i++) {
      if (shellBlocks[i] == loc) {
	underlying_shell = true;
	shell_block = i+2;
	break;
      }
    }

    if (underlying_shell) {
      // Get ids of shell elements at this location...
      element_map(shell_block, elem_sides);
      
      // Insert face_ordinal in between each entry in elem_sides...
      // Face will be 0 for all shells...
      elem_sides.resize(2*sideset_side_count_proc(id));
      ssize_t face_ordinal = 0;
      ssize_t i = 2* sideset_side_count_proc(id) - 1;
      ssize_t j =    sideset_side_count_proc(id) - 1;
      while (i >= 0) {
	elem_sides[i--] = face_ordinal;
	elem_sides[i--] = elem_sides[j--];
      }
    } else {
      element_surface_map(loc, elem_sides);
    }
  }

  void GeneratedMesh::set_rotation(const std::string &axis, double angle_degrees)
  {
    // PI / 180. Used in converting angle in degrees to radians
    static double degang = std::atan2(0.0, -1.0) / 180.0;

    doRotation = true;

    int n1 = -1;
    int n2 = -1;
    int n3 = -1;

    if (axis == "x" || axis == "X") {
      n1 = 1; n2 = 2; n3 = 0;
    } else if (axis == "y" || axis == "Y") {
      n1 = 2; n2 = 0; n3 = 1;
    } else if (axis == "z" || axis == "Z") {
      n1 = 0; n2 = 1; n3 = 2;
    } else {
      std::cerr << "\nInvalid axis specification '" << axis << "'. Valid options are 'x', 'y', or 'z'\n";
      return;
    }

    double ang = angle_degrees * degang; // Convert angle in degrees to radians
    double cosang = std::cos(ang);
    double sinang = std::sin(ang);

    assert(n1 >= 0 && n2 >= 0 && n3 >= 0);
    double by[3][3];
    by[n1][n1] =  cosang;
    by[n2][n1] = -sinang;
    by[n1][n3] = 0.0;
    by[n1][n2] =  sinang;
    by[n2][n2] =  cosang;
    by[n2][n3] = 0.0;
    by[n3][n1] = 0.0;
    by[n3][n2] = 0.0;
    by[n3][n3] = 1.0;

    double res[3][3];
    for (int i=0; i < 3; i++) {
      res[i][0] = rotmat[i][0]*by[0][0] + rotmat[i][1]*by[1][0] + rotmat[i][2]*by[2][0];
      res[i][1] = rotmat[i][0]*by[0][1] + rotmat[i][1]*by[1][1] + rotmat[i][2]*by[2][1];
      res[i][2] = rotmat[i][0]*by[0][2] + rotmat[i][1]*by[1][2] + rotmat[i][2]*by[2][2];
    }

#if 1
    std::memcpy(rotmat, res, 9*sizeof(double));
#else
    for (int i=0; i < 3; i++) {
      for (int j=0; j < 3; j++) {
	rotmat[i][j] = res[i][j];
      }
    }
#endif
  }

  int64_t DashSurfaceMesh::node_count() const
{
    return mCoordinates.size()/SPATIAL_DIMENSION;
}

int64_t DashSurfaceMesh::node_count_proc() const
{
    return node_count();
}

int64_t DashSurfaceMesh::element_count() const
{
    return (mQuadSurface1.size()+mQuadSurface2.size())/NUM_NODES_PER_QUAD_FACE;
}

int64_t DashSurfaceMesh::element_count(int64_t block_number) const
{
    if(block_number == 1)
    {
        return mQuadSurface1.size()/NUM_NODES_PER_QUAD_FACE;
    }
    else if(block_number == 2)
    {
        return mQuadSurface2.size()/NUM_NODES_PER_QUAD_FACE;
    }
    throw std::exception();

    return INVALID;
}

int64_t DashSurfaceMesh::block_count() const
{
    return 2;
}

int64_t DashSurfaceMesh::nodeset_count() const
{
    return 0;
}

int64_t DashSurfaceMesh::sideset_count() const
{
    return 2;
}

int64_t DashSurfaceMesh::element_count_proc() const
{
    return element_count();
}

int64_t DashSurfaceMesh::element_count_proc(int64_t block_number) const
{
    return element_count(block_number);
}

int64_t DashSurfaceMesh::nodeset_node_count_proc(int64_t id) const
{
    return 0;
}
int64_t DashSurfaceMesh::sideset_side_count_proc(int64_t id) const
{
    return element_count(id);
}
int64_t DashSurfaceMesh::communication_node_count_proc() const
{
    return 0;
}

void DashSurfaceMesh::coordinates(double *coord) const
{
    std::copy(mCoordinates.begin(),mCoordinates.end(), coord);
}

void DashSurfaceMesh::coordinates(std::vector<double> &coord) const
{
    throw std::exception();
}

void DashSurfaceMesh::coordinates(int component, std::vector<double> &xyz) const
{
    throw std::exception();
}

void DashSurfaceMesh::coordinates(std::vector<double> &x, std::vector<double> &y, std::vector<double> &z) const
{
    throw std::exception();
}

void DashSurfaceMesh::connectivity(int64_t block_number, int* connect) const
{
    switch(block_number)
    {
        case 1:
            std::copy(mQuadSurface1.begin(),mQuadSurface1.end(), connect);
            return;
        case 2:
            std::copy(mQuadSurface2.begin(),mQuadSurface2.end(), connect);
            return;
        default:
            throw std::exception();
    }
}

std::pair<std::string, int> DashSurfaceMesh::topology_type(int64_t block_number) const
{
    const int numNodesPerElement = 4;
    return std::make_pair(std::string("shell4"), numNodesPerElement);
}

void DashSurfaceMesh::sideset_elem_sides(int64_t setId, Int64Vector &elem_sides) const
{
    elem_sides.clear();
    size_t numElementsInSurface1 = mQuadSurface1.size()/NUM_NODES_PER_QUAD_FACE;
    switch(setId)
    {
        case 1:
            for(size_t i=0; i<numElementsInSurface1; ++i)
            {
                elem_sides.push_back(i+1);
                elem_sides.push_back(0);
            }
            return;
        case 2:
            for(size_t i=0; i<mQuadSurface2.size()/NUM_NODES_PER_QUAD_FACE; ++i)
            {
                elem_sides.push_back(numElementsInSurface1+i+1);
                elem_sides.push_back(0);
            }
            return;
        default:
            throw std::exception();
    }
}

void DashSurfaceMesh::nodeset_nodes(int64_t nset_id, Int64Vector &nodes) const
{
    return;
}

void DashSurfaceMesh::node_communication_map(MapVector &map, std::vector<int> &proc)
{
    return;
}

void DashSurfaceMesh::node_map(IntVector &map)
{
    map.resize(node_count());
    for(int i = 0; i < node_count(); i++)
    {
        map[i] = i + 1;
    }
}

void DashSurfaceMesh::node_map(MapVector &map)
{
    map.resize(node_count());
    for(int i = 0; i < node_count(); i++)
    {
        map[i] = i + 1;
    }
}

void DashSurfaceMesh::element_map(int block_number, IntVector &map) const
{
    size_t numElementsInSurface1 = mQuadSurface1.size() / NUM_NODES_PER_QUAD_FACE;
    size_t numElementsInSurface2 = mQuadSurface2.size() / NUM_NODES_PER_QUAD_FACE;
    switch(block_number)
    {
        case 1:
            for(size_t i = 0; i < numElementsInSurface1; ++i)
            {
                map[i] = i + 1;
            }
            return;
        case 2:
            for(size_t i = 0; i < numElementsInSurface2; ++i)
            {
                map[numElementsInSurface1 + i] = numElementsInSurface1 + i + 1;
            }
            return;
        default:
            throw std::exception();
    }
}

void DashSurfaceMesh::element_map(int64_t block_number, MapVector &map) const
{
    size_t numElementsInSurface1 = mQuadSurface1.size() / NUM_NODES_PER_QUAD_FACE;
    size_t numElementsInSurface2 = mQuadSurface2.size() / NUM_NODES_PER_QUAD_FACE;
    switch(block_number)
    {
        case 1:
            for(size_t i = 0; i < numElementsInSurface1; ++i)
            {
                map[i] = i + 1;
            }
            return;
        case 2:
            for(size_t i = 0; i < numElementsInSurface2; ++i)
            {
                map[numElementsInSurface1 + i] = numElementsInSurface1 + i + 1;
            }
            return;
        default:
            throw std::exception();
    }
}

void DashSurfaceMesh::element_map(MapVector &map) const
{
    int count = element_count_proc();
    map.resize(count);

    for(int i = 0; i < count; i++)
    {
        map[i] = i + 1;
    }
}

void DashSurfaceMesh::element_map(IntVector &map) const
{
    int count = element_count_proc();
    map.resize(count);

    for(int i = 0; i < count; i++)
    {
        map[i] = i + 1;
    }
}

}

