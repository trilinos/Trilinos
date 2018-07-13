/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/unit_test_support/GeneratedMesh.hpp>
#include <stk_util/util/tokenize.hpp>

#include <cmath>
#include <cstring>

#include <vector>
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <assert.h>

namespace stk_classic {
  namespace io {
    namespace util {
      GeneratedMesh::GeneratedMesh(int num_x, int num_y, int num_z,
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
        std::vector<std::string> groups;
        stk_classic::util::tokenize(parameters, "|+", groups);

        // First 'group' is the interval specification -- IxJxK
        std::vector<std::string> tokens;
        stk_classic::util::tokenize(groups[0], "x", tokens);
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
        assert(numZ >= processorCount);
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

      size_t GeneratedMesh::add_shell_block(ShellLocation loc)
      {
        shellBlocks.push_back(loc);
        return shellBlocks.size();
      }

      size_t GeneratedMesh::add_nodeset(ShellLocation loc)
      {
	nodesets.push_back(loc);
	return nodesets.size();
      }
      
      size_t GeneratedMesh::add_sideset(ShellLocation loc)
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
          stk_classic::util::tokenize(groups[i], ":", option);
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
            stk_classic::util::tokenize(option[1], ",", tokens);
            assert(tokens.size() == 3);
            sclX = std::strtod(tokens[0].c_str(), NULL);
            sclY = std::strtod(tokens[1].c_str(), NULL);
            sclZ = std::strtod(tokens[2].c_str(), NULL);
          }

          else if (option[0] == "offset") {
            // Option of the form  "offset:xo,yo,zo
            std::vector<std::string> tokens;
            stk_classic::util::tokenize(option[1], ",", tokens);
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
            stk_classic::util::tokenize(option[1], ",", tokens);
            assert(tokens.size() == processorCount);
            std::vector<size_t> Zs;
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
            stk_classic::util::tokenize(option[1], ",", tokens);
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
            stk_classic::util::tokenize(option[1], ",", tokens);
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

      size_t GeneratedMesh::node_count() const
      {
        return (numX+1) * (numY+1) * (numZ+1);
      }

      size_t GeneratedMesh::node_count_proc() const
      {
        return (numX+1) * (numY+1) * (myNumZ+1);
      }

      size_t GeneratedMesh::block_count() const
      {
        return shellBlocks.size() + 1;
      }

      size_t GeneratedMesh::nodeset_count() const
      {
	return nodesets.size();
      }
      
      size_t GeneratedMesh::sideset_count() const
      {
	return sidesets.size();
      }

      size_t GeneratedMesh::element_count() const
      {
        size_t count = element_count(1);
        for (size_t i=0; i < shellBlocks.size(); i++) {
          count += element_count(i+2);
        }
        return count;
      }

      size_t GeneratedMesh::element_count_proc() const
      {
        size_t count = 0;
        for (size_t i=0; i < block_count(); i++) {
          count += element_count_proc(i+1);
        }
        return count;
      }

      size_t GeneratedMesh::element_count(size_t block_number) const
      {
        assert(block_number <= block_count());

        if (block_number == 1) {
          return numX * numY * numZ;
        } else {
          ShellLocation loc = shellBlocks[block_number-2];
          return shell_element_count(loc);
        }
      }

      size_t GeneratedMesh::shell_element_count(ShellLocation loc) const
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

      size_t GeneratedMesh::element_count_proc(size_t block_number) const
      {
        assert(block_number <= block_count());

        if (block_number == 1) {
          return numX * numY * myNumZ;
        } else {
          ShellLocation loc = shellBlocks[block_number-2];
          return shell_element_count_proc(loc);
        }
      }

      size_t GeneratedMesh::shell_element_count_proc(ShellLocation loc) const
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

      size_t GeneratedMesh::nodeset_node_count(size_t id) const
      {
	// id is position in nodeset list + 1
	assert(id > 0 && id <= nodesets.size());
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

      size_t GeneratedMesh::nodeset_node_count_proc(size_t id) const
      {
	// id is position in nodeset list + 1
	assert(id > 0 && id <= nodesets.size());
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

      size_t GeneratedMesh::sideset_side_count(size_t id) const
      {
	// id is position in sideset list + 1
	assert(id > 0 && id <= sidesets.size());
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

      size_t GeneratedMesh::sideset_side_count_proc(size_t id) const
      {
	// id is position in sideset list + 1
	assert(id > 0 && id <= sidesets.size());
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

      std::pair<std::string,int> GeneratedMesh::topology_type(size_t block_number) const
      {
        assert(block_number <= block_count() && block_number > 0);

        if (block_number == 1) {
          return std::make_pair(std::string("hex8"), 8);
        } else {
          return std::make_pair(std::string("shell4"), 4);
        }
      }

      void GeneratedMesh::node_map(std::vector<int> &map)
      {
        size_t count = node_count_proc();
        map.resize(count);
        size_t offset = myStartZ * (numX+1) * (numY+1);
        for (size_t i=0; i < count; i++) {
          map[i] = static_cast<int>(offset + i + 1);
        }
      }

      size_t GeneratedMesh::communication_node_count_proc() const
      {
	size_t count = (numX+1) * (numY+1);
	if (myProcessor != 0 && myProcessor != processorCount-1)
	  count *= 2;
    
	return count;
      }

      void GeneratedMesh::node_communication_map(std::vector<int> &map, std::vector<int> &proc)
      {
        size_t count = (numX+1) * (numY+1);
        size_t slab = count;
        if (myProcessor != 0 && myProcessor != processorCount-1)
          count *= 2;

        map.resize(count);
        proc.resize(count);
        size_t j = 0;
        if (myProcessor != 0) {
          size_t offset = myStartZ * (numX+1) * (numY+1);
          for (size_t i=0; i < slab; i++) {
            map[j] = static_cast<int>(offset + i + 1);
            proc[j++] = static_cast<int>(myProcessor-1);
          }
        }
        if (myProcessor != processorCount-1) {
          size_t offset = (myStartZ + myNumZ) * (numX+1) * (numY+1);
          for (size_t i=0; i < slab; i++) {
            map[j] = static_cast<int>(offset + i + 1);
            proc[j++] = static_cast<int>(myProcessor+1);
          }
        }
      }

      void GeneratedMesh::element_map(size_t block_number, std::vector<int> &map) const
      {
        assert(block_number <= block_count() && block_number > 0);

        size_t count = element_count_proc(block_number);
        map.resize(count);

        if (block_number == 1) {
          // Hex block...
          count = element_count_proc(1);
          size_t offset = myStartZ * numX * numY;
          for (size_t i=0; i < count; i++) {
            map[i] = static_cast<int>(offset + i + 1);
          }
        } else {
          size_t start = element_count(1);

          // Shell blocks...
          for (size_t ib=0; ib < shellBlocks.size(); ib++) {
            count = element_count_proc(ib+2);
            if (static_cast<size_t>(block_number) == ib + 2) {
              size_t offset = 0;
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
              for (size_t i=0; i < count; i++) {
                map[i] = static_cast<int>(start + offset + i + 1);
              }
            } else {
              start += element_count(ib+2);
            }
          }
        }
      }

      void GeneratedMesh::element_map(std::vector<int> &map) const
      {
        size_t count = element_count_proc();
        map.resize(count);

        size_t k = 0;
        // Hex block...
        count = element_count_proc(1);
        size_t offset = myStartZ * numX * numY;
        for (size_t i=0; i < count; i++) {
          map[k++] = static_cast<int>(offset + i + 1);
        }

        size_t start = element_count(1);

        // Shell blocks...
        for (size_t ib=0; ib < shellBlocks.size(); ib++) {
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
          for (size_t i=0; i < count; i++) {
            map[k++] = static_cast<int>(start + offset + i + 1);
          }
          start += element_count(ib+2);
        }
      }

      void GeneratedMesh::element_surface_map(ShellLocation loc, std::vector<int> &map) const
      {
	size_t count = shell_element_count_proc(loc);
	map.resize(2*count);

	size_t index  = 0;
	size_t offset = 0;

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
        size_t count = node_count_proc();
        coord.resize(count * 3);

        size_t k = 0;
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
          for (size_t i=0; i < count*3; i+=3) {
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
        size_t count = node_count_proc();
        x.resize(count);
        y.resize(count);
        z.resize(count);

        size_t k = 0;
        for (size_t m=myStartZ; m < myStartZ+myNumZ+1; m++) {
          for (size_t i=0; i < numY+1; i++) {
            for (size_t j=0; j < numX+1; j++) {
              x[k] = sclX * static_cast<double>(j) + offX;
              y[k] = sclY * static_cast<double>(i) + offY;
              z[k] = sclZ * static_cast<double>(m) + offZ;
              ++k;
            }
          }
        }
        if (doRotation) {
          for (size_t i=0; i < count; i++) {
            double xn = x[i];
            double yn = y[i];
            double zn = z[i];
            x[i] = xn * rotmat[0][0] + yn * rotmat[1][0] + zn * rotmat[2][0];
            y[i] = xn * rotmat[0][1] + yn * rotmat[1][1] + zn * rotmat[2][1];
            z[i] = xn * rotmat[0][2] + yn * rotmat[1][2] + zn * rotmat[2][2];
          }
        }
      }

      void GeneratedMesh::connectivity(size_t block_number, std::vector<int> &connect) const
      {
	assert(block_number <= block_count());

	size_t xp1yp1 = (numX+1) * (numY+1);

	/* build connectivity array (node list) for mesh */
	if (block_number == 1) {  // HEX Element Block
	  connect.resize(element_count_proc(block_number)*8);

	  size_t cnt = 0;
	  for (size_t m=myStartZ; m < myNumZ+myStartZ; m++) {
	    for (size_t i=0, k=0; i < numY; i++) {
	      for (size_t j=0; j < numX; j++, k++) {
		size_t base = (m*xp1yp1) + k + i + 1;
		;
		connect[cnt++] = static_cast<int>(base);
		connect[cnt++] = static_cast<int>(base+1);
		connect[cnt++] = static_cast<int>(base+numX+2);
		connect[cnt++] = static_cast<int>(base+numX+1);

		connect[cnt++] = static_cast<int>(xp1yp1 + base);
		connect[cnt++] = static_cast<int>(xp1yp1 + base+1);
		connect[cnt++] = static_cast<int>(xp1yp1 + base+numX+2);
		connect[cnt++] = static_cast<int>(xp1yp1 + base+numX+1);
	      }
	    }
	  }
	} else { // Shell blocks....
	  ShellLocation loc = shellBlocks[block_number-2];
	  connect.resize(element_count_proc(block_number)*4);

	  size_t cnt = 0;
	  switch (loc) {
	  case MX:  // Minumum X Face
	    for (size_t i=0; i < myNumZ; i++) {
	      size_t layer_off = i * xp1yp1;
	      for (size_t j=0; j < numY; j++) {
		size_t base = layer_off + j * (numX+1) + 1 + myStartZ * xp1yp1;
		connect[cnt++] = static_cast<int>(base);
		connect[cnt++] = static_cast<int>(base + xp1yp1);
		connect[cnt++] = static_cast<int>(base + xp1yp1 + (numX+1));
		connect[cnt++] = static_cast<int>(base + (numX+1));
	      }
	    }
	    break;
	  case PX: // Maximum X Face
	    for (size_t i=0; i < myNumZ; i++) {
	      size_t layer_off = i * xp1yp1;
	      for (size_t j=0; j < numY; j++) {
		size_t base = layer_off + j * (numX+1) + numX + 1 + myStartZ * xp1yp1;
		connect[cnt++] = static_cast<int>(base);
		connect[cnt++] = static_cast<int>(base + (numX+1));
		connect[cnt++] = static_cast<int>(base + xp1yp1 + (numX+1));
		connect[cnt++] = static_cast<int>(base + xp1yp1);
	      }
	    }
	    break;
	  case MY: // Minumum Y Face
	    for (size_t i=0; i < myNumZ; i++) {
	      size_t layer_off = i * xp1yp1;
	      for (size_t j=0; j < numX; j++) {
		size_t base = layer_off + j + 1 + myStartZ * xp1yp1;
		connect[cnt++] = static_cast<int>(base);
		connect[cnt++] = static_cast<int>(base + 1);
		connect[cnt++] = static_cast<int>(base + xp1yp1 + 1);
		connect[cnt++] = static_cast<int>(base + xp1yp1);
	      }
	    }
	    break;
	  case PY: // Maximum Y Face
	    for (size_t i=0; i < myNumZ; i++) {
	      size_t layer_off = i * xp1yp1;
	      for (size_t j=0; j < numX; j++) {
		size_t base = layer_off + (numX+1)*(numY) + j + 1 + myStartZ * xp1yp1;
		connect[cnt++] = static_cast<int>(base);
		connect[cnt++] = static_cast<int>(base + xp1yp1);
		connect[cnt++] = static_cast<int>(base + xp1yp1 + 1);
		connect[cnt++] = static_cast<int>(base + 1);
	      }
	    }
	    break;
	  case MZ: // Minumum Z Face
	    if (myProcessor == 0) {
	      for (size_t i=0, k=0; i < numY; i++) {
		for (size_t j=0; j < numX; j++, k++) {
		  size_t base = i + k + 1 + myStartZ * xp1yp1;
		  connect[cnt++] = static_cast<int>(base);
		  connect[cnt++] = static_cast<int>(base+numX+1);
		  connect[cnt++] = static_cast<int>(base+numX+2);
		  connect[cnt++] = static_cast<int>(base+1);
		}
	      }
	    }
	    break;
	  case PZ: // Maximum Z Face
	    if (myProcessor == processorCount-1) {
	      for (size_t i=0, k=0; i < numY; i++) {
		for (size_t j=0; j < numX; j++, k++) {
		  size_t base = xp1yp1 * (numZ - myStartZ) + k + i + 1 + myStartZ * xp1yp1;
		  connect[cnt++] = static_cast<int>(base);
		  connect[cnt++] = static_cast<int>(base+1);
		  connect[cnt++] = static_cast<int>(base+numX+2);
		  connect[cnt++] = static_cast<int>(base+numX+1);
		}
	      }
	    }
	    break;
	  }
	  assert(cnt == 4 * element_count_proc(block_number));
	}
	return;
      }

      void GeneratedMesh::nodeset_nodes(size_t id, std::vector<int> &nodes) const
      {
	// id is position in nodeset list + 1
	assert(id > 0 && id <= nodesets.size());
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

      void GeneratedMesh::sideset_elem_sides(size_t id, std::vector<int> &elem_sides) const
      {
	// id is position in sideset list + 1
	assert(id > 0 && id <= sidesets.size());
	ShellLocation loc = sidesets[id-1];
    
	// If there is a shell block on this face, then the sideset is
	// applied to the shell block; if not, it is applied to the
	// underlying hex elements.
	bool underlying_shell = false;
	size_t shell_block = 0;
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
	  int face_ordinal = 0;
	  int i = 2* (int)sideset_side_count_proc(id) - 1;
	  int j =    (int)sideset_side_count_proc(id) - 1;
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
        for (size_t i=0; i < 3; i++) {
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
    }
  }
}

#if defined(DEBUG)
#include <sstream>
  std::string ToString(size_t t) {
    std::ostringstream os;
    os << t;
    return os.str();
  }

#include <exodusII.h>

  int main() {
    int num_processors = 8;
    for (int proc = 0; proc < num_processors; proc++) {

      stk_classic::GeneratedMesh mesh(100, 125, 10*num_processors, num_processors, proc);

      std::cerr << "Node Count (total)    = " << mesh.node_count() << "\n";
      std::cerr << "Node Count (proc)     = " << mesh.node_count_proc() << "\n";
      std::cerr << "Element Count (total) = " << mesh.element_count() << "\n";
      std::cerr << "Element Count (proc)  = " << mesh.element_count_proc() << "\n";
      std::cerr << "Block Count           = " << mesh.block_count() << "\n";

      int CPU_word_size = 8;                   /* sizeof(float) */
      int IO_word_size = 8;                    /* (4 bytes) */
      std::string name = "test-scale.e";
      if (num_processors > 1) {
        name += "." + ToString(num_processors) + "." + ToString(proc);
      }
      int exoid = ex_create (name.c_str(), EX_CLOBBER, &CPU_word_size, &IO_word_size);

      int num_nodes = mesh.node_count_proc();
      int num_elems = mesh.element_count_proc();
      int num_elem_blk = mesh.block_count();
      int error = ex_put_init (exoid, "title", 3, num_nodes, num_elems,
                               num_elem_blk, 0, 0);

      if (num_processors > 1) {
        std::vector<int> nodes;
        std::vector<int> procs;
        mesh.node_communication_map(nodes, procs);

        int node_map_ids[1] = {1};
        int node_map_node_cnts[1] = {procs.size()};
        ex_put_init_info(exoid, num_processors, 1, "p");
        ex_put_loadbal_param(exoid, 0, 0, 0, 0, 0, 1, 0, proc);
        ex_put_cmap_params(exoid, node_map_ids, node_map_node_cnts, 0, 0, proc);
        ex_put_node_cmap(exoid, 1, &nodes[0], &procs[0], proc);
      }

      for (int i=1; i < mesh.block_count(); i++) {
        std::cerr << "Block " << i+1 << " has " << mesh.element_count_proc(i+1) << " "
                  << mesh.topology_type(i+1).first <<" elements\n";

        std::string btype = mesh.topology_type(i+1).first;
        error = ex_put_elem_block(exoid, i+1, btype.c_str(),
                                  mesh.element_count_proc(i+1),
                                  mesh.topology_type(i+1).second, 0);
      }
      {
        std::cerr << "Block " << 1 << " has " << mesh.element_count_proc(1) << " "
                  << mesh.topology_type(1).first <<" elements\n";

        std::string btype = mesh.topology_type(1).first;
        error = ex_put_elem_block(exoid, 1, btype.c_str(),
                                  mesh.element_count_proc(1),
                                  mesh.topology_type(1).second, 0);
      }

      if (num_processors > 1) {
        std::vector<int> map;
        mesh.node_map(map);
        ex_put_id_map(exoid, EX_NODE_MAP, &map[0]);

        mesh.element_map(map);
        ex_put_id_map(exoid, EX_ELEM_MAP, &map[0]);
      }

      std::cerr << "Outputting connectivity...\n";
      for (int i=1; i < mesh.block_count(); i++) {
        if (mesh.element_count_proc(i+1) > 0) {
          std::vector<int> connectivity;
          mesh.connectivity(i+1, connectivity);
          ex_put_elem_conn(exoid, i+1, &connectivity[0]);
        }
      }
      {
        std::vector<int> connectivity;
        mesh.connectivity(1, connectivity);
        ex_put_elem_conn(exoid, 1, &connectivity[0]);
      }

      {
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
        mesh.coordinates(x, y, z);
        error = ex_put_coord (exoid, &x[0], &y[0], &z[0]);
      }

      ex_close(exoid);
    }
  }
#endif
