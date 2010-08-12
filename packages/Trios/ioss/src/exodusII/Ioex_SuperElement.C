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

#include <exodusII/Ioex_SuperElement.h>

#include <Ioss_Property.h>
#include <Ioss_Field.h>
#include <Ioss_Utils.h>
#include <string>

#include <assert.h>

#include <netcdf.h>

namespace {
  int nc_get_array(int ncid, const char *name, double *data)
  {
    // NOTE: size of data array is validated in
    // calling code.  Just read and store here.

    int varid = 0;
    int status = nc_inq_varid(ncid, name, &varid);
    if (status != NC_NOERR) {
      return status;
    }

    status = nc_get_var_double(ncid, varid, data);
    return status;
  }

  int nc_get_dimension(int ncid, const char* DIMENSION,
		       const char *label, size_t *count)
  {
    std::ostringstream errmsg;

    *count = 0;
    int dimid = -1;

    int status = nc_inq_dimid(ncid, DIMENSION, &dimid);
    if (status != NC_NOERR) {
      if (status == NC_EBADDIM) {
	// Value is zero if the dimension is not defined.
	*count = 0;
	return 0;
      } else {
	errmsg << "ERROR: Failed to locate number of " << label
	       << " in superelement file.";
	IOSS_ERROR(errmsg);
      }
      return -1;
    }

    status = nc_inq_dimlen (ncid, dimid, count);
    if (status != NC_NOERR) {
      errmsg << "ERROR: failed to get number of " << label
	     << " in superelement file.";
      IOSS_ERROR(errmsg);
      return -1;
    }
    return status;
  }

  const std::string SCALAR() { return std::string("scalar");}
}

Ioex::SuperElement::SuperElement(const std::string &filename,
				 const std::string &my_name)
  : Ioss::GroupingEntity(NULL, my_name), fileName(filename),
    numDOF(0), numEIG(0), filePtr(-1)
{

  // For now, we will open the raw netcdf file here and parse the
  // dimensions. This is probably not how this should be done long
  // term, but is better than putting netcdf calls in application...

  // Check that file specified by filename exists...
  // Add working directory if needed.
  std::string local_filename = fileName;

  int status = nc_open(local_filename.c_str(), NC_NOWRITE, &filePtr);
  if (status != NC_NOERR) {
    std::ostringstream errmsg;
    errmsg << "ERROR: Failed to open superelement file '"
	   << local_filename << "'.";
    IOSS_ERROR(errmsg);
  }

  // At this point have a valid netcdf file handle.
  // Read some dimensions to determine size of Mass and Stiffness
  // matrix.
  status = nc_get_dimension(filePtr, "NumDof",
			    "number of degrees of freedom",
			    &numDOF);
  
  status = nc_get_dimension(filePtr, "NumEig",
			    "number of eigenvalues",
			    &numEIG);

  size_t num_constraints = 0;
  status = nc_get_dimension(filePtr, "NumConstraints",
			    "number of interface dof",
			    &num_constraints);
  assert(num_constraints == numDOF - numEIG);
    
  // NumCols and NumDof are redundant dimensions in the netcdf file.
  // Eventually, NumCols will be removed and NumDof is the long-term value.
  // Verify that they match.
  size_t numCols = 0;
  status = nc_get_dimension(filePtr, "NumCols",
			    "number of matrix columns",
			    &numCols);
  assert(numCols == numDOF);
  

  // Add the standard properties...
  properties.add(Ioss::Property(this, "numDOF",
				Ioss::Property::INTEGER));
  properties.add(Ioss::Property(this, "numEIG",
				Ioss::Property::INTEGER));
  properties.add(Ioss::Property(this, "numConstraints",
				Ioss::Property::INTEGER));

  // Add the standard fields...
  fields.add(Ioss::Field("Kr", Ioss::Field::REAL, SCALAR(),
			 Ioss::Field::MESH, numDOF*numDOF));

  fields.add(Ioss::Field("Mr", Ioss::Field::REAL, SCALAR(),
			 Ioss::Field::MESH, numDOF*numDOF));


  // There are additional properties and fields on the netcdf file,
  // but for now we only need "Kr" and "Mr"
}

int Ioex::SuperElement::internal_get_field_data(const Ioss::Field& field,
				      void *data, size_t data_size) const
{
  size_t num_to_get = field.verify(data_size);
  assert(num_to_get == numDOF * numDOF);
  
  if (field.get_name() == "Kr") {
    int status = nc_get_array(filePtr, "Kr", (double*)data);
    if (status != 0) {
      std::ostringstream errmsg;
      errmsg << "ERROR: Could not load stiffness matrix field 'Kr' from file '"
	     << fileName << "'.";
      IOSS_ERROR(errmsg);
    }
  }
  else if (field.get_name() == "Mr") {
    int status = nc_get_array(filePtr, "Mr", (double*)data);
    if (status != 0) {
      std::ostringstream errmsg;
      errmsg << "ERROR: Could not load mass matrix field 'Mr' from file '"
	     << fileName << "'.";
      IOSS_ERROR(errmsg);
    }
  }
  else {
    std::cerr << "WARNING: " << type() << " '" << name()
	      << "'. Unknown input field '" << field.get_name() << "'";
    return -4;
  }
  return num_to_get;
}

Ioex::SuperElement::~SuperElement()
{
  if (filePtr) nc_close(filePtr);
}

int Ioex::SuperElement::internal_put_field_data(const Ioss::Field& /* field */,
						void* /* data */, size_t /* data_size */) const
{
  return -1;
}

Ioss::Property Ioex::SuperElement::get_implicit_property(const std::string& the_name) const
{
  if (Ioss::Utils::case_strcmp(the_name, "numDOF") == 0) {
    return Ioss::Property(the_name, (int)numDOF);
  }
  else if (Ioss::Utils::case_strcmp(the_name, "numEIG") == 0) {
    return Ioss::Property(the_name, (int)numEIG);
  }
  else if (Ioss::Utils::case_strcmp(the_name, "numConstraints") == 0) {
    return Ioss::Property(the_name, (int)numDOF-(int)numEIG);
  }
  else {
    return Ioss::GroupingEntity::get_implicit_property(the_name);
  }
}
