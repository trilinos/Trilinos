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

#include <heartbeat/Iohb_DatabaseIO.h>
#include <heartbeat/Iohb_Layout.h>

#include <Ioss_CodeTypes.h>
#include <Ioss_SubSystem.h>
#include <Ioss_Utils.h>
#include <string>

#include <fstream>
#include <time.h>

namespace {
  std::string time_stamp(const std::string &format)
  {
    if (format == "") {
      return std::string("");
    } else {
      const int length=256;
      static char time_string[length];

      time_t calendar_time = time(NULL);
      struct tm *local_time = localtime(&calendar_time);

      size_t error = strftime(time_string, length, format.c_str(), local_time);
      if (error != 0) {
	time_string[length-1] = '\0';
	return std::string(time_string);
      } else {
	return std::string("[ERROR]");
      }
    }
  }

  std::ostream *open_stream(const std::string &filename, bool *needs_delete) {
    // A little wierdness and ambiguity is possible here.  We want to
    // minimize the number of commands, but maximize the
    // functionality. For example, we want to be able to specify output
    // to existing streams (cout/stdout, cerr/stderr, outputP0), files,
    // sockets, XML-RPC?, SOAP?.  However, we want to be able to specify
    // this with a single command.
    // So..., we first check for some 'reserved' stream names.  These
    // are the 'cout, stdout, cerr, stderr, output, outputP0' type.
    // Note that this means the user can't specify a file name of
    // 'cerr', but that shouldn't be too much of a hardship.  [If it is,
    // we can devise a syntax as: 'FILE:cout' or do a './cout'

    std::ostream *log_stream = NULL;
    *needs_delete = false;
    if (filename == "cout" || filename == "stdout") {
      log_stream = &std::cout;
    } else if (filename == "cerr" || filename == "stderr") {
      log_stream = &std::cerr;  // This is also default if nothing specified
    } else if (filename == "output" || filename == "outputP0") {
      log_stream = &std::cout;  // This should be sierra::Env::outputP0(), but not
                                // available during transition to stk...
    } else if (filename == "clog" || filename == "log") {
      log_stream = &std::clog; // Same as cerr, but not flushed automatically.
    } else {
      // Open the file (on processor 0 only) Might need to do
      // something better here if we want to share streams among
      // different heartbeats or logging mechanisms.  Need perhaps a
      // 'logger' class which handles sharing and destruction...
      log_stream = new std::ofstream(filename.c_str());
      *needs_delete = true;
    }
    return log_stream;
  }
}

namespace Iohb {

  // ========================================================================
  const IOFactory* IOFactory::factory()
  {
    static IOFactory registerThis;
    return &registerThis;
  }

  IOFactory::IOFactory()
    : Ioss::IOFactory("heartbeat")
  {}

  Ioss::DatabaseIO* IOFactory::make_IO(const std::string& filename,
				       Ioss::DatabaseUsage db_usage,
				       MPI_Comm communicator) const
  { return new DatabaseIO(NULL, filename, db_usage, communicator); }

  // ========================================================================
  DatabaseIO::DatabaseIO(Ioss::Region *region, const std::string& filename,
			 Ioss::DatabaseUsage db_usage,
			 MPI_Comm communicator) :
    Ioss::DatabaseIO(region, filename, db_usage, communicator),
    logStream(NULL), layout_(NULL), legend_(NULL),
    tsFormat("[%H:%M:%S]"), precision_(5), showLabels(true), showLegend(false), initialized_(false),
    streamNeedsDelete(false), fileFormat(DEFAULT)
  { }

  DatabaseIO::~DatabaseIO()
  {
    delete layout_;
    delete legend_;
    if (streamNeedsDelete && logStream) {
      delete logStream;
    }
  }

  void DatabaseIO::initialize(const Ioss::Region *region) const
  {
    if (!initialized_) {
      assert(layout_ == NULL);
      assert(legend_ == NULL);

      DatabaseIO *new_this = const_cast<DatabaseIO*>(this);

      if (region->property_exists("file_format")) {
	std::string format = region->get_property("file_format").get_string();
	if (format == "spyhis")
	  new_this->fileFormat = SPYHIS;
      }

      if (util().parallel_rank() == 0) {
	new_this->logStream = open_stream(get_filename().c_str(),
					  &(new_this->streamNeedsDelete));
      } else {
	// All processors except processor 0
	new_this->logStream = NULL;
      }

      // Pull variables from the regions property data...
      if (region->property_exists("time_stamp_format")) {
	new_this->tsFormat = region->get_property("time_stamp_format").get_string();
      }

      if (region->property_exists("precision")) {
	new_this->precision_ = region->get_property("precision").get_int();
      }

      if (region->property_exists("show_labels")) {
	new_this->showLabels = (region->get_property("show_labels").get_int() == 1);
      }

      if (region->property_exists("show_legend")) {
	new_this->showLegend = (region->get_property("show_legend").get_int() == 1);
      }

      if (fileFormat == SPYHIS) {
	new_this->showLegend = true;
	new_this->showLabels = false;
	new_this->tsFormat = "";
      }
      
      if (showLegend) {
	new_this->legend_ = new Layout(false, precision_);
	if (tsFormat != "") {
	  new_this->legend_->add_literal("+");
	  new_this->legend_->add_literal(time_stamp(tsFormat));
	  new_this->legend_->add_literal(" ");
	}
	if (!fileFormat == SPYHIS)
	  new_this->legend_->add_literal("Legend: ");

	new_this->legend_->add_literal("TIME, ");
      }
      new_this->initialized_ = true;
    }
  }

  bool DatabaseIO::begin(Ioss::State /* state */)
  {
    return true;
  }

  bool   DatabaseIO::end(Ioss::State /* state */)
  {
    return true;
  }

  bool DatabaseIO::begin_state(Ioss::Region *region, int /* state */, double time )
  {
    // If this is the first time, open the output stream and see if user wants a legend
    initialize(region);

    layout_ = new Layout(showLabels, precision_);
    if (tsFormat != "") {
      layout_->add_literal("+");
      layout_->add_literal(time_stamp(tsFormat));
      layout_->add_literal(" ");
    }

    if (fileFormat == SPYHIS) {
      layout_->add("TIME", time);
    }

    return true;
  }

  bool   DatabaseIO::end_state(Ioss::Region */* region */, int /* state */, double /* time */)
  {
    if (legend_ != NULL) {
      if (fileFormat == SPYHIS) {
	time_t calendar_time = time(NULL);
	*logStream << "% Sierra SPYHIS Output " << ctime(&calendar_time);
	// This should be the comma-separated long names for the variables and the units.
	// For example Density (g/cc)
	// I think it will work with for now to label the plots with the "short name"
	*logStream << *legend_ << std::endl;
      }

      *logStream << *legend_ << std::endl;
      delete legend_;
      legend_ = NULL;
    }

    *logStream << *layout_ << std::endl;
    delete layout_;
    layout_ = NULL;
    return true;
  }

  int DatabaseIO::get_field_internal(const Ioss::Region* /* reg */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }

  int DatabaseIO::get_field_internal(const Ioss::NodeBlock* /* nb */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int DatabaseIO::get_field_internal(const Ioss::EdgeBlock* /* nb */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int DatabaseIO::get_field_internal(const Ioss::FaceBlock* /* nb */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int DatabaseIO::get_field_internal(const Ioss::ElementBlock* /* eb */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }

  int DatabaseIO::get_field_internal(const Ioss::NodeSet* /* ns */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int DatabaseIO::get_field_internal(const Ioss::EdgeSet* /* ns */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int DatabaseIO::get_field_internal(const Ioss::FaceSet* /* ns */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int DatabaseIO::get_field_internal(const Ioss::ElementSet* /* ns */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int DatabaseIO::get_field_internal(const Ioss::SideBlock* /* eb */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int DatabaseIO::get_field_internal(const Ioss::SideSet* /* fs */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int DatabaseIO::get_field_internal(const Ioss::CommSet* /* cs */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }

  int DatabaseIO::put_field_internal(const Ioss::Region* region, const Ioss::Field& field,
				     void *data, size_t data_size) const
  {
    initialize(region);
    Ioss::Field::RoleType role = field.get_role();
    int num_to_get = field.verify(data_size);

    if ((role == Ioss::Field::TRANSIENT || role == Ioss::Field::REDUCTION) &&
	num_to_get == 1) {

      int ncomp = field.transformed_storage()->component_count();

      if (legend_ != NULL && layout_ != NULL) {
	legend_->add(field.get_name(), field.get_name());
	if (ncomp > 1) {
	  legend_->add_literal("(");
	  legend_->add_literal(Ioss::Utils::to_string(ncomp));
	  legend_->add_literal(")");
	}
      }

      if (field.get_type() == Ioss::Field::STRING) {
	// Assume that if layout_ is NULL, then we want special one-line output.
	if (layout_ == NULL) {
	  Layout layout(false, 0);
	  layout.add_literal("-");
	  layout.add_literal(time_stamp(tsFormat));
	  layout.add_literal(" ");
	  layout.add_literal(*(std::string*)data);
	  if (logStream != NULL)
	    *logStream << layout << std::endl;
	} else {
	  layout_->add(field.get_name(), *(std::string*)data);
	}
      } else {
	if (layout_ == NULL) {
	  std::ostringstream errmsg;
	  errmsg << "INTERNAL ERROR: Unexpected NULL layout.\n";
	  IOSS_ERROR(errmsg);
	}
	if (field.get_type() == Ioss::Field::INTEGER) {
	  assert(field.transformed_count() == 1);

	  int *i_data = (int*)data;
	  std::vector<int> idata(ncomp);
	  for (int i=0; i < ncomp; i++) {
	    idata[i] = i_data[i];
	  }
	  layout_->add(field.get_name(), idata);
	} else {
	  std::vector<double> rdata(ncomp);
	  double *r_data = (double*)data;
	  for (int i=0; i < ncomp; i++) {
	    rdata[i] = r_data[i];
	  }
	  layout_->add(field.get_name(), rdata);
	}
      }
    } else {
      std::ostringstream errmsg;
      errmsg << "FATAL: Can not handle non-TRANSIENT or non-REDUCTION fields on regions.\n";
      IOSS_ERROR(errmsg);
    }
    return num_to_get;
  }

  int DatabaseIO::put_field_internal(const Ioss::ElementBlock* /* eb */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int DatabaseIO::put_field_internal(const Ioss::FaceBlock* /* nb */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int DatabaseIO::put_field_internal(const Ioss::EdgeBlock* /* nb */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int DatabaseIO::put_field_internal(const Ioss::NodeBlock* /* nb */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }

  int DatabaseIO::put_field_internal(const Ioss::NodeSet* /* ns */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int DatabaseIO::put_field_internal(const Ioss::EdgeSet* /* ns */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int DatabaseIO::put_field_internal(const Ioss::FaceSet* /* ns */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int DatabaseIO::put_field_internal(const Ioss::ElementSet* /* ns */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int DatabaseIO::put_field_internal(const Ioss::SideBlock* /* fb */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int DatabaseIO::put_field_internal(const Ioss::SideSet* /* fs */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }
  int DatabaseIO::put_field_internal(const Ioss::CommSet* /* cs */, const Ioss::Field& /* field */,
				     void */* data */, size_t /* data_size */) const
  {
    return -1;
  }

  unsigned DatabaseIO::entity_field_support() const
  {
    return Ioss::REGION;
  }
}

