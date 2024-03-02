// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_MeshSurface_h
#define Akri_MeshSurface_h

#include <Akri_Faceted_Surface.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_util/diag/Timer.hpp>

#include <stk_math/StkVector.hpp>
#include <memory>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Entity; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Selector; } }

namespace krino {

template <class FACET>
class MeshSurface : public Faceted_Surface<FACET> {
public:
  MeshSurface(const stk::mesh::MetaData & meta,
              const stk::mesh::Field<double>& coord_ref,
              const stk::mesh::Selector & surface_selector,
              const int sign);

  virtual ~MeshSurface() {}
  virtual void build_local_facets(const BoundingBox & proc_bbox) override;

protected:
  void add_facet2d(const stk::mesh::BulkData& mesh, stk::mesh::Entity side, unsigned p0, unsigned p1);
  void add_facet3d(const stk::mesh::BulkData& mesh, stk::mesh::Entity side, unsigned p0, unsigned p1, unsigned p2);
  void add_quad3d(const stk::mesh::BulkData& mesh, stk::mesh::Entity side, unsigned p0, unsigned p1, unsigned p2, unsigned p3);
  void add_facet2d(stk::math::Vector3d & p0, stk::math::Vector3d & p1);
  void add_facet3d(stk::math::Vector3d & p0, stk::math::Vector3d & p1, stk::math::Vector3d & p2);
  
private:
  int my_sign;
  const stk::mesh::MetaData & my_mesh_meta;
  const stk::mesh::Field<double>& my_coord_ref;
  const stk::mesh::Selector my_surface_selector;
};

std::unique_ptr<FacetedSurfaceBase> build_mesh_surface(const stk::mesh::MetaData & meta,
    const stk::mesh::Field<double>& coordsField,
    const stk::mesh::Selector & surfaceSelector,
    const int sign);

class Parallel_Facet_File_Reader {
public:
  Parallel_Facet_File_Reader(const std::string & in_filename) : my_filename(in_filename) {}

  static void get_batch_size(const int local_num_facets, int & batch_size, int & num_batches);

  std::ifstream & input() { return my_input; }
  const std::string & filename() const { return my_filename; }

  void open_file();
  void close_file();
  void read(const std::string & read_description, const std::function<void(void)> & read_function);

private:
  std::string my_filename;
  std::ifstream my_input;
};

class Faceted_Surface_From_File : public Faceted_Surface<Facet3d> {
public:
  Faceted_Surface_From_File(const std::string & surface_name, const stk::diag::Timer &parent_timer);
  virtual ~Faceted_Surface_From_File() {}
  virtual void build_local_facets(const BoundingBox & proc_bbox) override;
  virtual BoundingBox get_bounding_box() = 0;
  void insert_into(BoundingBox & bbox) const override;

protected:
  virtual void read_file(const std::vector<BoundingBox> & proc_bboxes) = 0;
private:
  stk::diag::Timer my_timer;
  bool my_built_local_facets;
};

class STLSurface : public Faceted_Surface_From_File {
public:
  STLSurface(const std::string & surface_name,
      const stk::diag::Timer &parent_timer,
      const std::string & filename,
      const int sign,
      const stk::math::Vector3d & scale);

  virtual ~STLSurface() {}
  virtual BoundingBox get_bounding_box() override;

private:
  virtual void read_file(const std::vector<BoundingBox> & proc_bboxes) override;
  void read_header();
  unsigned read_num_binary_facets();
  unsigned read_ascii_facets(const unsigned max_batch_size);
  unsigned read_binary_facets(const unsigned maxBatchSize, unsigned & numRemainingInFile);
  BoundingBox get_ascii_facet_bounding_box();
  BoundingBox get_binary_facet_bounding_box(const unsigned numFacets);

private:
  bool my_is_ascii;
  Parallel_Facet_File_Reader my_reader;
  int my_dist_sign;
  stk::math::Vector3d my_scale;
};

class FACSurface : public Faceted_Surface_From_File {
public:
  FACSurface(const std::string & surface_name,
      const stk::diag::Timer &parent_timer,
      const std::string & filename,
      const int sign,
      const stk::math::Vector3d & scale);

  virtual ~FACSurface() {}
  virtual BoundingBox get_bounding_box() override;

private:
  virtual void read_file(const std::vector<BoundingBox> & proc_bboxes) override;
  void read_points(std::vector<stk::math::Vector3d> & points);
  void read_facets(const int batch_size, const int num_facets, const std::vector<stk::math::Vector3d> & points);

private:
  Parallel_Facet_File_Reader my_reader;
  int my_dist_sign;
  stk::math::Vector3d my_scale;
};

class PLYSurface : public Faceted_Surface_From_File {
public:
  PLYSurface(const std::string & surface_name,
      const stk::diag::Timer &parent_timer,
      const std::string & filename,
      const int sign,
      const stk::math::Vector3d & scale);

  virtual ~PLYSurface() {}
  virtual BoundingBox get_bounding_box() override;

private:
  Parallel_Facet_File_Reader my_reader;
  virtual void read_file(const std::vector<BoundingBox> & proc_bboxes) override;
  void read_header(int & num_points, int & num_facets);
  void read_points(const int num_points, std::vector<stk::math::Vector3d> & points);
  void read_facets(const int batch_size, const std::vector<stk::math::Vector3d> & points);

private:
  int my_dist_sign;
  stk::math::Vector3d my_scale;
};

class EXOSurface : public Faceted_Surface_From_File {
public:
  EXOSurface(const std::string & surface_name,
      const stk::diag::Timer &parent_timer,
      const std::string & filename,
      const int sign,
      const stk::math::Vector3d & scale);

  virtual ~EXOSurface() {}
  virtual BoundingBox get_bounding_box() override;

private:
  virtual void read_file(const std::vector<BoundingBox> & proc_bboxes) override;

private:
  std::string my_filename;
  int my_dist_sign;
  stk::math::Vector3d my_scale;
};

} // namespace krino

#endif // Akri_MeshSurface_h
