// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_SPARSE_3_TENSOR_PARTITION_HPP
#define STOKHOS_SPARSE_3_TENSOR_PARTITION_HPP

#include "Stokhos_Sparse3Tensor.hpp"

#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Array.hpp"

namespace Stokhos {

  template <typename TupleType>
  class RCB {
  public:
    typedef typename TupleType::value_type coord_type;
    typedef typename Teuchos::ArrayView<TupleType>::size_type size_type;
    typedef typename TupleType::id_type id_type;

    struct CoordCompare {
      size_type d;
      CoordCompare(const size_type& d_) : d(d_) {}
      bool operator() (const TupleType& a, const TupleType& b) const {
        return a(d) < b(d);
      }
    };

    struct Box {
      coord_type xmin, ymin, zmin, xmax, ymax, zmax;
      coord_type delta_x, delta_y, delta_z;
      size_type split_dim;
      Teuchos::RCP<Box> left, right;
      Teuchos::Array<TupleType> coords;

      Box(const Teuchos::ArrayView<TupleType>& c) : coords(c.begin(), c.end()) {
        computeBoundingBox();
        computeSplittingDimension();
      }

      // Compute bounding box around points
      void computeBoundingBox(){
        xmin = coords[0](0); xmax = coords[0](0);
        ymin = coords[0](1); ymax = coords[0](1);
        zmin = coords[0](2); zmax = coords[0](2);
        for (size_type i=0; i<coords.size(); ++i) {
          coord_type x = coords[i](0);
          coord_type y = coords[i](1);
          coord_type z = coords[i](2);
          if (x < xmin) xmin = x;
          if (y < ymin) ymin = y;
          if (z < zmin) zmin = z;
          if (x > xmax) xmax = x;
          if (y > ymax) ymax = y;
          if (z > zmax) zmax = z;
        }
        delta_x = xmax - xmin + 1;
        delta_y = ymax - ymin + 1;
        delta_z = zmax - zmin + 1;

        // std::cout << "delta_x = " << delta_x
        //           << " delta_y = " << delta_y
        //           << " delta_z = " << delta_z << std::endl;
      }

      // Compute dimension to split
      void computeSplittingDimension() {
        split_dim = 0;
        if (delta_y >= delta_x && delta_y >= delta_z) split_dim = 1;
        if (delta_z >= delta_x && delta_z >= delta_y) split_dim = 2;

        //std::cout << "splitting dimension = " << split_dim << std::endl;
      }

      // Split box into two pieces with roughly equal numbers of points
      void split() {
        // Sort points based on splitting dimension
        CoordCompare cmp(split_dim);
        std::sort(coords.begin(), coords.end(), cmp);

        // Divide coords into two bins of roughly equal size, keeping
        // coords with equal values for split dimension together
        size_type n = coords.size();
        size_type s = n / 2;

        while (s < n-1 && coords[s](split_dim) == coords[s+1](split_dim)) ++s;
        //std::cout << "n = " << n << " s = " << s << std::endl;

        if (s > 0)
          left = Teuchos::rcp(new Box(coords.view(0, s)));
        if (s < n)
          right = Teuchos::rcp(new Box(coords.view(s, n-s)));

        // Clear my coordinate array since we aren't a leaf
        //Teuchos::Array<TupleType>().swap(coords);
        //coords.resize(0);
      }

    };

    //! Constructor
    RCB(const coord_type& max_length_,
        const size_type& max_parts_,
        const Teuchos::ArrayView<TupleType>& coords_) :
      max_length(max_length_),
      max_parts(max_parts_),
      coords(coords_.begin(), coords_.end()) {
      partition();
    }

    //! Destructor
    ~RCB() {}

    //! Get number of parts
    size_type get_num_parts() const { return num_parts; }

    //! Get root of partition
    Teuchos::RCP<Box> get_partition_root() const { return root; }

    //! Get parts array
    Teuchos::RCP< Teuchos::Array< Teuchos::RCP<Box> > > get_parts() const {
      return parts; }

    // Create list of part IDs for each tuple
    Teuchos::ArrayRCP<id_type> get_part_IDs() const {
      Teuchos::ArrayRCP<id_type> part_ids(coords.size());
      for (size_type part=0; part<num_parts; ++part) {
        Teuchos::RCP<Box> box = (*parts)[part];
        size_type n = box->coords.size();
        for (size_type i=0; i<n; ++i)
          part_ids[ box->coords[i].ID() ] = part;
      }
      return part_ids;
    }

  private:

    coord_type max_length;
    size_type max_parts, num_parts;
    Teuchos::Array<TupleType> coords;
    Teuchos::RCP<Box> root;
    Teuchos::RCP< Teuchos::Array< Teuchos::RCP<Box> > > parts;

    //! Partition
    void partition() {

      // Create root bounding box
      root = Teuchos::rcp(new Box(coords()));
      num_parts = 1;
      parts = Teuchos::rcp(new Teuchos::Array< Teuchos::RCP<Box> >);

      // Create array of boxes that are too big
      Teuchos::Array< Teuchos::RCP<Box> > boxes_to_split;
      if (root->delta_x > max_length ||
          root->delta_y > max_length ||
          root->delta_z > max_length)
        boxes_to_split.push_back(root);
      else
        parts->push_back(root);

      // Split each box until all boxes are less than tolerance
      while(boxes_to_split.size() > 0 && num_parts < max_parts) {
        Teuchos::RCP<Box> box = boxes_to_split.back();
        boxes_to_split.pop_back();
        box->split();
        ++num_parts;
        if (box->left != Teuchos::null) {
          if (box->left->delta_x > max_length ||
              box->left->delta_y > max_length ||
              box->left->delta_z > max_length)
            boxes_to_split.push_back(box->left);
          else
            parts->push_back(box->left);
        }
        if (box->right != Teuchos::null) {
          if (box->right->delta_x > max_length ||
              box->right->delta_y > max_length ||
              box->right->delta_z > max_length)
            boxes_to_split.push_back(box->right);
          else
            parts->push_back(box->right);
        }
      }

      TEUCHOS_ASSERT(parts->size() == num_parts);
    }

  };

  template <typename ordinal_type, typename scalar_type>
  struct CijkData {
    typedef ordinal_type value_type;
    typedef ordinal_type id_type;
    ordinal_type gid;
    ordinal_type i, j, k;
    scalar_type c;

    ordinal_type operator() (ordinal_type d) const {
      if (d == 0) return i;
      if (d == 1) return j;
      if (d == 2) return k;
      return -1;
    }

    ordinal_type ID() const { return gid; }
  };

  enum CijkSymmetryType {
    CIJK_NO_SYMMETRY,
    CIJK_TWO_WAY_SYMMETRY,
    CIJK_SIX_WAY_SYMMETRY
  };

  template <typename ordinal_type, typename scalar_type>
  Teuchos::ArrayRCP< CijkData<ordinal_type,scalar_type> >
  build_cijk_coordinate_list(
    const Sparse3Tensor<ordinal_type,scalar_type>& Cijk,
    CijkSymmetryType symmetry_type) {
    typedef Sparse3Tensor<ordinal_type,scalar_type> Cijk_type;
    typedef typename Cijk_type::k_iterator k_iterator;
    typedef typename Cijk_type::kj_iterator kj_iterator;
    typedef typename Cijk_type::kji_iterator kji_iterator;

    ordinal_type num_cijk = Cijk.num_entries();
    Teuchos::ArrayRCP< CijkData<ordinal_type,scalar_type> > coordinate_list(
      num_cijk);
    ordinal_type idx = 0;
    k_iterator k_begin = Cijk.k_begin();
    k_iterator k_end = Cijk.k_end();
    for (k_iterator k_it=k_begin; k_it!=k_end; ++k_it) {
      ordinal_type k = index(k_it);
      kj_iterator j_begin = Cijk.j_begin(k_it);
      kj_iterator j_end = Cijk.j_end(k_it);
      for (kj_iterator j_it = j_begin; j_it != j_end; ++j_it) {
        ordinal_type j = index(j_it);
        kji_iterator i_begin = Cijk.i_begin(j_it);
        kji_iterator i_end = Cijk.i_end(j_it);
        for (kji_iterator i_it = i_begin; i_it != i_end; ++i_it) {
          ordinal_type i = index(i_it);
          if (symmetry_type == CIJK_NO_SYMMETRY) {
            coordinate_list[idx].i = i;
            coordinate_list[idx].j = j;
            coordinate_list[idx].k = k;
            coordinate_list[idx].c = value(i_it);
            coordinate_list[idx].gid = idx;
            ++idx;
          }
          else if (symmetry_type == CIJK_TWO_WAY_SYMMETRY && j >= k) {
            coordinate_list[idx].i = i;
            coordinate_list[idx].j = j;
            coordinate_list[idx].k = k;
            if (j == k)
              coordinate_list[idx].c = 0.5*value(i_it);
            else
              coordinate_list[idx].c = value(i_it);
            coordinate_list[idx].gid = idx;
            ++idx;
          }
          else if (symmetry_type == CIJK_SIX_WAY_SYMMETRY && i >= j && j >= k) {
            coordinate_list[idx].i = i;
            coordinate_list[idx].j = j;
            coordinate_list[idx].k = k;
            if (i == j && j == k)
              coordinate_list[idx].c = (1.0/6.0)*value(i_it);
            else
              coordinate_list[idx].c = value(i_it);
            coordinate_list[idx].gid = idx;
            ++idx;
          }
        }
      }
    }
    coordinate_list.resize(idx);

    return coordinate_list;
  }

} // namespace Stokhos

#endif // STOKHOS_SPARSE_3_TENSOR_PARTITION_HPP
