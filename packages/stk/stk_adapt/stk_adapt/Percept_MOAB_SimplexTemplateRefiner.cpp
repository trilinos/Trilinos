
//p #include "EdgeSizeEvaluator.hpp"
//p #include "RefinerTagManager.hpp"

//p #include "moab/Interface.hpp"

#include <iostream>
#include <stack>
#include <string>
#include <stdexcept>

#include <boost/lexical_cast.hpp>

#include <stk_adapt/Percept_MOAB_SimplexTemplateRefiner.hpp>

// Static arrays holding parametric coordinates of element vertices
#if 0
static double MBVertexParametric[] = { 0., 0., 0. };
static double MBEdgeParametric[]   = { 0., 0., 0.,   1., 0., 0. };
static double MBTriParametric[]    = { 0., 0., 0.,   1., 0., 0.,   0., 1., 0. };
static double MBTetParametric[]    = { 0., 0., 0.,   1., 0., 0.,   0., 1., 0.,   0., 0., 1. };
#endif

#ifdef MB_DEBUG_TESSELLATOR
#  define MB_TESSELLATOR_INCR_CASE_COUNT(cs) this->case_counts[cs]++
#  define MB_TESSELLATOR_INCR_SUBCASE_COUNT(cs,sc) this->subcase_counts[cs][sc]++
#else // MB_DEBUG_TESSELLATOR
#  define MB_TESSELLATOR_INCR_CASE_COUNT(cs)
#  define MB_TESSELLATOR_INCR_SUBCASE_COUNT(cs,sc)
#endif // MB_DEBUG_TESSELLATOR

#define PERCEPT_DEBUG 1

namespace moab {

  static void error(int i)
  {
    std::string ii =  "Percept_MOAB_SimplexTemplateRefiner:: err # "+ boost::lexical_cast<std::string>(i);
    throw std::logic_error("test");
  }

  //p this is from some comments below, and coded in tet_edges.  It's only used for disambiguation of same-length edges.
  //p   * Edge 0-1, Edge 1-2, Edge 2-0, Edge 0-3, Edge 1-3, Edge 2-3,
  int tet_edges[6][2] = {{0,1}, {1,2}, {2,0}, {0,3}, {1,3}, {2,3}};

  double compute_edge_length_squared(double *c0, double *c1)
  {
    return 
      (c0[0]-c1[0])*(c0[0]-c1[0])+
      (c0[1]-c1[1])*(c0[1]-c1[1])+
      (c0[2]-c1[2])*(c0[2]-c1[2]);
  }

  /**\brief Refine a tetrahedron.

   */
  bool SimplexTemplateRefiner::refine_3_simplex(std::vector<TetTupleInt>& new_tets,
                                                unsigned edge_marks[6],
                                                int max_depth,
                                                double* v0, void* t0, EntityHandle h0,
                                                double* v1, void* t1, EntityHandle h1,
                                                double* v2, void* t2, EntityHandle h2,
                                                double* v3, void* t3, EntityHandle h3 )
  {
    int edge_code = 0;

    double* midpt0c=0;
    double* midpt1c=0;
    double* midpt2c=0;
    double* midpt3c=0;
    double* midpt4c=0;
    double* midpt5c=0;

    void* midpt0t=0;
    void* midpt1t=0;
    void* midpt2t=0;
    void* midpt3t=0;
    void* midpt4t=0;
    void* midpt5t=0;

    EntityHandle midpt0h=0;
    EntityHandle midpt1h=0;
    EntityHandle midpt2h=0;
    EntityHandle midpt3h=0;
    EntityHandle midpt4h=0;
    EntityHandle midpt5h=0;

    if ( max_depth-- > 0 )
      {
        midpt0c = this->heap_coord_storage();
        midpt1c = this->heap_coord_storage();
        midpt2c = this->heap_coord_storage();
        midpt3c = this->heap_coord_storage();
        midpt4c = this->heap_coord_storage();
        midpt5c = this->heap_coord_storage();

        midpt0t = this->heap_tag_storage();
        midpt1t = this->heap_tag_storage();
        midpt2t = this->heap_tag_storage();
        midpt3t = this->heap_tag_storage();
        midpt4t = this->heap_tag_storage();
        midpt5t = this->heap_tag_storage();

        //pfor ( int i = 0; i < 6; ++ i )
        for ( int i = 0; i < 3; ++ i )
          {
            midpt0c[i] = ( v0[i] + v1[i] ) * .5;
            midpt1c[i] = ( v1[i] + v2[i] ) * .5;
            midpt2c[i] = ( v2[i] + v0[i] ) * .5;
            midpt3c[i] = ( v0[i] + v3[i] ) * .5;
            midpt4c[i] = ( v1[i] + v3[i] ) * .5;
            midpt5c[i] = ( v2[i] + v3[i] ) * .5;
          }

        /*p
          (*this->tag_assigner)( v0, t0, h0, midpt0c, midpt0t, v1, t1, h1 );
          (*this->tag_assigner)( v1, t1, h1, midpt1c, midpt1t, v2, t2, h2 );
          (*this->tag_assigner)( v2, t2, h2, midpt2c, midpt2t, v0, t0, h0 );
          (*this->tag_assigner)( v0, t0, h0, midpt3c, midpt3t, v3, t3, h3 );
          (*this->tag_assigner)( v1, t1, h1, midpt4c, midpt4t, v3, t3, h3 );
          (*this->tag_assigner)( v2, t2, h2, midpt5c, midpt5t, v3, t3, h3 );
        */

        for (int ie=0; ie < 6; ie++)
          {
            if (edge_marks[ie]) edge_code |= (1 << ie);
          }
      }

    //p perform edge length calculations in a consistent way (see below, we do it on permuted values)
    double edge_length2[6];
    for ( int ei = 0; ei < 6; ++ ei )
      edge_length2[ei] = 0.;


    if (PERCEPT_DEBUG)
      std::cout << "tmp PM edge_code= " << edge_code << std::endl;

    if ( ! edge_code )
      {
        // No edges to subdivide

        return false;
      }

    double* facept0c;
    double* facept1c;
    double* facept2c;
    double* facept3c;

    facept0c = this->heap_coord_storage();
    facept1c = this->heap_coord_storage();
    facept2c = this->heap_coord_storage();
    facept3c = this->heap_coord_storage();

    double* vertex_coords[14] = {
      v0, v1, v2, v3, 
      midpt0c, midpt1c, midpt2c, 
      midpt3c, midpt4c, midpt5c,
      facept0c, facept1c, facept2c, facept3c
    };

    void* facept0t = this->heap_tag_storage();
    void* facept1t = this->heap_tag_storage();
    void* facept2t = this->heap_tag_storage();
    void* facept3t = this->heap_tag_storage();

    void* vertex_tags[14] = {
      t0, t1, t2, t3, 
      midpt0t, midpt1t, midpt2t, 
      midpt3t, midpt4t, midpt5t,
      facept0t, facept1t, facept2t, facept3t
    };

    EntityHandle vertex_hash[14] = {
      h0, h1, h2, h3,
      midpt0h, midpt1h, midpt2h,
      midpt3h, midpt4h, midpt5h,
      0, 0, 0, 0
    };

    // Generate tetrahedra that are compatible except when edge
    // lengths are equal on indeterminately subdivided faces.
    double* permuted_coords[14];
    void* permuted_tags[14];
    EntityHandle permuted_hash[14];
    int permuted_local_ids[14];   // this is just a copy of permutations_from_index, for readability
    double permlen[6]; // permuted edge lengths
    int C = SimplexTemplateRefiner::template_index[edge_code][0];
    int P = SimplexTemplateRefiner::template_index[edge_code][1];

    if (PERCEPT_DEBUG)
      std::cout << "tmp PM C,P= " << C << " " << P << std::endl;
  
    // 1. Permute the tetrahedron into our canonical configuration
    for ( int i = 0; i < 14; ++ i )
      {
        permuted_local_ids[i] = SimplexTemplateRefiner::permutations_from_index[P][i];
        permuted_coords[i] = vertex_coords[SimplexTemplateRefiner::permutations_from_index[P][i]];
        permuted_tags[i] = vertex_tags[SimplexTemplateRefiner::permutations_from_index[P][i]];
        permuted_hash[i] = vertex_hash[SimplexTemplateRefiner::permutations_from_index[P][i]];
      }

    /*p
    for ( int i = 4 ; i < 10; ++ i )
      {
        // permute edge lengths too
        permlen[i-4]  = edge_length2[SimplexTemplateRefiner::permutations_from_index[P][i] - 4];
      }
    */
    for (int i = 0; i < 6; i++)
      {

        double *c0 = permuted_coords[tet_edges[i][0]];
        double *c1 = permuted_coords[tet_edges[i][1]];
        bool reverse=false;
        if (permuted_hash[tet_edges[i][0]] > permuted_hash[tet_edges[i][1]])
          {
            c0 = permuted_coords[tet_edges[i][1]];
            c1 = permuted_coords[tet_edges[i][0]];
            reverse=true;
          }
        permlen[i] = compute_edge_length_squared(c0, c1);
      }

    int comparison_bits;
    std::stack<int*> output_tets;
    std::stack<int*> output_perm;
    std::stack<int>  output_sign;

    // cout << "Case " << C << "  Permutation " << P << endl;
    // 2. Generate tetrahedra based on the configuration.
    //    Note that case 0 is handled above (edgeCode == 0).
  
    switch ( C )
      {
      case 1: // Ruprecht-Müller Case 1
        MB_TESSELLATOR_INCR_CASE_COUNT(0);
        output_tets.push( SimplexTemplateRefiner::templates + 0 );
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );
        MB_TESSELLATOR_INCR_SUBCASE_COUNT(0,0);
        break;
      case 2: // Ruprecht-Müller Case 2a
        comparison_bits = 
          ( permlen[0] <= permlen[1] ? 1 : 0 ) | ( permlen[0] >= permlen[1] ? 2 : 0 ) |
          0;

        // a tie-breaker for a pair of edges that have the same length - add all the vertex handles
        //   and choose one case or the other depending on if the sum is odd or even

#define VH(iedge) (permuted_hash[tet_edges[iedge][0]] + permuted_hash[tet_edges[iedge][1]])
#define CMP_VH(ie,je) ((VH(ie)+VH(je)) % 2 == 0)

        if ( ( comparison_bits & 3 ) == 3 )
          {
#if PERCEPT_MOAB_ALLOW_FACE_DISAMB
            // Compute face point and tag
            for ( int i = 0; i < 6; ++ i )
              {
                permuted_coords[10][i] = ( permuted_coords[0][i] + permuted_coords[2][i] ) * .375 + permuted_coords[1][i] * .25;
              }
            (*this->tag_assigner)( permuted_tags[0], permuted_tags[2], permuted_tags[1], permuted_tags[10] );
            permuted_hash[10] = (*this->output_functor)(
                                                        permuted_hash[0], permuted_hash[2], permuted_hash[1], permuted_coords[10], permuted_tags[10] );
#else
            comparison_bits -=  3 ;
            if (CMP_VH(0,1))
              comparison_bits |= 1;
            else
              comparison_bits |= 2;
#endif
          }

        MB_TESSELLATOR_INCR_CASE_COUNT(1);
        output_tets.push( SimplexTemplateRefiner::templates + 9 );
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );
        MB_TESSELLATOR_INCR_SUBCASE_COUNT(1,0);
        switch ( comparison_bits )
          {
          case 2: // 0>1
            output_tets.push( SimplexTemplateRefiner::templates + 14 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(1,1);
            break;
          case 1: // 1>0
            output_tets.push( SimplexTemplateRefiner::templates + 14 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[13] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(1,2);
            break;
          case 3: // 0=1
            error(1);
            output_tets.push( SimplexTemplateRefiner::templates + 23 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(1,3);
            break;
          }
        break;
      case 3: // Ruprecht-Müller Case 2b
        MB_TESSELLATOR_INCR_CASE_COUNT(2);
        output_tets.push( SimplexTemplateRefiner::templates + 40 );
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );
        MB_TESSELLATOR_INCR_SUBCASE_COUNT(2,0);
        break;
      case 4: // Ruprecht-Müller Case 3a
        comparison_bits = 
          ( permlen[0] <= permlen[3] ? 1 : 0 ) | ( permlen[0] >= permlen[3] ? 2 : 0 ) |
          ( permlen[2] <= permlen[3] ? 4 : 0 ) | ( permlen[2] >= permlen[3] ? 8 : 0 ) |
          ( permlen[0] <= permlen[2] ? 16 : 0 ) | ( permlen[0] >= permlen[2] ? 32 : 0 ) |
          0;
        if ( ( comparison_bits & 3 ) == 3 )
          {
#if PERCEPT_MOAB_ALLOW_FACE_DISAMB
            // Compute face point and tag
            for ( int i = 0; i < 6; ++ i )
              {
                permuted_coords[11][i] = ( permuted_coords[1][i] + permuted_coords[3][i] ) * .375 + permuted_coords[0][i] * .25;
              }
            (*this->tag_assigner)( permuted_tags[1], permuted_tags[3], permuted_tags[0], permuted_tags[11] );
            permuted_hash[11] = (*this->output_functor)(
                                                        permuted_hash[1], permuted_hash[3], permuted_hash[0], permuted_coords[11], permuted_tags[11] );
#else
            comparison_bits -=  3 ;
            if (CMP_VH(0,3))
              comparison_bits |= 1;
            else
              comparison_bits |= 2;
#endif
          }
        if ( ( comparison_bits & 12 ) == 12 )
          {
#if PERCEPT_MOAB_ALLOW_FACE_DISAMB
            // Compute face point and tag
            for ( int i = 0; i < 6; ++ i )
              {
                permuted_coords[13][i] = ( permuted_coords[2][i] + permuted_coords[3][i] ) * .375 + permuted_coords[0][i] * .25;
              }
            (*this->tag_assigner)( permuted_tags[2], permuted_tags[3], permuted_tags[0], permuted_tags[13] );
            permuted_hash[13] = (*this->output_functor)(
                                                        permuted_hash[2], permuted_hash[3], permuted_hash[0], permuted_coords[13], permuted_tags[13] );
#else
            comparison_bits -=  12 ;
            if (CMP_VH(2,3))
              comparison_bits |= 4;
            else
              comparison_bits |= 8;
#endif
          }
        if ( ( comparison_bits & 48 ) == 48 )
          {
#if PERCEPT_MOAB_ALLOW_FACE_DISAMB
            // Compute face point and tag
            for ( int i = 0; i < 6; ++ i )
              {
                permuted_coords[10][i] = ( permuted_coords[1][i] + permuted_coords[2][i] ) * .375 + permuted_coords[0][i] * .25;
              }
            (*this->tag_assigner)( permuted_tags[1], permuted_tags[2], permuted_tags[0], permuted_tags[10] );
            permuted_hash[10] = (*this->output_functor)(
                                                        permuted_hash[1], permuted_hash[2], permuted_hash[0], permuted_coords[10] , permuted_tags[10] );
#else
            comparison_bits -=  48 ;
            if (CMP_VH(0,2))
              comparison_bits |= 16;
            else
              comparison_bits |= 32;
#endif
          }

        MB_TESSELLATOR_INCR_CASE_COUNT(3);
        output_tets.push( SimplexTemplateRefiner::templates + 57 );
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );
        MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,0);
        switch ( comparison_bits )
          {
          case 42: // 0>2>3<0
            output_tets.push( SimplexTemplateRefiner::templates + 62 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,1);
            break;
          case 25: // 2>3>0<2
            output_tets.push( SimplexTemplateRefiner::templates + 62 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[11] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,2);
            break;
          case 37: // 3>0>2<3
            output_tets.push( SimplexTemplateRefiner::templates + 62 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[3] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,3);
            break;
          case 21: // 3>2>0<3
            output_tets.push( SimplexTemplateRefiner::templates + 62 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[22] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,4);
            break;
          case 26: // 2>0>3<2
            output_tets.push( SimplexTemplateRefiner::templates + 62 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[12] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,5);
            break;
          case 38: // 0>3>2<0
            output_tets.push( SimplexTemplateRefiner::templates + 62 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[15] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,6);
            break;
          case 58: // 0=2>3<0
            error(2);
            output_tets.push( SimplexTemplateRefiner::templates + 75 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,7);
            break;
          case 29: // 2=3>0<2
            error(3);
            output_tets.push( SimplexTemplateRefiner::templates + 75 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[11] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,8);
            break;
          case 39: // 0=3>2<0
            error(4);
            output_tets.push( SimplexTemplateRefiner::templates + 75 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[3] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,9);
            break;
          case 53: // 3>0=2<3
            error(5);
            output_tets.push( SimplexTemplateRefiner::templates + 96 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,10);
            break;
          case 46: // 0>2=3<0
            error(6);
            output_tets.push( SimplexTemplateRefiner::templates + 96 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[11] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,11);
            break;
          case 27: // 2>0=3<2
            error(7);
            output_tets.push( SimplexTemplateRefiner::templates + 96 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[3] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,12);
            break;
          case 63: // 0=2=3=0
            error(8);
            output_tets.push( SimplexTemplateRefiner::templates + 117 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(3,13);
            break;
          }
        break;
      case 5: // Ruprecht-Müller Case 3b
        MB_TESSELLATOR_INCR_CASE_COUNT(4);
        output_tets.push( SimplexTemplateRefiner::templates + 162 );
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );
        MB_TESSELLATOR_INCR_SUBCASE_COUNT(4,0);
        break;
      case 6: // Ruprecht-Müller Case 3c
        comparison_bits = 
          ( permlen[0] <= permlen[1] ? 1 : 0 ) | ( permlen[0] >= permlen[1] ? 2 : 0 ) |
          ( permlen[0] <= permlen[3] ? 4 : 0 ) | ( permlen[0] >= permlen[3] ? 8 : 0 ) |
          0;
        if ( ( comparison_bits & 3 ) == 3 )
          {
#if PERCEPT_MOAB_ALLOW_FACE_DISAMB
            // Compute face point and tag
            for ( int i = 0; i < 6; ++ i )
              {
                permuted_coords[10][i] = ( permuted_coords[0][i] + permuted_coords[2][i] ) * .375 + permuted_coords[1][i] * .25;
              }
            (*this->tag_assigner)( permuted_tags[0], permuted_tags[2], permuted_tags[1], permuted_tags[10] );
            permuted_hash[10] = (*this->output_functor)(
                                                        permuted_hash[0], permuted_hash[2], permuted_hash[1], permuted_coords[10], permuted_tags[10] );
#else
            comparison_bits -=  3 ;
            if (CMP_VH(0,1))
              comparison_bits |= 1;
            else
              comparison_bits |= 2;
#endif

          }
        if ( ( comparison_bits & 12 ) == 12 )
          {
#if PERCEPT_MOAB_ALLOW_FACE_DISAMB
            // Compute face point and tag
            for ( int i = 0; i < 6; ++ i )
              {
                permuted_coords[11][i] = ( permuted_coords[1][i] + permuted_coords[3][i] ) * .375 + permuted_coords[0][i] * .25;
              }
            (*this->tag_assigner)( permuted_tags[1], permuted_tags[3], permuted_tags[0], permuted_tags[11] );
            permuted_hash[11] = (*this->output_functor)(
                                                        permuted_hash[1], permuted_hash[3], permuted_hash[0], permuted_coords[11], permuted_tags[11] );
#else
            comparison_bits -=  12 ;
            if (CMP_VH(0,3))
              comparison_bits |= 4;
            else
              comparison_bits |= 8;
#endif

          }
        MB_TESSELLATOR_INCR_CASE_COUNT(5);
        switch ( comparison_bits )
          {
          case 10: // 0>1,0>3
            output_tets.push( SimplexTemplateRefiner::templates + 179 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(5,0);
            break;
          case 5: // 1>0,3>0
            output_tets.push( SimplexTemplateRefiner::templates + 200 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(5,1);
            break;
          case 6: // 0>1,3>0
            output_tets.push( SimplexTemplateRefiner::templates + 221 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(5,2);
            break;
          case 9: // 1>0,0>3
            output_tets.push( SimplexTemplateRefiner::templates + 242 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(5,3);
            break;
          case 11: // 0=1,0>3
            error(9);
            output_tets.push( SimplexTemplateRefiner::templates + 263 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(5,4);
            break;
          case 14: // 0=3,0>1
            error(10);
            output_tets.push( SimplexTemplateRefiner::templates + 263 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[5] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(5,5);
            break;
          case 7: // 3>0,0=1
            error(11);
            output_tets.push( SimplexTemplateRefiner::templates + 292 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(5,6);
            break;
          case 13: // 1>0,0=3
            error(12);
            output_tets.push( SimplexTemplateRefiner::templates + 292 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[5] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(5,7);
            break;
          case 15: // 0=1,0=3
            error(13);
            output_tets.push( SimplexTemplateRefiner::templates + 321 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(5,8);
            break;
          }
        break;
      case 7: // Ruprecht-Müller Case 3d
        comparison_bits = 
          ( permlen[0] <= permlen[2] ? 1 : 0 ) | ( permlen[0] >= permlen[2] ? 2 : 0 ) |
          ( permlen[0] <= permlen[4] ? 4 : 0 ) | ( permlen[0] >= permlen[4] ? 8 : 0 ) |
          0;
        if ( ( comparison_bits & 3 ) == 3 )
          {
#if PERCEPT_MOAB_ALLOW_FACE_DISAMB
            // Compute face point and tag
            for ( int i = 0; i < 6; ++ i )
              {
                permuted_coords[10][i] = ( permuted_coords[1][i] + permuted_coords[2][i] ) * .375 + permuted_coords[0][i] * .25;
              }
            (*this->tag_assigner)( permuted_tags[1], permuted_tags[2], permuted_tags[0], permuted_tags[10] );
            permuted_hash[10] = (*this->output_functor)(
                                                        permuted_hash[1], permuted_hash[2], permuted_hash[0], permuted_coords[10], permuted_tags[10] );
#else
            comparison_bits -=  3 ;
            if (CMP_VH(0,2))
              comparison_bits |= 1;
            else
              comparison_bits |= 2;
#endif
          }
        if ( ( comparison_bits & 12 ) == 12 )
          {
#if PERCEPT_MOAB_ALLOW_FACE_DISAMB
            // Compute face point and tag
            for ( int i = 0; i < 6; ++ i )
              {
                permuted_coords[11][i] = ( permuted_coords[0][i] + permuted_coords[3][i] ) * .375 + permuted_coords[1][i] * .25;
              }
            (*this->tag_assigner)( permuted_tags[0], permuted_tags[3], permuted_tags[1], permuted_tags[11] );
            permuted_hash[11] = (*this->output_functor)(
                                                        permuted_hash[0], permuted_hash[3], permuted_hash[1], permuted_coords[11], permuted_tags[11] );
#else
            comparison_bits -=  12 ;
            if (CMP_VH(0,4))
              comparison_bits |= 4;
            else
              comparison_bits |= 8;
#endif
          }
        MB_TESSELLATOR_INCR_CASE_COUNT(6);
        switch ( comparison_bits )
          {
          case 10: // 0>4,0>2
            output_tets.push( SimplexTemplateRefiner::templates + 362 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(6,0);
            break;
          case 5: // 4>0,2>0
            output_tets.push( SimplexTemplateRefiner::templates + 383 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(6,1);
            break;
          case 9: // 0>4,2>0
            output_tets.push( SimplexTemplateRefiner::templates + 404 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(6,2);
            break;
          case 6: // 4>0,0>2
            output_tets.push( SimplexTemplateRefiner::templates + 425 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(6,3);
            break;
          case 14: // 0=4,0>2
            error(14);
            output_tets.push( SimplexTemplateRefiner::templates + 446 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(6,4);
            break;
          case 11: // 0=2,0>4
            error(15);
            output_tets.push( SimplexTemplateRefiner::templates + 446 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[5] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(6,5);
            break;
          case 13: // 2>0,0=4
            error(16);
            output_tets.push( SimplexTemplateRefiner::templates + 475 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(6,6);
            break;
          case 7: // 4>0,0=2
            error(17);
            output_tets.push( SimplexTemplateRefiner::templates + 475 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[5] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(6,7);
            break;
          case 15: // 0=4,0=2
            error(18);
            output_tets.push( SimplexTemplateRefiner::templates + 504 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(6,8);
            break;
          }
        break;
      case 8: // Ruprecht-Müller Case 4a
        comparison_bits = 
          ( permlen[4] <= permlen[5] ? 1 : 0 ) | ( permlen[4] >= permlen[5] ? 2 : 0 ) |
          ( permlen[3] <= permlen[4] ? 4 : 0 ) | ( permlen[3] >= permlen[4] ? 8 : 0 ) |
          0;
        if ( ( comparison_bits & 3 ) == 3 )
          {
#if PERCEPT_MOAB_ALLOW_FACE_DISAMB
            // Compute face point and tag
            for ( int i = 0; i < 6; ++ i )
              {
                permuted_coords[12][i] = ( permuted_coords[1][i] + permuted_coords[2][i] ) * .375 + permuted_coords[3][i] * .25;
              }
            (*this->tag_assigner)( permuted_tags[1], permuted_tags[2], permuted_tags[3], permuted_tags[12] );
            permuted_hash[12] = (*this->output_functor)(
                                                        permuted_hash[1], permuted_hash[2], permuted_hash[3], permuted_coords[12], permuted_tags[12] );
#else
            comparison_bits -=  3 ;
            if (CMP_VH(4,5))
              comparison_bits |= 1;
            else
              comparison_bits |= 2;
#endif

          }
        if ( ( comparison_bits & 12 ) == 12 )
          {
#if PERCEPT_MOAB_ALLOW_FACE_DISAMB
            // Compute face point and tag
            for ( int i = 0; i < 6; ++ i )
              {
                permuted_coords[11][i] = ( permuted_coords[0][i] + permuted_coords[1][i] ) * .375 + permuted_coords[3][i] * .25;
              }
            (*this->tag_assigner)( permuted_tags[0], permuted_tags[1], permuted_tags[3], permuted_tags[11] );
            permuted_hash[11] = (*this->output_functor)(
                                                        permuted_hash[0], permuted_hash[1], permuted_hash[3], permuted_coords[11], permuted_tags[11] );
#else
            comparison_bits -=  12 ;
            if (CMP_VH(3,4))
              comparison_bits |= 4;
            else
              comparison_bits |= 8;
#endif
          }
        if (PERCEPT_DEBUG) std::cout << "tmp case 4a comparison_bits= " << comparison_bits << std::endl;
        MB_TESSELLATOR_INCR_CASE_COUNT(7);
        output_tets.push( SimplexTemplateRefiner::templates + 545 );
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );
        MB_TESSELLATOR_INCR_SUBCASE_COUNT(7,0);
        switch ( comparison_bits )
          {
          case 5: // 5>4>3
            output_tets.push( SimplexTemplateRefiner::templates + 554 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(7,1);
            break;
          case 10: // 3>4>5
            output_tets.push( SimplexTemplateRefiner::templates + 554 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[13] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(7,2);
            break;
          case 6: // 3<4>5
            output_tets.push( SimplexTemplateRefiner::templates + 571 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(7,3);
            break;
          case 9: // 3>4<5
            output_tets.push( SimplexTemplateRefiner::templates + 588 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(7,4);
            break;
          case 14: // 3=4>5
            error(19);
            output_tets.push( SimplexTemplateRefiner::templates + 605 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(7,5);
            break;
          case 7: // 4=5,4>3
            error(20);
            output_tets.push( SimplexTemplateRefiner::templates + 605 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[13] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(7,6);
            break;
          case 13: // 5>4,3=4
            error(21);
            output_tets.push( SimplexTemplateRefiner::templates + 630 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(7,7);
            break;
          case 11: // 3>4=5
            error(22);
            output_tets.push( SimplexTemplateRefiner::templates + 630 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[13] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(7,8);
            break;
          case 15: // 3=4=5
            error(23);
            output_tets.push( SimplexTemplateRefiner::templates + 655 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(7,9);
            break;
          }
        break;
      case 9: // Ruprecht-Müller Case 4b
        comparison_bits = 
          ( permlen[1] <= permlen[2] ? 1 : 0 ) | ( permlen[1] >= permlen[2] ? 2 : 0 ) |
          ( permlen[2] <= permlen[3] ? 4 : 0 ) | ( permlen[2] >= permlen[3] ? 8 : 0 ) |
          ( permlen[3] <= permlen[4] ? 16 : 0 ) | ( permlen[3] >= permlen[4] ? 32 : 0 ) |
          ( permlen[1] <= permlen[4] ? 64 : 0 ) | ( permlen[1] >= permlen[4] ? 128 : 0 ) |
          0;
        if ( ( comparison_bits & 3 ) == 3 )
          {
#if PERCEPT_MOAB_ALLOW_FACE_DISAMB
            // Compute face point and tag
            for ( int i = 0; i < 6; ++ i )
              {
                permuted_coords[10][i] = ( permuted_coords[1][i] + permuted_coords[0][i] ) * .375 + permuted_coords[2][i] * .25;
              }
            (*this->tag_assigner)( permuted_tags[1], permuted_tags[0], permuted_tags[2], permuted_tags[10] );
            permuted_hash[10] = (*this->output_functor)(
                                                        permuted_hash[1], permuted_hash[0], permuted_hash[2], permuted_coords[10], permuted_tags[10] );
#else
            comparison_bits -=  3 ;
            if (CMP_VH(1,2))
              comparison_bits |= 1;
            else
              comparison_bits |= 2;
#endif

          }
        if ( ( comparison_bits & 12 ) == 12 )
          {
#if PERCEPT_MOAB_ALLOW_FACE_DISAMB
            // Compute face point and tag
            for ( int i = 0; i < 6; ++ i )
              {
                permuted_coords[13][i] = ( permuted_coords[2][i] + permuted_coords[3][i] ) * .375 + permuted_coords[0][i] * .25;
              }
            (*this->tag_assigner)( permuted_tags[2], permuted_tags[3], permuted_tags[0], permuted_tags[13] );
            permuted_hash[13] = (*this->output_functor)(
                                                        permuted_hash[2], permuted_hash[3], permuted_hash[0], permuted_coords[13], permuted_tags[13] );
#else
            comparison_bits -=  12 ;
            if (CMP_VH(2,3))
              comparison_bits |= 4;
            else
              comparison_bits |= 8;
#endif
          }
        if ( ( comparison_bits & 48 ) == 48 )
          {
#if PERCEPT_MOAB_ALLOW_FACE_DISAMB
            // Compute face point and tag
            for ( int i = 0; i < 6; ++ i )
              {
                permuted_coords[11][i] = ( permuted_coords[0][i] + permuted_coords[1][i] ) * .375 + permuted_coords[3][i] * .25;
              }
            (*this->tag_assigner)( permuted_tags[0], permuted_tags[1], permuted_tags[3], permuted_tags[11] );
            permuted_hash[11] = (*this->output_functor)(
                                                        permuted_hash[0], permuted_hash[1], permuted_hash[3], permuted_coords[11], permuted_tags[11] );
#else
            comparison_bits -=  48 ;
            if (CMP_VH(3,4))
              comparison_bits |= 16;
            else
              comparison_bits |= 32;
#endif
          }
        if ( ( comparison_bits & 192 ) == 192 )
          {
#if PERCEPT_MOAB_ALLOW_FACE_DISAMB
            // Compute face point and tag
            for ( int i = 0; i < 6; ++ i )
              {
                permuted_coords[12][i] = ( permuted_coords[2][i] + permuted_coords[3][i] ) * .375 + permuted_coords[1][i] * .25;
              }
            (*this->tag_assigner)( permuted_tags[2], permuted_tags[3], permuted_tags[1], permuted_tags[12] );
            permuted_hash[12] = (*this->output_functor)(
                                                        permuted_hash[2], permuted_hash[3], permuted_hash[1], permuted_coords[12], permuted_tags[12] );
#else
            comparison_bits -=  192 ;
            if (CMP_VH(1,4))
              comparison_bits |= 64;
            else
              comparison_bits |= 128;
#endif
          }
        MB_TESSELLATOR_INCR_CASE_COUNT(8);
        switch ( comparison_bits )
          {
          case 85: // 2>1,3>2,4>3,4>1
            output_tets.push( SimplexTemplateRefiner::templates + 688 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,0);
            break;
          case 102: // 1>2,3>2,3>4,4>1
            output_tets.push( SimplexTemplateRefiner::templates + 688 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[14] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,1);
            break;
          case 170: // 1>2,2>3,3>4,1>4
            output_tets.push( SimplexTemplateRefiner::templates + 688 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[15] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,2);
            break;
          case 153: // 2>1,2>3,4>3,1>4
            output_tets.push( SimplexTemplateRefiner::templates + 688 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[5] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,3);
            break;
          case 90: // 1>2,2>3,4>3,4>1
            output_tets.push( SimplexTemplateRefiner::templates + 688 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[9] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,4);
            break;
          case 105: // 2>1,2>3,3>4,4>1
            output_tets.push( SimplexTemplateRefiner::templates + 688 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[7] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,5);
            break;
          case 165: // 2>1,3>2,3>4,1>4
            output_tets.push( SimplexTemplateRefiner::templates + 688 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[19] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,6);
            break;
          case 150: // 1>2,3>2,4>3,1>4
            output_tets.push( SimplexTemplateRefiner::templates + 688 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[23] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,7);
            break;
          case 101: // 2>1,3>2,3>4,4>1
            {
              int alternates[] = { 713, 738, -1 };
              output_tets.push( SimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 0, 1 ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,8);
            break;
          case 86: // 1>2,3>2,4>3,4>1
            {
              int alternates[] = {713, 738, -1 };
              output_tets.push( SimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 14, -1 ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[14] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,9);
            break;
          case 154: // 1>2,2>3,4>3,1>4
            {
              int alternates[] = {713, 738, -1 };
              output_tets.push( SimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 5, 1 ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[5] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,10);
            break;
          case 169: // 2>1,2>3,3>4,1>4
            {
              int alternates[] = {713, 738, -1 };
              output_tets.push( SimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 15, -1 ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[15] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,11);
            break;
          case 89: // 2>1,2>3,4>3,4>1
            output_tets.push( SimplexTemplateRefiner::templates + 763 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,12);
            break;
          case 166: // 1>2,3>2,3>4,1>4
            output_tets.push( SimplexTemplateRefiner::templates + 763 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[15] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,13);
            break;
          case 103: // 1=2,3>2,3>4,4>1
            error(24);
            output_tets.push( SimplexTemplateRefiner::templates + 788 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,14);
            break;
          case 87: // 1=2,3>2,4>3,4>1
            error(25);
            output_tets.push( SimplexTemplateRefiner::templates + 788 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[14] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,15);
            break;
          case 185: // 2>1,2>3,3=4,1>4
            error(26);
            output_tets.push( SimplexTemplateRefiner::templates + 788 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[15] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,16);
            break;
          case 186: // 1>2,2>3,3=4,1>4
            error(27);
            output_tets.push( SimplexTemplateRefiner::templates + 788 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[5] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,17);
            break;
          case 158: // 1>2,2=3,4>3,1>4
            error(28);
            output_tets.push( SimplexTemplateRefiner::templates + 788 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[9] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,18);
            break;
          case 229: // 2>1,3>2,3>4,1=4
            error(29);
            output_tets.push( SimplexTemplateRefiner::templates + 788 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[7] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,19);
            break;
          case 233: // 2>1,2>3,3>4,1=4
            error(30);
            output_tets.push( SimplexTemplateRefiner::templates + 788 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[19] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,20);
            break;
          case 94: // 1>2,2=3,4>3,4>1
            error(31);
            output_tets.push( SimplexTemplateRefiner::templates + 788 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[23] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,21);
            break;
          case 155: // 1=2,2>3,4>3,1>4
            error(32);
            output_tets.push( SimplexTemplateRefiner::templates + 825 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,22);
            break;
          case 171: // 1=2,2>3,3>4,1>4
            error(33);
            output_tets.push( SimplexTemplateRefiner::templates + 825 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[14] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,23);
            break;
          case 118: // 1>2,3>2,3=4,4>1
            error(34);
            output_tets.push( SimplexTemplateRefiner::templates + 825 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[15] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,24);
            break;
          case 117: // 2>1,3>2,3=4,4>1
            error(35);
            output_tets.push( SimplexTemplateRefiner::templates + 825 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[5] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,25);
            break;
          case 109: // 2>1,2=3,3>4,4>1
            error(36);
            output_tets.push( SimplexTemplateRefiner::templates + 825 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[9] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,26);
            break;
          case 218: // 1>2,2>3,4>3,1=4
            error(37);
            output_tets.push( SimplexTemplateRefiner::templates + 825 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[7] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,27);
            break;
          case 214: // 1>2,3>2,4>3,1=4
            error(38);
            output_tets.push( SimplexTemplateRefiner::templates + 825 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[19] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,28);
            break;
          case 173: // 2>1,2=3,3>4,1>4
            error(39);
            output_tets.push( SimplexTemplateRefiner::templates + 825 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[23] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,29);
            break;
          case 91: // 1=2,2>3,4>3,4>1
            error(40);
            output_tets.push( SimplexTemplateRefiner::templates + 862 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,30);
            break;
          case 167: // 1=2,3>2,3>4,1>4
            error(41);
            output_tets.push( SimplexTemplateRefiner::templates + 862 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[14] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,31);
            break;
          case 182: // 1>2,3>2,3=4,1>4
            error(42);
            output_tets.push( SimplexTemplateRefiner::templates + 862 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[15] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,32);
            break;
          case 121: // 2>1,2>3,3=4,4>1
            error(43);
            output_tets.push( SimplexTemplateRefiner::templates + 862 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[5] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,33);
            break;
          case 93: // 2>1,2=3,4>3,4>1
            error(44);
            output_tets.push( SimplexTemplateRefiner::templates + 862 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[9] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,34);
            break;
          case 217: // 2>1,2>3,4>3,1=4
            error(45);
            output_tets.push( SimplexTemplateRefiner::templates + 862 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[7] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,35);
            break;
          case 230: // 1>2,3>2,3>4,1=4
            error(46);
            output_tets.push( SimplexTemplateRefiner::templates + 862 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[19] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,36);
            break;
          case 174: // 1>2,2=3,3>4,1>4
            error(47);
            output_tets.push( SimplexTemplateRefiner::templates + 862 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[23] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,37);
            break;
          case 119: // 1=2,3>2,3=4,4>1
            error(48);
            output_tets.push( SimplexTemplateRefiner::templates + 899 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,38);
            break;
          case 187: // 1=2>3=4,1>4
            error(49);
            output_tets.push( SimplexTemplateRefiner::templates + 899 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[15] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,39);
            break;
          case 222: // 1>2,2=3,4>3,1=4
            error(50);
            output_tets.push( SimplexTemplateRefiner::templates + 899 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[9] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,40);
            break;
          case 237: // 2>1,2=3,3>4,1=4
            error(51);
            output_tets.push( SimplexTemplateRefiner::templates + 899 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[7] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,41);
            break;
          case 95: // 4>1=2=3,4>3
            error(52);
            output_tets.push( SimplexTemplateRefiner::templates + 944 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,42);
            break;
          case 231: // 1=2,3>2,3>4,1=4
            error(53);
            output_tets.push( SimplexTemplateRefiner::templates + 944 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[14] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,43);
            break;
          case 190: // 1>2=3=4,1>4
            error(54);
            output_tets.push( SimplexTemplateRefiner::templates + 944 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[15] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,44);
            break;
          case 249: // 2>1,2>3,3=4,1=4
            error(55);
            output_tets.push( SimplexTemplateRefiner::templates + 944 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[5] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,45);
            break;
          case 175: // 1=2=3>4,1>4
            error(56);
            output_tets.push( SimplexTemplateRefiner::templates + 993 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,46);
            break;
          case 219: // 1=2>3,4>3,1=4
            error(57);
            output_tets.push( SimplexTemplateRefiner::templates + 993 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[14] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,47);
            break;
          case 125: // 2>1,2=3=4>1
            error(58);
            output_tets.push( SimplexTemplateRefiner::templates + 993 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[15] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,48);
            break;
          case 246: // 1>2,3>2,3=4=1
            error(59);
            output_tets.push( SimplexTemplateRefiner::templates + 993 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[5] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,49);
            break;
          case 255: // 1=2=3=4=1
            error(60);
            output_tets.push( SimplexTemplateRefiner::templates + 1042 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(8,50);
            break;
          }
        break;
      case 10: // Ruprecht-Müller Case 5
        comparison_bits = 
          ( permlen[1] <= permlen[2] ? 1 : 0 ) | ( permlen[1] >= permlen[2] ? 2 : 0 ) |
          ( permlen[3] <= permlen[4] ? 4 : 0 ) | ( permlen[3] >= permlen[4] ? 8 : 0 ) |
          0;
        if ( ( comparison_bits & 3 ) == 3 )
          {
#if PERCEPT_MOAB_ALLOW_FACE_DISAMB
            // Compute face point and tag
            for ( int i = 0; i < 6; ++ i )
              {
                permuted_coords[10][i] = ( permuted_coords[1][i] + permuted_coords[0][i] ) * .375 + permuted_coords[2][i] * .25;
              }
            (*this->tag_assigner)( permuted_tags[1], permuted_tags[0], permuted_tags[2], permuted_tags[10] );
            permuted_hash[10] = (*this->output_functor)(
                                                        permuted_hash[1], permuted_hash[0], permuted_hash[2], permuted_coords[10], permuted_tags[10] );
#else
            comparison_bits -=  3 ;
            if (CMP_VH(1,2))
              comparison_bits |= 1;
            else
              comparison_bits |= 2;
#endif
          }
        if ( ( comparison_bits & 12 ) == 12 )
          {
#if PERCEPT_MOAB_ALLOW_FACE_DISAMB
            // Compute face point and tag
            for ( int i = 0; i < 6; ++ i )
              {
                permuted_coords[11][i] = ( permuted_coords[0][i] + permuted_coords[1][i] ) * .375 + permuted_coords[3][i] * .25;
              }
            (*this->tag_assigner)( permuted_tags[0], permuted_tags[1], permuted_tags[3], permuted_tags[11] );
            permuted_hash[11] = (*this->output_functor)(
                                                        permuted_hash[0], permuted_hash[1], permuted_hash[3], permuted_coords[11], permuted_tags[11] );
#else
            comparison_bits -=  12 ;
            if (CMP_VH(3,4))
              comparison_bits |= 4;
            else
              comparison_bits |= 8;
#endif
          }
        MB_TESSELLATOR_INCR_CASE_COUNT(9);
        output_tets.push( SimplexTemplateRefiner::templates + 1107 );
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );
        MB_TESSELLATOR_INCR_SUBCASE_COUNT(9,0);
        switch ( comparison_bits )
          {
          case 10: // 1>2,3>4
            output_tets.push( SimplexTemplateRefiner::templates + 1116 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(9,1);
            break;
          case 5: // 2>1,4>3
            output_tets.push( SimplexTemplateRefiner::templates + 1116 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[14] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(9,2);
            break;
          case 6: // 1>2,4>3
            {
              int alternates[] = { 1137, 1158, -1 };
              output_tets.push( SimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 0, 1 ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(9,3);
            break;
          case 9: // 2>1,3>4
            {
              int alternates[] = {1137, 1158, -1 };
              output_tets.push( SimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 14, -1 ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[14] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(9,4);
            break;
          case 11: // 1=2,3>4
            error(61);
            {
              int alternates[] = { 1179, 1212, 1245, -1 };
              output_tets.push( SimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 0, 1 ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(9,5);
            break;
          case 7: // 1=2,4>3
            error(62);
            {
              int alternates[] = {1179, 1212, 1245, -1 };
              output_tets.push( SimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 14, -1 ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[14] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(9,6);
            break;
          case 14: // 3=4,1>2
            error(63);
            {
              int alternates[] = {1179, 1212, 1245, -1 };
              output_tets.push( SimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 5, 1 ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[5] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(9,7);
            break;
          case 13: // 3=4,2>1
            error(64);
            {
              int alternates[] = {1179, 1212, 1245, -1 };
              output_tets.push( SimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 15, -1 ) );
            }
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[15] );
            output_sign.push( -1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(9,8);
            break;
          case 15: // 1=2,3=4
            error(65);
            output_tets.push( SimplexTemplateRefiner::templates + 1278 );
            output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
            output_sign.push( 1 );
            MB_TESSELLATOR_INCR_SUBCASE_COUNT(9,9);
            break;
          }
        break;
      case 11: // Ruprecht-Müller Case 6
        MB_TESSELLATOR_INCR_CASE_COUNT(10);
        output_tets.push( SimplexTemplateRefiner::templates + 1319 );
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );
        MB_TESSELLATOR_INCR_SUBCASE_COUNT(10,0);
        {
          int alternates[] = { 1336, 1353, 1370, -1 };
          output_tets.push( SimplexTemplateRefiner::templates + this->best_tets( alternates, permuted_coords, 0, 1 ) );
        }
        output_perm.push( SimplexTemplateRefiner::permutations_from_index[0] );
        output_sign.push( 1 );
        MB_TESSELLATOR_INCR_SUBCASE_COUNT(10,1);
        break;
      }

    int* tets;
    int  ntets;
    int* perm;
    int  sgn;
#ifdef MB_DEBUG_TESSELLATOR
    if ( output_tets.empty() )
      {
        cout << "Argh! Case " << C << " Perm " << P << " has no output!" << endl;
      }
#endif // MB_DEBUG_TESSELLATOR
    while ( ! output_tets.empty() )
      {
        tets = output_tets.top();
        ntets = *tets;
        tets++;
        perm = output_perm.top();
        sgn = output_sign.top();

        output_tets.pop();
        output_perm.pop();
        output_sign.pop();

        if (PERCEPT_DEBUG)
          std::cout << "tmp PM ntets= " << ntets << std::endl;
        int t;
        if ( sgn > 0 )
          {
            for ( t = 0; t < ntets; ++t )
              {
#if 0
                this->refine_3_simplex( max_depth,
                                        permuted_coords[perm[tets[0]]], permuted_tags[perm[tets[0]]], permuted_hash[perm[tets[0]]],
                                        permuted_coords[perm[tets[1]]], permuted_tags[perm[tets[1]]], permuted_hash[perm[tets[1]]],
                                        permuted_coords[perm[tets[2]]], permuted_tags[perm[tets[2]]], permuted_hash[perm[tets[2]]],
                                        permuted_coords[perm[tets[3]]], permuted_tags[perm[tets[3]]], permuted_hash[perm[tets[3]]]
                                        );
#endif
                TetTupleInt nt = TetTupleInt(permuted_local_ids[perm[tets[0]]],
                                             permuted_local_ids[perm[tets[1]]],
                                             permuted_local_ids[perm[tets[2]]],
                                             permuted_local_ids[perm[tets[3]]]);

                new_tets.push_back(nt);

                if (PERCEPT_DEBUG)
                  std::cout << "tmp PM new tet= " << nt << std::endl;

                tets += 4;
              }
          }
        else
          {
            // we have an inverted tet... reverse the first 2 vertices
            // so the orientation is positive.
            for ( t = 0; t < ntets; ++t )
              {
#if 0
                this->refine_3_simplex( max_depth,
                                        permuted_coords[perm[tets[1]]], permuted_tags[perm[tets[1]]], permuted_hash[perm[tets[1]]],
                                        permuted_coords[perm[tets[0]]], permuted_tags[perm[tets[0]]], permuted_hash[perm[tets[0]]],
                                        permuted_coords[perm[tets[2]]], permuted_tags[perm[tets[2]]], permuted_hash[perm[tets[2]]],
                                        permuted_coords[perm[tets[3]]], permuted_tags[perm[tets[3]]], permuted_hash[perm[tets[3]]]
                                        );
#endif
                TetTupleInt nt = TetTupleInt(permuted_local_ids[perm[tets[1]]],
                                             permuted_local_ids[perm[tets[0]]],
                                             permuted_local_ids[perm[tets[2]]],
                                             permuted_local_ids[perm[tets[3]]]);

                new_tets.push_back(nt);

                if (PERCEPT_DEBUG)
                  std::cout << "tmp PM new tet= " << nt << std::endl;

                tets += 4;
              }
          }
      }

    return true;
  }

  /*
   * The array below is indexed by the edge code for a tetrahedron.
   * Looking up a row with a tet's edge code will return C and P.
   * C is a configuration number and P is a permutation index. 
   *
   * C is based on the case number from Ruprecht and
   * Müller's (1998) paper on adaptive tetrahedra. (The case
   * numbers are shown to the left of the row in the column
   * labeled case. The only difference is that we introduce
   * a case 3d which is part of case 3c in the paper.)
   *
   * P is an index into the permutations_from_index array below,
   * and is used to transform the current tetrahedron into
   * the canonical configuration associated with C.
   *
   * The 6-digit binary number to the left (which is shown in
   * the horribly UNconventional LSB->MSB order) is the edge
   * code for the row. The 6 digits correspond to the 6 edges
   * of the tetrahedron; a '0' implies no subdivision while
   * a '1' implies subdivision should occur. The ordering of
   * the bits is
   *
   * Edge 0-1, Edge 1-2, Edge 2-0, Edge 0-3, Edge 1-3, Edge 2-3,
   *
   * where the numbers are vertices of the tetrahedron 0-1-2-3.
   * Note that Tet 0-1-2-3 must be positive (i.e., the plane
   * specified by Triangle 0-1-2 must have a normal pointing
   * towards vertex 3, and Triangle 0-1-2's normal must be
   * calculated using the cross-product (Edge 0-1) x (Edge 0-2)).
   *
   * ===========
   * References:
   * (Ruprect and Müller, 1998) A Scheme for Edge-based Adaptive
   *   Tetrahedron Subdivision, Mathematical Visualization (eds.
   *   Hege and Polthier), pp. 61--70. Springer-Verlag. 1998.
   */
  int SimplexTemplateRefiner::template_index[64][2] =
    {
      /*      code case      C    P */
      /* 000000  0  0  */ {  0,   0 },
      /* 100000  1  1  */ {  1,   0 },
      /* 010000  2  1  */ {  1,   1 },
      /* 110000  3  2a */ {  2,   0 },
      /* 001000  4  1  */ {  1,   2 },
      /* 101000  5  2a */ {  2,   2 },
      /* 011000  6  2a */ {  2,   1 },
      /* 111000  7  3b */ {  5,  11 },
      /* 000100  8  1  */ {  1,  10 },
      /* 100100  9  2a */ {  2,   5 },
      /* 010100 10  2b */ {  3,   1 },
      /* 110100 11  3c */ {  6,   0 },
      /* 001100 12  2a */ {  2,  10 },
      /* 101100 13  3a */ {  4,   0 },
      /* 011100 14  3d */ {  7,   2 },
      /* 111100 15  4a */ {  8,   6 },
      /* 000010 16  1  */ {  1,   6 },
      /* 100010 17  2a */ {  2,   4 },
      /* 010010 18  2a */ {  2,   8 },
      /* 110010 19  3a */ {  4,   1 },
      /* 001010 20  2b */ {  3,   2 },
      /* 101010 21  3d */ {  7,   0 },
      /* 011010 22  3c */ {  6,   1 },
      /* 111010 23  4a */ {  8,   9 },
      /* 000110 24  2a */ {  2,   3 },
      /* 100110 25  3b */ {  5,   0 },
      /* 010110 26  3d */ {  7,   4 },
      /* 110110 27  4a */ {  8,  11 },
      /* 001110 28  3c */ {  6,  10 },
      /* 101110 29  4a */ {  8,   7 },
      /* 011110 30  4b */ {  9,   0 },
      /* 111110 31  5  */ { 10,   7 },
      /* 000001 32  1  */ {  1,   7 },
      /* 100001 33  2b */ {  3,   0 },
      /* 010001 34  2a */ {  2,   7 },
      /* 110001 35  3d */ {  7,   1 },
      /* 001001 36  2a */ {  2,  11 },
      /* 101001 37  3c */ {  6,   2 },
      /* 011001 38  3a */ {  4,   2 },
      /* 111001 39  4a */ {  8,   3 },
      /* 000101 40  2a */ {  2,   9 },
      /* 100101 41  3d */ {  7,  10 },
      /* 010101 42  3c */ {  6,   7 },
      /* 110101 43  4b */ {  9,   2 },
      /* 001101 44  3b */ {  5,   7 },
      /* 101101 45  4a */ {  8,   8 },
      /* 011101 46  4a */ {  8,   4 },
      /* 111101 47  5  */ { 10,   6 },
      /* 000011 48  2a */ {  2,   6 },
      /* 100011 49  3c */ {  6,   4 },
      /* 010011 50  3b */ {  5,   1 },
      /* 110011 51  4a */ {  8,  10 },
      /* 001011 52  3d */ {  7,   7 },
      /* 101011 53  4b */ {  9,   1 },
      /* 011011 54  4a */ {  8,   5 },
      /* 111011 55  5  */ { 10,  10 },
      /* 000111 56  3a */ {  4,  10 },
      /* 100111 57  4a */ {  8,   1 },
      /* 010111 58  4a */ {  8,   2 },
      /* 110111 59  5  */ { 10,   2 },
      /* 001111 60  4a */ {  8,   0 },
      /* 101111 61  5  */ { 10,   1 },
      /* 011111 62  5  */ { 10,   0 },
      /* 111111 63  6  */ { 11,   0 },
    };


  /* Does this mean anything? If so, then you are either 
   * superstitious or much more clever than I (or both?).
   */
  /* permutation index, P:  0  1  2  3  4  5  6  7  8  9 10 11 */
  /* number of references: 12  9  9  3  4  2  5  6  2  3  7  2 */


  /*
   * The array below is a list of all the _positive_
   * permutations of Tetrahedron 0-1-2-3. Given a
   * permutation index, it returns a row of 14 values:
   * these are the vertex numbers of the permuted
   * tetrahedron. The first 4 values are the permuted
   * corner indices, the next 6 values are the
   * permuted edge midpoint indices, and the final
   * entries reference mid-face points inserted
   * to maintain a compatible tetrahedralization.
   *
   * There are 24 entries, 6 for each of the 4 faces of
   * the tetrahedron.
   */
  int SimplexTemplateRefiner::permutations_from_index[24][14] =
    {
      /* corners      midpoints          face points   */
      /* POSITIVE ARRANGEMENTS                         */
      { 0, 1, 2, 3,   4, 5, 6, 7, 8, 9,  10, 11, 12, 13 }, /* Face 0-1-2 */
      { 1, 2, 0, 3,   5, 6, 4, 8, 9, 7,  10, 12, 13, 11 },
      { 2, 0, 1, 3,   6, 4, 5, 9, 7, 8,  10, 13, 11, 12 },

      { 0, 3, 1, 2,   7, 8, 4, 6, 9, 5,  11, 13, 12, 10 }, /* Face 0-3-1 */
      { 3, 1, 0, 2,   8, 4, 7, 9, 5, 6,  11, 12, 10, 13 },
      { 1, 0, 3, 2,   4, 7, 8, 5, 6, 9,  11, 10, 13, 12 },

      { 1, 3, 2, 0,   8, 9, 5, 4, 7, 6,  12, 11, 13, 10 }, /* Face 1-3-2 */
      { 3, 2, 1, 0,   9, 5, 8, 7, 6, 4,  12, 13, 10, 11 },
      { 2, 1, 3, 0,   5, 8, 9, 6, 4, 7,  12, 10, 11, 13 },

      { 2, 3, 0, 1,   9, 7, 6, 5, 8, 4,  13, 12, 11, 10 }, /* Face 2-3-0 */
      { 3, 0, 2, 1,   7, 6, 9, 8, 4, 5,  13, 11, 10, 12 },
      { 0, 2, 3, 1,   6, 9, 7, 4, 5, 8,  13, 10, 12, 11 },

      /* NEGATIVE ARRANGEMENTS                         */
      { 0, 2, 1, 3,   6, 5, 4, 7, 9, 8,  10, 13, 12, 11 }, /* Face 0-1-2 */
      { 2, 1, 0, 3,   5, 4, 6, 9, 8, 7,  10, 12, 11, 13 },
      { 1, 0, 2, 3,   4, 6, 5, 8, 7, 9,  10, 11, 13, 12 },

      { 0, 1, 3, 2,   4, 8, 7, 6, 5, 9,  11, 10, 12, 13 }, /* Face 0-3-1 */
      { 1, 3, 0, 2,   8, 7, 4, 5, 9, 6,  11, 12, 13, 10 },
      { 3, 0, 1, 2,   7, 4, 8, 9, 6, 5,  11, 13, 10, 12 },

      { 1, 2, 3, 0,   5, 9, 8, 4, 6, 7,  12, 10, 13, 11 }, /* Face 1-3-2 */
      { 2, 3, 1, 0,   9, 8, 5, 6, 7, 4,  12, 13, 11, 10 },
      { 3, 1, 2, 0,   8, 5, 9, 7, 4, 6,  12, 11, 10, 13 },

      { 2, 0, 3, 1,   6, 7, 9, 5, 4, 8,  13, 10, 11, 12 }, /* Face 2-3-0 */
      { 0, 3, 2, 1,   7, 9, 6, 4, 8, 5,  13, 11, 12, 10 },
      { 3, 2, 0, 1,   9, 6, 7, 8, 5, 4,  13, 12, 10, 11 }
    };

  /*
   * Below is a list of output tetrahedra. The array is
   * generated by TessellatorGenerator.py
   * which also generates the code that references it.
   * Each set of tetrahedra begins with a single integer
   * that is the number of tetrahedra for that particular
   * case. It is followed by 5 integers for each output
   * tetrahedron; the first four numbers on each row are
   * indices of the output tetrahedron. The final number
   * is a bit vector specifying which edges of the
   * tetrahedron are internal to the parent tetrahedron
   * being decomposed.
   *
   * Multiple lists of output tetrahedra may be
   * combined to create the tessellation of a single
   * input tetrahedron.
   */

  int SimplexTemplateRefiner::templates[] = 
    {
      // case 1_0
      2,
      0,  4,  2,  3,
      4,  1,  2,  3,

      // case 2a_0
      1,
      3,  4,  5,  1,

      // case 2a, 0>1
      2,
      0,  4,  2,  3,
      4,  5,  2,  3,

      // case 2a, 0=1
      4,
      10,  3,  0,  4,
      10,  3,  4,  5,
      10,  3,  5,  2,
      10,  3,  2,  0,

      // case 2b_0
      4,
      0,  4,  9,  3,
      4,  1,  9,  3,
      0,  4,  2,  9,
      4,  1,  2,  9,

      // case 3a_0
      1,
      4,  7,  6,  0,

      // case 3a, 0>2>3<0
      3,
      1,  3,  2,  4,
      4,  6,  3,  2,
      4,  6,  7,  3,

      // case 3a, 0=2>3<0
      5,
      4,  6,  7,  3,
      10,  1,  2,  3,
      10,  2,  6,  3,
      10,  6,  4,  3,
      10,  4,  1,  3,

      // case 3a, 3>0=2<3
      5,
      1,  3,  2,  7,
      10,  1,  2,  7,
      10,  2,  6,  7,
      10,  6,  4,  7,
      10,  4,  1,  7,

      // case 3a, 0=2=3=0
      11,
      2,  6, 10, 13,
      3,  7, 13, 11,
      4,  1, 10, 11,
      11,  6, 10,  4,
      11,  6, 13, 10,
      11,  6,  7, 13,
      11,  6,  4,  7,
      2, 10, 11, 13,
      1, 10, 11,  2,
      2, 11,  3, 13,
      3,  2,  1, 11,

      // case 3b_0
      4,
      0,  7,  4,  2,
      4,  7,  8,  2,
      4,  8,  1,  2,
      7,  3,  8,  2,

      // case 3c, 0>1,0>3
      5,
      4,  2,  7,  5,
      4,  2,  0,  7,
      4,  3,  1,  5,
      4,  3,  5,  7,
      3,  5,  7,  2,

      // case 3c, 1>0,3>0
      5,
      0,  5,  2,  7,
      0,  5,  7,  4,
      7,  1,  4,  5,
      7,  1,  5,  3,
      3,  5,  7,  2,

      // case 3c, 0>1,3>0
      5,
      4,  2,  7,  5,
      4,  2,  0,  7,
      7,  1,  4,  5,
      7,  1,  5,  3,
      3,  5,  7,  2,

      // case 3c, 1>0,0>3
      5,
      0,  5,  2,  7,
      0,  5,  7,  4,
      4,  3,  1,  5,
      4,  3,  5,  7,
      3,  5,  7,  2,

      // case 3c, 0=1,0>3
      7,
      4,  1,  5,  3,
      10,  0,  4,  7,
      10,  2,  0,  7,
      10,  7,  4,  3,
      10,  2,  7,  3,
      10,  5,  2,  3,
      10,  4,  5,  3,

      // case 3c, 3>0,0=1
      7,
      7,  1,  5,  3,
      7,  5,  2,  3,
      10,  0,  4,  7,
      10,  2,  0,  7,
      10,  5,  2,  7,
      10,  4,  5,  7,
      1,  5,  4,  7,

      // case 3c, 0=1,0=3
      10,
      4,  1,  5, 11,
      11,  1,  5,  3,
      10,  0,  4,  7,
      10,  2,  0,  7,
      10,  5,  2,  3,
      10,  2,  7,  3,
      10,  7,  4, 11,
      10,  7, 11,  3,
      10,  4,  5, 11,
      10, 11,  5,  3,

      // case 3d, 0>4,0>2
      5,
      4,  3,  6,  0,
      4,  3,  8,  6,
      4,  2,  8,  1,
      4,  2,  6,  8,
      2,  3,  6,  8,

      // case 3d, 4>0,2>0
      5,
      8,  0,  6,  4,
      8,  0,  3,  6,
      6,  1,  8,  4,
      6,  1,  2,  8,
      2,  3,  6,  8,

      // case 3d, 0>4,2>0
      5,
      4,  3,  6,  0,
      4,  3,  8,  6,
      6,  1,  8,  4,
      6,  1,  2,  8,
      2,  3,  6,  8,

      // case 3d, 4>0,0>2
      5,
      8,  0,  6,  4,
      8,  0,  3,  6,
      4,  2,  8,  1,
      4,  2,  6,  8,
      2,  3,  6,  8,

      // case 3d, 0=4,0>2
      7,
      4,  1,  2,  8,
      11,  4,  0,  6,
      11,  0,  3,  6,
      11,  2,  4,  6,
      11,  3,  2,  6,
      11,  3,  8,  2,
      11,  8,  4,  2,

      // case 3d, 2>0,0=4
      7,
      6,  2,  8,  1,
      6,  8,  2,  3,
      11,  4,  0,  6,
      11,  0,  3,  6,
      8, 11,  3,  6,
      8,  4, 11,  6,
      1,  6,  4,  8,

      // case 3d, 0=4,0=2
      10,
      4,  1, 10,  8,
      10,  2,  8,  1,
      11,  4,  0,  6,
      11,  0,  3,  6,
      11,  3,  8,  2,
      11,  3,  2,  6,
      11, 10,  4,  6,
      11, 10,  6,  2,
      8,  4, 11, 10,
      11, 10,  2,  8,

      // case 4a_0
      2,
      7,  8,  9,  3,
      7,  9,  8,  6,

      // case 4a, 5>4>3
      4,
      8,  0,  6,  1,
      8,  0,  7,  6,
      9,  1,  6,  2,
      9,  1,  8,  6,

      // case 4a, 3<4>5
      4,
      8,  0,  6,  1,
      8,  0,  7,  6,
      8,  2,  6,  9,
      8,  2,  1,  6,

      // case 4a, 3>4<5
      4,
      6,  9,  8,  1,
      6,  9,  1,  2,
      6,  7,  0,  1,
      6,  7,  1,  8,

      // case 4a, 3=4>5
      6,
      6,  7,  0, 11,
      6,  0,  1, 11,
      6,  7, 11,  8,
      6, 11,  1,  8,
      1,  2,  6,  8,
      2,  6,  8,  9,

      // case 4a, 5>4,3=4
      6,
      6,  7,  0, 11,
      6,  0,  1, 11,
      6,  7, 11,  8,
      6, 11,  1,  8,
      1,  2,  6,  9,
      1,  6,  8,  9,

      // case 4a, 3=4=5
      8,
      6,  7,  0, 11,
      6,  0,  1, 11,
      6,  7, 11,  8,
      6, 11,  1,  8,
      6,  1,  2, 12,
      6,  2,  9, 12,
      6,  9,  8, 12,
      6,  8,  1, 12,

      // case 4b, 2>1,3>2,4>3,4>1
      6,
      6,  8,  1,  5,
      6,  8,  0,  1,
      6,  8,  7,  0,
      6,  8,  2,  7,
      7,  8,  2,  3,
      6,  8,  5,  2,

      // case 4b, 2>1,3>2,3>4,4>1
      6,
      6,  8,  1,  5,
      6,  8,  7,  1,
      6,  7,  0,  1,
      8,  7,  3,  2,
      6,  8,  5,  2,
      6,  8,  2,  7,

      // case 4b, 2>1,3>2,3>4,4>1, a
      6,
      7,  8,  1,  5,
      6,  5,  7,  1,
      6,  7,  0,  1,
      8,  7,  3,  2,
      7,  8,  5,  2,
      6,  5,  2,  7,

      // case 4b, 2>1,2>3,4>3,4>1
      6,
      6,  8,  5,  2,
      6,  8,  2,  3,
      6,  8,  3,  7,
      6,  8,  7,  0,
      6,  8,  0,  1,
      6,  8,  1,  5,

      // case 4b, 1=2,3>2,3>4,4>1
      9,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1,  7,
      10,  7,  1,  8,
      6,  7, 10,  8,
      6, 10,  5,  8,
      6,  2,  7,  8,
      6,  5,  2,  8,
      7,  8,  2,  3,

      // case 4b, 1=2,2>3,4>3,1>4
      9,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1,  8,
      10,  7,  0,  8,
      6,  7, 10,  8,
      6, 10,  5,  8,
      6,  3,  7,  8,
      6,  5,  3,  8,
      6,  5,  2,  3,

      // case 4b, 1=2,2>3,4>3,4>1
      9,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1,  8,
      10,  7,  0,  8,
      6,  7, 10,  8,
      6, 10,  5,  8,
      6,  3,  7,  8,
      6,  5,  2,  8,
      6,  2,  3,  8,

      // case 4b, 1=2,3>2,3=4,4>1
      11,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1, 11,
      10, 11,  1,  8,
      10,  0, 11,  7,
      10,  7, 11,  8,
      6,  7, 10,  8,
      6, 10,  5,  8,
      6,  2,  7,  8,
      6,  5,  2,  8,
      7,  8,  2,  3,

      // case 4b, 4>1=2=3,4>3
      12,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1,  8,
      10,  7,  0,  8,
      13,  6,  2,  5,
      13,  3,  7,  8,
      13,  2,  3,  8,
      13,  2,  8,  5,
      6,  7, 10,  8,
      6, 10,  5,  8,
      6, 13,  7,  8,
      6,  5, 13,  8,

      // case 4b, 1=2=3>4,1>4
      12,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1,  7,
      10,  7,  1,  8,
      13,  6,  2,  5,
      13,  3,  7,  8,
      13,  2,  3,  5,
      13,  3,  8,  5,
      6,  7, 10,  8,
      6, 10,  5,  8,
      6, 13,  7,  8,
      6,  5, 13,  8,

      // case 4b, 1=2=3=4=1
      16,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1, 11,
      10, 11,  1,  8,
      10,  0, 11,  7,
      10,  7, 11,  8,
      13,  6,  2,  5,
      13,  3,  7,  8,
      13,  2,  3, 12,
      13,  2, 12,  5,
      13, 12,  3,  8,
      13, 12,  5,  8,
      6,  7, 10,  8,
      6, 10,  5,  8,
      6,  5, 13,  8,
      6, 13,  7,  8,

      // case 5_0
      2,
      7,  8,  9,  3,
      6,  5,  2,  9,

      // case 5, 1>2,3>4
      5,
      5,  7,  1,  8,
      5,  7,  0,  1,
      5,  7,  6,  0,
      5,  7,  9,  6,
      5,  7,  8,  9,

      // case 5, 1>2,4>3
      5,
      0,  5,  6,  7,
      0,  5,  7,  8,
      0,  5,  8,  1,
      5,  7,  9,  6,
      5,  7,  8,  9,

      // case 5, 1>2,4>3, a
      5,
      0,  5,  6,  8,
      0,  6,  7,  8,
      0,  5,  8,  1,
      5,  8,  9,  6,
      6,  7,  8,  9,

      // case 5, 1=2,3>4
      8,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1,  7,
      10,  7,  1,  8,
      10,  8,  5,  9,
      10,  6,  7,  9,
      10,  7,  8,  9,
      10,  5,  6,  9,

      // case 5, 1=2,3>4, a
      8,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1,  7,
      10,  7,  1,  8,
      7,  8,  5,  9,
      10,  6,  7,  5,
      10,  7,  8,  5,
      5,  9,  6,  7,

      // case 5, 1=2,3>4, b
      8,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1,  7,
      10,  7,  1,  8,
      6,  8,  5,  9,
      10,  6,  7,  8,
      10,  6,  8,  5,
      8,  9,  6,  7,

      // case 5, 1=2,3=4
      10,
      10,  6,  0,  7,
      10,  1,  5,  8,
      10,  0,  1, 11,
      10, 11,  1,  8,
      10,  0, 11,  7,
      10,  7, 11,  8,
      10,  8,  5,  9,
      10,  6,  7,  9,
      10,  7,  8,  9,
      10,  5,  6,  9,

      // case 6_0
      4,
      7,  8,  9,  3,
      6,  5,  2,  9,
      4,  1,  5,  8,
      0,  4,  6,  7,

      // case 6_1
      4,
      6,  4,  5,  8,
      6,  5,  9,  8,
      6,  9,  7,  8,
      6,  7,  4,  8,

      // case 6_1, a
      4,
      5,  8,  9,  7,
      5,  9,  6,  7,
      5,  6,  4,  7,
      5,  4,  8,  7,

      // case 6_1, b
      4,
      4,  5,  6,  9,
      4,  6,  7,  9,
      4,  7,  8,  9,
      4,  8,  5,  9,

    };

} // namespace moab 

