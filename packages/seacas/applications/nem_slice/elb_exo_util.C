/*
 * Copyright(C) 1999-2024 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include <exodusII.h> // for ex_close, ex_inquire, etc

#include <algorithm>
#include <cassert>
#include <copy_string_cpp.h>
#include <cstddef> // for size_t
#include <cstdio>  // for nullptr
#include <cstdlib> // for malloc, free, calloc
#include <cstring> // for strlen
#include <fmt/ostream.h>
#include <string>
#include <vector>

#include "elb.h"      // for Weight_Description, etc
#include "elb_elem.h" // for get_elem_type, E_Type, etc
#include "elb_err.h"  // for Gen_Error, MAX_ERR_MSG
#include "elb_exo.h"
#include "elb_groups.h" // for parse_groups
#include "elb_util.h"   // for in_list, roundfloat
#include "vector_data.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function read_exo_weights() begins:
 *----------------------------------------------------------------------------
 * This function reads the nodal or elemental values from an ExodusII file
 * which will be used by Chaco for weighting of the graph.
 *****************************************************************************/
template int read_exo_weights(Problem_Description *prob, Weight_Description *weight, int dummy);
template int read_exo_weights(Problem_Description *prob, Weight_Description *weight, int64_t dummy);

template <typename INT>
int read_exo_weights(Problem_Description *prob, Weight_Description *weight, INT /*dummy*/)
{
  int exoid;
  /*---------------------------Execution Begins--------------------------------*/

  /* Open the ExodusII file containing the weights */
  int   mode   = EX_READ | prob->int64api;
  int   cpu_ws = 0;
  int   io_ws  = 0;
  float version;
  if ((exoid = ex_open(weight->exo_filename.c_str(), mode, &cpu_ws, &io_ws, &version)) < 0) {
    std::string ctemp = fmt::format("fatal: could not open ExodusII file {}", weight->exo_filename);
    Gen_Error(0, ctemp);
    return 0;
  }

  std::vector<float> values(weight->nvals);
  if (prob->type == NODAL) {
    size_t tmp_nodes = ex_inquire_int(exoid, EX_INQ_NODES);
    /* check to make sure the sizes agree */
    if ((size_t)weight->nvals != tmp_nodes) {
      Gen_Error(0, "fatal: different number of nodes in mesh and weight files");
      ex_close(exoid);
      return 0;
    }

    weight->ow.resize(weight->nvals);
    /* Read in the nodal values */
    if (ex_get_var(exoid, weight->exo_tindx, EX_NODAL, weight->exo_vindx, 1, weight->nvals,
                   Data(values)) < 0) {
      Gen_Error(0, "fatal: unable to read nodal values");
      ex_close(exoid);
      return 0;
    }
  }
  else {
    size_t tmp_elem = ex_inquire_int(exoid, EX_INQ_ELEM);
    /* check to make sure the sizes agree */
    if ((size_t)weight->nvals != tmp_elem) {
      Gen_Error(0, "fatal: different number of elems in mesh and weight files");
      ex_close(exoid);
      return 0;
    }

    /* Get the number of element blocks */
    int              neblks = ex_inquire_int(exoid, EX_INQ_ELEM_BLK);
    std::vector<INT> eblk_ids(neblks);
    std::vector<INT> eblk_ecnts(neblks);

    if (ex_get_ids(exoid, EX_ELEM_BLOCK, Data(eblk_ids)) < 0) {
      Gen_Error(0, "fatal: unable to get element block IDs");
      ex_close(exoid);
      return 0;
    }

    /* Get the count of elements in each element block */
    for (int cnt = 0; cnt < neblks; cnt++) {
      INT  dum1;
      INT  dum2;
      char elem_type[MAX_STR_LENGTH + 1];
      if (ex_get_block(exoid, EX_ELEM_BLOCK, eblk_ids[cnt], elem_type, &(eblk_ecnts[cnt]), &dum1,
                       nullptr, nullptr, &dum2) < 0) {
        Gen_Error(0, "fatal: unable to get element block");
        ex_close(exoid);
        return 0;
      }
    }

    /* Get the element variables */
    size_t offset = 0;
    for (int cnt = 0; cnt < neblks; cnt++) {
      if (ex_get_var(exoid, weight->exo_tindx, EX_ELEM_BLOCK, weight->exo_vindx, eblk_ids[cnt],
                     eblk_ecnts[cnt], &(values[offset])) < 0) {
        Gen_Error(0, "fatal: unable to get element variable");
        ex_close(exoid);
        return 0;
      }
      offset += eblk_ecnts[cnt];
    }
  }

  /* Close the ExodusII weighting file */
  if (ex_close(exoid) < 0) {
    std::string ctemp =
        fmt::format("warning: failed to close ExodusII file {}", weight->exo_filename);
    Gen_Error(0, ctemp);
  }

  /* now I need to translate the values to positive integers */

  /* first find the minimum value */
  float minval = *std::min_element(values.begin(), values.end());

  /* now translate the values to be greater than 1 and convert to ints */
  for (int cnt = 0; cnt < weight->nvals; cnt++) {
    values[cnt] += 1.0f - minval;
    weight->vertices[cnt] = roundfloat(values[cnt]);
  }
  return 1;
} /*------------------------End read_exo_weights()----------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function read_mesh_params() begins:
 *----------------------------------------------------------------------------
 * This function reads in information about the finite element mesh.
 *****************************************************************************/
template int read_mesh_params(const std::string &exo_file, Problem_Description *problem,
                              Mesh_Description<int> *mesh, Sphere_Info *sphere);
template int read_mesh_params(const std::string &exo_file, Problem_Description *problem,
                              Mesh_Description<int64_t> *mesh, Sphere_Info *sphere);

template <typename INT>
int read_mesh_params(const std::string &exo_file, Problem_Description *problem,
                     Mesh_Description<INT> *mesh, Sphere_Info *sphere)
{
  int   exoid;
  int   cpu_ws = 0;
  int   io_ws  = 0;
  float version;
  char  elem_type[MAX_STR_LENGTH + 1];
  /*---------------------------Execution Begins--------------------------------*/

  /* Open the ExodusII geometry file */
  int mode = EX_READ | problem->int64api;
  if ((exoid = ex_open(exo_file.c_str(), mode, &cpu_ws, &io_ws, &version)) < 0) {
    Gen_Error(0, "fatal: unable to open ExodusII file for mesh params");
    return 0;
  }

  /* Get the init info */
  ex_init_params exo{};
  if (ex_get_init_ext(exoid, &exo)) {
    Gen_Error(0, "fatal: unable to get init info from ExodusII file");
    ex_close(exoid);
    return 0;
  }
  copy_string(mesh->title, exo.title);
  mesh->num_dims      = exo.num_dim;
  mesh->num_nodes     = exo.num_nodes;
  mesh->num_elems     = exo.num_elem;
  mesh->num_el_blks   = exo.num_elem_blk;
  mesh->num_node_sets = exo.num_node_sets;
  mesh->num_side_sets = exo.num_side_sets;

  /* Read the element block IDs */
  mesh->eb_ids.resize(mesh->num_el_blks);
  mesh->eb_cnts.resize(mesh->num_el_blks);
  mesh->eb_npe.resize(mesh->num_el_blks);
  mesh->eb_type.resize(mesh->num_el_blks);

  if (ex_get_ids(exoid, EX_ELEM_BLOCK, Data(mesh->eb_ids)) < 0) {
    Gen_Error(0, "fatal: unable to get element block IDs");
    ex_close(exoid);
    return 0;
  }

  /* Get the length of the concatenated node set node list */
  if (mesh->num_node_sets > 0) {
    mesh->ns_list_len = ex_inquire_int(exoid, EX_INQ_NS_NODE_LEN);
  }
  else {
    mesh->ns_list_len = 0;
  }

  /* Allocate and initialize memory for the sphere adjustment */
  sphere->adjust.resize(mesh->num_el_blks);
  sphere->begin.resize(mesh->num_el_blks);
  sphere->end.resize(mesh->num_el_blks);

  /* Determine the maximum number of nodes per element */
  mesh->max_np_elem = 0;
  for (size_t cnt = 0; cnt < mesh->num_el_blks; cnt++) {
    INT nodes_in_elem;

    if (ex_get_block(exoid, EX_ELEM_BLOCK, mesh->eb_ids[cnt], elem_type, &(mesh->eb_cnts[cnt]),
                     &nodes_in_elem, nullptr, nullptr, nullptr) < 0) {
      Gen_Error(0, "fatal: unable to get element block");
      ex_close(exoid);
      return 0;
    }

    if (mesh->eb_cnts[cnt] == 0) {
      continue;
    }

    mesh->eb_npe[cnt]  = nodes_in_elem;
    mesh->eb_type[cnt] = get_elem_type(elem_type, nodes_in_elem, mesh->num_dims);

    if (cnt == 0) {
      sphere->end[0] = mesh->eb_cnts[cnt];
    }

    if (mesh->eb_type[cnt] == SPHERE && problem->no_sph != 1) {
      sphere->num += mesh->eb_cnts[cnt];
      sphere->adjust[cnt] = 0;
    }
    else {
      sphere->adjust[cnt] = sphere->num;
    }

    if (cnt != 0) {
      sphere->begin[cnt] = sphere->end[cnt - 1];
      sphere->end[cnt]   = sphere->begin[cnt] + mesh->eb_cnts[cnt];
    }

    mesh->max_np_elem = MAX(mesh->max_np_elem, (size_t)nodes_in_elem);
  }

  /* Close the ExodusII file */
  if (ex_close(exoid) < 0) {
    Gen_Error(1, "warning: unable to close ExodusII file");
  }

  fmt::print("ExodusII mesh information\n");
  if (strlen(mesh->title) > 0) {
    fmt::print("\ttitle: {}\n", mesh->title);
  }
  fmt::print("\tgeometry dimension: {}\n", mesh->num_dims);
  fmt::print("\tnumber of nodes: {}\tnumber of elements: {}\n", fmt::group_digits(mesh->num_nodes),
             fmt::group_digits(mesh->num_elems));
  fmt::print("\tnumber of element blocks: {}\n", mesh->num_el_blks);
  fmt::print("\tnumber of node sets: {}\tnumber of side sets: {}\n", mesh->num_node_sets,
             mesh->num_side_sets);

  return 1;

} /*--------------------------End read_mesh_params()-------------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function read_mesh_params() begins:
 *----------------------------------------------------------------------------
 * This function reads in the finite element mesh.
 *****************************************************************************/
template int read_mesh(const std::string &exo_file, Problem_Description *problem,
                       Mesh_Description<int> *mesh, Weight_Description *weight);
template int read_mesh(const std::string &exo_file, Problem_Description *problem,
                       Mesh_Description<int64_t> *mesh, Weight_Description *weight);

template <typename INT>
int read_mesh(const std::string &exo_file, Problem_Description *problem,
              Mesh_Description<INT> *mesh, Weight_Description *weight)
{
  float  version;
  float *xptr;
  float *yptr;
  float *zptr;
  /*---------------------------Execution Begins--------------------------------*/

  /* Open the ExodusII file */
  int exoid;
  int cpu_ws = 0;
  int io_ws  = 0;
  int mode   = EX_READ | problem->int64api;
  if ((exoid = ex_open(exo_file.c_str(), mode, &cpu_ws, &io_ws, &version)) < 0) {
    Gen_Error(0, "fatal: unable to open ExodusII mesh file");
    return 0;
  }

  /* Read the coordinates, if desired */
  xptr = yptr = zptr = nullptr;

  if (problem->read_coords == ELB_TRUE) {
    switch (mesh->num_dims) {
    case 3: zptr = Data(mesh->coords) + 2 * (mesh->num_nodes); FALL_THROUGH;
    case 2: yptr = Data(mesh->coords) + (mesh->num_nodes); FALL_THROUGH;
    case 1: xptr = Data(mesh->coords);
    }

    if (ex_get_coord(exoid, xptr, yptr, zptr) < 0) {
      Gen_Error(0, "fatal: unable to read coordinate values for mesh");
      return 0;
    }

  } /* End "if(problem->read_coords == ELB_TRUE)" */

  /* Read the element connectivity */
  size_t gelem_cnt = 0;
  for (size_t cnt = 0; cnt < mesh->num_el_blks; cnt++) {
    if (mesh->eb_cnts[cnt] == 0) {
      continue;
    }

    std::vector<INT> blk_connect(mesh->eb_cnts[cnt] * mesh->eb_npe[cnt]);

    /* Get the connectivity for this element block */
    if (ex_get_conn(exoid, EX_ELEM_BLOCK, mesh->eb_ids[cnt], Data(blk_connect), nullptr, nullptr) <
        0) {
      Gen_Error(0, "fatal: failed to get element connectivity");
      return 0;
    }

    /* find out if this element block is weighted */
    int wgt = -1;
    if (weight->type & EL_BLK) {
      wgt = in_list(mesh->eb_ids[cnt], weight->elemblk);
    }

    /* Fill the 2D global connectivity array */
    if (((problem->type == ELEMENTAL) && (weight->type & EL_BLK)) ||
        ((problem->type == NODAL) && (weight->type & EL_BLK))) {

      for (int64_t cnt2 = 0; cnt2 < mesh->eb_cnts[cnt]; cnt2++) {
        mesh->elem_type[gelem_cnt] = mesh->eb_type[cnt];

        /* while going through the blocks, take care of the weighting */
        if ((problem->type == ELEMENTAL) && (weight->type & EL_BLK)) {
          /* is this block weighted */
          if (wgt >= 0) {
            /* check if there is a read value */
            if (weight->vertices[gelem_cnt] >= 1) {
              /* and if it should be overwritten */
              if (weight->ow_read) {
                weight->vertices[gelem_cnt] = weight->elemblk_wgt[wgt];
              }
            }
            else {
              weight->vertices[gelem_cnt] = weight->elemblk_wgt[wgt];
            }
          }
          else {
            /* now check if this weight has been initialized */
            if (weight->vertices[gelem_cnt] < 1) {
              weight->vertices[gelem_cnt] = 1;
            }
          }
        }

        for (int64_t cnt3 = 0; cnt3 < mesh->eb_npe[cnt]; cnt3++) {
          INT node = blk_connect[cnt3 + cnt2 * mesh->eb_npe[cnt]] - 1;
          assert(node >= 0);
          mesh->connect[gelem_cnt][cnt3] = node;

          /* deal with the weighting if necessary */
          if ((problem->type == NODAL) && (weight->type & EL_BLK)) {
            /* is this block weighted */
            if (wgt >= 0) {
              /* check if I read an exodus file */
              if (weight->type & READ_EXO) {
                /* check if it can be overwritten */
                if (weight->ow_read) {
                  /* check if it has been overwritten already */
                  if (weight->ow[node]) {
                    weight->vertices[node] = MAX(weight->vertices[node], weight->elemblk_wgt[wgt]);
                  }
                  else {
                    weight->vertices[node] = weight->elemblk_wgt[wgt];
                    weight->ow[node]       = 1; /* read value has been overwritten */
                  }
                }
              }
              else {
                weight->vertices[node] = MAX(weight->vertices[node], weight->elemblk_wgt[wgt]);
              }
            }
            else {
              /* now check if this weight has been initialized */
              if (weight->vertices[node] < 1) {
                weight->vertices[node] = 1;
              }
            }
          }
        }
        gelem_cnt++;
      }
    }
    else {
      // No weights...
      for (int64_t cnt2 = 0; cnt2 < mesh->eb_cnts[cnt]; cnt2++) {
        mesh->elem_type[gelem_cnt] = mesh->eb_type[cnt];

        for (int64_t cnt3 = 0; cnt3 < mesh->eb_npe[cnt]; cnt3++) {
          INT node = blk_connect[cnt2 * mesh->eb_npe[cnt] + cnt3] - 1;
          assert(node >= 0);
          mesh->connect[gelem_cnt][cnt3] = node;
        }
        gelem_cnt++;
      }
    }
  } /* End "for(cnt=0; cnt < mesh->num_el_blks; cnt++)" */

  /* if there is a group designator, then parse it here */
  if (problem->groups != nullptr) {
    if (!parse_groups(mesh, problem)) {
      Gen_Error(0, "fatal: unable to parse group designator");
      ex_close(exoid);
      return 0;
    }
  }
  else {
    problem->num_groups = 1; /* there is always one group */
  }

  /* Close the ExodusII file */
  if (ex_close(exoid) < 0) {
    Gen_Error(0, "warning: failed to close ExodusII mesh file");
  }

  return 1;

} /*---------------------------End read_mesh()-------------------------------*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Function init_weight_struct() begins:
 *----------------------------------------------------------------------------
 * This function initializes the weight structure given the current mesh.
 *****************************************************************************/
template int init_weight_struct(Problem_Description *problem, Mesh_Description<int> *mesh,
                                Weight_Description *weight);
template int init_weight_struct(Problem_Description *problem, Mesh_Description<int64_t> *mesh,
                                Weight_Description *weight);

template <typename INT>
int init_weight_struct(Problem_Description *problem, Mesh_Description<INT> *mesh,
                       Weight_Description *weight)
{
  if (problem->type == NODAL) {
    weight->nvals = mesh->num_nodes;
  }
  else {
    weight->nvals = mesh->num_elems;
  }

  /* Allocate memory */
  weight->vertices.resize(weight->nvals);
  return 1;
} /*-----------------------End init_weight_struct()--------------------------*/
