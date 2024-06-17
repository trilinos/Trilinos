/*
 * Copyright(C) 1999-2024 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#undef NDEBUG
#include "exodusII.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define STRINGIFY(x) #x
#define TOSTRING(x)  STRINGIFY(x)

static char *my_strdup(const char *s)
{
  size_t slen   = strlen(s);
  char  *result = malloc(slen + 1);
  if (result == NULL) {
    return NULL;
  }

  memcpy(result, s, slen + 1);
  return result;
}

static char *my_strsep(char **stringp, const char *delim)
{
  assert(delim != NULL);
  char *rv = *stringp;
  if (rv) {
    assert(stringp != NULL);
    assert(*stringp != NULL);
    *stringp += strcspn(*stringp, delim);
    if (**stringp)
      *(*stringp)++ = '\0';
    else
      *stringp = NULL;
  }
  return rv;
}

#define EXCHECK(funcall)                                                                           \
  do {                                                                                             \
    int error = (funcall);                                                                         \
    printf("after %s, error = %d\n", TOSTRING(funcall), error);                                    \
    if (error != EX_NOERR) {                                                                       \
      fprintf(stderr, "Error calling %s\n", TOSTRING(funcall));                                    \
      ex_close(exoid);                                                                             \
      exit(-1);                                                                                    \
    }                                                                                              \
  } while (0)

static char *get_type_name(char *type_name, size_t which)
{
  if (type_name != NULL && type_name[0] != '\0') {
    char *string = my_strdup(type_name);
    char *tofree = string;
    char *token  = my_strsep(&string, ",");
    for (int i = 0; i < which; i++) {
      token = my_strsep(&string, ",");
    }
    if (token != NULL) {
      static char tmp_type_name[256 + 1];
      ex_copy_string(tmp_type_name, token, 256);
      free(tofree);
      return tmp_type_name;
    }
    free(tofree);
  }
  return NULL;
}

static void get_field_cardinality(ex_field *field, ex_basis *basis, int bas_cnt,
                                  ex_quadrature *quad, int quad_cnt)
{
  for (int j = 0; j < field->nesting; j++) {
    if (field->cardinality[j] == 0) {
      if (field->type[j] == EX_BASIS) {
        int   found      = 0;
        char *basis_type = get_type_name(field->type_name, j);
        if (basis != NULL) {
          for (int k = 0; k < bas_cnt; k++) {
            if (strcmp(basis[k].name, basis_type) == 0) {
              field->cardinality[j] = basis[k].cardinality;
              found                 = 1;
              break;
            }
          }
        }
        if (found == 0) {
          fprintf(stderr, "ERROR: Could not find basis named `%s` for field `%s`\n", basis_type,
                  field->name);
        }
      }
      else if (field->type[j] == EX_QUADRATURE) {
        int   found     = 0;
        char *quad_type = get_type_name(field->type_name, j);
        if (quad != NULL) {
          for (int k = 0; k < quad_cnt; k++) {
            if (strcmp(quad[k].name, quad_type) == 0) {
              field->cardinality[j] = quad[k].cardinality;
              found                 = 1;
              break;
            }
          }
        }
        if (found == 0) {
          fprintf(stderr, "ERROR: Could not find quadrature named `%s` for field `%s`\n", quad_type,
                  field->name);
        }
      }
      else {
        field->cardinality[j] = ex_field_cardinality(field->type[j]);
      }
    }
  }
}

static void print_basis_metadata(ex_basis *basis, size_t num_basis)
{
  for (size_t j = 0; j < num_basis; j++) {
    printf("\nBasis Metadata: Name: `%s`, Cardinality: %d\n", basis[j].name, basis[j].cardinality);
    printf("ordinal,\t subc:  _dim\t_ordinal\t_dof_ordinal\t_num_dof\t xi      eta     zeta\n");
    for (int i = 0; i < basis[j].cardinality; i++) {
      double xi   = basis[j].xi != NULL ? basis[j].xi[i] : 0.0;
      double eta  = basis[j].eta != NULL ? basis[j].eta[i] : 0.0;
      double zeta = basis[j].zeta != NULL ? basis[j].zeta[i] : 0.0;
      printf("%8d\t%8d\t%8d\t%8d\t%8d\t%6.3f\t%6.3f\t%6.3f\n", i, basis[j].subc_dim[i],
             basis[j].subc_ordinal[i], basis[j].subc_dof_ordinal[i], basis[j].subc_num_dof[i], xi,
             eta, zeta);
    }
  }
}

static void print_quad_metadata(ex_quadrature *quad, size_t num_quad)
{
  for (size_t j = 0; j < num_quad; j++) {
    printf("\nQuadrature Metadata: Name: `%s`, Cardinality: %d\n", quad[j].name,
           quad[j].cardinality);
    printf("ordinal,\t  xi      eta     zeta    weight\n");
    for (int i = 0; i < quad[j].cardinality; i++) {
      double xi   = quad[j].xi != NULL ? quad[j].xi[i] : 0.0;
      double eta  = quad[j].eta != NULL ? quad[j].eta[i] : 0.0;
      double zeta = quad[j].zeta != NULL ? quad[j].zeta[i] : 0.0;
      double wgt  = quad[j].weight != NULL ? quad[j].weight[i] : 1.0;
      printf("%8d\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n", i, xi, eta, zeta, wgt);
    }
  }
}

static void print_field_metadata(ex_field *field)
{
  printf("\n");
  printf("Field Metadata: Name: `%s`, Nesting: %d\n", field->name, field->nesting);
  for (int j = 0; j < field->nesting; j++) {
    char sep = field->component_separator[j] == 0 ? ' ' : field->component_separator[j];
    if (field->type[j] == EX_BASIS) {
      char *basis_type = get_type_name(field->type_name, j);
      printf("\tNesting level: %d, Type: %s (%s), Cardinality: %d, Separator: \"%c\"\n", j,
             ex_field_type_enum_to_string(field->type[j]), basis_type, field->cardinality[j], sep);
    }
    else if (field->type[j] == EX_QUADRATURE) {
      char *quad_type = get_type_name(field->type_name, j);
      printf("\tNesting level: %d, Type: %s (%s), Cardinality: %d, Separator: \"%c\"\n", j,
             ex_field_type_enum_to_string(field->type[j]), quad_type, field->cardinality[j], sep);
    }
    else {
      printf("\tNesting level: %d, Type: %s, Cardinality: %d, Separator: \"%c\"\n", j,
             ex_field_type_enum_to_string(field->type[j]), field->cardinality[j], sep);
    }
    if (field->type[0] == EX_FIELD_TYPE_USER_DEFINED) {
      printf("\tUser-defined suffices: %s\n", field->suffices);
    }
  }
}

static void print_full_field_names(ex_field *field)
{
  int component[EX_MAX_FIELD_NESTING];
  if (field->nesting == 1) {
    for (int jj = 1; jj <= field->cardinality[0]; jj++) {
      component[0]     = jj;
      const char *name = ex_component_field_name(field, component);
      printf("\t\tComponent %d, Full name = %s\n", jj, name);
    }
  }
  else if (field->nesting == 2) {
    for (int kk = 1; kk <= field->cardinality[1]; kk++) {
      for (int jj = 1; jj <= field->cardinality[0]; jj++) {
        component[0]     = jj;
        component[1]     = kk;
        const char *name = ex_component_field_name(field, component);
        printf("\t\tComponent %d %d, Full name = %s\n", jj, kk, name);
      }
    }
  }
#if EX_MAX_FIELD_NESTING > 2
  else if (field->nesting == 3) {
    for (int ii = 1; ii <= field->cardinality[2]; ii++) {
      for (int kk = 1; kk <= field->cardinality[1]; kk++) {
        for (int jj = 1; jj <= field->cardinality[0]; jj++) {
          component[0]     = jj;
          component[1]     = kk;
          component[2]     = ii;
          const char *name = ex_component_field_name(field, component);
          printf("\t\tComponent %d %d %d, Full name = %s\n", jj, kk, ii, name);
        }
      }
    }
  }
#endif
#if EX_MAX_FIELD_NESTING > 3
  else if (field->nesting == 4) {
    for (int mm = 1; mm <= field->cardinality[3]; mm++) {
      for (int ii = 1; ii <= field->cardinality[2]; ii++) {
        for (int kk = 1; kk <= field->cardinality[1]; kk++) {
          for (int jj = 1; jj <= field->cardinality[0]; jj++) {
            component[0]     = jj;
            component[1]     = kk;
            component[2]     = ii;
            component[3]     = mm;
            const char *name = ex_component_field_name(field, component);
            printf("\t\tComponent %d %d %d %d, Full name = %s\n", jj, kk, ii, mm, name);
          }
        }
      }
    }
  }
#endif
}

int main(int argc, char **argv)
{
  ex_opts(EX_VERBOSE | EX_ABORT);

  /* open EXODUS II files */
  float version;
  int   CPU_word_size = 0; /* sizeof(float) */
  int   IO_word_size  = 0; /* use what is stored in file */

  int exoid = ex_open("test-field-metadata.exo", EX_READ, &CPU_word_size, &IO_word_size, &version);

  printf("\nafter ex_open\n");
  if (exoid < 0) {
    exit(1);
  }

  // ------------------------------------------------------------------------
  int            quad_cnt = 0;
  ex_quadrature *quad     = NULL;
  EXCHECK(ex_get_quadrature(exoid, &quad, &quad_cnt));
  print_quad_metadata(quad, quad_cnt);

  // ------------------------------------------------------------------------
  int       bas_cnt = 0;
  ex_basis *basis   = NULL;
  EXCHECK(ex_get_basis(exoid, &basis, &bas_cnt));
  print_basis_metadata(basis, bas_cnt);

  /* Check for nodal fields... */
  {
    int fld_cnt = ex_get_field_metadata_count(exoid, EX_NODAL, 0);
    assert(fld_cnt == 2);
    ex_field fields[2] = {{.entity_type = EX_NODAL}, {.entity_type = EX_NODAL}};
    EXCHECK(ex_get_field_metadata(exoid, fields));

    for (int i = 0; i < fld_cnt; i++) {
      get_field_cardinality(&fields[i], NULL, 0, NULL, 0);
      print_field_metadata(&fields[i]);
      print_full_field_names(&fields[i]);
    }
  }

  {
    int fld_cnt = ex_get_field_metadata_count(exoid, EX_ELEM_BLOCK, 10);
    assert(fld_cnt == 2);
    ex_field fields[2] = {{.entity_id = 10, .entity_type = EX_ELEM_BLOCK},
                          {.entity_id = 10, .entity_type = EX_ELEM_BLOCK}};
    EXCHECK(ex_get_field_metadata(exoid, fields));

    for (int i = 0; i < fld_cnt; i++) {
      get_field_cardinality(&fields[i], NULL, 0, NULL, 0);
      print_field_metadata(&fields[i]);
      print_full_field_names(&fields[i]);
    }
  }
  {
    int fld_cnt = ex_get_field_metadata_count(exoid, EX_ELEM_BLOCK, 11);
    assert(fld_cnt == 3);
    ex_field fields[3] = {{.entity_id = 11, .entity_type = EX_ELEM_BLOCK},
                          {.entity_id = 11, .entity_type = EX_ELEM_BLOCK},
                          {.entity_id = 11, .entity_type = EX_ELEM_BLOCK}};
    EXCHECK(ex_get_field_metadata(exoid, fields));

    for (int i = 0; i < fld_cnt; i++) {
      get_field_cardinality(&fields[i], basis, bas_cnt, quad, quad_cnt);
      print_field_metadata(&fields[i]);
      print_full_field_names(&fields[i]);
    }
  }

  // Now, deallocate any memory allocated on the `basis` struct.
  EXCHECK(ex_initialize_basis_struct(basis, bas_cnt, -1));
  EXCHECK(ex_initialize_quadrature_struct(quad, quad_cnt, -1));
  free(basis);
  free(quad);

  int fld_cnt = ex_get_field_metadata_count(exoid, EX_ELEM_BLOCK, 12);
  assert(fld_cnt == 0);

  int error = ex_close(exoid);
  printf("\nafter ex_close, error = %3d\n", error);
  return 0;
}
