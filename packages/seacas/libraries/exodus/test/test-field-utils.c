#include <exodusII.h>
#include <stdio.h>
#undef NDEBUG
#include <assert.h>

#define SIZE(X) sizeof(X) / sizeof(X[0])

int main(int argc, char **argv)
{
  ex_field_type types[] = {EX_FIELD_TYPE_USER_DEFINED,
                           EX_FIELD_TYPE_SEQUENCE,
                           EX_QUADRATURE,
                           EX_BASIS,
                           EX_SCALAR,
                           EX_VECTOR_1D,
                           EX_VECTOR_2D,
                           EX_VECTOR_3D,
                           EX_QUATERNION_2D,
                           EX_QUATERNION_3D,
                           EX_FULL_TENSOR_12,
                           EX_FULL_TENSOR_16,
                           EX_FULL_TENSOR_22,
                           EX_FULL_TENSOR_32,
                           EX_FULL_TENSOR_36,
                           EX_SYM_TENSOR_10,
                           EX_SYM_TENSOR_11,
                           EX_SYM_TENSOR_13,
                           EX_SYM_TENSOR_21,
                           EX_SYM_TENSOR_31,
                           EX_SYM_TENSOR_33,
                           EX_ASYM_TENSOR_01,
                           EX_ASYM_TENSOR_02,
                           EX_ASYM_TENSOR_03,
                           EX_MATRIX_2X2,
                           EX_MATRIX_3X3,
                           EX_FIELD_TYPE_INVALID};

  int num_enum = SIZE(types);
  for (int i = 0; i < num_enum; i++) {
    const char   *string_type = ex_field_type_enum_to_string(types[i]);
    ex_field_type return_type = ex_string_to_field_type_enum(string_type);
    assert(return_type == types[i]);
  }

  for (int i = 0; i < num_enum; i++) {
    struct ex_field field = (ex_field){.entity_type            = EX_ELEM_BLOCK,
                                       .entity_id              = 100,
                                       .name                   = "Testing",
                                       .type                   = {types[i]},
                                       .nesting                = 1,
                                       .component_separator[0] = '_'};

    int cardinality = ex_field_cardinality(types[i]);
    // Make sure the numeric suffices are formatted correctly "00", "01", ... "10" -- all same
    // width.
    if (cardinality == -1 && (types[i] == EX_FIELD_TYPE_SEQUENCE || types[i] == EX_BASIS)) {
      field.cardinality[0] = 10;
      cardinality          = 10;
    }
    if (cardinality > 0 && types[i] != EX_SCALAR) {
      fprintf(stderr, "%s (%s) components: ", ex_field_type_enum_to_string(field.type[0]),
              ex_field_type_name(field.type[0]));
      for (int j = 1; j <= cardinality; j++) {
        fprintf(stderr, "%s%c", ex_field_component_suffix(&field, 0, j),
                j == cardinality ? ' ' : ',');
      }
      fprintf(stderr, "\n");
    }
    else {
      fprintf(stderr, "%s (%s)\n", ex_field_type_enum_to_string(field.type[0]),
              ex_field_type_name(field.type[0]));
    }
  }
}
