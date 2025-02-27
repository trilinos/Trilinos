# Allows panzer to read an exodus mesh that has not been sliced. It
# will be sliced at runtime. Source this to avoid needin nem_slice.
export IOSS_PROPERTIES="DECOMPOSITION_METHOD=rib"
