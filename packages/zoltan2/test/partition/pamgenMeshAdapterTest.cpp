/**********************************************************************************/
/******************************** MAIN ********************************************/
/**********************************************************************************/

int main(int argc, char *argv[]) {

  // Generate mesh with Pamgen
  long long maxInt = 9223372036854775807LL;
  Create_Pamgen_Mesh(meshInput.c_str(), dim, rank, numProcs, maxInt);

  // delete mesh
  Delete_Pamgen_Mesh();

  return 0;
}
