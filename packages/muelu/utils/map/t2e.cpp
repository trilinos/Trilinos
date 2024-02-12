#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
int main(int argc, char* argv[]) {
  if (argc < 3) {
    cout << "Usage: ./a.out <in> <out>" << endl;
    return 1;
  }
  ifstream ifs(argv[1]);
  ofstream ofs(argv[2]);

  std::cout << "Reading Tpetra map \"" << argv[1] << "\"" << std::endl;
  // Skip % MatrixMarket header
  char line[256];
  while (true) {
    ifs.getline(line, 256);
    if (line[0] != '%')
      break;
  }

  int n;
  sscanf(line, "%d", &n);
  std::cout << "n = " << n << std::endl;

  std::vector<std::vector<int> > gids;
  for (int i = 0; i < n; i += 2) {
    int id, procid;
    ifs >> id >> procid;

    if (procid >= gids.size())
      gids.resize(procid + 1);
    gids[procid].push_back(id);
  }

  ofs << "%%MatrixMarket matrix array integer general" << std::endl;
  ofs << "%Format Version:" << std::endl;
  ofs << "% 2 " << std::endl;
  ofs << "%NumProc: Number of processors:" << std::endl;
  ofs << "% " << gids.size() << " " << std::endl;
  ofs << "%MaxElementSize: Maximum element size:" << std::endl;
  ofs << "% 1 " << std::endl;
  ofs << "%MinElementSize: Minimum element size:" << std::endl;
  ofs << "% 1 " << std::endl;
  ofs << "%IndexBase: Index base of map:" << std::endl;
  ofs << "% 0 " << std::endl;
  ofs << "%NumGlobalElements: Total number of GIDs in map:" << std::endl;
  ofs << "% " << n / 2 << " " << std::endl;
  ofs << "%NumMyElements: BlockMap lengths per processor:" << std::endl;
  for (int k = 0; k < gids.size(); k++)
    ofs << "% " << gids[k].size() << std::endl;
  ofs << n / 2 << " 1" << std::endl;

  for (int k = 0; k < gids.size(); k++)
    for (int i = 0; i < gids[k].size(); i++)
      ofs << gids[k][i] << std::endl;

  return 0;
}
