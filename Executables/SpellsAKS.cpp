#include "Elliptic/EllipticSolver/SolveEllipticEquations.hpp"
#include "Utils/ErrorHandling/BasicMpi.hpp"
#include "Utils/LowLevelUtils/InfoAtRuntime.hpp"
#include "Utils/IO/H5DatWriter.hpp"

int main(int argc, char **argv) {
  MpiInit(&argc,&argv);
  coutProc0 << InfoAtRuntime(argc,argv);

  SolveEllipticEquations();

  FlushAllH5DatWriters();
  MpiFinalize();
  return 0;
}
