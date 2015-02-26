#include <iostream>
#include "capd/capdlib.h"
//#include "capd/dynsys/DiscreteDynSys.h"
#include "time.h"

using std::cout;

using namespace capd;
using namespace matrixAlgorithms;
using namespace dynsys;

const double accuracy = 1e-12;           // accuracy for nonrigorous numerics (i.e. approximation of the slow manifold)

#include "polybd.hpp" 
#include "proof.hpp"

// ---------------------------------------------------------------------------------
// ----------------------------------- MAIN ----------------------------------------
// ---------------------------------------------------------------------------------




int main(){

  time_t start1,end1;

  time (&start1);
  cout.precision(15);

  time (&end1);
  return 0;
} 

