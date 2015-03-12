#include <iostream>
#include "capd/capdlib.h"
#include "time.h"

using std::cout;
using namespace capd;

#include "polybd.hpp" 
#include "proof.hpp"



// ---------------------------------------------------------------------------------
// ----------------------------------- MAIN ----------------------------------------
// ---------------------------------------------------------------------------------




int main(){

  /*
  time_t start1,end1;
  time (&start1);
  */

  cout.precision(15);
  bool verbose(0);
  
  interval sigma("3.","3.");

  interval beta1("1.5","1.5");
  interval eps1("-0.05","0.05");
  int M1(6);
  int min_smoothness1(6);

  interval beta2("1.75","1.75");
  interval eps2("-0.1","0.1");
  int M2(6);              
  int min_smoothness2(6);
 
  interval beta3("2.5","2.5");
  interval eps3("-0.3","0.3");
  int M3(6);              
  int min_smoothness3(6);

  IVector forcingA = IVector({1.});
  IVector forcingB = IVector({1.,1.,1.,1.});
 
  bsqVerifyExistenceOfPeriodicOrbit( beta1, sigma, eps1, forcingA, M1, min_smoothness1, verbose );
  bsqVerifyExistenceOfPeriodicOrbit( beta1, sigma, eps1, forcingB, M1, min_smoothness1, verbose );

  bsqVerifyExistenceOfPeriodicOrbit( beta2, sigma, eps2, forcingA, M2, min_smoothness2, verbose );
  bsqVerifyExistenceOfPeriodicOrbit( beta2, sigma, eps2, forcingB, M2, min_smoothness2, verbose );


  bsqVerifyExistenceOfPeriodicOrbit( beta3, sigma, eps3, forcingA, M3, min_smoothness3, verbose );
  bsqVerifyExistenceOfPeriodicOrbit( beta3, sigma, eps3, forcingB, M3, min_smoothness3, verbose );

  /*
  time (&end1);
  double elapsTime = difftime(end1, start1);
  cout << "Computation time: " << elapsTime << " seconds.\n"; // will show 0 seconds
  */

  return 0;
} 

