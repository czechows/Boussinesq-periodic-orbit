
/* -----------------------------------------------------------------------------------------
 * This is a header file to bsq.cpp providing a rigorous procedure of verification 
 * of existence of a periodic orbit
 * in the FitzHugh-Nagumo system for given parameters beta, sigma. 
 * * ----------------------------------------------------------------------------------------*/

/* ------------------------------------------------------------------------------------ */
/* ----- VERIFICATION OF EXISTENCE OF PERIODIC ORBITS FOR GIVEN PARAMETER VALUES ------ */
/* ------------------------------------------------------------------------------------ */


void bsqVerifyExistenceOfPeriodicOrbit( interval _beta, interval _sigma, bool _verbose = 0 ) 
{
  try                   // we check negations of all assumptions to throw exceptions, if no exception is thrown existence of the orbit is verified
  {
    PolyBd bsqBd( _beta, _sigma );
  }
  catch(const char* Message)
  {
    cout << Message << "EXISTENCE OF PERIODIC ORBIT FOR PARAMETER VALUES BETA=" << _beta << " AND SIGMA=" << _sigma << " NOT VERIFIED! \n";
  }
};


