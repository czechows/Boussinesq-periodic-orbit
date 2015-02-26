/* -----------------------------------------------------------------------------------------
 * This is a header file to bsq.cpp providing a class for self-consistent bounds
 * for the Boussinesq system
 * * ----------------------------------------------------------------------------------------*/



/* ---------------------------------------------------------------------------------------------- */
/* ---------------------------- SELF CONSISTENT BOUNDS  ----------------------------------------- */
/* ---------------------------------------------------------------------------------------------- */



class PolyBd                 // class for self-consistent bounds
{
public:
  int m;
  int M;

  PolyBd( interval _beta, interval _sigma, int _m=3, int _M=10 )
    : m(_m),
      M(_M)
  {
    cout << "Class constructed! \n";
  }

  interval D( int k ) // computes D from the nonlinearity bound D/k^s as in the paper
  {
    return interval(5.);
  }
};


