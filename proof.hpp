
/* -----------------------------------------------------------------------------------------
 * This is a header file to bsq.cpp providing a rigorous procedure of verification 
 * of existence of a periodic orbit
 * in the FitzHugh-Nagumo system for given parameters beta, sigma. 
 * * ----------------------------------------------------------------------------------------*/

/* ------------------------------------------------------------------------------------ */
/* ----- VERIFICATION OF EXISTENCE OF PERIODIC ORBITS FOR GIVEN PARAMETER VALUES ------ */
/* ------------------------------------------------------------------------------------ */

// function yields ulr k=1..M and C in that order (as an DVector of size M+1) only purpose for that is to export it later to a .tex table
// all intervals are symmetric wr to 0 (can be verified by setting verbose=1) so we only give a "rightBound" for the output

DVector bsqVerifyExistenceOfPeriodicOrbit( interval _beta, interval _sigma, interval _eps, IVector _f, int M, int min_smoothness = 6, bool _verbose=0 ) 
{
  DVector resultVect( M+1 );

  try                   // we check negations of all assumptions to throw exceptions, if no exception is thrown existence of the orbit is verified
  {
    interval C_init = interval(10.)*( _eps.rightBound() );
    int s_init(4); 
    int m(0);

    IVector ulr_init( 10 );

    ulr_init(1)=interval(-0.5, 0.5);          // we took these bounds as slightly modified bounds from ZM - seem to work after rescaling! if the dimension is greater than M 
    ulr_init(2)=interval(-0.2, 0.2);          // then the leftovers are discarded for C/k^s
    ulr_init(3)=interval(-0.157, 0.157); 
    ulr_init(4)=interval(-0.046,0.046); 
    ulr_init(5)=interval(-0.018,0.018); 
    ulr_init(6)=interval(-0.0087,0.0087); 
    ulr_init(7)=interval(-0.0046,0.0046); 
    ulr_init(8)=interval(-0.0027,0.0027); 
    ulr_init(9)=interval(-0.0016,0.0016); 
    ulr_init(10)=interval(-0.0011,0.0011);

    ulr_init = ulr_init*( _eps );

    PolyBd bsqBd( _beta, _sigma, C_init, s_init, ulr_init, m, M, _eps, _f, _verbose );

    int min_no_iters( min_smoothness - s_init + 1 );  // (n+1) iterates are needed to refine the bounds n times
    int max_no_iters( min_no_iters + 5 );              // it is easy to fall out of convergence domain of self-consistent bounds so we limit max number of refinements

    for( int i=1; i<=max_no_iters+1; i++ )   
    {
      bool result = bsqBd.isolationTest();
      
      if( _verbose )
        cout << "\n -------------------------- \n";
      
      if( result == 1 && i >= min_no_iters )
      {
        cout << "PROOF OF EXISTENCE OF A PERIODIC ORBIT FOR PARAMETER VALUES BETA = " << _beta << " AND SIGMA = " << _sigma << " COMPLETED! \n";
        cout << "The forcing term for the potential: " << _eps << "*2*f(t)*(";
        for( unsigned int k = 1; k <= _f.dimension() - 1; k++ )
          cout << _f(k) << "*cos(" << k << "x) +";
        cout << _f( _f.dimension() ) << "*cos(" << _f.dimension() << "x))\n";
        cout << "where f is smooth, 2pi periodic, |f(t)| <= 1 for all t (eg. sint)\n \n";

        cout << "C=" << bsqBd.C << ", s=" << bsqBd.s << ", M=" << bsqBd.M << "\n";
        cout << "Refinement iterates needed: " << i-1 << "\n";
        cout << "Bound for the L2 norm of solution (i.e. u): " << bsqBd.computeL2norm() << "\n";
        cout << "Bound for the C0 norm of solution (i.e. u): " << bsqBd.computeC0norm() << "\n";
        cout << "Bound for the L2 norm of time derivative of solution (i.e. du/dt): " << bsqBd.computeL2DerNorm() << "\n";
        cout << "Bound for the C0 norm of time derivative of solution (i.e. du/dt): " << bsqBd.computeC0DerNorm() << "\n";
        cout << "Smoothness of solution (i.e. u, also s in 2C/k^s): C^" << bsqBd.s << "\n";
        cout << "\n -------------------------- \n \n \n";
  
        for( int i = 1; i <= M; i++ )
          resultVect(i) = bsqBd.ulr(i).rightBound();

        resultVect(M+1) = bsqBd.C.rightBound();
        break;
      }
      else
      {
        if( i < max_no_iters+1 )
        {
          if( _verbose )
            cout << "\n \n REFINING THE BOUNDS ITERATE " << i << " !\n \n";

          bsqBd.refineBounds();
        }
        else
        {
          if( i == max_no_iters + 1 )
            throw "No convergence of self-consistent bounds! \n";
        }
      }
    }
  }
  catch(const char* Message)
  {
    cout << Message << "EXISTENCE OF PERIODIC ORBIT FOR PARAMETER VALUES BETA=" << _beta << " AND SIGMA=" << _sigma << " NOT VERIFIED! \n";
    cout << "The forcing term for the potential: " << _eps << "*2*f(t)*(";
    for( unsigned int k = 1; k <= _f.dimension() - 1; k++ )
      cout << _f(k) << "*cos(" << k << "x) +";
    cout << _f( _f.dimension() ) << "*cos(" << _f.dimension() << "x))\n";
    cout << "where f is smooth, 2pi periodic, |f(t)| <= 1 for all t (eg. sint)\n";
    cout << "\n -------------------------- \n \n \n";
  }


  return resultVect;
};


