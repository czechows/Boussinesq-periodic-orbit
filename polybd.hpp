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
  interval C;
  int s;
  IVector ulr;         // |u_k| <= 2|ulr_k|
  interval beta;
  interval sigma;
  interval eps;
  IVector f;   // we assume perturbation of the form epsilon*f where f has a finite Fourier expansion stored in IVector f
  bool verbose;

  PolyBd( interval _beta, interval _sigma, interval _C, int _s, IVector _ulr, int _m, int _M, interval _eps = interval(0.), IVector _f = IVector({0.}), bool _verbose=0 )
    : m(_m),
      M( _M ),
      C(_C),
      s(_s),
      ulr(_ulr),
      beta( _beta ),
      sigma( _sigma ),
      eps( _eps ),
      f( _f ),
      verbose( _verbose )
  {
    unsigned int M_un( M );   
    if( f.dimension() > M_un )
      throw "Perturbations exceeding Fourier dimension M not implemented yet";

    if( !( beta*(m+1)*(m+1) > 1) )
      throw "Parameter m too small. Increase m and repeat.";

    if( !( beta > 1 ) )
      throw "Low mode isolation not implemented yet.";

    // below: we need to resize f so it is of dimension M for some evaluations
    {
      IVector tempf(f);
      f.resize(M);

      for( unsigned int k = 1; k <= tempf.dimension(); k++ )
        f(k) = tempf(k);

      for( int k = tempf.dimension() + 1; k <=M; k++ )
        f( k ) = 0.;
    }

    // below: we need to resize ulr so it is of dimension 2*M for some evaluations
    {                               
      IVector temp_ulr(_ulr);
      ulr.resize(2*M);

      for( int k = 1; k <= M; k++ )
        ulr(k) = temp_ulr(k);
   
      for( int k = M + 1; k <=2*M; k++ )
        ulr(k) = interval( (-C/pow(k,s)).leftBound(), (C/pow(k,s)).rightBound() );
    }
  }


  interval convolution( int k, int llimit, int rlimit ) // computes convolution \sum |u_k| * |u_k+k1|, k1 = llimit...rlimit
                                                        // for rlimit~infty this is ~IS(k)
  {
    interval result(0.);

    for( int k1=llimit; k1 <= rlimit; k1++ )
      result = result + 4*abs(ulr( k1 ))*abs(ulr( k1+k ));

    return result;
  }

  interval uSum( int llimit, int rlimit ) // computes \sum |u_k1|, k1 = llimit..rlimit 
  {
    interval result(0.);
 
    for( int k1 = llimit; k1 <= rlimit; k1++ )
      result = result + 2*abs(ulr(k1));

    return result;
  }

  interval IS_bound( int k ) // computes bound for IS as in the paper
  {
    interval result(0.);        

    if( k > M )
    {
      result = ( 2*C/( pow(k,s-1) * (M+1) ) ) * ( 2*C/( pow(M+1,s-1)*(s-1) ) + uSum(1,M) ) * interval(-1.,1.);
      return result;
    }
    else
    {
      interval tempSum(0.);

      for( int k1 = M-k+1; k1 <= M; k1++ )
       tempSum = tempSum + ( ( 2*abs(ulr(k1)) )/pow(k+k1,s) )*interval(-1.,1.);
      
      result = convolution( k, 1, M-k )*interval(-1,1) + 2*C*tempSum + ( 4*C*C/( pow(k+M+1,s) * (s-1) * pow(M, s-1) ) )*interval(-1,1);
      return result;
    }
  }

  interval FS( int k ) // computes bound for FS(k) explicitly, maybe we overestimate here and should compute FS_left(k), FS_right(k) but we don't bother (yet)
  {
    if( k > 2*M )                                             // some checks we are not out of range with our loop - we had a bug before
      throw "Procedure FS(k) error: k too large! \n";
    if( k < 1 )
      throw "Procedure FS(k) error: k too small! \n";

    interval result(0.);

    for( int k1 = 1; k1 <= k-1; k1++ )
      result = result + 4*abs(ulr( k1 ))*abs(ulr( k-k1 ));
    
    result = result*interval(-1,1);
    return result;
  }

  interval FS_bound( int k ) // FS bound for k > 2M
  {
    interval result(0.);
    result = ( 2*C/pow(k,s-1) )*( ( pow(2,s+1)*uSum(1,M) )/(2*M + 1) + ( C*pow(2,2*s+1) )/pow(2*M+1,s+1) + ( C*pow(2,s+1) )/( (s-1)*pow(M,s) ) );
    result = result*interval(-1,1);
    return result;
  }

  interval D1() // computes D1 from the nonlinearity bound D1/k^(s-1) as in the paper
  {
    interval result(0.);

    for( int k = M+1; k <= 2*M; k++ )
      result = max( result, pow(k,s-1)*abs(FS(k)) );

    result = max( result, abs( FS_bound(1) ) ).rightBound(); // FS_bound(1) is precisely FS_bound(k)*k^(s-1) (that expression is independent of k)
    return result;
  }

  interval D2()
  {
    return ( abs( IS_bound(M+1) )*pow(M+1,s-1) ).rightBound(); // this is shorter than retyping all IS(k)
  }

  interval D()
  {
    return ( D1()+2*D2() ).rightBound();
  }

  void refineBounds() 
  {
    s=s+1;
    C = ( ( sigma*D() ) / ( 2*( beta - pow( (M+1), -2 ) ) ) ).rightBound() ;

    for( int k = m+1; k <= M; k++ )
    {
      interval G( 2*IS_bound( k ) + FS( k ) );
      ulr( k ) = ( sigma*G - ( eps* f( k )/ (k*k) ) )/( 2*(beta*k*k - 1) );   
    }

  }

  bool checkIsolation( int k )
  {
    interval G( 2*IS_bound( k ) + FS( k ) );

    if( subsetInterior( ( sigma*G - ( eps*f( k )/ (k*k) ) )/( 2*(beta*k*k - 1) ), ulr(k) ) ) 
      return 1;
    else
      return 0;  
  }

  bool checkFarTail()
  {
    if( C.leftBound() >  ( ( 1./(M+1) ) * ( sigma*D() / (2*(beta - pow( (M+1) ,-2)) ) ) ).rightBound() )
      return 1;
    else
      return 0;
  }

  bool verifyBounds() // verifies self-consistent bounds for k >= m
  {
    bool result( checkFarTail() );

    for( int k = m+1; k<=M; k++ )
      result = result*checkIsolation( k );
    
    return result;
  }

  bool isolationTest()  // here we unnecesarily repeat some evaluations in verbose part to have more 
                        // self-documented code - to improve efficiency, eg. for rigorous integration I should clean this
  {
    bool result(1);

    if( verbose )
    {
      cout << "Current C is: " << C << "\n";
      cout << "Current s is: " << s << "\n";
      cout << "Current status of far tail isolation is: " << checkFarTail() << "\n";

      cout << "D1 is " << D1() << "\n";
      cout << "D2 is " << D2() << "\n";
      cout << "D is " << D() << "\n";
      cout << "For far tail isolation C needs to be greater than: " << ( ( sigma*D()/(M+1) ) * ( 1. / (2*(beta - pow(M+1,-2)) ) ) ).rightBound() << "\n";
      
      cout << "\n \n";
    }

    result = result*checkFarTail(); 

    for( int k = m+1; k<=M; k++ )
    {
      if( verbose )
      {
        cout << "Current ulr(" << k << ") is: " << ulr(k) << "\n";
        cout << "Current status of the " << k << "-th mode isolation is: " << checkIsolation(k) << "\n";
    
        interval G( 2*IS_bound( k ) + FS( k ) );
        cout << "For isolation one needs the above to be an overset of " << ( sigma*G - ( eps*f( k )/ (k*k) ) )/( 2*(beta*k*k - 1) ) << "\n";
      }

      result = result*checkIsolation(k);
    }

    return result;
  }

  /* computing the norms we remember that |uk| = 2*|ulr|*/

  interval computeC0norm()
  {
    interval result(0.);
    result = result + uSum(1,2*M);
    result = result + 2*C/( (s-1)*pow(2*M,s-1) );            // integral estimate for the tail terms int_(2*M+1)^infty dx/x^s - see Lemma 3.1 in ZM
  
    return result.rightBound();
  }

  interval computeL2norm()
  {
    interval result(0.);
 
    for( int k1 = 1; k1 <= 2*M; k1++ )
      result = result + 4*abs(ulr(k1))*abs(ulr(k1));      

    result = result + (4*C*C) / ( (2*s-1)* pow(2*M, 2*s-1) );   // ak^2 \leq C^2/k^{2s}
    result = sqrt(2*interval::pi()*abs(result)); // we use Parsevals theorem \sum a_k^2 = L2 norm/2*pi


    return result.rightBound();
  }

  interval computeC0DerNorm()
  {
    interval result(0.);
 
    for( int k1 = 1; k1 <= 2*M; k1++ )
      result = result + 2*sqrt( k1*k1*(beta*k1*k1-1) )*abs(ulr(k1));

    result = result + 2*C*sqrt(beta)/( (s-3)*pow(2*M,s-3) );            // integral estimate for the tail terms int_(2*M+1)^infty dx/x^(s-2) - see Lemma 3.1 in ZM
  
    return result.rightBound();
  }

  interval computeL2DerNorm()
  {
    interval result(0.);
 
    for( int k1 = 1; k1 <= 2*M; k1++ )
      result = result + 4*sqrt( k1*k1*(beta*k1*k1-1) )*abs(ulr(k1))*abs(ulr(k1));

    result = result + (4*C*C*beta) / ( (2*s-5)* pow(2*M, 2*s-5) );
    result = sqrt(2*interval::pi()*abs(result));

    return result.rightBound();
  }
};


