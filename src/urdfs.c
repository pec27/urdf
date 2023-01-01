/*
  Module for calculating Unitary BRDF terms

  Peter Creasey 2015-2023
 */

#define PI 3.14159265359
// Limit for n-choose-r in uint32_t
#define MAX_K 30 
#include <math.h>
#include <complex.h>

static double term_k(const unsigned int k, const double complex z,
		     const double half_one_minus_x, const double two_inv_x1,
		     const double half_one_minus_y, const double two_inv_y1)
{

  double inv_k_2[MAX_K*2+1];      // inv_k_2[m] = 1 / (2+k+m)
  double z0_term_re[MAX_K*2+1]; // Re[z^m     / (1-z)^(m+2) / (3k choose m)]
  double z1_term_re[MAX_K*2+1]; // Re[z^(m+1) / (1-z)^(m+2) / (3k choose m)]
  double z2_term_re[MAX_K*2+1]; // Re[z^(m+2) / (1-z)^(m+2) / (3k choose m)]  

  // Initialise inv_k_2, z?_term_re for m=0...2k
  {
    const double one_minus_z_re = 1 - creal(z);
    const double z_im2 = cimag(z)*cimag(z);
    const double one_minus_z_re2 = one_minus_z_re * one_minus_z_re;
    // |z|^2
    const double z_mag2 = creal(z)*creal(z) + z_im2;
    // 1 / |1-z|^2
    const double inv_mag2_omz = 1 / (one_minus_z_re2+z_im2);
    // e_k (1-z)^{-2} / (k+1)		
    const double gamma = inv_mag2_omz*inv_mag2_omz / (1+k);
    
    // Initialise v to e_k (1-z)^{-2} / (k+1)	
    double v_re = gamma*(one_minus_z_re2 - z_im2);
    double v_im = 2*gamma*cimag(z)*one_minus_z_re;
    
    // Setup the values (2+k+m)^-1, Re[v_m], Re[v_m z] and Re[v_m z^2] for
    // v_m = e_k z^m / (1-z)^{m+2} (k!) / (k+m+2)!	
    for (unsigned int m=0,m_2=2, m_end=2*k; m<=m_end;++m, ++m_2)
    {
      
      const double inv = 1.0 / (m_2+k);
      inv_k_2[m] = inv;
      
      // Now we set v to 
      // v = v[m] = e_k z^m / (1-z)^{m+2} (k!) / (k+m+2)!	
      v_re *= inv;
      v_im *= inv;
      
      // Re[v]
      z0_term_re[m] = v_re;	  
      
      // v*z
      const double vz_re = v_re*creal(z) - v_im*cimag(z);
      const double vz_im = v_re*cimag(z) + v_im*creal(z);	  
      // Re[v*z]
      z1_term_re[m] = vz_re;
      // Re[v*z*z]
      z2_term_re[m] = vz_re*creal(z) - vz_im*cimag(z);
      
      // Update for next loop iteration
      const double tau = m_2*inv_mag2_omz;
      // v -> v (1+m) z/(1-z)
      v_re = (vz_re - v_re*z_mag2) * tau;
      v_im = (vz_im - v_im*z_mag2) * tau;
    }
  }

  // Accumulate the result (loop over a and b)
  double res = 0;

  // Initial values 
  double BBB0=1;
  double AAA=1;
  double AA=0; 
  
  const unsigned int k1k2 = (k+1) * (k+2);
  // Loop over a=k...0
  for (unsigned int a=k+1, ka=2*k, kCa=1, kma=0; a--; ka--)
  {
    // Next Jacobi polys a+1-> a
    AA  = ((2+k)*AAA - kma*AA) * inv_k_2[a]; // P^(a,1-a)_k (x)

    // (k choose a) P^(a,-a) from P^(a+1,-a) and P^(a+2,-a)
    const double A  = ((a+1 + 2*(k+1)*half_one_minus_x)*AA -(2+k)*AAA*half_one_minus_x)/(ka+1); 
    
    // Coefficient of 2nd and 3rd term
    const double rho1 = kCa*(2*k+1)*A;  
    const double rho2 = kCa*2*(k+1)*AA; 	  
    const double rho3 = kCa*k1k2*(2*k+1)*AAA; 
    
    // k,a,b=k term (written BB and B in terms of BBB)
    double res_k_a = BBB0*((ka+2)*((ka+1)*rho1*z0_term_re[ka] + rho2*z1_term_re[ka]) - 
			   rho3*z2_term_re[ka]);
    
    // Temporary variables for the inner loop
    double BBB =BBB0;
    BBB0 = (ka+2) * BBB0 * inv_k_2[ka]; // update for next a
    double BB = BBB0;
    
    // loop over b=k-1...0, (k choose b)
    for (unsigned int b=k,kCb=k,kmb=1;b--;)
    {
      const unsigned int ab = a+b;
      // Eqn. 38: P^(a+b+2,-b) from P^(a+b+3,-b-1) and P^(a+b+2,-b-1)
      BBB = (BB - half_one_minus_y*BBB)*two_inv_y1;
      
      const double tau = (2+ka)*BBB - kmb*BB;
      const double mu = z0_term_re[ab] * rho1;
      // Eqn. 39: P^(a+1+b,-b) from P^(a+b+2,-b-1) and P^(a+b+2,-b)
      BB  = tau* inv_k_2[ab];
      // Eqn. 40: (1+a+b+k)* P^(a+b,-b)_k (y) from P^(a+1+b,-b) and P^(a+2+b,-b)
      
      // k,a,b term
      res_k_a += (tau*(rho2*z1_term_re[ab]+ mu * ((ab+1+half_one_minus_y*(2+k+ka)))) - 
		  BBB*(rho3*z2_term_re[ab]+(2+ka)*(ka+b+2)*mu*half_one_minus_y))*kCb;
      
      // Update (k choose b) for next b
      kCb = (kCb * b)/(++kmb);
    }
    
    // Next Jacobi polys a-> a-1
    AAA = (AA - half_one_minus_x*AAA)*two_inv_x1; // P^(a+1,1-a) from P^(a+2,-a) and P^(a+1,-a)

    // Next (k choose a)
    kCa = (kCa * a) / (++kma);
    
    // Accumulate
    res += res_k_a;
  }
  
  return res;
}

static double iso_jacobi_sum(const double nu, const double complex z0, const double x, const double y,
			     const int kmax)
{
  /*
    Compute the isotropic BRDF terms using expansion of the Zernike 
    eigenfunctions in terms of the Jacobi polynomials

    
    i.e the sum

    P(z0,x,y) := 1/pi * Sum_{k=0}^kmax exp(-8 sig^2 k(k+1)) Sum_{a,b=0}^k (k choose a) (k choose b) *
                              (1+a+b) * z^(a+b) / (1-z)^(2+a+b) *    
			      real[ (2k+1) P^(a,-a)_k(x) P^(a+b,-b)_k(y) / (k+a+b choose k)
		                    +2z P^(1+a,-a)_k(x) P^(1+a+b,-b)_k(y) / (k+a+b+1 choose k+1)
				    -(2k+1) z^2 P^(2+a,-a)_k(x) P^(2+a+b,-b)_k(y) / (k+a+b+2 choose k+2)]

				    s.t. z = (z0 * exp(-4(2k+1)sig^2))

    where we are given
    z0   := (xr+I*yr)*(xs-I*ys) where xr,yr,xs,ys are the projected coordinates
           of the rays v_r,v_s
    x    := 1 - 2*(xr*xr+yr*yr)
    y    := 1 - 2*(xs*xs+ys*ys)
    nu   := exp(-4 sig^2) determines the falloff
    kmax := largest k term to include

    returns P
   */

  // prefactor used in all the z values
  if (kmax>MAX_K)
    return -1; 

  const double nu2 = nu*nu;// exp(-8 sig^2)
  const double nu4 = nu2*nu2;
  const double two_inv_y1 = 2.0 / (1+y);
  const double two_inv_x1 = 2.0 / (1+x);  
  const double half_one_minus_y = 0.5-0.5*y;
  const double half_one_minus_x = 0.5-0.5*x;
  
  double eta_nu4k = 1.5 * (x+1)*(y+1);
  double e_k=1.0/PI, res=0.0;
  double complex z = z0*nu;
  
  for (int k=0,half_3k4_3k5=10,k1_k3=3;k<=kmax;k++)
  {
    res += e_k * term_k(k, z, half_one_minus_x, two_inv_x1, half_one_minus_y, two_inv_y1);
      
    // new k (k->k+1)
    z        *= nu2; // next exp(-4(2k+1)sig^2)
    eta_nu4k *= nu4;
    e_k      *= half_3k4_3k5 * eta_nu4k / k1_k3;
    
    half_3k4_3k5 += 9*(k+2); // (3k+4)(3k+5)/2
    k1_k3        += 2*k+5;   // (k+1)(k+3)
  }

  return res;
}

void brdf_iso_jacobi(const double *xrxs, const double *nu, const int num_pts, const int kmax, double *out)
{
  /* 
     Calculate the isotropic BRDF using the sum of Jacobi polynomials

  */

  for (int i=0;i<num_pts;i++)
    {
      const double xr = xrxs[i<<2],yr=xrxs[i<<2|1], xs=xrxs[i<<2|2], ys=xrxs[i<<2|3];
      const double complex z0 = (xr + I*yr)*(xs - I*ys);
      const double x = 1-2*(xr*xr+yr*yr), y = 1-2*(xs*xs+ys*ys);

      out[i] = iso_jacobi_sum(nu[i],z0,x,y, kmax);
    }
}

