// 
// Kinetic Alfvén Wave Spectrum and Particle Interaction
//
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <complex.h>
#include <omp.h>
#include <sys/time.h>
//#include <mpi.h>

#define number_of_particles 10000
#define time_steps 200000
//~1 (temperature of electrons : temperature of protons)
#define TeTp 0.5
//slow solar wind (temperature of alphas : temperature of protons)
#define TaTp 4.0
//proton beta ~ (v_th / v_A)^2
#define beta 1.0
//fluids gamma
#define gamma 1.0
//mildly different definition of beta by Hollweg 1999
#define beta_joe (beta * gamma * 0.5 * (1. + TeTp))
//(mass of electron : mass of proton)
#define MeMp (1./1836.)
//normalised magnetic field fluctuation
//INITIAL:dB = 0.007049//give dvrho = 0.15 v_th//amplitude dB/B_0 of fluctuation at gyro-scale (kperp rho = 1) for waves with omega > 0
#define dB 0.15
//one side of the volume of the space centred around origin
#define BOX 100.0
//Alfvén velocity over speed of light
#define vAc 3e-4
#define PI 4*atan(1)
//minimum spatial frequency for AWs
#define kperp_min (1./3.793667895)
//maximum spatial frequency for KAWs -- do not increase this further as the model breaks down and results will be non-physical
#define kperp_max 3.793667895
//PROTONS -- we work with an isotropic Maxwellian distribution (i.e., a Gaussian where the standard deviation is modified by the plamsa beta)
#define mass (1.0)
#define charge (1.0)
#define stdev (sqrt(0.5 * beta))
//ALPHAS
//#define mass (4.0)
//#define charge (2.0)
//#define stdev (sqrt(0.5 * beta * TaTp / mass))
//ELECTRONS
//#define mass (1./1836.)
//#define charge (-1.)
//#define stdev (sqrt(0.5 * beta * TeTp / MeMp))



static inline double WTime(void){
  //timing the simulation for efficiency
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec / 1e6;
}

void Boris(double x[3], double v[3], double E[3], double B[3], double dt){
  //Boris Rotation for non-relativistic particle tracing -- "It really is state of the art" - D. Verscharen
  double v_minus[3];
  double v_prime[3];
  double v_plus[3];
  double t[3];
  double s[3];
  double t_mag2;

  //t vector
  t[0] = charge / mass * B[0] * 0.5 * dt;
  t[1] = charge / mass * B[1] * 0.5 * dt;
  t[2] = charge / mass * B[2] * 0.5 * dt;

  //|t|^2
  t_mag2 = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];

  //s vector
  s[0] = 2.0 * t[0] / (1.0 + t_mag2);
  s[1] = 2.0 * t[1] / (1.0 + t_mag2);
  s[2] = 2.0 * t[2] / (1.0 + t_mag2);

  //v minus
  v_minus[0] = v[0] + charge / (mass * vAc) * E[0] * 0.5 * dt;
  v_minus[1] = v[1] + charge / (mass * vAc) * E[1] * 0.5 * dt;
  v_minus[2] = v[2] + charge / (mass * vAc) * E[2] * 0.5 * dt;

  //v prime
  v_prime[0] = v_minus[0] + ( v_minus[1] * t[2] - v_minus[2] * t[1]);
  v_prime[1] = v_minus[1] + (-v_minus[0] * t[2] + v_minus[2] * t[0]);
  v_prime[2] = v_minus[2] + ( v_minus[0] * t[1] - v_minus[1] * t[0]);

  //v plus:
  v_plus[0] = v_minus[0] + ( v_prime[1] * s[2] - v_prime[2] * s[1]);
  v_plus[1] = v_minus[1] + (-v_prime[0] * s[2] + v_prime[2] * s[0]);
  v_plus[2] = v_minus[2] + ( v_prime[0] * s[1] - v_prime[1] * s[0]);

  //final v_n+1/2
  v[0] = v_plus[0] + charge / (mass * vAc) * E[0] * 0.5 * dt;
  v[1] = v_plus[1] + charge / (mass * vAc) * E[1] * 0.5 * dt;
  v[2] = v_plus[2] + charge / (mass * vAc) * E[2] * 0.5 * dt;

  //pusher
  x[0] += v[0] * dt;
  x[1] += v[1] * dt;
  x[2] += v[2] * dt;
}

void reset_totalEandBsq(double totalEandBsq[2]){
  //reset the memory of the array for each particle
  //I just use a two-element array with index (0,0) as the total electric field and (1,0) as the magnetic field squared
  totalEandBsq[0] = 0.;
  totalEandBsq[1] = 0.;
}

double gaussian_number(void){
  //obtain a gaussian number from two random numbers on a range of [0,1] inclusive.
  //this procedure CAN produce two gaussian numbers, but we only need one...c'est la vie
  //also, as we multiply by a standard deviation that is beta-dependent, this is really Maxwellian, aber es ist egal...
  double number1 = (double)rand() / (double)RAND_MAX ;
  double number2 = (double)rand() / (double)RAND_MAX ;
  return stdev * sqrt(-2. * log(number1)) * sin(2. * PI * number2);
}

void get_fields(double E[3], double B[3], double x[3]){
  //set the background/ambient field per the Kraichnan hypothesis. WLOG the magnetic field is orientated on the z-axis
  //below is a dipole code but should not be used in stochastic heating -- just particle tracing (simply turn off the KAW spectrum)
  E[0] = 0.;
  E[1] = 0.;
  E[2] = 0.;
  ///*
  B[0] = 0.;
  B[1] = 0.;
  B[2] = 1.;
  //*/
  //Magnetic dipole:
  /*
  double position_mag = sqrt(dot_product(x,x));
  const double m = 1000.; //magnetic moment
  B[0] = (3 * m * x[0] * x[2]) / (position_mag * position_mag * position_mag * position_mag * position_mag);
  B[1] = (3 * m * x[1] * x[2]) / (position_mag * position_mag * position_mag * position_mag * position_mag);
  B[2] = ((m) / (position_mag * position_mag * position_mag)) * (( (3 * x[2] * x[2]) / (position_mag * position_mag) ) - 1);
  */
}

void get_bhat(double B[3], double bhat[3]){
  //Normalise the Magnetic Field
  //b_hat = B / sqrt(B dot B) (i.e., the magnitude)
  double magB = (sqrt((B[0] * B[0]) + (B[1] * B[1]) + (B[2] * B[2])));
  bhat[0] = B[0] / magB;
  bhat[1] = B[1] / magB;
  bhat[2] = B[2] / magB;
}

void KAW(double x[3], double E[3], double B[3], double dt, int ts, double random_phase[162], double anisotropy, double totalEandBsq[2]){
  /*
  This calculates the dispersion relation for the KAW spectrum and obtains the electromagnetic field from an given particle's position
  */
  //frequency
  double omega;
  //k parallel
  double kparallel;
  //k perpendicular
  double kperp;
  //azimuthal coordinate
  double phi;
  //rotated k vector
  double kx,ky;
  //distance between k's
  double deltak;
  //components of the fields
  double complex Ex,Ey,Ez,Bx,By,Bz;
  //ratios of field components
  double complex BzBx,ByBz,ExEy,EzEy,EyBx;
  //rotation components for azimuthal sampling
  double complex Ezrot,Exrot,Eyrot,Bxrot,Byrot,Bzrot;
  //variations in the componets
  double Exvar,Eyvar,Ezvar,Bxvar,Byvar,Bzvar;
  //time
  double t = ts * dt;
  //amplidtude of B
  double B_amplitude;
  //Phases
  double phase,tot_phase;
  //number of positions in kperp
  int nkperp = 9;
  //number of positions in phi
  int nphi = 9;
  //upper power index for kperp rho < 1 (per Goldreich-Sridhar 1995)
  double pinup = -5./3.;
  //lower power index for kperp rho > 1 (per Cho-Lazarian 2004)
  double pindown = -7./3.;
  //ratio of max kparallel to kperp at kperp = ikperp_rhoi //range in kperp
  //double kperp_range[2] = {1./3.793667895,3.793667895};
  double kperp_range[2] = {kperp_min, kperp_max};

  get_fields(E,B,x);
  reset_totalEandBsq(totalEandBsq);

  //We sample 9 points in kperp at 9 points in the azimuthal direction, each having two suboptions for direction of alfven wave
  for (int i = 1; i <= nkperp; i++){
    //Obtain normalised kperp and kparallel:
    kperp = kperp_range[0] * pow(pow((kperp_range[1] / kperp_range[0]), 1./((double)nkperp -1.)), (i-1));

    //Distance in k
    deltak = pow((kperp_range[1] / kperp_range[0]), 1./((double)nkperp -1.));

    //determine the amplitude of the waves based off of our power indices
    if (kperp <= 1.){
      B_amplitude = dB * pow(kperp, ((pinup + 1.) / 2.));
      kparallel = anisotropy * pow(kperp, 2./3.);
    } else {
      B_amplitude = dB * pow(kperp, ((pindown + 1.) / 2.));
      kparallel = anisotropy * pow(kperp, 1./3.);
    }

    //normalise to units of inverse proton inertial length Omega_p / v_A
    kperp /= sqrt(beta);
    kparallel /= sqrt(beta);

    //Alfven waves are up or down, so we loop over both
    for (int k = 0; k < 2; k++){
      //Obtain omega by sqrt(Joe's eq 41)
      omega = 1. + (2. * beta_joe) + (kperp * kperp *beta_joe) + sqrt(pow((1. + kperp * kperp * beta_joe),2.) + 4. * kperp * kperp * beta_joe * (beta_joe - MeMp));
      omega = sqrt(omega / (2. * (1. + beta_joe + kperp * kperp * MeMp)));

      //account for up or down propagation
      if (k == 0){
        omega *= kparallel;
      } else {
        omega *= -kparallel;
      }

      //We use now ratios to determine components by Joe's equation 50 for up and down waves
      BzBx = I * ((omega * omega / (kparallel * kparallel))-1.) * kperp * beta_joe / ((omega / kparallel) * (beta_joe - MeMp * omega * omega / (kparallel * kparallel)));
      BzBx /= (omega * omega - kperp * kperp - kparallel * kparallel);
      ByBz =  -kparallel / kperp;

      //Likewise for electric field eq (44)
      ExEy = I * beta_joe * omega * (omega * omega / (kparallel * kparallel) - 1.) / (kperp * kperp + kparallel * kparallel - omega * omega);
      ExEy /= ((omega * omega / (kparallel * kparallel)) * (MeMp - beta * gamma * 0.5) - 0.5 * beta * gamma * TeTp);

      //and eq (47)
      EzEy = kparallel * (0.5 * gamma * beta * TeTp - MeMp * omega * omega / (kparallel * kparallel)) * (omega * omega / (kparallel * kparallel) -1.);
      EzEy /= (kperp * ((omega * omega / (kparallel * kparallel)) * (MeMp - 0.5 * beta * gamma) - 0.5 * beta * gamma * TeTp));
      EyBx = omega * vAc / (kperp * EzEy - kparallel);

      //thusly, our components fall out
      Bx = B_amplitude * sqrt(log(deltak) * 2. / ((double)nphi * (1. + fabs(BzBx * BzBx))));
      By = BzBx * ByBz * Bx;
      Bz = BzBx * Bx;

      Ex = ExEy * EyBx * Bx;
      Ey = EyBx * Bx;
      Ez = EzEy * Ey;

      //In terms of dvrho and dB_rms, we are interested in the fields only at these positions
      if (i == 4 || i == 5 || i == 6){
        totalEandBsq[0] += fabs(Ey) * fabs(Ey); //i.e. total E
        totalEandBsq[1] += fabs(Bx) * fabs(Bx) + fabs(Bz) + fabs(Bz); //i.e. total Bsq
      }

      //do the points in the azimuthal direction
      for (int j = 0; j < nphi; j++){
        //define the angle
        phi = j * 2. * PI / (double)nphi - PI/2.;

        //rotate k
        kx = -kperp * sin(phi);
        ky = kperp * cos(phi);

        //rotate amplitidues (cylindrically)
        Exrot = Ex * cos(phi) - Ey * sin(phi);
        Eyrot = Ex * sin(phi) + Ey * cos(phi);
        Ezrot = Ez;

        Bxrot = Bx * cos(phi) - By * sin(phi);
        Byrot = Bx * sin(phi) + By * cos(phi);
        Bzrot = Bz;

        //random phase distribution
        phase = random_phase[i + j * nkperp + k * nphi * nkperp];
        //vary periodically
        tot_phase = kx *  x[0] + ky * x[1] + kparallel * x[2] - omega * t - phase;

        //calculate variations in the fields
        Exvar = 1. * (creal(Exrot) * cos(tot_phase) - cimag(Exrot) * sin(tot_phase));
        Eyvar = 1. * (creal(Eyrot) * cos(tot_phase) - cimag(Eyrot) * sin(tot_phase));
        Ezvar = 1. * (creal(Ezrot) * cos(tot_phase) - cimag(Ezrot) * sin(tot_phase));

        Bxvar = 1. * (creal(Bxrot) * cos(tot_phase) - cimag(Bxrot) * sin(tot_phase));
        Byvar = 1. * (creal(Byrot) * cos(tot_phase) - cimag(Byrot) * sin(tot_phase));
        Bzvar = 1. * (creal(Bzrot) * cos(tot_phase) - cimag(Bzrot) * sin(tot_phase));

        //Add these values to the electric and magnetic fields:
        E[0] += Exvar;
        E[1] += Eyvar;
        E[2] += Ezvar;

        B[0] += Bxvar;
        B[1] += Byvar;
        B[2] += Bzvar;
      }
    }
  }

  ///*
  //Parallel electric field CORRECTION
  double bhat[3];
  get_bhat(B,bhat);
  double dummE = E[2];
  double bhat_dot_E = (bhat[0] * E[0] + bhat[1] * E[1] + bhat[2] * E[2]);
  E[0] += bhat[0] * (dummE - bhat_dot_E);
  E[1] += bhat[1] * (dummE - bhat_dot_E);
  E[2] += bhat[2] * (dummE - bhat_dot_E);
  //*/

  //total E
  totalEandBsq[0] *= 0.5 * (double)nphi;
  //total Bsq
  totalEandBsq[1] *= 0.5 * (double)nphi;
}

double get_dvrho(double totalEandBsq[2]){
  //returns dvrho from array
  return sqrt(totalEandBsq[0])/vAc;
}

double get_dB_rms(double totalEandBsq[2]){
  //returns dBrho from array
  return sqrt(totalEandBsq[1]);
}

double get_aniostropy(double x[3], double E[3], double B[3], double dt, int ts, double random_phase[162], double totalEandBsq[2]){
	//Use Newton-Secant method to obtain the anisotropy of the system.
  double dvrho;
  //number of iterations
  int iterations = 50;
  //iteration number
	int iteration = 0;
  //jump distance
	double jump;
  //we are considering the fourth k per 2010 paper
	int i = 5;
  //current root of the iteration
	double root;
  //previous root of the iteration
	double root_previous;
  //use this anisotropy...just a starting value...will reduce to true value
  double anisotropy = 2. * dB;
  //define an appropriate shift
  double previous_anisotropy = anisotropy - 0.00001;
  //define per Holwegs 1999 paper
  double omega;
  //k parallel and perpendicular directions
  double kparallel, kperp;
  //nine samplings
  int nkperp = 9;
  //range of k perp: 1d array with 2 values
  //double kperp_range[2] = {1./3.793667895,3.793667895};
  double kperp_range[2] = {kperp_min,kperp_max};

  //obtain dvrho
  KAW(x,E,B,dt,ts,random_phase,anisotropy,totalEandBsq);
  dvrho = get_dvrho(totalEandBsq);
	//printf("first: %g\n",dvrho);

  //k perpendicular (constant)
  kperp = kperp_range[0] * pow(pow((kperp_range[1] / kperp_range[0]), 1./(1.*nkperp -1.)), (i-1));

  //k parallel
  if (kperp <= 1){
    kparallel = previous_anisotropy * pow(kperp, 2./3.);
  } else {
    kparallel = previous_anisotropy * pow(kperp, 1./3.);
  }

  //Normalise
	kperp /= sqrt(beta);
  kparallel /= sqrt(beta);

  //Omega (multiply by kparallel here)
  omega = 1. + 2. * beta_joe + kperp * kperp * beta_joe + sqrt(pow((1. + kperp * kperp * beta_joe),2.) + 4. * kperp * kperp * beta_joe * (beta_joe - MeMp));
  omega = sqrt(omega / (2. * (1. + beta_joe + kperp * kperp * MeMp))) * kparallel;

  //define previous root
  root_previous = omega - kperp * dvrho;

  //boolean to determine accurate root to machine precision
  int DOIT = 1;

  //finding root
  while (iteration < iterations && DOIT == 1){
    iteration += 1;
    //k perpendicular (constant)
    kperp = kperp_range[0] * pow(pow((kperp_range[1] / kperp_range[0]), 1./(1.*nkperp -1.)), (i-1));
    //k parallel
    if (kperp <= 1){
      kparallel = anisotropy * pow(kperp, 2./3.);
    } else {
      kparallel = anisotropy * pow(kperp, 1./3.);
    }
    //Normalise
    kperp /= sqrt(beta);
    kparallel /= sqrt(beta);

    //Omega
    omega = 1. + 2. * beta_joe + kperp * kperp * beta_joe + sqrt(pow((1. + kperp * kperp * beta_joe),2.) + 4. * kperp * kperp * beta_joe * (beta_joe - MeMp));
    omega = sqrt(omega / (2. * (1. + beta_joe + kperp * kperp * MeMp))) * kparallel;

    //obtain a new dvrho
    KAW(x,E,B,dt,ts,random_phase,anisotropy,totalEandBsq);
    dvrho = get_dvrho(totalEandBsq);
    	
    //define current root
    root = omega - kperp * dvrho;

    //determine root to machine precicion or keep going using boolean
    //if not close enough, jump by the amount of the shifted anisotropy and continue
    if (fabs(root - root_previous) < 1e-80 || fabs(root) < 1e-20){
      jump = 0.;
      DOIT = 0;
    } else {
      jump = root * (anisotropy - previous_anisotropy) / (root - root_previous);
      //shift values to go to next time step
      previous_anisotropy = anisotropy;
      anisotropy -= jump;
      root_previous = root;
    }
  }
  //return the anisotropy
  return anisotropy;
}

void fill_vector(double X[3][number_of_particles], double x[3], int particle_number){
  //used for shared-memory OpenMP caculations
  x[0] = X[0][particle_number];
  x[1] = X[1][particle_number];
  x[2] = X[2][particle_number];
}

void update_matrix(double X[3][number_of_particles], double x[3], int particle_number){
  //used for shared memory OpenMP calculations
	X[0][particle_number] = x[0];
  X[1][particle_number] = x[1];
  X[2][particle_number] = x[2];
}

void fill_matrices(double X[3][number_of_particles], double V[3][number_of_particles]){
  //initialise the matrices for position and velocity
  //the if statement simply allows one to input specific values for tracing a single particle
  if (number_of_particles > 1){
    for (int j = 0; j< 3; j++){
      for (int i = 0; i < number_of_particles; i++) {
        V[j][i] = gaussian_number();
        X[j][i] = BOX * ((double)rand() / (double)RAND_MAX) - (BOX / 2.);
      }
    }
  } else {
    V[0][0] = 0.;
    V[1][0] = 1.;
    V[2][0] = 0.;
    X[0][0] = -1.;
    X[1][0] = 0.;
    X[2][0] = 0.;
  }
}

void get_rms(double V[3][number_of_particles], double rms[2]){
  //we wish to determine the rms perpendicular (to B) and parallel velocities over all of the particles
  double vavg[3];
  int i,j;
  for (i = 0; i < 3; i++){
    vavg[i] = 0;
    for (j = 0; j < number_of_particles; j++){
      vavg[i] += V[i][j];
    }
    vavg[i] /= (double)number_of_particles;
  }
  for (j = 0; j < number_of_particles; j++){
    rms[0] += (V[0][j] - vavg[0]) * (V[0][j] - vavg[0]) + (V[1][j] - vavg[1]) * (V[1][j] - vavg[1]);
    rms[1] += (V[2][j] - vavg[2]) * (V[2][j] - vavg[2]);
  }

  //NB: we actually want the rms-square velocities, but one can get both readily :)
  //perp rms
  //rms[0] = sqrt(rms[0] / (double)number_of_particles);
  rms[0] /= (double)number_of_particles;
  //parallel rms
  //rms[1] = sqrt(rms[1] / (double)number_of_particles);
  rms[1] /= (double)number_of_particles;
}


int main(){
  //time step size
  double dt = 1e-2;
  //time step number
  int ts = 0;
  //initialise particle count
  int particle;
  //define anisotropy value (obtained later) //anisotropy = 0.0091767732351740461;
	double anisotropy;

	//initialise position and velocity matrices
  double X[3][number_of_particles];
  double V[3][number_of_particles];

  //FILL MATRICES:
  fill_matrices(X,V);

  //INITALISE WAVE PARAMETERS: Determine random phases:
  double random_phase[162];
  for (int dim = 0; dim < 162; dim++){
    random_phase[dim] = 2 * PI * (double)rand() / (double)RAND_MAX;
  }
  //GET ANISOTROPY:
  double E[3];
  double B[3];
  double x[3];
  double totalEandBsq[2];
  fill_vector(X,x,0);
  anisotropy = get_aniostropy(x,E,B,dt,ts,random_phase,totalEandBsq);
  //Openfile
	FILE *fout1 = NULL;
	fout1 = fopen("rms_square.csv", "w");
	//FILE *fout2 = NULL;
	//fout2 = fopen("distribution.csv", "w");
	//TIMEING
  double time1 = WTime();
  //TIME LOOP:
  for (ts = 0; ts < time_steps; ts++){
    #pragma omp parallel for
    //PARTICLE LOOP:
    for (particle = 0; particle < number_of_particles; particle++){
      double E[3];
      double B[3];
      double x[3];
      double v[3];
      double totalEandBsq[2];
      fill_vector(X,x,particle);
      fill_vector(V,v,particle);
      KAW(x,E,B,dt,ts,random_phase,anisotropy,totalEandBsq);
      Boris(x, v, E, B, dt);
      update_matrix(X,x,particle);
      update_matrix(V,v,particle);
      /* //This prints all partices in phase space -- a massive amount of data !!
      if (ts%2==0){
          fprintf(fout1,"%d, %g, %g, %g, %g, %g, %g, %g\n",particle,ts*dt,x[0],x[1],x[2],v[0],v[1],v[2]);
      }
      //*/
      /*
      // Velocity Distribution
      if (ts == 10000 || ts == 199999){
        double vperp = sqrt(v[0] * v[0] + v[1] * v[1]);
        double vpara = v[2];
        fprintf(fout2, "%d, %g, %g\n", particle, vperp, vpara);
      }
      //*/
    }
    ///*
    //MEAN VELOCITY SQUARE PER TIME STEP: (nb if (ts * dt > 100 && ts%2 == 0))
    //if (ts%2 == 0){
        double rms[2];
        get_rms(V,rms);
        fprintf(fout1,"%g, %g, %g\n",ts*dt,rms[0],rms[1]);
    //}
    //*/
  }
  fclose(fout1);
  //fclose(fout2);
  double time2 = WTime();
  printf("TIME = %g\n",time2-time1);
  return 0;
} 
