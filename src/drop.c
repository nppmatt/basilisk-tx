/**

A Newtonian/viscoelastic drop impacting a solid surface, with or 
without gravity. The equilibrium shape should be controlled by 
surface tension, in the absence of gravity, defined by the 
"contact angle" $\theta$ between the solid surface and the 
interface. 

A drop is initialized as a disk (2D) with an initial impacting 
velocity and the contact angle is $\theta$. The contact line 
position is tracked. If no breakup or rebound, drop oscillates and 
eventually relaxes to its equilibrium position. This equilibrium 
should be exact to within machine accuracy. The curvature along 
interface should be constant at equilibrium (this is with no
gravity present).

Note that small contact angles are not accessible yet (/src/contact.h).

*/

//#include "grid/multigrid.h"         // For Cartesian mesh
//#include "grid/octree.h"            // For 3D mesh adaptivity
#include "navier-stokes/centered.h"
#define FILTERED 1                    // Smooth out density and viscosity jumps (and relaxation time?)
#include "two-phase.h"
#include "curvature.h"
#include "reduced.h"                  // Using reduced gravity approach
#include "contact.h"
#include "vof.h"
#include "tension.h"
#include "log-conform.h"

/*
 * Viscoelastic properties used are the ratio of the solvent 
 * to the total viscoelastic viscosity (polymeric
 * plus solvent), BETA, and the relaxation time LAM.
 * BETA is non-dimensional.
 * LAM is the relaxation time in seconds.
 * */
#define BETA 1.0
#define LAM 0.0

/*
 * We initialize the maximum and minimum levels of refinement.
 * Best results are given by max values of 9-10, min may be 4 under.
 * */
#define LEVEL 10
const int maxlevel = LEVEL;
const int minlevel = LEVEL - 4;

/* Drop initial radius in meters. */
const double R0 = 0.0025;


scalar lambdav[], mupv[];


/*
 * The drop moves from top to bottom.
 * We allow the fluid to get through top boundary unimpeded.
 * */
u.n[top] = neumann(0);
p[top]   = dirichlet(0);


/**
The wall is at the bottom side. We apply a no-slip boundary condition. */

u.t[bottom] = dirichlet(0);          // Comment out for imposing slip boundary condition
//u.t[bottom] = neumann(0);          // This slip boundary condition
//u.n[bottom]  = dirichlet(0.);      // This is imposed by default
//tau_qq[bottom] = dirichlet(0);     // This is the i=3, j=3 component of viscoelastic stress in axi-symmetric case
//f[bottom] = neumann(0);            // This is for imposing fully non-wetting condition, i.e. theta0 = 180

/*
 * To set the contact angle, we allocate a [height-function field](/src/heights.h)
 * and set the contact angle boundary condition on its tangential component.
 * 53 deg corresponds to wax on plaskolite
 * */
vector h[];
double theta0 = 53.0;
h.t[bottom] = contact_angle (theta0*pi/180.0);

/* Signature: name, velocity */
double velocity = 1.0;
int main(int argc, char** argv)
{
    /* Domain size in meters. */
    L0 = R0 * 10;


  /**
  We initialize the physical properties of the
  two-phase system and the gravity value, all in SI. */

  rho1 = 850.0;
  rho2 = 1.225;
  mu1 = BETA*0.1;
  mu2 = 1.0e-5;      
  f.sigma = 0.03;                   // These are wax physical properties
  G.y = -9.807;                      // Gravitational acceleration
  velocity = atof(argv[1]);

  /**
  We initialize viscoelastic fields below. */
  
  mup = mupv;
  lambda = lambdav;

  /**
  We must associate the height function field with the VOF tracer, so
  that it is used by the relevant functions (curvature calculation in
  particular). */

  f.height = h;

  /**
  We set a maximum timestep, if needed for stability. */
  
  DT = 1.0e-4;
  //DT = HUGE;                      // For dimensionless time

  /**
  We run for the range of contact angles.
  (e.g. Basilisk Sessile Drop demo) */
  
  init_grid (1 << LEVEL);
  //for (theta0 = 15; theta0 <= 105; theta0 += 15) {
    run();
  //}
}

/**
The initial drop is a semi circle. */

event init (t = 0)
{
  /**
  At a wall of normal $\mathbf{n}$ the component of the viscoelastic
  stress tensor $tau_p_{nn}$ is zero. Since the bottom boundary is a wall, we
  set $tau_p_{yy}$ equal to zero at that boundary. */
  
  scalar s = tau_p.y.y;
  s[bottom] = dirichlet(0.0);

  /**
  * The drop is centered on (0,0.3*L0) and has a radius of R0.
  * 0.94864 of a 25cm domain corresponds to experiment of wax on plaskolite
  */

  fraction (f, - (sq(x) + sq(y - 0.15*L0) - sq(R0)));

  /**
  The initial velocity of the droplet is -1.0 (m/s) */
  foreach()
    u.y[] = -f[] * velocity;
}

/**
We add the acceleration of gravity, if not using reduced method.
event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] -= 9.81;
}
*/

event properties (i++) {
  foreach() {
    mupv[] = 0.1*(1.0 - BETA)*clamp(f[],0,1);  // Polymeric viscosity of the drop
    lambdav[] = LAM*clamp(f[],0,1);           // Relaxation time of the drop
  }
}

#if 0
event logfile (i++)
{
  fprintf (fout, "%g %g\n", t, normf(u.x).max);
}

event snapshots (t += 1)
{
  p.nodump = false;
  dump();
}
#endif

/**
We refine the region around the interface of the droplet. */

#if TREE
event adapt (i++) {
  adapt_wavelet ({f,u.x,u.y}, (double []){1.0e-2, 1.0e-3, 1.0e-3}, maxlevel, minlevel);
}
#endif

/**
We track the normalized spreading diameter of the droplet. */

event logfile (i += 1; t <= 0.6) {
  scalar pos[];
  position (f, pos, {1,0});
  fprintf ( stderr, "%.15f,%.15f,%.15f\n", t, statsf(pos).max, (statsf(pos).max)/R0 );
}

/* Movies: in 3D, these are in a z=0 cross-section. */
/*
event output_interface (i += 50; t <= 1) {
 
  {
    char *outfile2 = NULL;
    outfile2 = (char *) malloc(sizeof(char) * 256);
    sprintf(outfile2, "interface/%d.out", i);
    FILE * fp_interface = fopen (outfile2, "w");
    output_facets (f, fp_interface);
    fclose(fp_interface);
  }

  {
    static FILE * fp = popen ("ppm2mp4 f.mp4", "w");
    output_ppm (f, fp, min = 0, max = 1, n = 512);
  }

  //{
    //scalar l[];
    //foreach()
    //l[] = level;
    //static FILE * fp2 = popen ("ppm2mp4 level.mp4", "w");
    //output_ppm (l, fp2, min = 5, max = 7, n = 512);
  //}
}
*/
/**
At equilibrium (t = 10 seems sufficient), we output the interface
shape and compute the (constant) curvature. */

/*
event end (t = end)
{
  output_facets (f, stdout);
  
  scalar kappa[];
  curvature (f, kappa);
  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 2.*statsf(f).sum;
  fprintf (stderr, "%d %g %.5g %.3g\n", N, theta0, R/sqrt(V/pi), s.stddev);
}
*/

