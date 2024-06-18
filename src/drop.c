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

/* For 2D adaptive mesh. */
//#include "grid/multigrid.h"

/* For 3D adaptive mesh. */
//#include "grid/octree.h"

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "curvature.h"

/* Enables reduced-gravity approach. */
#include "reduced.h"

#include "contact.h"
#include "vof.h"
#include "tension.h"
#include "log-conform.h"

/* For reading configuration file. */
#include "include/toml.h"
#include <errno.h>

/* Directory creation. */
#include <sys/stat.h>
#include <sys/types.h>

/* Smooth out jumps in density and viscosity. */
#define FILTERED 1

/* Define default viscoelastic fields. */
scalar lambdav[], mupv[];

/* The drop moves from top to bottom.
 * We allow the fluid to get through top boundary unimpeded.
 */
u.n[top] = neumann(0);
p[top]   = dirichlet(0);

/**
The wall is at the bottom side. We apply a no-slip boundary condition. */
u.t[bottom] = dirichlet(0);          // Comment out for imposing slip boundary condition
//u.t[bottom] = neumann(0);          // This slip boundary condition
//u.n[bottom]  = dirichlet(0.);      // This is imposed by default
//tau_qq[bottom] = dirichlet(0);     // This is the i=3, j=3 component of viscoelastic stress in axi-symmetric case
//f[bottom] = neumann(0);            // This is for imposing fully non-wetting condition, i.e. theta0 = 180

/* Drop experiment parameters. */
double R0;
double initHeight;
double velocity;

/* To set the contact angle, we allocate a [height-function field](/src/heights.h)
 * and set the contact angle boundary condition on its tangential component.
 * */
vector h[];
double theta0;

/* Viscoelastic properties used are the ratio of the solvent 
 * to the total viscoelastic viscosity (polymeric
 * plus solvent), BETA, and the relaxation time LAM.
 * BETA is non-dimensional.
 * LAM is the relaxation time in seconds.
 */
//double BETA;
//double LAM;
#define BETA 1.0
#define LAM 0.0

/* We initialize the maximum and minimum levels of refinement.
 * Best results are given by max values of 9-10, min may be 4 under.
 */
int maxLevel;
int minLevel;

/* How long the simulation will last. */
double simDuration;

/* TODO: Break config parsing and globals to utility header */
toml_table_t* liquid;
toml_table_t* gas;
toml_table_t* drop;
toml_table_t* sim;

/* File naming parameters, for data export. */
char *name = NULL;
char *out_dir = NULL;

static void tomlError(const char* msg, const char* msg1)
{
    fprintf(stderr, "TOML PARSING ERROR: %s%s\n", msg, msg1?msg1:"");
    exit(1);
}

int main(int argc, char** argv)
{
    FILE* configFile;
    char errbuf[256];
    /* For now, execution ought to be something like:
     * ./bin/drop config/wax.toml
     */
    configFile = fopen(argv[1], "r");
    if (!configFile) {
        tomlError("Failed to open configuration - ", strerror(errno));
    }

    toml_table_t* config = toml_parse_file(configFile, errbuf, sizeof(errbuf));
    if (!config) {
        tomlError("Cannot parse loaded config - ", errbuf);
    }
    fclose(configFile);

    /* Config loaded and parsed, start loading and assigning parameters. */
    liquid = toml_table_in(config, "liquid");
    gas = toml_table_in(config, "gas");
    drop = toml_table_in(config, "drop");
    sim = toml_table_in(config, "sim");

    /* Specify name of experiment and the directory it will go in. */
    name = (char *) malloc(sizeof(char) * 64);
    out_dir = (char *) malloc(sizeof(char) * 64);

    name = (toml_string_in(config, "name")).u.s;
    sprintf(out_dir, "out/%s", name);
    mkdir(out_dir, 0755);

    /**
    We initialize the physical properties of the
    two-phase system and the gravity value, all in SI units.
    */
    rho1 = (double)(toml_double_in(liquid, "rho")).u.d;
    mu1 = BETA * (double)(toml_double_in(liquid, "mu")).u.d;
    f.sigma = (double)(toml_double_in(liquid, "sigma")).u.d;

    rho2 = (double)(toml_double_in(gas, "rho")).u.d;
    mu2 = (double)(toml_double_in(gas, "mu")).u.d;

    /* Acceleration of gravity */
    G.y = -9.80665;

    /* Set experiment parameters,
     * beginning with drop radius and the domain size dependent on it.
     */
    R0 = (double)(toml_double_in(drop, "radius")).u.d;
    velocity = (double)(toml_double_in(drop, "velocity")).u.d;
    theta0 = (double)(toml_double_in(drop, "contact-angle")).u.d;
    h.t[bottom] = contact_angle(theta0*pi/180.0);
    L0 = R0 * 10;

    /* We initialize viscoelastic fields below. */
    mup = mupv;
    lambda = lambdav;

    /**
    We must associate the height function field with the VOF tracer, so
    that it is used by the relevant functions (curvature calculation in
    particular). */
    f.height = h;

    /* Set viscoelastic polymer properties. */
    //BETA = (double)(toml_double_in(sim, "beta")).u.d;
    //LAM = (double)(toml_double_in(sim, "lam")).u.d;

    /* Set grid level. */
    maxLevel = (int)(toml_int_in(sim, "level")).u.i;
    minLevel = maxLevel - 4;

    /* We set a maximum timestep, if needed for stability.
     * DT is an interpreted constant. Hence the indirect assignment.
     */
    double config_dt = (double)(toml_double_in(sim, "max-timestep")).u.d;
    DT = config_dt;
    simDuration = (double)(toml_double_in(sim, "duration")).u.d;

    /*
     * We run for the range of contact angles.
     * (e.g. Basilisk Sessile Drop demo)
    init_grid (1 << LEVEL);
    for (theta0 = 15; theta0 <= 105; theta0 += 15) {
        run();
    }
    */

    init_grid (1 << maxLevel);
    run();
    printf("end main\n");
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
    * The drop is centered on a offset (eq. of circle) and has a radius of R0.
    */
    initHeight = (double)(toml_double_in(drop, "init-height")).u.d;
    fraction (f, - (sq(x) + sq(y - initHeight*L0) - sq(R0)));

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
    adapt_wavelet ({f,u.x,u.y}, (double []){1.0e-2, 1.0e-3, 1.0e-3}, maxLevel, minLevel);
}
#endif

/* We track the normalized spreading diameter of the droplet.*/
/*
event logfile (i += 1; t <= simDuration) {
    scalar pos[];
    position (f, pos, {1,0});
    fprintf ( stdout, "%.15f,%.15f,%.15f\n", t, statsf(pos).max, (statsf(pos).max)/R0 );
}
*/

/* TODO 2024-05-20: Move rendering code to separate section.*/
/* Movies: in 3D, these are in a z=0 cross-section. */
char *interfaceFile = NULL;
char *movieFile = NULL;
event define_output (t = 0) {
    interfaceFile = (char *) malloc(sizeof(char) * 64);
    movieFile = (char *) malloc(sizeof(char) * 64);
}
event output_interface (i += 50; t <= simDuration) {
    {
        sprintf(interfaceFile, "out/%s/%d.out", name, i);
        FILE * fp_interface = fopen (interfaceFile, "w");
        output_facets (f, fp_interface);
        fclose(fp_interface);
    }

    {
        sprintf(movieFile, "ppm2mp4 %s.mp4", name);
        static FILE * fp = popen (movieFile, "w");
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

