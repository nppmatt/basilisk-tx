/* Adapted from pw_fast_fine variant */

/* #include "axi.h" */
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "log-conform.h"
#include "curvature.h"
#include "view.h"
#include "contact.h"
#include "vof.h"

#include "tension.h"

#define RHO_r 0.01
#define MU_r 0.01
#define RE 5.0
#define FR 2.26
#define LEVEL 8
#define BETA 0.1
#define WI 1.0 // Weissenberg number

scalar lambdav[], mupv[];
u.n[top] = neumann(0);
p[top] = dirichlet(0);

scalar tau_qq[];
tau_qq[bottom] = dirichlet(0);
u.t[bottom] = dirichlet(0);

/* The so-called "wetting condition"(?), which is the fraction of fluid enforced by VoF method. */
//f[bottom] = 1.0;
//f[bottom] = neumann(0);


/* Contact angle specification (removed axial symmetry) */
vector h[];
double theta0 = 90.0;
h.t[bottom] = contact_angle(theta0 * pi / 180.0);

/* Surface tension energy */
scalar sigma[] =  3.0;

int main() {
	size(3.2);
	init_grid(1 << LEVEL);

	rho1 = 1.0;
	rho2 = RHO_r;
	mu1 = BETA / RE;
	mu2 = MU_r / RE;

	mup = mupv;
	lambda = lambdav;
	DT = 2e-3;

	/* Contact angle height fcn field applied to VOF */
	f.height = h;
	run();
}

event init (t = 0) {
	/* Viscoelastic stress tensor at wall with normal y. */
	scalar s = tau_p.y.y;
	s[bottom] = dirichlet(0.0);
	fraction(f, -sq(x - 1.6) - sq(y - 2.0) + sq(0.5));

	foreach()
		u.y[] = -f[] * 2.0;
}

/*
event acceleration (i++) {
	face vector av = a;
	foreach_face(x)
		av.x[] -= 0.0 / sq(FR);
}
*/

event properties (i++) {
	foreach() {
		mupv[] = (1.0 - BETA) * clamp(f[], 0, 1) / RE;
		lambdav[] = WI * clamp(f[], 0, 1);
	}
}

#if TREE
event adapt (i++) {
	adapt_wavelet ({f, u.x, u.y}, (double[]){1e-2, 5e-3, 5e-3},
			maxlevel = LEVEL, minlevel = LEVEL - 2);
}
#endif

/* Logging spread diameter */
event logfile (i += 20; t <= 10) {
	scalar pos[];
	position (f, pos, {0,1});
	fprintf(stderr, "%g %g\n", t, 2.0 * statsf(pos).max);
}

/* Refer to draw.h documentation. 
 * tx/ty are fractional and direction is negative.
 * */
event viewing (i += 10) {
  view (width = 1024, height = 1024, fov = 20, tx = -0.5, ty = -0.5,
	quat = {0, 0, 0, 0});

  clear();
  draw_vof ("f", lw = 2);
  squares ("u.y", linear = true);
  box (notics = true);
  /*
  mirror ({0,1}) {
    draw_vof ("f", lw = 2);
    squares ("u.y", linear = true);
    box (notics = true);
  }
  */

  save ("movie/tension_test_90.mp4");
#if 0
  static FILE * fp = popen ("bppm","w");
  save (fp = fp);
#endif
}

