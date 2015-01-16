#include <math.h>
//*** Corr: Patricia Martinez-Garzon/Grzegorz Kwiatek 30.08.2012
// Wrong constant was used.
//#define TODEG 57.2577951
#define TODEG 57.29577951
//*** END: Correction: Patricia Martinez-Garzon/Grzegorz Kwiatek 30.08.2012

void dirplg(double e, double n, double u, double *pdir, double *pplg)

//void dirplg(e, n, u, pdir, pplg)
/* to find direction and plunge of a vector */
//double e, n, u; /* the vector in east,north,up coordinates */
//double *pdir, *pplg; /* pointers to the direction in east of north */
/* and the plunge down the direction */
{
  double z; /* dummy variable */
  z = e * e + n * n;
  z = sqrt(z);

  //*** Corr: Patricia Martinez-Garzon/Grzegorz Kwiatek 30.08.2012
  // Protection against atan2 domain exception.
  //*pplg=atan2(-u,z)*TODEG;
  if (u == 0.0f && z == 0.0f)
    *pplg = 0.0;
  else
    *pplg = atan2(-u, z) * TODEG;
  //*** END Corr: Patricia Martinez-Garzon/Grzegorz Kwiatek 30.08.2012

  if (*pplg < 0) {
    *pplg = -*pplg;
    e = -e;
    n = -n;
  }

  //*** Corr: Patricia Martinez-Garzon/Grzegorz Kwiatek 30.08.2012
  // Protection against atan2 domain exception.
  //*pdir=atan2(e,n)*TODEG;
  if (e == 0.0f && n == 0.0f)
    *pdir = 0.0;
  else
    *pdir = atan2(e, n) * TODEG;
  //*** END Corr: Patricia Martinez-Garzon/Grzegorz Kwiatek 30.08.2012
}
