#include <math.h>
#define TORAD  57.29577951

/* finds the strike and dip of a plane given its normal */
/* vector, output is in degrees north of east and then  */
/* uses a right hand rule for the dip of the plane */
void stridip(n, e, u, strike, dip)
  double n, e, u;double *strike, *dip; {
  double x;
  if (u < 0.) {
    n = -n;
    e = -e;
    u = -u;
  }
  //*** Corr: Patricia Martinez-Garzon/Grzegorz Kwiatek 30.08.2012
  // does not really make a sense, but at least atan2 domain exception is
  // not thrown.
  //*strike=atan2(e,n)*TORAD;
  if (e == 0.0 && n == 0.0)
    *strike = 0.0;
  else
    *strike = atan2(e, n) * TORAD;
  //*** ENDCorr: Patricia Martinez-Garzon/Grzegorz Kwiatek 30.08.2012

  *strike = *strike - 90.;
  if (*strike < 0.)
    *strike += 360.;
  if (*strike > 360.)
    *strike -= 360.;
  x = sqrt(n * n + e * e); /* x is the horizontal magnitude */
  //*** Corr: Patricia Martinez-Garzon/Grzegorz Kwiatek 30.08.2012
  // Does not make any sense with dip=0, but at least program should not
  // crash due to the atan2 domain error.
  //*dip=atan2(x,u)*TORAD;
  if (x == 0.0 && u == 0.0)
    *dip = 0.0;
  else
    *dip = atan2(x, u) * TORAD;
  //*** ENDCorr: Patricia Martinez-Garzon/Grzegorz Kwiatek 30.08.2012
  return;
}
