//-------------------------------------------------------------------------------------------------
#include <math.h>
//-------------------------------------------------------------------------------------------------
/* to find ddir, dip, rake of other plane given same of first */

#define TODEG 57.29577951

//-------------------------------------------------------------------------------------------------
void stridip(double n, double e, double u, double *strike, double *dip);

//-------------------------------------------------------------------------------------------------
void ranger(z)
  float *z;
/* makes z in 0 to 360 */
  {
  while (*z >= 360)
    *z -= 360;
  while (*z < 0)
    *z += 360;
  }

//-------------------------------------------------------------------------------------------------
void switcher(float ddir1, float dip1, float rake1, float *ddir2, float *dip2,
    float *rake2)
  {
  double z, z2, z3, s1, s2, s3;
  double n1, n2, h1, h2;
  //int j;

  z = ddir1 / TODEG;
  if (dip1 == 90)
    dip1 = 89.99999;
  z2 = dip1 / TODEG;
  z3 = rake1 / TODEG;
  /* slick vector in plane 1 */
  s1 = -cos(z3) * cos(z) - sin(z3) * sin(z) * cos(z2);
  s2 = cos(z3) * sin(z) - sin(z3) * cos(z) * cos(z2);
  s3 = sin(z3) * sin(z2);
  n1 = sin(z) * sin(z2); /* normal vector to plane 1 */
  n2 = cos(z) * sin(z2);
  //n3 = cos(z2);
  h1 = -s2; /* strike vector of plane 2 */
  h2 = s1;
  /* note h3=0 always so we leave it out */
  stridip(s2, s1, s3, &z, &z2);
  z += 90.;
  *ddir2 = z;
  ranger(ddir2);
  *dip2 = z2;
  z = h1 * n1 + h2 * n2;
  z /= sqrt(h1 * h1 + h2 * h2);
  z = acos(z);
  //return;
  if (s3 >= 0)
    *rake2 = (z * TODEG);
  else
    *rake2 = (-z * TODEG);
  }

