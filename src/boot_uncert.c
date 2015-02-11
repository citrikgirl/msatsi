#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define pi 3.14159265
#define torad pi/180
#define todeg 180/pi

main(argc, argv)
  int argc;char **argv; {
  char line[100]; /* character line */
  float tr[3], pl[3]; /* trend and plunge */
  float t1, t2, t3, d1, d2, d3; /* trend and plunge of best-fit solution */
  float trad[3], drad[3], trrad, plrad; /* trend and plunge, radians */
  FILE *fpin, *fpout; /* input and output file pointers */
  char namein[20], nameout[20]; /* input and output file names */
  float best[3][3], bestmag; /* best-fit stress tensor and it's length */
  float stress[3][3], mag; /* stress tensor and it's length */
  float dot[2000], bdot, level; /* tensor dot products */
  float conf, level95; /* confidence level, corresponding dot product */
  float phib, phi, phimin, phimax; /* best-fit stress ratio, min and max in conf region */
  float tr_min[3], tr_max[3], pl_min[3], pl_max[3]; /* min and max trend and plunge in conf region */
  float axis[3], bestax[9]; /* stress axes, best-fit stress axes */
  float z, pl_off, tr_off; /* dummy variables */
  int i, j, k, ic; /* more dummy variables */

  phimin = 1;
  phimax = 0;
  for (ic = 0; ic < 3; ic++) {
    pl_min[ic] = 0;
    pl_max[ic] = 0;
    tr_min[ic] = 0;
    tr_max[ic] = 0;
  }

  if (argc != 4) {
    printf("usage: boot_uncert *.slboot output_file confidence_level\n");
    return;
  }

  ++argv;
  sscanf(*argv, "%s", namein);
  fpin = fopen(namein, "r");
  if (fpin == NULL) {
    printf("unable to open %s\n", namein);
    return;
  }
  ++argv;
  sscanf(*argv, "%s", nameout);
  fpout = fopen(nameout, "a");
  if (fpout == NULL) {
    printf("unable to open %s\n", nameout);
    return;
  }
  ++argv;
  sscanf(*argv, "%f", &conf);

  /* first go through file to find confidence level */
  fgets(line, 100, fpin);
  sscanf(line, "%f %f %f %f %f %f", &best[0][0], &best[0][1], &best[0][2],
      &best[1][1], &best[1][2], &best[2][2]);
  best[1][0] = best[0][1];
  best[2][0] = best[0][2];
  best[2][1] = best[1][2];
  tenmag(best, &bestmag);
  fgets(line, 100, fpin);
  sscanf(line, "%f %f %f %f %f %f %f", &phib, &t1, &d1, &t2, &d2, &t3, &d3);
  trad[0] = t1 * torad;
  drad[0] = d1 * torad;
  trad[1] = t2 * torad;
  drad[1] = d2 * torad;
  trad[2] = t3 * torad;
  drad[2] = d3 * torad;
  vectorR(t1, d1, &bestax[0]);
  vectorR(t2, d2, &bestax[3]);
  vectorR(t3, d3, &bestax[6]);
  i = 0;
  while (fgets(line, 100, fpin) != NULL) {
    sscanf(line, "%f %f %f %f %f %f", &stress[0][0], &stress[0][1],
        &stress[0][2], &stress[1][1], &stress[1][2], &stress[2][2]);
    stress[1][0] = stress[0][1];
    stress[2][0] = stress[0][2];
    stress[2][1] = stress[1][2];
    fgets(line, 100, fpin);
    tenmag(stress, &mag);
    tendot(best, stress, bestmag, mag, &dot[i]);
    i++;
  }
  sort(dot, i);
  j = i * ((100. - conf) / 100.);
  level95 = dot[j];
  fclose(fpin);

  /* find the range of phi, trend, and plunge values within the uncertainty */
  fpin = fopen(namein, "r");
  fgets(line, 100, fpin);
  fgets(line, 100, fpin);
  k = 0;
  for (j = 0; j < i; ++j) {
    fgets(line, 100, fpin);
    sscanf(line, "%f %f %f %f %f %f", &stress[0][0], &stress[0][1],
        &stress[0][2], &stress[1][1], &stress[1][2], &stress[2][2]);
    stress[1][0] = stress[0][1];
    stress[2][0] = stress[0][2];
    stress[2][1] = stress[1][2];
    fgets(line, 100, fpin);
    sscanf(line, "%f %f %f %f %f %f %f", &phi, &tr[0], &pl[0], &tr[1], &pl[1],
        &tr[2], &pl[2]);
    tenmag(stress, &mag);
    tendot(best, stress, bestmag, mag, &level);
    if (level >= level95) {
      for (ic = 0; ic < 3; ic++) {
        vectorR(tr[ic], pl[ic], &axis);
        trrad = tr[ic] * torad;
        plrad = pl[ic] * torad;
        bdot = axis[0] * bestax[ic * 3] + axis[1] * bestax[ic * 3 + 1]
            + axis[2] * bestax[ic * 3 + 2];
        if (bdot < 0) {
          trrad -= pi;
          plrad = -plrad;
        }
        while ((trrad - trad[ic]) > pi)
          trrad -= (2 * pi);
        while ((trrad - trad[ic]) < -pi)
          trrad += (2 * pi);
        pl_off = todeg * (plrad - drad[ic]);
        tr_off = todeg * (trrad - trad[ic]);
        if (pl_off < pl_min[ic]) pl_min[ic] = pl_off;
        if (pl_off > pl_max[ic]) pl_max[ic] = pl_off;
        if (tr_off < tr_min[ic]) tr_min[ic] = tr_off;
        if (tr_off > tr_max[ic]) tr_max[ic] = tr_off;
      }
      k++;
      if (phi < phimin) phimin = phi;
      if (phi > phimax) phimax = phi;
    }
  }
  for (ic = 0; ic < 3; ic++) {
    pl_min[ic] = todeg * drad[ic] + pl_min[ic];
    pl_max[ic] = todeg * drad[ic] + pl_max[ic];
    tr_min[ic] = todeg * trad[ic] + tr_min[ic];
    tr_max[ic] = todeg * trad[ic] + tr_max[ic];
  }
  fprintf(fpout,
      "phi= %4.2f %4.2f %4.2f tr1= %7.2f %7.2f %7.2f pl1= %6.2f %6.2f %6.2f tr2= %7.2f %7.2f %7.2f ",
      phib, phimin, phimax, t1, tr_min[0], tr_max[0], d1, pl_min[0], pl_max[0],
      t2, tr_min[1], tr_max[1]);
  fprintf(fpout,
      "pl2= %6.2f %6.2f %6.2f tr3= %7.2f %7.2f %7.2f pl3= %6.2f %6.2f %6.2f\n",
      d2, pl_min[1], pl_max[1], t3, tr_min[2], tr_max[2], d3, pl_min[2],
      pl_max[2]);
}

tendot(ten1, ten2, mag1, mag2, pdot)
  float ten1[3][3];float ten2[3][3];float mag1, mag2;float *pdot; {
  int i, j;

  *pdot = 0;
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      *pdot += ten1[i][j] * ten2[i][j];
  *pdot /= mag1 * mag2;
  return;
}

tenmag(ten, pmag)
  float ten[3][3];float *pmag; {
  int i, j;
  double z;

  z = 0;
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      z += ten[i][j] * ten[i][j];
  z = sqrt(z);
  *pmag = z;
  return;
}

vectorR(trend, plunge, vec)
  float trend, plunge, vec[3]; {
  vec[0] = cos(plunge * torad) * cos(trend * torad);
  vec[1] = cos(plunge * torad) * sin(trend * torad);
  vec[2] = sin(plunge * torad);
  return;
}
