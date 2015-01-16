//-------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
//-------------------------------------------------------------------------------------------------

/* a program to run the bootstrap operation for satsi */

//-------------------------------------------------------------------------------------------------
// BOOTMECH_2D
//
// Original code:
//   Jeanne Hardebeck <jhardebeck@usgs.gov>
//   available at: http://earthquake.usgs.gov/research/software/
// 
// Corrections to the original code:
//   Grzegorz Kwiatek [GK] <kwiatek@gfz-potsdam.de> <http://www.sejsmologia-gornicza.pl/about>
//   Patricia Martinez-Garzon [PM] <patricia@gfz-potsdam.de>
// 
//   Code updated to C99 standard. 
//
// $Last revision: 1.0 $  $Date: 2012/07/11  $  
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// [GK 2013.05.15] Original satsi setup.
//#define MAXDATA 7000
// [GK 2013.05.15] Standard version.
#define MAXDATA 70000
//-------------------------------------------------------------------------------------------------

// GK: 2013.03.03 Additional declarations to suppress warning messages;
int slfast_2D(char name_in[], int nx_in[], int ny_in[], float dipf[],
    float ddirf[], float rakef[], int nobs_t, float cwt);
void switcher(float ddir1, float dip1, float rake1, float *ddir2, float *dip2,
    float *rake2);

//-------------------------------------------------------------------------------------------------
int main(int argc, char **argv) {
  float *dip1, *ddir1, *rake1, *dip2, *ddir2, *rake2; /* focal mechanism data, input + auxiliary plane */
  float *dipf, *ddirf, *rakef; /* focal mechanism data, for inversion */
  int *x, *y; /* focal mechanism bins, input */
  int *nx, *ny; /* focal mechanism bins, for inversion */
  float frac; /* fraction correct plane picked */
  double seed(), myrand(), rand; /* random number generator */
  int nobs, ntries; /* number of observations, number bootstrap trials */
  float cwt; /* damping parameter */
  FILE *fpin; /* input file pointer */
  char name[40], headline[200]; /* character line */
  int n, i, j; /* dummy variables */
  int returnvalue = 0;

  rand = seed();

  dip1 = (float *) malloc(MAXDATA * sizeof(float));
  ddir1 = (float *) malloc(MAXDATA * sizeof(float));
  rake1 = (float *) malloc(MAXDATA * sizeof(float));
  dip2 = (float *) malloc(MAXDATA * sizeof(float));
  ddir2 = (float *) malloc(MAXDATA * sizeof(float));
  rake2 = (float *) malloc(MAXDATA * sizeof(float));
  dipf = (float *) malloc(MAXDATA * sizeof(float));
  ddirf = (float *) malloc(MAXDATA * sizeof(float));
  rakef = (float *) malloc(MAXDATA * sizeof(float));
  x = (int *) malloc(MAXDATA * sizeof(int));
  y = (int *) malloc(MAXDATA * sizeof(int));
  nx = (int *) malloc(MAXDATA * sizeof(int));
  ny = (int *) malloc(MAXDATA * sizeof(int));

  /* read in data */
  --argc; /* [PM 11.04.2013] Added for consistency with other programs*/
  ++argv; /* [PM 11.04.2013] Added for consistency with other programs*/
  if (argc != 4) /* [PM 11.04.2013] changed from 5 to 4 */
  {
    fprintf(stderr,"usage: bootmech_2D.exe file ntries frac damping\n");
    return -5001; // [GK 2013.03.08] Corrected output value, 11.04.2013 PM: GK - 5000.
      }

    /* ++argv; 11.04.2013: suppressed acordingly I hope :D*/
  sscanf(*argv, "%s", name);
  ++argv;
  sscanf(*argv, "%d", &ntries);
  ++argv;
  sscanf(*argv, "%f", &frac);
  ++argv;
  sscanf(*argv, "%f", &cwt);

  fpin = fopen(name, "r");
  if (fpin == NULL ) {
    fprintf(stderr,"unable to open %s.\n",name);
    return -5002; // [GK 2013.03.08] Corrected output value., 2013.04.11 PM : GK - 5000
      }

  fgets(headline, 200, fpin);
  nobs = 0;
  while (fscanf(fpin, "%d %d %f %f %f", &x[nobs], &y[nobs], &ddir1[nobs],
      &dip1[nobs], &rake1[nobs]) != EOF) {
    switcher(ddir1[nobs], dip1[nobs], rake1[nobs], &ddir2[nobs], &dip2[nobs],
        &rake2[nobs]);
    nobs++;
  }

  for (i = 0; i < ntries; ++i) {

    for (n = 0; n < nobs; ++n) {
      rand = myrand(&rand);
      j = (int) (rand * (double) nobs); /* eh mod 091991*/
      if (j == nobs)
        j = nobs - 1;
      if (myrand(&rand) >= frac) {
        nx[n] = x[j];
        ny[n] = y[j];
        ddirf[n] = ddir1[j];
        dipf[n] = dip1[j];
        rakef[n] = rake1[j];
      } else {
        nx[n] = x[j];
        ny[n] = y[j];
        ddirf[n] = ddir2[j];
        dipf[n] = dip2[j];
        rakef[n] = rake2[j];
      }
      if (dipf[n] < 0) {
        dipf[n] = -dipf[n];
        ddirf[n] = ddirf[n] + 180;
        rakef[n] = -rakef[n];
      }
      if (dipf[n] > 90) {
        dipf[n] = 180 - dipf[n];
        ddirf[n] = ddirf[n] + 180;
        rakef[n] = -rakef[n];
      }
      if (ddirf[n] > 360)
        ddirf[n] -= 360;
      if (ddirf[n] < 0)
        ddirf[n] += 360;
      if (rakef[n] > 360)
        rakef[n] -= 360;
      if (rakef[n] < 0)
        rakef[n] += 360;
    }
    printf("%d\n", i);
    returnvalue = slfast_2D(name, nx, ny, dipf, ddirf, rakef, nobs, cwt);
    if (returnvalue < 0)
      break;

    printf("%d\n", i);
  }
  fclose(fpin);

  free(dip1);
  free(ddir1);
  free(rake1);
  free(dip2);
  free(ddir2);
  free(rake2);
  free(dipf);
  free(ddirf);
  free(rakef);
  free(nx);
  free(ny);

  return returnvalue; // [GK 2013.03.08] Fixed output value for the purpose of MSATSI.
}
