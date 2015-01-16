//-------------------------------------------------------------------------------------------------
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//-------------------------------------------------------------------------------------------------

/* COORDINATES ARE EAST,NORTH,UP */

//-------------------------------------------------------------------------------------------------
// SLFAST_2D
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
#define TODEG 57.29577951
//-------------------------------------------------------------------------------------------------
// [GK 2013.05.15] Original SATSI setup.
//#define MAXDATA 7000
//#define MAXX 56
//#define MAXY 56
//#define MAXBOX 900
//#define SRMAXBOX 30
// [GK 2013.05.15] Standard version.
//#define MAXDATA 70000
//#define MAXX 50
//#define MAXY 50
//#define MAXBOX 10000
//#define SRMAXBOX 100
// [GK 2013.05.15] Extended SATSI setup.
#define MAXDATA 70000
#define MAXX 200
#define MAXY 200
#define MAXBOX 10000
#define SRMAXBOX 100
//-------------------------------------------------------------------------------------------------

// [GK 2013.03.03] Additional declarations to suppress warning messages.
void dirplg(double e, double n, double u, double *pdir, double *pplg);
void eigen(double a[3][3], double lam[3], double q[3][3]);
void leasq_sparse(int a_ija[], double a_sa[], int d_ija[], double d_sa[], int m,
    int n, int p, double x[], double b[]);

/*slfast(name_in) slickenside inversion program */
//-------------------------------------------------------------------------------------------------
int slfast_2D(char name_in[], int nx_in[], int ny_in[], float dipf[],
    float ddirf[], float rakef[], int nobs_t, float cwt) {
  int x, y; /* focal mechanism bins */
  int nobs, nloc, nrows, nlocfill; /* number of observations, bins, rows */
  double *diag_sa, *amat_sa, *d_sa; /* inversion matrices in sparse matrix, */
  int *diag_ija, *amat_ija, *d_ija; /*      row-indexed form */
  int *loclist, nx, ny, index; /* book-keeping, which bins where in matrix*/
  double stress[5 * MAXBOX]; /* stress tensor in vector form xx,xy,xz,yy,yz,zz */
  double strten[3][3]; /* stress tensor in tensor form */
  double *slick; /* slickenside vector elements vector */
  double n1, n2, n3; /* normal vector elements */
  double lam[3]; /* eigenvalues */
  double vecs[3][3]; /* eigenvectors */
  double dev_stress; /* deviatoric stress mag */
  float phi; /* stress ratio */
  char name[40]; /* output file name */
  FILE *fpout; /* output file pointer */
  int i, j, k, k2, m, n, p; /* dummy variables */
  double z, z2, z3, temp[5]; /* more dummy variables */

  n = 3 * MAXDATA + 25 * MAXBOX;
  diag_sa = (double *) malloc(n * sizeof(double));
  diag_ija = (int *) malloc(n * sizeof(int));
  n = 18 * MAXDATA + 5 * MAXBOX + 1;
  amat_sa = (double *) malloc(n * sizeof(double));
  amat_ija = (int *) malloc(n * sizeof(int));
  n = 30 * (MAXBOX - SRMAXBOX) + 1;
  d_sa = (double *) malloc(n * sizeof(double));
  d_ija = (int *) malloc(n * sizeof(int));
  n = 3 * MAXDATA;
  slick = (double *) malloc(n * sizeof(double));
  loclist = (int *) malloc(MAXX * MAXY * sizeof(int));

  sprintf(name, "%s.slboot", name_in);
  fpout = fopen(name, "a");
  if (fpout == NULL ) {
    printf("unable to open %s.\n", name);
    return -7001 /* 11.04.2013 PM: From -1 to -7001*/;
  }

  for (i = 0; i < 3 * MAXDATA; i++)
    slick[i] = 0.0;
  for (i = 0; i < 3 * MAXDATA + 25 * MAXBOX; i++) {
    diag_sa[i] = 0;
    diag_ija[i] = 0;
  }
  for (i = 0; i < MAXX; i++)
    for (j = 0; j < MAXY; j++)
      loclist[MAXY * i + j] = -999;

  /* loop to get data and make up equation */
  nloc = 0;
  index = 0;
  for (nobs = 0; nobs < nobs_t; nobs++) {
    x = nx_in[nobs];
    y = ny_in[nobs];
    i = nobs;
    j = 3 * nobs;
    z = ddirf[nobs] / TODEG;
    z2 = dipf[nobs] / TODEG;
    z3 = rakef[nobs] / TODEG;

    n1 = sin(z) * sin(z2); /* normal vector to fault plane */
    n2 = cos(z) * sin(z2);
    n3 = cos(z2);

    /* slickenside vector calculation */
    slick[j] = -cos(z3) * cos(z) - sin(z3) * sin(z) * cos(z2);
    slick[j + 1] = cos(z3) * sin(z) - sin(z3) * cos(z) * cos(z2);
    slick[j + 2] = sin(z3) * sin(z2);

    /* find the matrix elements */
    if (loclist[MAXY * x + y] > -999)
      k = 5 * loclist[MAXY * x + y];
    else {
      loclist[MAXY * x + y] = nloc;
      k = 5 * nloc;
      nloc++;
    }

    temp[0] = n1 - n1 * n1 * n1 + n1 * n3 * n3;
    temp[1] = n2 - 2. * n1 * n1 * n2;
    temp[2] = n3 - 2. * n1 * n1 * n3;
    temp[3] = -n1 * n2 * n2 + n1 * n3 * n3;
    temp[4] = -2. * n1 * n2 * n3;
    diag_ija[j] = index;
    for (i = 0; i < 5; i++) {
      if ((k + i) == j)
        diag_sa[j] = temp[i];
      else {
        amat_ija[index] = k + i;
        amat_sa[index] = temp[i];
        index++;
      }
    }

    temp[0] = -n2 * n1 * n1 + n2 * n3 * n3;
    temp[1] = n1 - 2. * n1 * n2 * n2;
    temp[2] = -2. * n1 * n2 * n3;
    temp[3] = n2 - n2 * n2 * n2 + n2 * n3 * n3;
    temp[4] = n3 - 2. * n2 * n2 * n3;
    diag_ija[j + 1] = index;
    for (i = 0; i < 5; i++) {
      if ((k + i) == (j + 1))
        diag_sa[j + 1] = temp[i];
      else {
        amat_ija[index] = k + i;
        amat_sa[index] = temp[i];
        index++;
      }
    }

    temp[0] = -n3 * n1 * n1 - n3 + n3 * n3 * n3;
    temp[1] = -2. * n1 * n2 * n3;
    temp[2] = n1 - 2. * n1 * n3 * n3;
    temp[3] = -n3 * n2 * n2 - n3 + n3 * n3 * n3;
    temp[4] = n2 - 2. * n2 * n3 * n3;
    diag_ija[j + 2] = index;
    for (i = 0; i < 5; i++) {
      if ((k + i) == (j + 2))
        diag_sa[j + 2] = temp[i];
      else {
        amat_ija[index] = k + i;
        amat_sa[index] = temp[i];
        index++;
      }
    }
    /* check to see if all possible data has been read */
    if (nobs == MAXDATA) {
      printf("NOT ALL DATA COULD BE READ.\n");
      break;
    }
  } /* end of data read loop */
  /* fill in holes in grid */
  nlocfill = nloc;
  for (nx = 0; nx < MAXX; nx++) {
    i = 0;
    while ((i < MAXY) && (loclist[MAXY * nx + i] == -999))
      i++;
    j = MAXY - 1;
    while ((j > -1) && (loclist[MAXY * nx + j] == -999))
      j--;
    if (i < j)
      for (ny = i + 1; ny < j; ny++)
        if (loclist[MAXY * nx + ny] == -999) {
          loclist[MAXY * nx + ny] = nloc;
          nloc++;
        }
  }
  for (ny = 0; ny < MAXY; ny++) {
    i = 0;
    while ((i < MAXX) && (loclist[MAXY * i + ny] == -999))
      i++;
    j = MAXX - 1;
    while ((j > -1) && (loclist[MAXY * j + ny] == -999))
      j--;
    if (i < j)
      for (nx = i + 1; nx < j; nx++)
        if (loclist[MAXY * nx + ny] == -999) {
          loclist[MAXY * nx + ny] = nloc;
          nloc++;
        }
  }
  /* fill in diagonal */
  nrows = 3 * nobs;
  m = 5 * nloc;
  if (nrows < m) {
    for (i = nrows; i < m; i++)
      diag_ija[i] = index;
    nrows = m;
  }
  for (i = index - 1; i >= 0; i--) {
    amat_ija[i + nrows + 1] = amat_ija[i];
    amat_sa[i + nrows + 1] = amat_sa[i];
  }
  for (i = 0; i < nrows; i++) {
    amat_ija[i] = diag_ija[i] + nrows + 1;
    amat_sa[i] = diag_sa[i];
  }
  amat_ija[nrows] = index + nrows + 1;
  amat_sa[nrows] = 0;
  n = nrows;

  /* set up smoothing constraints */
  for (i = 0; i < 3 * MAXDATA + 25 * MAXBOX; i++) {
    diag_sa[i] = 0;
    diag_ija[i] = 0;
  }
  index = 0;
  j = 0;
  for (nx = 0; nx < MAXX; nx++)
    for (ny = 0; ny < MAXY; ny++)
      if (loclist[MAXY * nx + ny] > -999) {
        k = 5 * loclist[MAXY * nx + ny];
        if ((nx < MAXX - 1) && (loclist[MAXY * (nx + 1) + ny] > -999)) {
          k2 = 5 * loclist[MAXY * (nx + 1) + ny];
          for (i = 0; i < 5; i++) {
            diag_ija[j + i] = index;
            if ((k + i) == (j + i))
              diag_sa[j + i] = 1;
            else {
              d_ija[index] = k + i;
              d_sa[index] = 1;
              index++;
            }
            if ((k2 + i) == (j + i))
              diag_sa[j + i] = -1;
            else {
              d_ija[index] = k2 + i;
              d_sa[index] = -1;
              index++;
            }
          }
          j += 5;
        }
        if ((ny < MAXY - 1) && (loclist[MAXY * nx + ny + 1] > -999)) {
          k2 = 5 * loclist[MAXY * nx + ny + 1];
          for (i = 0; i < 5; i++) {
            diag_ija[j + i] = index;
            if ((k + i) == (j + i))
              diag_sa[j + i] = 1;
            else {
              d_ija[index] = k + i;
              d_sa[index] = 1;
              index++;
            }
            if ((k2 + i) == (j + i))
              diag_sa[j + i] = -1;
            else {
              d_ija[index] = k2 + i;
              d_sa[index] = -1;
              index++;
            }
          }
          j += 5;
        }
      }
  nrows = j;
  if (5 * nloc > nrows) {
    for (i = nrows; i < 5 * nloc; i++)
      diag_ija[i] = index;
    nrows = 5 * nloc;
  }
  for (i = index - 1; i >= 0; i--) {
    d_ija[i + nrows + 1] = d_ija[i];
    d_sa[i + nrows + 1] = d_sa[i];
  }
  for (i = 0; i < nrows; i++) {
    d_ija[i] = diag_ija[i] + nrows + 1;
    d_sa[i] = diag_sa[i];
  }
  d_ija[nrows] = index + nrows + 1;
  d_sa[nrows] = 0;
  for (i = 0; i < d_ija[d_ija[0] - 1]; i++)
    d_sa[i] = cwt * cwt * d_sa[i];
  p = nrows;

  /* solve equations via linear least squares */
  leasq_sparse(amat_ija, amat_sa, d_ija, d_sa, m, n, p, stress, slick);

  for (nx = 0; nx < MAXX; nx++)
    for (ny = 0; ny < MAXY; ny++)
      if ((loclist[MAXY * nx + ny] > -999)
          && (loclist[MAXY * nx + ny] < nlocfill)) {
        k = 5 * loclist[MAXY * nx + ny];
        strten[0][0] = stress[k];
        strten[0][1] = stress[k + 1];
        strten[1][0] = stress[k + 1];
        strten[0][2] = stress[k + 2];
        strten[2][0] = stress[k + 2];
        strten[1][1] = stress[k + 3];
        strten[1][2] = stress[k + 4];
        strten[2][1] = stress[k + 4];
        strten[2][2] = -(stress[k] + stress[k + 3]);
        eigen(strten, lam, vecs);
        i = 1;
        while (i) {
          i = 0;
          for (j = 0; j < 2; ++j) {
            if (lam[j] > lam[j + 1]) {
              z = lam[j];
              lam[j] = lam[j + 1];
              lam[j + 1] = z;
              z = vecs[0][j];
              vecs[0][j] = vecs[0][j + 1];
              vecs[0][j + 1] = z;
              z = vecs[1][j];
              vecs[1][j] = vecs[1][j + 1];
              vecs[1][j + 1] = z;
              z = vecs[2][j];
              vecs[2][j] = vecs[2][j + 1];
              vecs[2][j + 1] = z;
              i = 1;
            }
          }
        }
        dev_stress = lam[2] - lam[0];
        for (i = 0; i < 3; i++)
          for (j = 0; j < 3; j++)
            strten[i][j] /= dev_stress;
        fprintf(fpout, "%d %d %g %g %g %g %g %g\n", nx, ny, strten[0][0],
            strten[0][1], strten[0][2], strten[1][1], strten[1][2],
            strten[2][2]);
        if (lam[0] != lam[2]) {
          phi = (lam[1] - lam[2]) / (lam[0] - lam[2]);
          fprintf(fpout, "%d %d %g ", nx, ny, phi);
        } else
          fprintf(fpout, "2. "); /* error flag */
        for (i = 0; i < 3; ++i) {
          dirplg(vecs[0][i], vecs[1][i], vecs[2][i], &z, &z2);
          fprintf(fpout, "%5.1f  %5.1f  ", z, z2);
        }
        fprintf(fpout, "\n");
      }
  free(diag_sa);
  free(diag_ija);
  free(amat_sa);
  free(amat_ija);
  free(d_sa);
  free(d_ija);
  free(slick);
  free(loclist);

  fclose(fpout);

  return 0;
}
