//-------------------------------------------------------------------------------------------------
// SATSI_4D_TRADEOFF
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define TODEG 57.29577951
/* [GK 2013.05.15] Original satsi setup.
 #define MAXDATA 7200
 #define MAXX 14
 #define MAXY 14
 #define MAXZ 14
 #define MAXT 14
 #define MAXBOX 1764
 #define SRMAXBOX 42
 */
// [GK 2013.05.15] Standard version.
#define MAXDATA 70000
#define MAXX 50
#define MAXY 50
#define MAXZ 50
#define MAXT 50
#define MAXBOX 10000
#define SRMAXBOX 100
/* COORDINATES ARE EAST,NORTH,UP */

// [GK 2013.03.03] Additional declarations to suppress warning messages;
void sprsax(double sa[], int ija[], double x[], double b[], int m, int n);
void leasq_sparse(int a_ija[], double a_sa[], int d_ija[], double d_sa[], int m,
    int n, int p, double x[], double b[]);

int main(argc, argv)
  /* slickenside inversion program */
  int argc; /* argument count */
  char **argv; /* argument string */
{
  double ddir, dip, rake; /* focal mechanism data */
  int x, y, dep, time; /* focal mechanism bins */
  int nobs, nloc, nrows; /* number of observations, bins, rows */
  double cwt, twt; /* damping parameter, time/space damping ratio */
  double *diag_sa, *amat_sa, *d_sa, *dtemp_sa; /* inversion matrices in sparse matrix, */
  int *diag_ija, *amat_ija, *d_ija; /*      row-indexed form */
  int *lindex, nx, ny, nz, nt, index; /* book-keeping, which bins where in matrix */
  double stress[5 * MAXBOX]; /* stress tensor in vector form xx,xy,xz,yy,yz,zz */
  double *slick; /* slickenside vector elements vector */
  double n1, n2, n3; /* normal vector elements */
  //double norm[MAXDATA][3]; /* storage of n1,n2,n3 */
  double mech_misfit, mvar; /* mechanism misfit, model variance */
  double *stress_len; /* stress field model length (vector) */
  double *slick_pre; /* predicted slip vector */
  char line[80]; /* character line */
  FILE *fpin; /* input file pointer */
  FILE *fpout; /* output file pointer */
  int i, j, k, k2, ip, m, n, p, i2; /* dummy variables */
  double z, z2, z3, temp[5]; /* more dummy variables */

  n = 3 * MAXDATA + 25 * MAXBOX;
  diag_sa = (double *) malloc(n * sizeof(double));
  diag_ija = (int *) malloc(n * sizeof(int));
  n = 18 * MAXDATA + 5 * MAXBOX + 1;
  amat_sa = (double *) malloc(n * sizeof(double));
  amat_ija = (int *) malloc(n * sizeof(int));
  n = 60 * (MAXBOX - SRMAXBOX) + 1;
  d_sa = (double *) malloc(n * sizeof(double));
  d_ija = (int *) malloc(n * sizeof(int));
  dtemp_sa = (double *) malloc(n * sizeof(double));
  n = 20 * MAXBOX;
  stress_len = (double *) malloc(n * sizeof(double));
  n = 3 * MAXDATA;
  slick = (double *) malloc(n * sizeof(double));
  slick_pre = (double *) malloc(n * sizeof(double));
  lindex = (int *) malloc(MAXBOX * sizeof(int));

  /* get file pointers */
  --argc;
  ++argv;
  if (argc != 4) /* [PM 11.04.2013] Corrected from 3 to 4, changed from < to != */
  {
    printf(
        "usage: satsi_4D_tradeoff.exe data_file outfile damping time/space_damping\n");
    return -2001; /* [PM 11.04.2013] Changed from -1 to -2001 */
  }
  fpin = fopen(*argv, "r");
  if (fpin == NULL) {
    printf("unable to open %s.\n", *argv);
    return -2002; /* [PM 11.04.2013] Changed from -2 to -2002 */
  }
  ++argv;
  fpout = fopen(*argv, "a");
  if (fpout == NULL) {
    printf("unable to open %s.\n", *argv);
    return -2003; /* [PM 11.04.2013] Changed from -2 to -2002 */
  }
  ++argv;
  sscanf(*argv, "%lf", &cwt);
  ++argv;
  sscanf(*argv, "%lf", &twt);

  fgets(line, 80, fpin);

  for (i = 0; i < 3 * MAXDATA; i++)
    slick[i] = 0;
  for (i = 0; i < 3 * MAXDATA + 25 * MAXBOX; i++) {
    diag_sa[i] = 0;
    diag_ija[i] = 0;
  }

  /* loop to get data and make up equation */
  nobs = 0;
  nloc = 0;
  index = 0;
  while (fscanf(fpin, "%d %d %d %d %lf %lf %lf", &x, &y, &dep, &time, &ddir,
      &dip, &rake) != EOF) {
    j = 3 * nobs;
    z = ddir / TODEG;
    z2 = dip / TODEG;
    z3 = rake / TODEG;

    n1 = sin(z) * sin(z2); /* normal vector to fault plane */
    n2 = cos(z) * sin(z2);
    n3 = cos(z2);

    // [GK 2013.03.03] Commented out as it is of no use.
    //norm[nobs][0] = n1;
    //norm[nobs][1] = n2;
    //norm[nobs][2] = n3;

    /* slickenside vector calculation */
    slick[j] = -cos(z3) * cos(z) - sin(z3) * sin(z) * cos(z2);
    slick[j + 1] = cos(z3) * sin(z) - sin(z3) * cos(z) * cos(z2);
    slick[j + 2] = sin(z3) * sin(z2);

    /* find the matrix elements */
    k = -1;
    for (i = 0; i < nloc; i++)
      if (lindex[i]
          == MAXY * MAXZ * MAXT * x + MAXZ * MAXT * y + MAXT * dep + time)
        k = 5 * i;
    if (k == -1) {
      lindex[nloc] = MAXY * MAXZ * MAXT * x + MAXZ * MAXT * y + MAXT * dep
          + time;
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

    ++nobs;
    /* check to see if all possible data has been read */
    if (nobs == MAXDATA) {
      printf("NOT ALL DATA COULD BE READ.\n");
      break;
    }
  } /* end of data read loop */

  /* fill in holes in grid */
  for (nz = 0; nz < MAXZ; nz++)
    for (nt = 0; nt < MAXT; nt++)
      for (nx = 0; nx < MAXX; nx++) {
        i = MAXY;
        j = 0;
        for (ny = 0; ny < MAXY; ny++) {
          for (i2 = 0; i2 < nloc; i2++) {
            if ((lindex[i2]
                == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt)
                && (ny < i)) i = ny;
            if ((lindex[i2]
                == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt)
                && (ny > j)) j = ny;
          }
        }
        if (i < j) for (ny = i + 1; ny < j; ny++) {
          ip = 0;
          for (i2 = 0; i2 < nloc; i2++)
            if (lindex[i2] == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz
        + nt)
        ip = 1;
          if (ip == 0) {
            lindex[nloc] = MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny
                + MAXT * nz + nt;
            nloc++;
          }
        }
      }
  for (nz = 0; nz < MAXZ; nz++)
    for (nt = 0; nt < MAXT; nt++)
      for (ny = 0; ny < MAXY; ny++) {
        i = MAXX;
        j = 0;
        for (nx = 0; nx < MAXX; nx++) {
          for (i2 = 0; i2 < nloc; i2++) {
            if ((lindex[i2]
                == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt)
                && (nx < i)) i = nx;
            if ((lindex[i2]
                == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt)
                && (nx > j)) j = nx;
          }
        }
        if (i < j) for (nx = i + 1; nx < j; nx++) {
          ip = 0;
          for (i2 = 0; i2 < nloc; i2++)
            if (lindex[i2] == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz
        + nt)
        ip = 1;
          if (ip == 0) {
            lindex[nloc] = MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny
                + MAXT * nz + nt;
            nloc++;
          }
        }
      }
  for (nx = 0; nx < MAXX; nx++)
    for (ny = 0; ny < MAXY; ny++)
      for (nt = 0; nt < MAXT; nt++) {
        i = MAXZ;
        j = 0;
        for (nz = 0; nz < MAXZ; nz++) {
          for (i2 = 0; i2 < nloc; i2++) {
            if ((lindex[i2]
                == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt)
                && (nz < i)) i = nz;
            if ((lindex[i2]
                == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt)
                && (nz > j)) j = nz;
          }
        }
        if (i < j) for (nz = i + 1; nz < j; nz++) {
          ip = 0;
          for (i2 = 0; i2 < nloc; i2++)
            if (lindex[i2] == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz
        + nt)
        ip = 1;
          if (ip == 0) {
            lindex[nloc] = MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny
                + MAXT * nz + nt;
            nloc++;
          }
        }
      }
  for (nx = 0; nx < MAXX; nx++)
    for (ny = 0; ny < MAXY; ny++)
      for (nz = 0; nz < MAXZ; nz++) {
        i = MAXT;
        j = 0;
        for (nt = 0; nt < MAXT; nt++) {
          for (i2 = 0; i2 < nloc; i2++) {
            if ((lindex[i2]
                == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt)
                && (nt < i)) i = nt;
            if ((lindex[i2]
                == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt)
                && (nt > j)) j = nt;
          }
        }
        if (i < j) for (nt = i + 1; nt < j; nt++) {
          ip = 0;
          for (i2 = 0; i2 < nloc; i2++)
            if (lindex[i2] == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz
        + nt)
        ip = 1;
          if (ip == 0) {
            lindex[nloc] = MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny
                + MAXT * nz + nt;
            nloc++;
          }
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
      for (nz = 0; nz < MAXZ; nz++)
        for (nt = 0; nt < MAXT; nt++) {
          k = -1;
          for (i2 = 0; i2 < nloc; i2++)
            if (lindex[i2]
                == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt)
              k = 5 * i2;
          if (k > -1) {
            k2 = -1;
            for (i2 = 0; i2 < nloc; i2++)
              if (lindex[i2]
                  == MAXY * MAXZ * MAXT * (nx + 1) + MAXZ * MAXT * ny
                      + MAXT * nz + nt) k2 = 5 * i2;
            if ((nx < MAXX - 1) && (k2 > -1)) {
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
            k2 = -1;
            for (i2 = 0; i2 < nloc; i2++)
              if (lindex[i2]
                  == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * (ny + 1)
                      + MAXT * nz + nt) k2 = 5 * i2;
            if ((ny < MAXY - 1) && (k2 > -1)) {
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
            k2 = -1;
            for (i2 = 0; i2 < nloc; i2++)
              if (lindex[i2]
                  == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny
                      + MAXT * (nz + 1) + nt) k2 = 5 * i2;
            if ((nz < MAXZ - 1) && (k2 > -1)) {
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
            k2 = -1;
            for (i2 = 0; i2 < nloc; i2++)
              if (lindex[i2]
                  == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt
                      + 1) k2 = 5 * i2;
            if ((nt < MAXT - 1) && (k2 > -1)) {
              for (i = 0; i < 5; i++) {
                diag_ija[j + i] = index;
                if ((k + i) == (j + i))
                  diag_sa[j + i] = twt;
                else {
                  d_ija[index] = k + i;
                  d_sa[index] = twt;
                  index++;
                }
                if ((k2 + i) == (j + i))
                  diag_sa[j + i] = -twt;
                else {
                  d_ija[index] = k2 + i;
                  d_sa[index] = -twt;
                  index++;
                }
              }
              j += 5;
            }
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
  for (i = 0; i < d_ija[d_ija[0] - 1]; i++) {
    dtemp_sa[i] = d_sa[i];
    d_sa[i] = cwt * cwt * d_sa[i];
  }
  p = nrows;

  /* solve equations via linear least squares */
  leasq_sparse(amat_ija, amat_sa, d_ija, d_sa, m, n, p, stress, slick);

  /* compute RMS mechanism misfit */
  sprsax(amat_sa, amat_ija, stress, slick_pre, m, n);
  mech_misfit = 0;
  for (i = 0; i < n; i++) {
    mech_misfit += (slick_pre[i] - slick[i]) * (slick_pre[i] - slick[i]);
  }
  mech_misfit /= ((double) n);
  mech_misfit = sqrt(mech_misfit);

  /* compute stress field model length (variance) */
  sprsax(dtemp_sa, d_ija, stress, stress_len, m, p);
  mvar = 0;
  for (i = 0; i < p; i++) {
    mvar += (stress_len[i]) * (stress_len[i]);
  }
  mvar /= ((double) p);
  mvar = sqrt(mvar);

  fprintf(fpout, "%g %g\n", mech_misfit, mvar);

  free(diag_sa);
  free(diag_ija);
  free(amat_sa);
  free(amat_ija);
  free(d_sa);
  free(d_ija);
  free(dtemp_sa);
  free(stress_len);
  free(slick);
  free(slick_pre);
  free(lindex);

  return 0; // [GK 2013.03.08] Default return value.
}
