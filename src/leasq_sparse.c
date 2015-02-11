#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPS 1.0e-14
#define sub(I,J) (J+I*mp)
#define SWAP(a,b) itempri=(a);(a)=(b);(b)=itempri;

void leasq_sparse(a_ija, a_sa, d_ija, d_sa, m, n, p, x, b)
  /* finds the least squares solution of ax=b with damping matrix d */
  int a_ija[]; /* the coefficients matrix of size m by n*/
  double a_sa[]; /* in row-indexed form */
  int d_ija[]; /* the weighting matrix of size p by m*/
  double d_sa[]; /* in row-indexed form */
  double x[]; /* the solution vector of length m */
  double b[]; /* the constant vector of length n */
  int m, n, p; /* matrix sizes */

/* steps 0: atrans = a transpose */
/*          dtrans = d transpose */
/*       1: a2= atrans*a + dtrans*d */
/*       2: c= atrans*b       PM:  b?????? */
/*       3: solve a2x=c by conjagate gradient */

{
  int itol, itmax, *iter;
  double tol, *err;
  int nlen;
  int *a1_ija, *a2_ija, *a3_ija; /* square matrices of size m for internal use */
  double *a1_sa, *a2_sa, *a3_sa; /* in row-indexed form */
  int *atrans_ija; /* matrix of size m by n for internal use */
  double *atrans_sa; /* in row-indexed form */
  int *dtrans_ija; /* matrix of size m by p for internal use */
  double *dtrans_sa; /* in row-indexed form */
  double *c; /* vector of length m for internal use */

  c = (double *) malloc(m * sizeof(double));
  nlen = a_ija[a_ija[0] - 1];
  atrans_sa = (double *) malloc(nlen * sizeof(double));
  atrans_ija = (int *) malloc(nlen * sizeof(int));
  nlen = d_ija[d_ija[0] - 1];
  dtrans_sa = (double *) malloc(nlen * sizeof(double));
  dtrans_ija = (int *) malloc(nlen * sizeof(int));
  nlen = m * m;
  a1_sa = (double *) malloc(nlen * sizeof(double));
  a1_ija = (int *) malloc(nlen * sizeof(int));
  a2_sa = (double *) malloc(nlen * sizeof(double));
  a2_ija = (int *) malloc(nlen * sizeof(int));
  a3_sa = (double *) malloc(nlen * sizeof(double));
  a3_ija = (int *) malloc(nlen * sizeof(int));

  printf("m,n,p= %d %d %d\n", m, n, p);

  sprstp(a_sa, a_ija, atrans_sa, atrans_ija);
  sprstp(d_sa, d_ija, dtrans_sa, dtrans_ija);

  sprstm(atrans_sa, atrans_ija, atrans_sa, atrans_ija, 0.00001, nlen, a3_sa,
      a3_ija, n, m);
  sprstm(dtrans_sa, dtrans_ija, dtrans_sa, dtrans_ija, 0.00001, nlen, a1_sa,
      a1_ija, p, m);
  sprsax(atrans_sa, atrans_ija, b, c, n, m);
  sprsad(a1_sa, a1_ija, a3_sa, a3_ija, a2_sa, a2_ija, nlen);

  itol = 1;
  itmax = 100000;
  tol = 1.0e-6;
  iter = (int *) malloc(sizeof(int));
  err = (double *) malloc(sizeof(double));

  linbcg(a2_ija, a2_sa, m, c, x, itol, tol, itmax, iter, err);

  free(c);
  free(atrans_sa);
  free(atrans_ija);
  free(dtrans_sa);
  free(dtrans_ija);
  free(a1_sa);
  free(a1_ija);
  free(a2_sa);
  free(a2_ija);
  free(a3_sa);
  free(a3_ija);

  return;
}

void sprsax(double sa[], int ija[], double x[], double b[], int m, int n) {
  int i, k;
  for (i = 0; i < n; i++) {
    if (i < m)
      b[i] = sa[i] * x[i];
    else
      b[i] = 0;
    for (k = ija[i]; k < ija[i + 1]; k++)
      b[i] += sa[k] * x[ija[k]];
  }
  return;
}

void sprsad(sa, ija, sb, ijb, sc, ijc, nmax)
  double sa[], sb[], sc[];int ija[], ijb[], ijc[];int nmax; // GK: 03.03.2013 Missing argument.
{
  int i, na, nb, index;

  index = ija[0];
  for (i = 0; i < ija[0] - 1; i++) {
    sc[i] = sa[i] + sb[i];
    ijc[i] = index;
    na = ija[i];
    nb = ijb[i];
    while ((na < ija[i + 1]) || (nb < ijb[i + 1])) {
      if (na >= ija[i + 1])
        for (; nb < ijb[i + 1]; nb++, index++) {
          ijc[index] = ijb[nb];
          sc[index] = sb[nb];
        }
      else if (nb >= ijb[i + 1])
        for (; na < ija[i + 1]; na++, index++) {
          ijc[index] = ija[na];
          sc[index] = sa[na];
        }
      else if (ija[na] == ijb[nb]) {
        ijc[index] = ija[na];
        sc[index] = sa[na] + sb[nb];
        na++;
        nb++;
        index++;
      }
      else if (ija[na] > ijb[nb]) {
        ijc[index] = ijb[nb];
        sc[index] = sb[nb];
        nb++;
        index++;
      }
      else {
        ijc[index] = ija[na];
        sc[index] = sa[na];
        na++;
        index++;
      }
      if (index >= nmax) {
        printf("WARNING, sprsad exceeds vector length\n");
        break;
      }
    }
    ijc[ija[0] - 1] = index;
    sc[ija[0] - 1] = 0;
  }
}

void sprstm(sa, ija, sb, ijb, thresh, nmax, sc, ijc, m, n)
  double sa[], sb[], sc[], thresh;int ija[], ijb[], ijc[], nmax;int m, n; {
  int i, ijma, ijmb, j, k, ma, mb, mbb;
  double sum;

  ijc[0] = n + 1;
  k = n + 1;
  for (i = 0; i < ija[0] - 1; i++) {
    for (j = 0; j < ijb[0] - 1; j++) {
      if (i == j)
        sum = sa[i] * sb[j];
      else
        sum = 0.0;
      mb = ijb[j];
      for (ma = ija[i]; ma < ija[i + 1]; ma++) {
        ijma = ija[ma];
        if (ijma == j)
          sum += sa[ma] * sb[j];
        else {
          while (mb < ijb[j + 1]) {
            ijmb = ijb[mb];
            if (ijmb == i) {
              sum += sa[i] * sb[mb++];
              continue;
            }
            else if (ijmb < ijma) {
              mb++;
              continue;
            }
            else if (ijmb == ijma) {
              sum += sa[ma] * sb[mb++];
              continue;
            }
            break;
          }
        }
      }
      for (mbb = mb; mbb < ijb[j + 1]; mbb++) {
        if (ijb[mbb] == i) sum += sa[i] * sb[mbb];
      }
      if ((i == j) && (i < n))
        sc[i] = sum;
      else if (fabs(sum) > thresh) {
        if (k > nmax) {
          printf("WARNING, sprstm exceeds vector length\n");
          return;
        }
        sc[k] = sum;
        ijc[k++] = j;
      }
    }
    if (i < n) {
      ijc[i + 1] = k;
    }
  }
}

void sprstp(sa, ija, sb, ijb)
  double sa[], sb[];int ija[], ijb[]; {
  void indexx(int n, int arr[], int indx[]);
  int j, jl, jm, jp, ju, k, m, n2, noff, inc, iv;
  double v;

  n2 = ija[0];
  for (j = 0; j < n2 - 1; j++)
    sb[j] = sa[j];

  indexx(ija[n2 - 1] - ija[0], &ija[n2], &ijb[n2]);

  jp = -1;
  for (k = ija[0]; k < ija[n2 - 1]; k++) {
    m = ijb[k] + n2;
    sb[k] = sa[m];
    for (j = jp + 1; j <= ija[m]; j++)
      ijb[j] = k;
    jp = ija[m];
    jl = 0;
    ju = n2 - 1;
    while (ju - jl > 1) {
      jm = (ju + jl) / 2;
      if (ija[jm] > m)
        ju = jm;
      else
        jl = jm;
    }
    ijb[k] = jl;
  }
  for (j = jp + 1; j < n2; j++)
    ijb[j] = ija[n2 - 1];
  for (j = 0; j < n2 - 1; j++) {
    jl = ijb[j + 1] - ijb[j];
    noff = ijb[j] - 1;
    inc = 1;
    do {
      inc *= 3;
      inc++;
    } while (inc <= jl);
    do {
      inc /= 3;
      for (k = noff + inc + 1; k <= noff + jl; k++) {
        iv = ijb[k];
        v = sb[k];
        m = k;
        while (ijb[m - inc] > iv) {
          ijb[m] = ijb[m - inc];
          sb[m] = sb[m - inc];
          m -= inc;
          if (m - noff <= inc) break;
        }
        ijb[m] = iv;
        sb[m] = v;
      }
    } while (inc > 1);
  }
}

void indexx(int n, int arr[], int indx[]) {
  int i, indxt, ir, itempri, j, k, l, jstack, *istack, a;
  ir = n - 1;
  l = 0;
  jstack = 0;
  istack = (int *) malloc(n * sizeof(int));
  for (j = 0; j < n; j++)
    indx[j] = j;
  for (;;) {
    if (ir - l < 7) {
      for (j = l + 1; j <= ir; j++) {
        indxt = indx[j];
        a = arr[indxt];
        for (i = j - 1; i >= l; i--) {
          if (arr[indx[i]] <= a) break;
          indx[i + 1] = indx[i];
        }
        indx[i + 1] = indxt;
      }
      if (jstack == 0) break;
      ir = istack[jstack--];
      l = istack[jstack--];
    }
    else {
      k = (l + ir) >> 1;
      SWAP(indx[k], indx[l + 1])
      if (arr[indx[l]] > arr[indx[ir]]) {
        SWAP(indx[l], indx[ir])
      }
      if (arr[indx[l + 1]] > arr[indx[ir]]) {
        SWAP(indx[l + 1], indx[ir])
      }
      if (arr[indx[l]] > arr[indx[l + 1]]) {
        SWAP(indx[l], indx[l + 1])
      }
      i = l + 1;
      j = ir;
      indxt = indx[l + 1];
      a = arr[indxt];
      for (;;) {
        do
          i++;
        while (arr[indx[i]] < a);
        do
          j--;
        while (arr[indx[j]] > a);
        if (j < i) break;
        SWAP(indx[i], indx[j])
      }
      indx[l + 1] = indx[j];
      indx[j] = indxt;
      jstack += 2;
      if (ir - i + 1 >= j - l) {
        istack[jstack] = ir;
        istack[jstack - 1] = i;
        ir = j - 1;
      }
      else {
        istack[jstack] = j - 1;
        istack[jstack - 1] = l;
        l = i;
      }
    }
  }
  free(istack);
}

void linbcg(a_ija, a_sa, n, b, x, itol, tol, itmax, iter, err)
  int a_ija[], n, itol, itmax, *iter;double a_sa[], b[], x[], tol, *err; {
  void asolve(int n, double b[], double x[], double a_sa[]);
  double snrm(int n, double sx[], int itol);
  int j;
  double ak, akden, bk, bkden, bknum, bnrm = 0;
  double *p, *pp, *r, *rr, *z, *zz;

  p = (double *) malloc(n * sizeof(double));
  pp = (double *) malloc(n * sizeof(double));
  r = (double *) malloc(n * sizeof(double));
  rr = (double *) malloc(n * sizeof(double));
  z = (double *) malloc(n * sizeof(double));
  zz = (double *) malloc(n * sizeof(double));
  for (j = 0; j < n; j++)
    x[j] = 0.5;
  *iter = 0;
  sprsax(a_sa, a_ija, x, r, n, n);
  for (j = 0; j < n; j++) {
    r[j] = b[j] - r[j];
    rr[j] = r[j];
  }
  if (itol == 1) {
    bnrm = snrm(n, b, itol);
    asolve(n, r, z, a_sa);
  }
  else if (itol == 2) {
    asolve(n, b, z, a_sa);
    bnrm = snrm(n, z, itol);
    asolve(n, r, z, a_sa);
  }
  while (*iter <= itmax) {
    ++(*iter);
    asolve(n, rr, zz, a_sa);
    for (bknum = 0.0, j = 0; j < n; j++)
      bknum += z[j] * rr[j];
    if (*iter == 1) {
      for (j = 0; j < n; j++) {
        p[j] = z[j];
        pp[j] = zz[j];
      }
    }
    else {
      bk = bknum / bkden;
      for (j = 0; j < n; j++) {
        p[j] = bk * p[j] + z[j];
        pp[j] = bk * pp[j] + zz[j];
      }
    }
    bkden = bknum;
    sprsax(a_sa, a_ija, p, z, n, n);
    for (akden = 0.0, j = 0; j < n; j++)
      akden += z[j] * pp[j];
    ak = bknum / akden;
    sprsax(a_sa, a_ija, pp, zz, n, n);
    for (j = 0; j < n; j++) {
      x[j] += ak * p[j];
      r[j] -= ak * z[j];
      rr[j] -= ak * zz[j];
    }
    asolve(n, r, z, a_sa);
    if (itol == 1)
      *err = snrm(n, r, itol) / bnrm;
    else if (itol == 2) *err = snrm(n, z, itol) / bnrm;
    //printf("iter=%4d err=%12.6f\n", *iter, *err); // 2013.06.17 GK Suppressed output.
    if (*err <= tol) break;
  }
  free(p);
  free(pp);
  free(r);
  free(rr);
  free(z);
  free(zz);
}

void asolve(int n, double b[], double x[], double sa[]) {
  int i;
  for (i = 0; i < n; i++)
    x[i] = (sa[i] != 0.0 ? b[i] / sa[i] : b[i]);
}

double snrm(int n, double sx[], int itol) {
  int i;
  double ans;

  ans = 0.0;
  for (i = 0; i < n; i++)
    ans += sx[i] * sx[i];
  return sqrt(ans);
}

