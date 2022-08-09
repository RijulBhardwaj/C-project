#include <stdio.h>
#include <math.h>
#define MMAX 2001
#define M1MX 1001
#define NMAX   19 

void main()
/* This program solves the time dependent temperature field
   around a nuclear waste rod in a two-dimensional model.
    */
{
int i,j,m,n;
double ht,tc,xk,ra,rb,h,h2,t0,s0,cs,a0,b0,c0,r,t;
double a[NMAX],b[NMAX],c[NMAX],y[NMAX],g[NMAX],w[NMAX],v[NMAX],u[NMAX];
double s[MMAX][NMAX],x[MMAX][NMAX+1];


m  = MMAX;
n  = NMAX;
ht = 1.0/(M1MX-1);
tc = 1;
xk = 3153.6;
ra = 25;
rb = 100;
h  = rb/(n+1);
h2 = h*h;
t0 = 10;
s0 = ht*xk*t0/(ra*ra);
cs = ht*xk/h2;

for (i = 0; i < m; ++i)
  {
  for (j = 0; j <= n; ++j)
    {
    x[i][j] = 0;
    }
  }

for (i = 0; i < n; ++i)
  {
  r = h*(i+1);
  a[i] =  2*(1+cs)*r;
  b[i] = -(1+0.5/(i+1))*cs*r;
  c[i] = -(1-0.5/(i+1))*cs*r;
  if (r < ra)
    {
    s[0][i] = s0*r;
    }
  else
    {
    s[0][i] = 0;
    }
  }

/* Assign the source of the radiation heat */
for (i = 1; i < m; ++i)
  {
  t = ht*i;
  for (j = 0; j < n; ++j)
    {
    r = h*(j+1);
    if (r < ra)
      {
      s[i][j] = r*s0/exp(t/tc);
      }
    else
      {
      s[i][j] = 0;
      }
    }
  }
for (i = 1; i < m; ++i)
  {
/* Find the values from the last time step */

  a0 = 2*(1-cs)*h;
  b0 = (1+0.5)*cs*h;
  c0 = (1-0.5)*cs*h;
  g[0] = a0*x[i-1][0]+b0*x[i-1][1]
    +c0*(4*x[i-1][0]-x[i-1][1])/3+s[i-1][0]+s[i][0];
  for (j = 1; j < n; ++j)
    {
    r = h*(j+1);
    a0 = 2*(1-cs)*r;
    b0 = (1+0.5/(j+1))*cs*r;
    c0 = (1-0.5/(j+1))*cs*r;
    g[j] = a0*x[i-1][j]+b0*x[i-1][j+1]+c0*x[i-1][j-1]+s[i-1][j]+s[i][j];
    }
/* Find the elements in L and U */

  w[0] = a[0];
  v[0] = c[0];
  u[0] = b[0]/w[0];
  for (j = 1; j < n; ++j)
    {
    v[j] = c[j];
    w[j] = a[j]-v[j-1]*u[j-1];
    u[j] = b[j]/w[j];
    }
/* Find the solution of the temperature */

  y[0] = g[0]/w[0];
  for (j = 1; j < n; ++j)
    {
    y[j] = (g[j]-v[j-1]*y[j-1])/w[j];
    }
  x[i][n-1] = y[n-1];
  for (j = n-2; j >= 0; j = j-1)
    {
    x[i][j] = y[j]-u[j]*x[i][j+1];
    }
  }
i = m-1;
for (j = 0; j <= n; ++j)
  {
  r = (j+1)*h;
  printf("%16.8lf %16.8lf\n", r,x[i][j]);
  }
}