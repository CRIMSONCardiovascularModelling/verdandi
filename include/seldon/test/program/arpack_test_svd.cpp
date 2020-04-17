// g++ -o arpack_test_svd arpack_test_svd.cpp -larpack -lblas -I../../

#include <cmath>

#define SELDON_WITH_BLAS
#define SELDON_WITH_ARPACK
#define ARPACK_WITH_ABORT
#include "Seldon.hxx"
using namespace Seldon;


// Compute y <- A' * w
void atv (int m, int n, float* w, float* y)
{  
  // w[m], y[n]
  int i, j;
  float h, k, s, t;
  float one = 1., zero = 0.;

  h = one / float(m + 1);
  k = one / float(n + 1);

  for (int i = 0; i < n; i++)
    y[i] = zero;
  
  t = zero;

  for (int j = 0; j < n; j++)
    {
      t = t + k;
      s = zero;
      for (int i = 0; i < j; i++)
	{
	  s = s + h;
	  y[j] = y[j] + k * s * (t - one) * w[i];
	}
      for (int i = j; i < m; i++)
	{
	  s = s + h;
	  y[j] = y[j] + k * t * (s - one) * w[i];
	}
    }
}


// Computes  w <- A*x
void av(int m, int n, float* x, float* w)
{
  // x[n], w[m]
  int i, j;
  float one = 1., zero = 0.;
  float h, k, s, t;

  h = one / float(m + 1);
  k = one / float(n + 1);

  for (int i = 0; i < m; i++)
    w[i] = zero;
  t = zero;
  
  for (int j = 0; j < n; j++)
    {
      t = t + k;
      s = zero;

      for (int i = 0; i < j; i++)
	{
	  s = s + h;
	  w[i] += k * s * (t - one) * x[j];
	}
      for (int i = j; i < m; i++)
	{
	  s = s + h;
	  w[i] += k * t * (s - one) * x[j];
	}
    }
}


int main()
{
  // The following sets dimensions for this problem.
  int m = 500;
  int n = 100;
  int nev = 4;
  int ncv = 20;
  int maxit = n;
  float tol = 0.;
  string solver_type = "symmetric";
  int mode = 1;
  string which = "LM";
  char bmat = 'I';
  char HowMny = 'A';
  bool with_arpack_verbose = false;

  int ldu = m;
  float* ax = new float[m];
  float* s = new float[ncv];
  float* d = new float[ncv];
  float* u = new float[ldu * nev]; // 2D

  ArpackSolver<float, float> S;
  S.Init(n, nev, ncv, maxit, tol, solver_type, mode,
	 which, bmat, HowMny, with_arpack_verbose);

  while (S.Continue())
    {
      int ido = S.GetReverseCommunicationFlag();
      if (ido == 1 || ido == -1)
	{
	  av(m, n, S.GetFirstWorkVector(), ax);
	  atv(m, n, ax, S.GetSecondWorkVector());
	}
    }

  bool success = S.Finish();

  if (success)
    {
      d = new float[S.GetConvergedNumber()];
      for (int j = 0; j < S.GetConvergedNumber(); j++)
	{
	  av(m, n, S.GetEigenVector(j), ax);
	  cblas_scopy(m, ax, 1, u + j * ldu, 1);
	  float temp = 1. / cblas_snrm2(m, u + j * ldu, 1);
	  cblas_sscal(m, temp, u + j * ldu, 1);
	  
	  s[j] = sqrt(S.GetEigenValue(j));
	  cblas_saxpy(m, -s[j], u + j * ldu, 1, ax, 1);
	  d[j] = cblas_snrm2(m, ax, 1);
	}

      cout << "Ritz Values" << "\t\t" << "Relative Residuals" << endl;
      for (int i = 0; i < S.GetConvergedNumber(); i++)
	cout << s[i] << "\t\t\t" << d[i] << endl;
    }
   
  delete[] u;
  delete[] s;
  delete[] d;
  delete[] ax;

  return 0;
}
