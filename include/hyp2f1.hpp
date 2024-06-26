
// Standard Template Library
#include <complex>
#include <cstdlib>
#include <iostream>
#include <cmath>

using namespace std;

#define SIGN(a) (((a) < 0) ? (-1) : (1))

#include "AEAE/complex_functions.H"
#include "AEAE/hyp_2F1.cpp"



double hyp2f1a(double beta, double z1, double z2)
{
  double x;
  if (beta >= 1)
  {
    return (hyp_2F1(std::complex<double>(1.0,                      0.0),
                    std::complex<double>(1.0 / beta,               0.0),
                    std::complex<double>(1.0 + (1.0 / beta),       0.0),
                    std::complex<double>(-std::pow(z1 / z2, beta), 0.0) )).real();
  }
  else if (beta > 0)
  {
    x = (hyp_2F1(std::complex<double>(1.0,                      0.0),
                 std::complex<double>(1.0 / beta,               0.0),
                 std::complex<double>(1.0 + (1.0 / beta),       0.0),
                 std::complex<double>(-std::pow(z1, beta) / z2, 0.0) )).real();
    if (isnan(x))
    {
      return 1.0 / (1.0 + 1.0 / z2);
    }
    else
    {
      return x;
    }
  }
  else
  {
    return 1.0 / (1.0 + 1.0 / z2);
  }
}

double hyp2f1b(double beta, double z1, double z2)
{
  double x;
  if (beta >= 1)
  {
    return (hyp_2F1(std::complex<double>(1.0,                      0.0),
                    std::complex<double>(2.0 / beta,               0.0),
                    std::complex<double>(1.0 + (2.0 / beta),       0.0),
                    std::complex<double>(-std::pow(z1 / z2, beta), 0.0) )).real();
  }
  else if (beta > 0)
  {
    x = (hyp_2F1(std::complex<double>(1.0,                      0.0),
                 std::complex<double>(2.0 / beta,               0.0),
                 std::complex<double>(1.0 + (2.0 / beta),       0.0),
                 std::complex<double>(-std::pow(z1, beta) / z2, 0.0) )).real();
    if (isnan(x))
    {
      return 1.0 / (1.0 + 1.0 / z2);
    }
    else
    {
      return x;
    }
  }
  else
  {
    return 1.0 / (1.0 + 1.0 / z2);
  }
}
