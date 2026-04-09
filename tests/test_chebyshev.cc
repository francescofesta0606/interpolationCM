#include "Interpolation/interpolation.hh"
using namespace Interpolation;

int test_exp_function()
{
   const auto testfunction = [](double x) {
      return exp(2.0 * x);
   };
   const auto testfunction_d = [](double x) {
      return 2.0 * exp(2.0 * x);
   };

   const double p = 20;
   Chebyshev::StandardGrid grid(p);

   int errcode = 1;

   if (grid._p != p) return errcode;
   if (grid._betaj.size() != p + 1) return errcode;
   if (grid._tj.size() != p + 1) return errcode;
   if (grid._Dij.size() != p + 1) return errcode;
   for (size_t i = 0; i <= p; i++) {
      if (grid._Dij[i].size() != p + 1) return errcode;
   }

   for (size_t i = 0; i <= p; i++) {
      if (std::abs(grid._tj[i] - cos(i * M_PI / ((double)p))) > 1.0e-14) return errcode;
   }

   auto v = grid.discretize(testfunction);

   double xmin = -1;
   double xmax = 1;
   size_t n    = 1e4;
   double dx   = (xmax - xmin) / ((double)n - 1);

   double m = 0, m_d = 0.;

   std::FILE *fptr = std::fopen("StandardGrid_interpolation_chebyshev_exp.dat", "w");
   for (size_t i = 0; i < n; i++) {
      const double x     = xmin + i * dx;
      const double exact = testfunction(x);
      const double inter = grid.interpolate(x, v, 0, grid._p);

      const double exact_d = testfunction_d(x);
      const double inter_d = grid.interpolate_der(x, v, 0, grid._p);

      m   = std::max(m, std::abs(exact - inter));
      m_d = std::max(m_d, std::abs(exact_d - inter_d));

      std::fprintf(fptr, "%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", x, exact, inter, exact_d, inter_d);
   }
   for (size_t i = 0; i <= grid._p; i++) {
      const double x     = grid._tj[i];
      const double exact = testfunction(x);
      const double inter = grid.interpolate(x, v, 0, grid._p);

      const double exact_d = testfunction_d(x);
      const double inter_d = grid.interpolate_der(x, v, 0, grid._p);

      m   = std::max(m, std::abs(exact - inter));
      m_d = std::max(m_d, std::abs(exact_d - inter_d));
   }
   std::fclose(fptr);

   std::printf("Max difference in interpolation: %.6e\n", m);
   std::printf("Max difference in interpolation derivatives: %.6e\n", m_d);
   if (m > 1.0e-12) return errcode;
   if (m_d > 1.0e-11) return errcode;
   return 0;
}

int test_runge_function()
{
   const auto testfunction = [](double x) {
      return 1.0 / (1. + 25 * x * x);
   };
   const auto testfunction_d = [](double x) {
      return -50. * x / std::pow(1. + 25 * x * x, 2);
   };

   const double p = 105;
   Chebyshev::StandardGrid grid(p);

   int errcode = 2;

   auto v = grid.discretize(testfunction);

   double xmin = -1;
   double xmax = 1;
   size_t n    = 1e4;
   double dx   = (xmax - xmin) / ((double)n - 1);

   double m = 0, m_d = 0.;

   std::FILE *fptr = std::fopen("StandardGrid_interpolation_chebyshev_runge.dat", "w");
   for (size_t i = 0; i < n; i++) {
      const double x     = xmin + i * dx;
      const double exact = testfunction(x);
      const double inter = grid.interpolate(x, v, 0, grid._p);

      const double exact_d = testfunction_d(x);
      const double inter_d = grid.interpolate_der(x, v, 0, grid._p);

      m   = std::max(m, std::abs(exact - inter));
      m_d = std::max(m_d, std::abs(exact_d - inter_d));

      std::fprintf(fptr, "%.16e\t%.16e\t%.16e\t%.16e\t%.16e\n", x, exact, inter, exact_d, inter_d);
   }
   for (size_t i = 0; i <= grid._p; i++) {
      const double x     = grid._tj[i];
      const double exact = testfunction(x);
      const double inter = grid.interpolate(x, v, 0, grid._p);

      const double exact_d = testfunction_d(x);
      const double inter_d = grid.interpolate_der(x, v, 0, grid._p);

      m   = std::max(m, std::abs(exact - inter));
      m_d = std::max(m_d, std::abs(exact_d - inter_d));
   }
   std::fclose(fptr);

   std::printf("Max difference in interpolation: %.6e\n", m);
   std::printf("Max difference in interpolation derivatives: %.6e\n", m_d);
   if (m > 1.0e-8) return errcode;
   if (m_d > 1.0e-6) return errcode;
   return 0;
}

int main()
{
   int i;
   i = test_exp_function();
   if (i != 0) return i;
   i = test_runge_function();
   if (i != 0) return i;
   return 0;
}
