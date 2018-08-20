#include <cmath>
#include <ostream>
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "fundem/lifetable.hpp"


using namespace fundem;

double epsilon = 1e-9;

// Simple test, does not use gmock
TEST(Dummy, foobar)
{
EXPECT_EQ(1, 1);
}


/*! Tests survival calculation for constant mortality rate.
 *
 */
TEST(FM_SURVIVAL_CONSTANT_MORTALITY, blah)
{
    std::vector<double> mx(23);
    std::vector<double> ax(23);
    std::vector<double> nx(23);
    std::vector<double> expected(23);

    double mortality_rate = 0.1;
    double base_interval = 5;

    for (auto i=0; i<23; i++) {
        nx[i] = base_interval * (i+1)/23;
        auto mean_age = (1 / mortality_rate) - nx[i] / (std::exp(mortality_rate * nx[i]) - 1);
        mx[i] = mortality_rate;
        ax[i] = mean_age;
        expected[i] = std::exp(-nx[i] * mortality_rate);
    }
    std::vector<double> answer(mx.size());
    FirstMomentSurvival(&mx[0], &ax[0], &nx[0], &answer[0], mx.size());
    EXPECT_EQ(answer.size(), mx.size());

    for (auto check_idx=0; check_idx<23; check_idx++) {
        EXPECT_LT(std::abs(answer[check_idx] - expected[check_idx]), epsilon);
    }
}


/*! Tests survival calculation for constant mortality rate.
 *  Given mu(t)=\mu_c t so that \int \mu_c t = \mu_c\left(\frac{t^2}{2}-\frac{x^2}{2}\right).
 *  S(t)=\exp(-\frac{1}{2} \mu_c n(n+2 t)) which leads to
 *  ax(t)= lots of terms.
 */
TEST(FM_SURVIVAL_RISING_MORTALITY, blah)
{
    std::vector<double> mx(23);
    std::vector<double> ax(23);
    std::vector<double> nx(23);
    std::vector<double> expected(23);

    const double m{0.01};
    const double base_interval{5};

    auto Pi = 4 * std::atan(1);
    double x{0};
    for (auto i=0; i<23; i++) {
        // Pick any interval
        auto n = base_interval * (i+1)/23;
        auto mortality_rate = ((-1 + std::exp((m*n*(n + 2*x))/2.))*std::sqrt(m)*std::sqrt(2/Pi))/
                 (std::exp((m*std::pow(n + x,2))/2.)*(-std::erf((std::sqrt(m)*x)/std::sqrt(2)) +
                                                  std::erf((std::sqrt(m)*(n + x))/std::sqrt(2))));
        auto mean_age = (-2*std::sqrt(m)*n + std::exp((m*std::pow(n + x,2))/2.)*std::sqrt(2*Pi)*
                                        (-std::erf((std::sqrt(m)*x)/std::sqrt(2)) + std::erf((std::sqrt(m)*(n + x))/std::sqrt(2))))/
                        (2.*(-1 + std::exp((m*n*(n + 2*x))/2.))*std::sqrt(m));
        auto S = std::exp(-0.5 * m * n * (n + 2 * x));

        // The mortality rate should be close to m * t at the half-interval.
        auto mt = mortality_rate / (m * (x+0.5*n));
        // The mean age should be near half of the interval size.

        mx[i] = mortality_rate;
        ax[i] = mean_age;
        nx[i] = n;
        expected[i] = S;
        x += nx[i];
    }
    std::vector<double> answer(mx.size());
    FirstMomentSurvival(&mx[0], &ax[0], &nx[0], &answer[0], mx.size());
    EXPECT_EQ(answer.size(), mx.size());

    for (auto check_idx=0; check_idx<23; check_idx++) {
        EXPECT_LT(std::abs(answer[check_idx] - expected[check_idx]), epsilon);
    }
}

