#include <cmath>
#include <ostream>
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "fundem/hazards.hpp"
#include "fundem/lifetable.hpp"


using namespace fundem;

double epsilon = 1e-9;




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
    FirstMomentSurvival(&mx[0], &ax[0], &nx[0], &answer[0], mx.size(), 1);
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
    FirstMomentSurvival(&mx[0], &ax[0], &nx[0], &answer[0], mx.size(), 1);
    EXPECT_EQ(answer.size(), mx.size());

    for (auto check_idx=0; check_idx<23; check_idx++) {
        EXPECT_LT(std::abs(answer[check_idx] - expected[check_idx]), epsilon);
    }
}




/*! See if the graduation method converges for a known distribution.
 *  Uses Siler distribution for the mx.
 */
 TEST(GRADUATION_SILER, graduation_siler)
{
     const int age_cnt = 20;
     std::vector<double> mx(age_cnt);
     std::vector<double> nx(age_cnt);

     for (int nx_init=0; nx_init < age_cnt; nx_init++) {
         nx[nx_init] = 5.0;
     }

     double x = 0;
     mx[0] = siler_default(0.5, 0.0);
     for (int mx_init=1; mx_init < age_cnt; mx_init++) {
         mx[mx_init] = siler_default(x + 0.5 * nx[mx_init], 0.0);
         x += nx[mx_init];
     }

     std::vector<double> ax(age_cnt);
     GraduationMethod(&mx[0], &nx[0], &ax[0], age_cnt, 1);

     for (int check_idx=0; check_idx < age_cnt; check_idx++) {
         EXPECT_LT(ax[check_idx], nx[check_idx]);
         EXPECT_LT(0, ax[check_idx]);
     }
}




/*! Graduation with this method should fail for different nx.
 */
TEST(GRADUATION_SILER_AGE_INCONSISTENT, graduation_siler_age_inconsistent)
{
    std::vector<double> mx(23);
    std::vector<double> nx(23);

    nx[0] = 7.0/365.0;
    nx[1] = 28.0/365.0;
    nx[2] = (365 - 7 - 28) / 365.0;
    nx[3] = 4;
    for (int nx_init=4; nx_init < 23; nx_init++) {
        nx[nx_init] = 5.0;
    }

    double x = 0;
    for (int mx_init=0; mx_init < 23; mx_init++) {
        mx[mx_init] = siler_default(x + 0.5 * nx[mx_init], 0.0);
        x += nx[mx_init];
    }

    std::vector<double> ax(23);
    ASSERT_THROW(GraduationMethod(&mx[0], &nx[0], &ax[0], 23, 1),
            std::runtime_error);
}
