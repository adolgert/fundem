#include <cmath>
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
    auto answer = FirstMomentSurvival(mx, ax, nx);
    EXPECT_EQ(answer.size(), mx.size());

    for (auto check_idx=0; check_idx<23; check_idx++) {
        EXPECT_LT(std::abs(answer[check_idx] - expected[check_idx]), epsilon);
    }
}
