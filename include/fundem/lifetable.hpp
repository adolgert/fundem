#ifndef FUNDEM_LIFETABLE_HPP
#define FUNDEM_LIFETABLE_HPP

#include <cmath>


namespace fundem {

template<typename REAL>
void FirstMomentSurvival(
        const REAL *const mx, const REAL *const ax, const REAL *const nx,
        REAL *const survival, size_t age_cnt, size_t N)
{
    auto mxi = mx;
    auto axi = ax;
    auto si = survival;
    for (size_t pop_idx=0; pop_idx < N; pop_idx++)
    {
        auto nxi = nx;
        auto mx_end = mxi + age_cnt;
        for (; mxi != mx_end; mxi++, axi++, nxi++, si++)
        {
            *si = (1.0 - *mxi * *axi) / (1.0 + *mxi * (*nxi - *axi));
        }
    }
}


template<typename REAL>
void FirstMomentPopulation(
        const REAL *const mx, const REAL *const ax, const REAL *const nx,
        REAL *const lx, REAL* dx, size_t age_cnt, size_t N)
{
    auto mxi = mx;
    auto axi = ax;
    auto lxi = lx;
    auto dxi = dx;
    for (size_t pop_idx=0; pop_idx < N; pop_idx++)
    {
        auto nxi = nx;
        auto mx_end = mxi + age_cnt;
        REAL l = 1;
        for (; mxi != mx_end; mxi++, axi++, nxi++, lxi++, dxi++)
        {
            *lxi = l;
            auto px = (1.0 - *mxi * *axi) / (1.0 + *mxi * (*nxi - *axi));
            *dxi = l * px;
            l *= px;
        }
    }
}


template<typename REAL>
void FirstMomentPeriodLifeExpectancy(
        const REAL *const mx, const REAL *const ax, const REAL *const nx,
        REAL *const le, size_t age_cnt, size_t N)
{
    for (size_t pop_idx=0; pop_idx < N; pop_idx++)
    {
        le[pop_idx * age_cnt + age_cnt - 1] = mx[pop_idx * age_cnt + age_cnt - 1];
        for (int age_idx = age_cnt - 2; age_idx >= 0; age_idx--)
        {
            size_t now_idx = pop_idx * age_cnt + age_idx;
            size_t prev_idx = pop_idx * age_cnt + age_idx + 1;
            le[now_idx] = (nx[age_idx] + (1 - mx[now_idx] * ax[now_idx]) * le[prev_idx]) /
                    (1 + mx[now_idx] * (nx[age_idx] - ax[now_idx]));
        }
    }
}


template<typename REAL>
void ConstantMortalityMeanAge(
        const REAL *const mxi,  const REAL *const nxi, REAL *const ax,
        size_t age_cnt, size_t N)
{
    const double taylor_a = 1e-4;
    const double taylor_b = 1e-3;

    for (size_t pop_idx=0; pop_idx < N; pop_idx++)
    {
        for (size_t age_idx = 0; age_idx < age_cnt; age_idx++)
        {
            auto mx = mxi[pop_idx * age_cnt + age_idx];
            auto nx = nxi[age_idx];

            if (mx <= taylor_a) {
                *ax = nx * (0.5 + nx * (mx/12 + std::pow(mx, 3) * std::pow(nx, 2) / 720));
            } else if (mx >= taylor_b) {
                auto expx = std::exp(-mx * nx);
                *ax = 1 / mx - nx * expx / (1 - expx);
            } else {
                auto expx = std::exp(-mx * nx);
                *ax = (nx * (0.5 + nx * (mx / 12 + std::pow(mx, 3) * std::pow(nx, 2) / 720))) *
                        (mx - taylor_a) / (taylor_b - taylor_a);
                *ax += ((1 / mx) - (nx * expx) / (1 - expx)) *
                        (taylor_b - mx) / (taylor_b - taylor_a);
            }
        }
    }
}


}

#endif //FUNDEM_LIFETABLE_HPP
