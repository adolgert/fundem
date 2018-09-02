#ifndef FUNDEM_LIFETABLE_HPP
#define FUNDEM_LIFETABLE_HPP

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>


namespace fundem {

template<typename REAL>
void FirstMomentSurvival(
        const REAL *const mx, const REAL *const ax, const REAL *const nx,
        REAL *const survival, int age_cnt, size_t N)
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
        REAL *const lx, REAL* dx, int age_cnt, size_t N)
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
            *dxi = l * (1 - px);
            l *= px;
        }
    }
}


template<typename REAL>
void FirstMomentPeriodLifeExpectancy(
        const REAL *const mx, const REAL *const ax, const REAL *const nx,
        REAL *const le, int age_cnt, size_t N)
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
        int age_cnt, size_t N)
{
    const REAL taylor_a = 1e-4;
    const REAL taylor_b = 1e-3;

    for (size_t pop_idx=0; pop_idx < N; pop_idx++)
    {
        for (int age_idx = 0; age_idx < age_cnt; age_idx++)
        {
            auto idx = pop_idx * age_cnt + age_idx;
            auto mx = mxi[idx];
            auto nx = nxi[age_idx];

            if (mx <= taylor_a) {
                ax[idx] = nx * (0.5 + nx * (mx/12 + std::pow(mx, 3) * std::pow(nx, 2) / 720));
            } else if (mx >= taylor_b) {
                auto expx = std::exp(-mx * nx);
                ax[idx] = 1 / mx - nx * expx / (1 - expx);
            } else {
                auto expx = std::exp(-mx * nx);
                ax[idx] = (nx * (0.5 + nx * (mx / 12 + std::pow(mx, 3) * std::pow(nx, 2) / 720))) *
                        (mx - taylor_a) / (taylor_b - taylor_a);
                ax[idx] += ((1 / mx) - (nx * expx) / (1 - expx)) *
                        (taylor_b - mx) / (taylor_b - taylor_a);
            }
        }
    }
}


template<typename REAL>
void GraduationMethod(
        const REAL *const mxi, const REAL *const nx, REAL *const axi,
        int age_cnt, size_t N)
{
    const REAL max_difference = 1e-5;
    const int look_back = 6;
    const int max_iterations = 20;

    for (int check_n_idx=0; check_n_idx < age_cnt; check_n_idx++) {
        if (nx[check_n_idx] != nx[0]) {
            throw std::runtime_error(
                    "All age intervals must match for graduation.");
        }
    }
    const REAL n = nx[0];

    // The constant mortality result is the starting approximation and
    // the fallback answer if graduation fails. This initial condition
    // is well-behaved in that the ax is >0 and <= n/2.
    ConstantMortalityMeanAge(mxi, nx, axi, age_cnt, N);

    auto working = std::vector<REAL>(2 * age_cnt);
    auto differences = std::vector<REAL>(max_iterations + look_back);
    auto dx = std::vector<REAL>(age_cnt);

    for (size_t pop_idx = 0; pop_idx < N; pop_idx++) {
        std::copy(axi + pop_idx * age_cnt, axi + (pop_idx + 1) * age_cnt,
                working.begin());
        std::copy(axi + pop_idx * age_cnt, axi + (pop_idx + 1) * age_cnt,
                working.begin() + age_cnt);
        auto mx = mxi + pop_idx * age_cnt;
        for (int look_init_idx = 0; look_init_idx < look_back; look_init_idx++) {
            differences[look_init_idx] = n;
        }

        REAL* answer = nullptr;
        for (int it_idx=0; it_idx < max_iterations; it_idx++) {
            REAL * ax = &working[(it_idx % 2) * age_cnt];
            const REAL * last_ax = &working[((it_idx + 1) % 2) * age_cnt];
            // Compute dx, but look for dx<0 b/c it indicates not converging.
            REAL l = 1;
            for (int dx_idx = 0; dx_idx < age_cnt; dx_idx++) {
                auto px = (1.0 - mx[dx_idx] * last_ax[dx_idx]) /
                        (1.0 + mx[dx_idx] * (n - last_ax[dx_idx]));
                dx[dx_idx] = l * (1 - px);
                l *= px;
            }

            // Compute the new ax.
            for (int shift_idx=1; shift_idx < age_cnt - 1; shift_idx++) {
                REAL axx = -n;
                if (dx[shift_idx] > 1e-8) {
                    // dx drops exponentially, so this equation gives wild
                    // values for large ages.
                    axx = n * (0.5 +
                            (dx[shift_idx + 1] - dx[shift_idx - 1])
                            / (24 * dx[shift_idx]));
                }
                if (axx < 0 || axx > n) { // Fall back to constant mortality.
                    axx = axi[pop_idx * age_cnt + shift_idx];
                }
                ax[shift_idx] = axx;
            }

            // Compute difference between this and last iteration.
            REAL iter_difference = 0;
            for (int diff_idx=0; diff_idx < age_cnt; diff_idx++) {
                iter_difference = std::max(iter_difference,
                        std::abs(ax[diff_idx] - last_ax[diff_idx]));
            }
            differences[it_idx + look_back] = iter_difference;
            if (iter_difference > differences[it_idx]) {
                break;
            } else if (iter_difference < max_difference) {
                answer = ax;
                break;
            } // else keep going.
        }
        if (nullptr != answer) {
            std::copy(answer, answer + age_cnt, axi + pop_idx * age_cnt);
        } // else it's already initialized with the constant-mortality answer.
    }
}

}

#endif //FUNDEM_LIFETABLE_HPP
