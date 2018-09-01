#ifndef FUNDEM_LIFETABLE_HPP
#define FUNDEM_LIFETABLE_HPP

#include <vector>


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

}

#endif //FUNDEM_LIFETABLE_HPP
