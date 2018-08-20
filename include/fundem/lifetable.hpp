#ifndef FUNDEM_LIFETABLE_HPP
#define FUNDEM_LIFETABLE_HPP

#include <vector>


namespace fundem {

template<typename REAL>
void FirstMomentSurvival(
        const REAL *const mx, const REAL *const ax, const REAL *const nx,
        REAL *const survival, size_t N)
{
    auto mxi = mx;
    auto axi = ax;
    auto nxi = nx;
    auto si = survival;
    auto end = mx + N;
    for (; mxi != end; mxi++, axi++, nxi++, si++)
    {
        *si = (1.0 - *mxi * *axi) / (1.0 + *mxi * (*nxi - *axi));
    }
}

}

#endif //FUNDEM_LIFETABLE_HPP
