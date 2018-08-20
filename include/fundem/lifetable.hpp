#ifndef FUNDEM_LIFETABLE_HPP
#define FUNDEM_LIFETABLE_HPP

#include <vector>


namespace fundem {

template<typename REAL>
void FirstMomentSurvival(
        const std::vector<REAL>& mx, const std::vector<REAL>& ax,
        const std::vector<REAL>& nx, std::vector<REAL>& survival)
{
    auto survival_iter = survival.begin();

    for (auto mxi = mx.begin(), axi = ax.begin(), nxi = nx.begin();
            mxi != mx.end();
            mxi++, axi++, nxi++, survival_iter++
            )
    {
        *survival_iter = (1.0 - *mxi * *axi) / (1.0 + *mxi * (*nxi - *axi));
    }
}

}

#endif //FUNDEM_LIFETABLE_HPP
