#ifndef FUNDEM_LIFETABLE_HPP
#define FUNDEM_LIFETABLE_HPP

namespace fundem {

template<REAL> std::vector<REAL> FirstMomentSurvival(
        const std::vector<REAL>& mx, const std::vector<REAL>& ax, const std::vector<REAL>& nx)
{
    std::vector<REAL> survival(mx.size());

    for (
            const auto mxi = mx.begin(),
            const auto axi = ax.begin(),
            const auto nxi = nx.begin(),
            auto survival_iter = survival.begin();
            mxi != mx.end();
            mxi++, axi++, nxi++, survival_iter++
            ) {
        *survival_iter = (1.0 - *mxi * *axi) / (1.0 + *mxi * (*nxi - *axi))
    }
    return survival;
}

}

#endif //FUNDEM_LIFETABLE_HPP
