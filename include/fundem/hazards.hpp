//
// Defines various hazard rates for total mortality.
//

#ifndef FUNDEM_HAZARDS_HPP
#define FUNDEM_HAZARDS_HPP

#include <cmath>


/*! Siler hazard rate with default values
 *
 * This Siler distribution is a good approximation to what a real total
 *  mortality rate looks like. Both the equations and the parameters come
 *  from a paper [1] where they were fit to a Scandinavian country.
 *  We will use this as the one true mortality rate for this session.
 *  [1] V. Canudas-Romo and R. Schoen, “Age-specific contributions to
 *  changes in the period and cohort life expectancy,” Demogr. Res.,
 *  vol. 13, pp. 63–82, 2005.
 *
 * @tparam REAL The type of floating point.
 * @param x Age group.
 * @param t Year from the starting year.
 * @return
 */
template<typename REAL>
REAL siler_default(REAL x, REAL t) {
    const REAL a1{0.2};
    const REAL a2{0.0002};
    const REAL a3{0.003};
    const REAL b1{1};
    const REAL b2{0.1};
    const REAL c1{0.015};
    const REAL c2{0.01};

    return a1 * std::exp(-b1 * x - c1 * t) + a2 * std::exp(b2 * x - c2 * t)
        + a3 * std::exp(-c2 * t);
}


#endif //FUNDEM_HAZARDS_HPP
