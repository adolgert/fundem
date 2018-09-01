#include <exception>
#include "Python.h"
#include "fundem/lifetable.hpp"

using namespace fundem;

#ifdef FUNDEM_DLL
#ifdef FUNDEM_EXPORTS
#define FUNDEM_API __declspec(dllexport)
#else
#define FUNDEM_API __declspec(dllimport)
#endif /* FUNDEM_EXPORTS */
#else
#define FUNDEM_API
#endif /* FUNDEM_DLL */

#ifdef __cplusplus
extern "C" {
#endif

FUNDEM_API void first_moment_survival(
        double* mx, double* ax, double* nx, double* s, size_t age_cnt, size_t N)
{
    try {
        FirstMomentSurvival(mx, ax, nx, s, age_cnt, N);
    } catch (std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
    }
}


#ifdef __cplusplus
}
#endif
