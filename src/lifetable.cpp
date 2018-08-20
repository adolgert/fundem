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
        double* mx, double* ax, double* nx, double* s, int N)
{
    FirstMomentSurvival(mx, ax, nx, s, N);
}


#ifdef __cplusplus
}
#endif
