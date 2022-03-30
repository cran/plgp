#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void calc_alcs_R(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void calc_ecis_R(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void calc_eis_R(void *, void *, void *, void *, void *, void *);
extern void calc_ents_R(void *, void *, void *, void *);
extern void calc_ieci_R(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void calc_iecis_R(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void calc_ktKikx_R(void *, void *, void *, void *, void *, void *, void *);
extern void calc2_ktKikx_R(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void covar_sep_R(void *, void *, void *, void *, void *, void *, void *, void *);
extern void covar_sep_symm_R(void *, void *, void *, void *, void *, void *);
extern void covar_sim_R(void *, void *, void *, void *, void *, void *, void *, void *);
extern void covar_sim_symm_R(void *, void *, void *, void *, void *, void *);
extern void dist2covar_R(void *, void *, void *, void *, void *, void *);
extern void dist2covar_symm_R(void *, void *, void *, void *, void *);
extern void distance_R(void *, void *, void *, void *, void *, void *);
extern void distance_symm_R(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"calc_alcs_R",       (DL_FUNC) &calc_alcs_R,       20},
    {"calc_ecis_R",       (DL_FUNC) &calc_ecis_R,       18},
    {"calc_eis_R",        (DL_FUNC) &calc_eis_R,         6},
    {"calc_ents_R",       (DL_FUNC) &calc_ents_R,        4},
    {"calc_ieci_R",       (DL_FUNC) &calc_ieci_R,       20},
    {"calc_iecis_R",      (DL_FUNC) &calc_iecis_R,      23},
    {"calc_ktKikx_R",     (DL_FUNC) &calc_ktKikx_R,      7},
    {"calc2_ktKikx_R",    (DL_FUNC) &calc2_ktKikx_R,    12},
    {"covar_sep_R",       (DL_FUNC) &covar_sep_R,        8},
    {"covar_sep_symm_R",  (DL_FUNC) &covar_sep_symm_R,   6},
    {"covar_sim_R",       (DL_FUNC) &covar_sim_R,        8},
    {"covar_sim_symm_R",  (DL_FUNC) &covar_sim_symm_R,   6},
    {"dist2covar_R",      (DL_FUNC) &dist2covar_R,       6},
    {"dist2covar_symm_R", (DL_FUNC) &dist2covar_symm_R,  5},
    {"distance_R",        (DL_FUNC) &distance_R,         6},
    {"distance_symm_R",   (DL_FUNC) &distance_symm_R,    4},
    {NULL, NULL, 0}
};

void R_init_plgp(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
