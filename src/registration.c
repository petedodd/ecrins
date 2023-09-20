#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void tbmod0_initmod_desolve(void *);
extern void tbmod0_rhs_dde(void *);
extern void tbmod0_rhs_desolve(void *);

/* .Call calls */
extern SEXP tbmod0_contents(void *);
extern SEXP tbmod0_create(void *);
extern SEXP tbmod0_initial_conditions(void *, void *);
extern SEXP tbmod0_metadata(void *);
extern SEXP tbmod0_rhs_r(void *, void *, void *);
extern SEXP tbmod0_set_initial(void *, void *, void *, void *);
extern SEXP tbmod0_set_user(void *, void *);

static const R_CMethodDef CEntries[] = {
    {"tbmod0_initmod_desolve", (DL_FUNC) &tbmod0_initmod_desolve, 1},
    {"tbmod0_rhs_dde",         (DL_FUNC) &tbmod0_rhs_dde,         1},
    {"tbmod0_rhs_desolve",     (DL_FUNC) &tbmod0_rhs_desolve,     1},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"tbmod0_contents",           (DL_FUNC) &tbmod0_contents,           1},
    {"tbmod0_create",             (DL_FUNC) &tbmod0_create,             1},
    {"tbmod0_initial_conditions", (DL_FUNC) &tbmod0_initial_conditions, 2},
    {"tbmod0_metadata",           (DL_FUNC) &tbmod0_metadata,           1},
    {"tbmod0_rhs_r",              (DL_FUNC) &tbmod0_rhs_r,              3},
    {"tbmod0_set_initial",        (DL_FUNC) &tbmod0_set_initial,        4},
    {"tbmod0_set_user",           (DL_FUNC) &tbmod0_set_user,           2},
    {NULL, NULL, 0}
};

void R_init_ecrins(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
