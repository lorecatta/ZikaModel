#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void odin_model_determ_rhs_dde(void *);
extern void odin_model_stoch_rhs_dde(void *);

/* .Call calls */
extern SEXP odin_model_determ_contents(SEXP);
extern SEXP odin_model_determ_create(SEXP);
extern SEXP odin_model_determ_initial_conditions(SEXP, SEXP);
extern SEXP odin_model_determ_metadata(SEXP);
extern SEXP odin_model_determ_rhs_r(SEXP, SEXP, SEXP);
extern SEXP odin_model_determ_set_initial(SEXP, SEXP, SEXP, SEXP);
extern SEXP odin_model_determ_set_user(SEXP, SEXP);
extern SEXP odin_model_stoch_contents(SEXP);
extern SEXP odin_model_stoch_create(SEXP);
extern SEXP odin_model_stoch_initial_conditions(SEXP, SEXP);
extern SEXP odin_model_stoch_metadata(SEXP);
extern SEXP odin_model_stoch_rhs_r(SEXP, SEXP, SEXP);
extern SEXP odin_model_stoch_set_initial(SEXP, SEXP, SEXP, SEXP);
extern SEXP odin_model_stoch_set_user(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"odin_model_determ_rhs_dde", (DL_FUNC) &odin_model_determ_rhs_dde, 1},
    {"odin_model_stoch_rhs_dde",  (DL_FUNC) &odin_model_stoch_rhs_dde,  1},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"odin_model_determ_contents",           (DL_FUNC) &odin_model_determ_contents,           1},
    {"odin_model_determ_create",             (DL_FUNC) &odin_model_determ_create,             1},
    {"odin_model_determ_initial_conditions", (DL_FUNC) &odin_model_determ_initial_conditions, 2},
    {"odin_model_determ_metadata",           (DL_FUNC) &odin_model_determ_metadata,           1},
    {"odin_model_determ_rhs_r",              (DL_FUNC) &odin_model_determ_rhs_r,              3},
    {"odin_model_determ_set_initial",        (DL_FUNC) &odin_model_determ_set_initial,        4},
    {"odin_model_determ_set_user",           (DL_FUNC) &odin_model_determ_set_user,           2},
    {"odin_model_stoch_contents",            (DL_FUNC) &odin_model_stoch_contents,            1},
    {"odin_model_stoch_create",              (DL_FUNC) &odin_model_stoch_create,              1},
    {"odin_model_stoch_initial_conditions",  (DL_FUNC) &odin_model_stoch_initial_conditions,  2},
    {"odin_model_stoch_metadata",            (DL_FUNC) &odin_model_stoch_metadata,            1},
    {"odin_model_stoch_rhs_r",               (DL_FUNC) &odin_model_stoch_rhs_r,               3},
    {"odin_model_stoch_set_initial",         (DL_FUNC) &odin_model_stoch_set_initial,         4},
    {"odin_model_stoch_set_user",            (DL_FUNC) &odin_model_stoch_set_user,            2},
    {NULL, NULL, 0}
};

void R_init_ZikaModel(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
