#ifndef BLASFUNCTIONS_H
#define BLASFUNCTIONS_H
#ifdef __cplusplus
extern "C" {
#endif
#if defined _WIN32
#define BLASWRAPPER_PUBLIC __declspec(dllexport) // Note: actually gcc seems to also supports this syntax.
#define BLASWRAPPER_LOCAL __declspec(dllimport)
#else
#define BLASWRAPPER_PUBLIC __attribute__ ((visibility("default")))
#define BLASWRAPPER_LOCAL  __attribute__ ((visibility("hidden")))
#endif

BLASWRAPPER_PUBLIC
void blasw_clear_table();

BLASWRAPPER_PUBLIC
int blasw_load_functions(void *dll_p);

BLASWRAPPER_PUBLIC
void *blasw_load_dll(const char *dllname, const char **error);

#ifdef __cplusplus
}
#endif
#endif
