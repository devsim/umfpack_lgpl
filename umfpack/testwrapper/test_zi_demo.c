#include <stdio.h>
#include "blasfunctions.h"
extern void umfpack_zi_demo();
int main (int argc, char **argv)
{
  const char *msg = NULL;
  void *dll = blasw_load_dll("libopenblas.so", &msg);
  if (!dll)
  {
    printf("%s\n", msg);
    return -1;
  }

  int i = blasw_load_functions(dll);
  if (i != 0)
  {

    printf("Missing %i symbols\n", i);
    return -1;
  }

  umfpack_zi_demo();
  return 0;
}

