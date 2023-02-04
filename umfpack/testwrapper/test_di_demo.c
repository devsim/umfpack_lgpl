#include "blasfunctions.h"
#include <stdio.h>

static int pprefix = 1;
static void myrp(const char *msg)
{
  if (pprefix)
  {
    printf("UMF: ");
  }
  int pc = printf("UMF: %s", msg);
  pprefix = (pc > 0) && (msg[pc] == '\n');
}

extern void umfpack_di_demo();
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

  blasw_set_printer_callback(myrp);
  umfpack_di_demo();
  return 0;
}

