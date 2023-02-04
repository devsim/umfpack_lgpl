
#include "blasfunctions.h"
#include <stdio.h>

void testload(const char *dll_path)
{

  printf("Loading %s\n", dll_path);

  const char *msg = NULL;
  void *handle = blasw_load_dll(dll_path, &msg);
  if (msg)
  {
    printf("ERROR %s\n", msg);
  }
  printf("DLL ADDRESS %p\n\n", handle);

  if (handle)
  {
    int i = blasw_load_functions(handle);
    printf("MISSING %d functions\n", i);
  }
}

int main()
{
  blasw_clear_table();
  testload("libfoo.so");
  testload("libopenblas.so");

}

