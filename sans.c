#include <stdio.h>
#include <stdlib.h>
#include <x86intrin.h>

int main () {
  float V[1000];
  float C[1000];
  srand(_rdtsc());
  for(register unsigned int i = 0; i < 1000; i++) {
    V[i] = rand();
  }

  unsigned long long start = _rdtsc();

  for(register unsigned int i = 0; i < 1000; i++) {
    C[i] = V[i];
  }

  printf("%lld\n", _rdtsc() - start);
}