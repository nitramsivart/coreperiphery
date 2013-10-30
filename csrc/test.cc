#include <time.h>
#include <stdio.h>
int add(int *i) {
  return *i+1;
}

int main(int argc, char* argv[]) {
  int runs = 10000000000;
  clock_t t;
  t = clock();
  for(int i = 0; i < runs; i++)
    i += 1;
  t = clock() - t;
  printf("Clicks: %d\n", t);
  t = clock();
  for(int i = 0; i < runs; i++)
    add(&i);
  t = clock() - t;
  printf("Clicks: %d\n", t);
  return 0;
}
