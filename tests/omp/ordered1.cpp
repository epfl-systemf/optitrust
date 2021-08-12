#include <stdio.h>

void work(int k)
{
  printf(" %d\n", k);
}

void ordered_example(int lb, int ub, int stride)
{
  for (int i=lb; i<ub; i+=stride)
    work(i);
}

int main()
{
  ordered_example(0, 100, 5);
  return 0;
}
