#include <stdio.h>
#include <omp.h>

void bar(float *a, int i, int j, int k);
int kl, ku, ks, jl, ju, js, il, iu,is;
void sub1(float *a)
{
  int i, j, k;
  for (k=kl; k<=ku; k+=ks)
    for (j=jl; j<=ju; j+=js)
      for (i=il; i<=iu; i+=is)
        bar(a,i,j,k);
}

void test()
{
  int j, k, jlast, klast;
  {
    for (k=1; k<=2; k++)
    for (j=1; j<=3; j++)
  {
    jlast=j;
    klast=k;
  }
    printf("%d %d\n", klast, jlast);
  }
}

void work(int a, int j, int k);
void sub2()
{
   int j, k, a;
   {
      for (k=1; k<=3; k++)
         for (j=1; j<=2; j++)
         {
            printf("%d %d %d\n", omp_get_thread_num(), k, j);
            /* end ordered */
            work(a,j,k);
         }
   }
}

