#include <omp.h>
#include <stdio.h>

int main() {
#pragma omp parallel sections firstprivate(section_count)
  int section_count = 0;
  omp_set_dynamic(0);
  omp_set_num_threads(4);
  {
#pragma omp section
    {
      section_count++;
      printf("section_count %d\n", section_count);
    }
#pragma omp section
    {
      section_count++;
      printf("section_count %d\n", section_count);
    }
  }
  return 0;
}
