#include <omp.h>
int main() {
  int x = 3;
  for (int i = 0; i < 3; i++) {
    if (true) {
      x++;
    } else {
      i++;
    }
  }

  x++;
  return 0;
}
