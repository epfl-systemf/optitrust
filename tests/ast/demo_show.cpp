#include "../../include/optitrust.h"

int main() {
  int a, b;
  int x = 3;
  x--;
  for (int i = 0; i < 3; i++) {
    __strict();
    __smodifies("x ~> Cell");
    x++;
  }
}
