void f(int j) {
  int s = 0;
  s += (2 * j);
  s -= j;
}

int main() {
  int i = 1;
  /*@body*/ {
    int s = 0;
    s += (2 * i);
    s -= i;
  } /*body@*/
  return 0;
}
