int *t;

int *u;

int n;

int main() {
  for (int i = 1; i < n; i++) {
    int a = i;
    t[i] += a;
  }
  for (int i = 1; i < n; i++) {
    int b = i;
    u[i] += b;
  }
  return 0;
}
