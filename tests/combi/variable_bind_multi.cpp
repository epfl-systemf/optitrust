
int main() {

  int s = 0;
  for (int i = 0; i < 3; i++) {
    s = 2 + 3;
  }

  int t = 0;
  for (int j = 0; j < 3; j++) {
    t += 2 + s;
    t += 2 + s;
  }

  return 0;
}

