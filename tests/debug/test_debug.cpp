int main(){

  int p = 10;
  int q = p;
  q = p;

  int &x = p;
  x = x + 1;

  int *y = &p;
  *y = *y + 1;

  return 0;
}