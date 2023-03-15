#include "../../include/optitrust.h"

int f(int x) {
  return x;
}

int g(int x) {
  return x+1;
}

void simpl_in_depth () {
  int x = g(1 + 2 + f(2 + 4)) + g(5) + g(5);
  int y = g(1 + 2 + f(2 + 4)) + g(5) + g(5);
}

int main()
{
  double a; double  b; double  c; double d; double  e;
  double r; double  s; double t; double  u; double  v; double w; double  x; double  y; double z;

  // test parsing and normalizing
  x = 2 * b;
  x = a + b * 2 - (3 * c) - d;
  x = a + ((a + 2*b) + 1);
  x = 3 + 5.0 + a + (a + b*2 - 3*c);
  x = 3 * 5.0 * a * b * b;
  x = a / (b * b* c) / d;

  // test gather sum
  y = 1 + 2; // = 3
  y = 1 + 3.2 + 4 + 5.2; // = 5 + 3.2 + 5.2
  y = a + a + a - a; // = 2*a
  y = 2*a - 3*a + b - 5*b; // = - a -4*b
  y = 3*a + 2*b - b + c - 4*a ; // = -a + b + c

  // test gather prod
  z = a * a * a / a; // = a * a
  z = a / b / a * b / b; // = 1 / b

  // test gather nested
  t = (a + (3 + 2)*a); // = 6 * a
  t = a + b * a / b; // = 2 * a
  t = (a + (3 + 2)*a) / a * b; // = 6 * b
  t = a / (b * a * b) / (1 / b); // = 1 / b

  // test expand
  u = a * (b + c); // = a * b + a * c
  u = (a + b) * (c + d); // = a * c + a * d + b * c + b * d
  u = 5*a * (2*b + 3*c + 4*d);

  // test expand from pic demo
  v = (a / b * c / (a * d * d / b * e)); // = c / (d * d) / e
  v = (a * b * c / a); // = b * c
  v = ((a / (b / c)) + d * c / b) * (b / c); // = a + d
  v = ((a * ((b / c) * (d / ((((b * a) * a) / c) / e)))) * (a / e)); // = d
  v = (((a * d) / (a / e)) * (1. / e)); // = d
  v = (((a / (1. / e)) + (b * (c / (b / e )))) * (1. / e)); // = a + c


  // test for recursion in atoms (LATER)
  w = a + a * 3 + f(b * f(c + c) / b); // = 4*a + f(f(2*c))

  // test for loop shift
  int ls;
  ls = 0 + x - (-2); // = x + 2
  ls = 0 + x; // = x
  for (int ls2 = 0 + 2; ls2 < 10 + 2; ls2++) {
    ls = 0 + 12 + (-2); // = 10
  }

  // test exact division
  int ed1, ed2;
  exact_div(ed1, ed2);
  exact_div(5 + 5, 2); // = 5

  return 0;
}

