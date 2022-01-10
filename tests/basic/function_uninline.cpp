

void g(int x, int y) {}

void gtwice(int x) {
  g(x, x);
}

void test_trivial() {
  gtwice_body: { g(3, 3); }
}

void f(int x) {
  int a = x+1;
  g(a, x);
}

void test_basic() {
  int r = 5;
  fbody:{
    int b = (r+2)+1;
    g(b, r+2);
  }
}

void iter_nat_for(int n, void body(int)) {
  for (int i = 0; i < n; i++) {
    body(i);
  }
}

void test_ho() {
  int s = 0;
  int m = 3;
  hobody: {
    for (int j = 0; j < m; j++) {
      {
        s += 2*j;
        s -= j;
      }
    }
  }
}

#include <stdio.h>
#include <stdlib.h>

typedef struct { } particle;
typedef struct { } bag;
typedef struct { } bag_iter;
bag_iter* bag_iter_begin(bag* b);
particle* bag_iter_get(bag_iter* it);
particle* bag_iter_next(bag_iter* it, bool destructive);

void iter_bag(bag* b, void body(particle*)) {
  bag_iter* const iter = bag_iter_begin(b);
  for (particle* p = bag_iter_get(iter); p != NULL; p = bag_iter_next(iter, true)) {
    body(p);
  }
  free(iter);
}

void test_bag() {
  int x = 0;
  bag* mybag;
  bagbody: {
    bag_iter* const myit = bag_iter_begin(mybag);
    for (particle* p = bag_iter_get(myit); p != NULL; p = bag_iter_next(myit, true)) {
      {
         if (p = p) { x++; }
      }
    }
    free(myit);
  }
}

/* TODO add support for this nicer version

typedef struct { } particle;
typedef struct { } bag;
typedef struct { } bag_iter;
particle* bag_iter_begin(bag_iter* it, bag* b);
particle* bag_iter_next(bag_iter* it, bool destructive);

void iter_bag(bag* b, void body(particle*)) {
  bag_iter it;
  for (particle* p = bag_iter_begin(&it, b); p != NULL; p = bag_iter_next(&it, true)) {
    body(p);
  }
}

void test_bag() {
  bag* mybag;
  bagbody: {
    // This is the code pattern to use in pic_demo
    bag_iter myit;
    for (particle* p = bag_iter_begin(&myit, mybag); p != NULL; p = bag_iter_next(&myit, true)) {
      {
         if (p = p) { return; }
      }
    }
  }
}

*/