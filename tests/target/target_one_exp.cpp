/*@2_0, 1_0*/ /*@5_0*/ typedef struct {
  int x;
  int y;
} vect; /*5_0@*/

typedef struct {
  vect pos;
  vect speed;
} particle;

/*@6_0*/ typedef int* intstar; /*6_0@*/

/*@24_0*/ int f(int n) { /*@20_0*/
  return n;              /*20_0@*/
} /*24_0@*/

/*@27_0, 26_0, 25_0*/ void g(int t[2], vect* varg) {
  /*@33_0*/ int b = t[0];    /*33_0@*/
  /*@33_1*/ int a = varg->x; /*33_1@*/
} /*27_0, 26_0, 25_0@*/

/*@23_0*/ int main() {
  /*@15_0*/ for (int i = 0; i < 10; i++) {
    /*@17_0, 16_0*/ for (int j = 0; j < 5; j += 1) {
      /*@36_0, 32_0*/ i++; /*36_0, 32_0@*/
      /*@18_0*/ break;     /*18_0@*/
    }                      /*17_0, 16_0@*/
    /*@19_0*/ continue;    /*19_0@*/
  }                        /*15_0@*/
  int val = /*@4_0*/ 13 /*4_0@*/;
  /*@41_0, 40_0, 38_0, 37_0, 33_2, 10_0*/ int r =
      /*@14_0, 13_0, 12_0*/ 3 /*14_0, 13_0, 12_0@*/; /*41_0, 40_0, 38_0, 37_0,
                                                        33_2, 10_0@*/
  r = r + 1;
  /*@31_0, 28_0*/ r += 2; /*31_0, 28_0@*/
  r++;
  /*@42_0, 33_3*/ int s = f(r); /*42_0, 33_3@*/
  /*@41_1, 33_4*/ int n = 3;    /*41_1, 33_4@*/
  int* m = &n;
  *m = *m + 4;
  int t[2] = {5, 6};
  /*@33_5*/ int u = t[0]; /*33_5@*/
  t[1] = /*@7_0*/ u /*7_0@*/ + 2;
  vect v = {5, 6};
  /*@33_6*/ int a = v.x; /*33_6@*/
  v.y = a + 2;
  /*@34_0*/ vect v2 = v; /*34_0@*/
  v2 = v;
  /*@39_0, 34_1*/ particle p1 = {v, v}; /*39_0, 34_1@*/
  /*@43_0, 39_1, 34_2*/ particle p2 = {
      v, {7, /*@3_0*/ 8 /*3_0@*/}}; /*43_0, 39_1, 34_2@*/
  /*@44_0*/ p1.pos.x /*44_0@*/ = p2.pos.y;
  {
    /*@41_2, 38_1, 34_3*/ int r1 = 1;                    /*41_2, 38_1, 34_3@*/
    /*@41_3, 38_2, 34_4*/ int r2 = 2;                    /*41_3, 38_2, 34_4@*/
    /*@41_4, 38_3, 34_5*/ int r3 = /*@8_0*/ r2 /*8_0@*/; /*41_4, 38_3, 34_5@*/
  }
  /*@33_7*/ int y = /*@22_0*/ f(2) /*22_0@*/; /*33_7@*/
/*@21_0*/ lbl2:
  g(t, &v); /*21_0@*/
  int k;
  for (k = 0; k < 10; k++) {
    y = k;
  }
  /*@20_1*/ return 0; /*20_1@*/
} /*23_0@*/           /*2_0, 1_0@*/
