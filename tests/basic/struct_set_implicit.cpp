
typedef struct {
    int x;
    int y; }
  vect;

typedef struct {
  int weight;
  vect pos;
  vect speed;
} obj;

vect f() { 
  return {1,1};
}
int main() {
  vect p = {0,0};
  vect b;
  {
    b.x = p.x;
    b.y = p.y;
  }

  vect e;
  e = f();
  obj a = {0,{0,0},0};
  a.pos = p;
}
