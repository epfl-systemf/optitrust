typedef struct {
  int x;
  int y;
} vect;

typedef struct {
  int weight;
  vect pos;
  vect speed;
} obj;

vect f() { return {1, 1}; }

int main() {
  vect p = {0, 0};
  vect b = p;
  vect d;
  d.x = 1;
  d.y = 2;
  vect e;
  e = f();
  obj a = {0, {0, 0}, 0};
  a.pos = p;
}
