typedef struct {
  int x;
  int y;
} vect;

typedef struct {
  vect pos;
  vect speed;
} particle;

int main() {
  
  vect v = {0,0};
  particle p = {{0,0},{0,0}};

  vect u = p.pos;
  int x = p.pos.x;

  particle p2 = p;

  return 0;
}

