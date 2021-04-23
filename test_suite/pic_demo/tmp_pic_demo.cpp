const int nbSteps = 100;

const double step_duration = 0.2;

const int gridSize = 64;

const int nbCells = ((gridSize * gridSize) * gridSize);

const int bagCapacity = 100;

typedef struct {
  double x;
  double y;
  double z;
} vect;

vect vect_add(vect v1, vect v2) {
  return {(v1.x + v2.x), (v1.y + v2.y), (v1.z + v2.z)};
}

vect vect_mul(double d, vect v) { return {(d * v.x), (d * v.y), (d * v.z)}; }

typedef struct {
  vect pos;
  vect speed;
} particle;

typedef struct {
  int nb;
  particle items[bagCapacity];
} bag;

void bag_push(bag *b, particle p) {
  (b->items)[(b->nb)] = p;
  (b->nb)++;
}

void bag_clear(bag *b) { (b->nb) = 0; }

void bag_transfer(bag *b1, bag *b2) {
  for (int i = 0; (i < (b2->nb)); i++) {
    bag_push(b1, (b2->items)[i]);
  }
  bag_clear(b2);
}

bag bagsCur[nbCells];

bag bagsNext[nbCells];

vect fields[nbCells];

double nextCharge[nbCells];

void updateFieldsUsingNextCharge();

int idCellOfPos(vect pos);

int main() {
  for (int step = 0; (step < nbSteps); step++) {
    for (int idCell = 0; (idCell < nbCells); idCell++) {
      vect field = fields[idCell];
      bag *b = (&bagsCur[idCell]);
      int nb = (b->nb);
      for (int idParticle = 0; (idParticle < nb); idParticle++) {
        particle p = (b->items)[idParticle];
        vect speed2;
        speed2 = vect_add(p.speed, field);
        vect pos2 = vect_add(p.pos, vect_mul(step_duration, speed2));
        int idCell2 = idCellOfPos(pos2);
        nextCharge[idCell2] += 1.;
        particle p2 = {speed2, pos2};
        bag *b2 = (&bagsNext[idCell2]);
        bag_push(b2, p2);
      }
      bag_clear((&bagsCur[idCell]));
    }
    updateFieldsUsingNextCharge();
    for (int idCell = 0; (idCell < nbCells); idCell++) {
      bag_transfer((&bagsCur[idCell]), (&bagsNext[idCell]));
    }
  }
}
