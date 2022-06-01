
int main() {}

class test_static_class {

private:

  // "static" is a trm_annot
  static int foo(int x) {
     return x;
  }

public:

  static int bar(int x) {
     return x;
  }

};


class test_class {

private:

  int x;

public:

  // for each function, mark it as "public" or "private" as a trm_annot
  void move(int d) {
     x += d;
  }

  /* encoded as:
  void move(test_class* this, int d) {
     this->x += d;
  }
  */

  bool test_this() {
    return this.x == x;
  }

};


// we will consider templated functions and templated classes,
// only with simple typename arguments

// the function/class binds a list of type variables
template<typename T>
bool test_poly(T* x, T* y) { // occurence of T is a Typ_var
  return x == y;
}


#include <vector>
#include <algorithm>

// foo::bar is a "qualified variable",
// same as modules in ocaml M.N.x
// type qvar =  { qvar_var = "x"; qvar_path = ["M";"N"]; qvar_str = "M.N.x" in OCaml or "M::N::x" in C }

void test_vector() {
  std::vector<int> v;   // vector<int> is a typ_constr (typ_constrid, typ_int)
  v.push_back(3);  // encoded as push_back(v,3)
  int a = v[0];
}

int test_iterator(std::vector<int> v) {
  int r = 0;
  for (auto it = std::begin(v); it != std::end(v); it++) { // auto type needs to be supported
      r += *it;
  }
  return r;
}

int test_lambda(std::vector<int> v) {
  int r = 0;
  std::for_each(std::begin(v), std::end(v), [&](int const& x) { // we only support arguments by references [&]
     r += x; }); // trm_fun
  return r;
}
