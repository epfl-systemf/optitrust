#ifndef PIC_VERT_DIAGNOSTICS
#define PIC_VERT_DIAGNOSTICS

#include "meshes.h"   // types  cartesian_mesh_1d, cartesian_mesh_2d, cartesian_mesh_3d
#include "variadic.h" // macros VARIADIC, NUMARG2

/*****************************************************************************
 *                        Theoretical Damping Values                         *
 *****************************************************************************/

// Linear Landau Damping : theoretical electric field values.
// Ex(x, t) = 4 * alpha * er * exp(−omega_imag * t) * sin(kmode * x) * cos(omega_real * t − psi)
// Sonnendrücker, "Numerical Methods for the Vlasov-Maxwell equations", 2016
// 4.4.2. Landau damping, p. 54
typedef struct damping_values damping_values;
struct damping_values {
    double er;
    double psi;
    double omega_real; // plasma frequency
    double omega_imag; // damping coefficient
};

// Hashtable data structure, slightly modified.
// Kernighan & Ritchie, "The C programming language", 1988
// 6.6, p. 143
typedef struct hashmap hashmap;
struct hashmap {
    struct hashmap *next;
    double kmode;
    damping_values values;
};

#define MAP_SIZE 100

hashmap* hashtable[MAP_SIZE];

/*
 * Form hash value in [|0 ; MAP_SIZE|[ for a given kmode.
 *
 * @param  kmode
 * @return the hash value of kmode
 */
unsigned int hash(double kmode);

/*
 * Look for a given kmode in the hashmap
 *
 * @param  kmode
 * @return the damping values when found, NULL if not found.
 */
damping_values* lookup(double kmode);

/*
 * If the kmode is already in the hashtable, do nothing.
 * Else, add it at its hash position, first in the list.
 *
 * @param kmode
 * @param values
 */
void add_pair(double kmode, damping_values values);

/*
 * Put the theoretical values inside the hashtable.
 *
 * er, psi values for 0.2, 0.3, 0.4 & 0.5 from
 *   Sonnendrücker, "Numerical Methods for the Vlasov-Maxwell equations", 2016
 *   4.4.2. Landau damping, p. 54
 * omega_real, omega_imag for 2 * PI / 22 from Zeal / Maple.
 * omega_real, omega_imag for 0.25 to 2.0 by 0.05 & 2.3088 from
 *   Canosa, "Numerical Solution of Landau's Dispersion Equation", 1973
 *   p. 159
 */
void init_damping_table();

// When calling with 1 argument, the macro will set compare_values to 0.
#define get_damping_values_1(a   ) a, 0
#define get_damping_values_2(a, b) a, b
#define get_damping_values(...) VARIADIC(get_damping_values, NUMARG2(__VA_ARGS__), __VA_ARGS__)

/*
 * Look for a given kmode in the hashmap.
 * If not found, applies a formula.
 *
 * @param  kmode
 * @param  compare_values set to non-0 if you want to print the comparison of formulae. Defaults to 0.
 * @return the damping values.
 */
damping_values* get_damping_values(double kmode, int compare_values);



/*****************************************************************************
 *                           L2-norm computation                             *
 *****************************************************************************/

/*
 * @param[in] mesh, the mesh on which we're working.
 * @param[in] e[ncx+1] the electric field on any axis.
 * @return    the electric energy of e.
 */
double normL2_field_1d(cartesian_mesh_1d mesh1d, double* e);

/*
 * @param[in] mesh, the mesh on which we're working.
 * @param[in] e[ncx+1][ncy+1] the electric field on any axis.
 * @return    the electric energy of e.
 */
double normL2_field_2d(cartesian_mesh_2d mesh2d, double** e);

/*
 * @param[in] mesh2d the spatial mesh.
 * @param[in] mode_x the mode to be computed.
 * @param[in] ex the electric field on the x-axis.
 * @return    the Fourier mode of the electric field
 */
double fourier_mode_1d(cartesian_mesh_2d mesh2d, int mode_x, double** ex);

/*
 * @param[in] mesh2d the spatial mesh.
 * @param[in] mode_x the mode to be computed.
 * @param[in] ex the electric field on the x-axis.
 * @return    the energy of the Fourier mode of the electric fields.
 */
double energy_fourier_mode_1d(cartesian_mesh_2d mesh2d, int mode_x, double** ex);

/*
 * @param[in] mesh, the mesh on which we're working.
 * @param[in] e[ncx+1][ncy+1][ncz+1] the electric field on any axis.
 * @return    the electric energy of e.
 */
double normL2_field_3d(cartesian_mesh_3d mesh, double*** e);

#endif // ifndef PIC_VERT_DIAGNOSTICS
