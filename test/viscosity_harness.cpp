// Exposes calc_visco_ponto of the C++ Mandyoc to test_rheology.jl (which compiles this file).
#include CALC_VISC
int rheol, geoq_on, tcont; double visc_MAX, visc_MIN, visc_MAX_comp, visc_MIN_comp, visco_r, Delta_T, h_air;
PetscInt WITH_NON_LINEAR, PLASTICITY, pressure_in_rheol; PetscReal pressure_const, visc0_scaled = 1, h0_scaled = 1, strain_rate0_scaled = 1, pressure0_scaled = 1;
PetscScalar weakening_min, weakening_max, *weakening_seed, *cohesion_min, *cohesion_max, *friction_angle_min, *friction_angle_max;
PetscBool seed_layer_set, weakening_from_interfaces_file;
static double c0, c1, f0, f1;
extern "C" void set_params(int r, int g, int nl, int pl, int pr, double *d) {
  rheol = r; geoq_on = g; WITH_NON_LINEAR = nl; PLASTICITY = pl; pressure_in_rheol = pr;
  visc_MAX = d[0]; visc_MIN = d[1]; visc_MAX_comp = d[2]; visc_MIN_comp = d[3]; visco_r = d[4]; Delta_T = d[5]; h_air = d[6];
  pressure_const = d[7]; weakening_min = d[8]; weakening_max = d[9]; c0 = d[10]; c1 = d[11]; f0 = d[12]; f1 = d[13];
  cohesion_min = &c0; cohesion_max = &c1; friction_angle_min = &f0; friction_angle_max = &f1; }
extern "C" double visc(double T, double P, double z, double g, double e2, double s, double A, double n, double Q, double V) {
  return calc_visco_ponto(T, P, 0.0, z, g, e2, s, A, n, Q, V, 0); }
