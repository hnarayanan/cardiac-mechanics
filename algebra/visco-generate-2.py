from sympy import *
from sympy.printing import print_ccode

# Declare useful symbols
I1_bar, I2_bar, J, I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar = \
    symbols("I1_bar, I2_bar, J, I4_f_bar, I4_s_bar, I8_fs_bar, I8_fn_bar")

a, b, a_f, b_f, a_s, b_s, a_fs, b_fs, kappa, beta = \
symbols("a, b, a_f, b_f, a_s, b_s, a_fs, b_fs, kappa, beta")

material_parameters = {a: 0.500, b: 8.023, a_f: 16.472, b_f: 16.026, \
                       a_s: 2.481, b_s: 11.120, a_fs: 0.356, b_fs: 11.436, \
                       kappa: 2.0e6, beta:9.0}

def print_header():
    print ("""
#include<iostream>
#include<math.h>

using namespace std;

int main(void)
{""")

def print_footer():
    print("""
    cout << S_iso_12 + S_vol_12 << endl;
    return 0;
}""")


print_header()

for parameter in material_parameters.keys():
    print("double %s\t= %g;" % (parameter.name, material_parameters[parameter]))

f11, f12, f13, f21, f22, f23, f31, f32, f33 = \
symbols("f11, f12, f13, f21, f22, f23, f31, f32, f33")

F = Matrix([[f11, f12, f13],
            [f21, f22, f23],
            [f31, f32, f33]])

C = F.T*F

for i in range(9):
    print("double %s = %g; " % (F[i], 0.0))

print_ccode(F.det(), assign_to="double J")

# Strain energy function in terms of the invariants of the right
# Cauchy-Green tensor
psi_iso =  a/(2*b)*exp(b*(I1_bar - 3)) \
         + a_f/(2*b_f)*(exp(b_f*(I4_f_bar - 1)**2) - 1) \
         + a_s/(2*b_s)*(exp(b_s*(I4_s_bar - 1)**2) - 1) \
         + a_fs/(2*b_fs)*(exp(b_fs*I8_fs_bar**2) - 1)
psi_vol = kappa*(1/(beta**2)*(beta*ln(J) + 1/(J**beta) - 1))

# Identity Matrix
I = eye(3)

# Reference fibre, sheet and sheet-normal directions
f0 = Matrix([1, 0, 0])
s0 = Matrix([0, 1, 0])
n0 = Matrix([0, 0, 1])

# Compute the elastic components of the stress given the deformation
# gradient
def elastic_stresses(C):

    # Modified right Cauchy-Green tensor
    C_bar = J**(-Rational(2, 3))*(C)

    # Define the second Piola-Kirchhoff stress in terms of the invariants
    S_bar =   2*(diff(psi_iso, I1_bar) + diff(psi_iso, I2_bar))*I \
            - 2*diff(psi_iso, I2_bar)*C_bar \
            + 2*diff(psi_iso, I4_f_bar)*(f0*f0.T) \
            + 2*diff(psi_iso, I4_s_bar)*(s0*s0.T) \
            + diff(psi_iso, I8_fs_bar)*(f0*s0.T + s0*f0.T) \
            + diff(psi_iso, I8_fn_bar)*(f0*n0.T + n0*f0.T)
    S_bar_contract_C = sum([sum([S_bar[a, b]*C[a, b]
                                 for a in range(3)]) for b in range(3)])
    Dev_S_bar = S_bar - Rational(1, 3)*(S_bar_contract_C)*C.inv()

    S_iso_inf = J**(-Rational(2, 3))*Matrix(Dev_S_bar)
    S_vol_inf = J*diff(psi_vol, J)*C.inv()

    # Substitute the current values of the invariants
    substitutions = {I8_fn_bar: (f0.T*C_bar*n0)[0]}
    print_ccode(C_bar.trace(), assign_to="double I1_bar")
    print_ccode(Rational(1, 2)*(I1_bar*I1_bar - (C_bar*C_bar).trace()), assign_to="double I2_bar")
    print_ccode((f0.T*C_bar*f0)[0], assign_to="double I4_f_bar")
    print_ccode((s0.T*C_bar*s0)[0], assign_to="double I4_s_bar")
    print_ccode((f0.T*C_bar*s0)[0], assign_to="double I8_fs_bar")
    print_ccode((f0.T*C_bar*n0)[0], assign_to="double I8_fn_bar")

    return S_vol_inf, S_iso_inf

S_vol_inf, S_iso_inf = elastic_stresses(C)


for i in range(3):
    for j in range(3):
        print_ccode(S_iso_inf[i, j], assign_to="double S_iso_%i%i" % (i + 1, j + 1))

for i in range(3):
    for j in range(3):
        print_ccode(S_vol_inf[i, j], assign_to="double S_vol_%i%i" % (i + 1, j + 1))


print_footer()
