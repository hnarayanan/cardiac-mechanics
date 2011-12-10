from sympy import *

C_bar_11, C_bar_12, C_bar_13 = symbols('C_bar[0][0] C_bar[0][1] C_bar[0][2]')
C_bar_21, C_bar_22, C_bar_23 = symbols('C_bar[1][0] C_bar[1][1] C_bar[1][2]')
C_bar_31, C_bar_32, C_bar_33 = symbols('C_bar[2][0] C_bar[2][1] C_bar[2][2]')

C_bar_inv_11, C_bar_inv_12, C_bar_inv_13 = symbols('C_bar_inv[0][0] C_bar_inv[0][1] C_bar_inv[0][2]')
C_bar_inv_21, C_bar_inv_22, C_bar_inv_23 = symbols('C_bar_inv[1][0] C_bar_inv[1][1] C_bar_inv[1][2]')
C_bar_inv_31, C_bar_inv_32, C_bar_inv_33 = symbols('C_bar_inv[2][0] C_bar_inv[2][1] C_bar_inv[2][2]')

S_bar_11, S_bar_12, S_bar_13 = symbols('S_bar[0][0] S_bar[0][1] S_bar[0][2]')
S_bar_21, S_bar_22, S_bar_23 = symbols('S_bar[1][0] S_bar[1][1] S_bar[1][2]')
S_bar_31, S_bar_32, S_bar_33 = symbols('S_bar[2][0] S_bar[2][1] S_bar[2][2]')

def delta(a, b):
    if a == b:
        return 1.0
    else:
        return 0.0

I4 = [[[[0.5*(delta(a, c)*delta(b, d) + delta(a, d)*delta(b, c))
         for a in range(3)] for b in range(3)] for c in range(3)] for d in range(3)]

C_bar = [[C_bar_11, C_bar_12, C_bar_13],
         [C_bar_21, C_bar_22, C_bar_23],
         [C_bar_31, C_bar_32, C_bar_33]]

C_bar_inv = [[C_bar_inv_11, C_bar_inv_12, C_bar_inv_13],
             [C_bar_inv_21, C_bar_inv_22, C_bar_inv_23],
             [C_bar_inv_31, C_bar_inv_32, C_bar_inv_33]]

C_bar_inv_outer_C_bar = [[[[C_bar_inv[a][b]*C_bar[c][d]
                            for a in range(3)] for b in range(3)] for c in range(3)] for d in range(3)]

Proj = [[[[I4[a][b][c][d] - (1.0/3.0)*C_bar_inv_outer_C_bar[a][b][c][d]
           for a in range(3)] for b in range(3)] for c in range(3)] for d in range(3)]

S_bar = [[S_bar_11, S_bar_12, S_bar_13],
         [S_bar_21, S_bar_22, S_bar_23],
         [S_bar_31, S_bar_32, S_bar_33]]

Proj_contract_S_bar = [[sum([sum([Proj[a][b][c][d]*S_bar[c][d] for c in range(3)]) for d in range(3)])
                        for a in range(3)] for b in range(3)]

print(Proj_contract_S_bar)
