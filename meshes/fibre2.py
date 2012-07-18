from sympy import *

# Perhaps x => y and y => x

x, y, z = symbols('coordsX, coordsY, coordsZ')
phi, theta = symbols('phi, theta')


#K_T, K_N, K_B = symbols('K_T, K_N, K_B')
E_T = Matrix([cos(phi)*cos(theta), cos(phi)*sin(theta), sin(phi)])
E_N = Matrix([-sin(theta), cos(theta), 0])
E_B = Matrix([-sin(phi)*cos(theta), -sin(phi)*sin(theta), cos(phi)])

K_B = -0.17
K_N =  0.003
K_T = -0.005

substitutions = {phi: 0, theta: atan((K_T*x + K_N*y)/(1 + K_N*x - K_T*y)) + K_B*z}

E_T = E_T.subs(substitutions)
E_T_string = "%s*iHat + %s*jHat + %s*kHat" % (E_T[0], E_T[1], E_T[2])
print E_T_string

print atan((K_T*x + K_N*y)/(1 + K_N*x - K_T*y)) + K_B*z



#expression = "%s*iHat + %s*jHat %s*kHat" % (diff(theta, x), diff(theta, y), diff(theta, z))

