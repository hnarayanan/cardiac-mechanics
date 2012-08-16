from sympy import *

x0, dSLdt, phi, XBPreR, XBPostR, xXBPreR, xXBPostR, XBDutyFractPreR, XBDutyFractPostR, fappT, hbT, hft = symbols('x0, dSLdt, phi, XBPreR, XBPostR, xXBPreR, xXBPostR, XBDutyFractPreR, XBDutyFractPostR, fappT, hbT, hft')

eq1 = Rational(1, 2)*dSLdt + phi/XBDutyFractPreR*(fappT*(-xXBPreR) + hbT*(xXBPostR - x0 - xXBPreR))
eq2 = Rational(1, 2)*dSLdt + phi/XBDutyFractPostR*(hft*(xXBPreR + x0 - xXBPostR))

print(latex(solve([eq1, eq2], set([xXBPreR, xXBPostR]))))
