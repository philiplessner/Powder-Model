from __future__ import print_function, unicode_literals, division
from pprint import pprint
from sympy import *


init_printing()
r, r0, L, t0, d, Ds, rhoTa, z = symbols('r, r0, L, t0, d, Ds, rhoTa, z')

print('\n', 'Neck diameter (d) as a function of particle radius (r)',
      '\n', 'particle thickness (t0), and neck length (L)', '\n')
d = (t0 - 2 * r) * (L / t0) + 2 * r
pprint(d)

print('\n', 'f(r)=0 and r as a function of L, d, t0,',
      '\n', 'and initial particle radius (r0)', '\n')
eq1 = sqrt(r0 ** 2 - (L * d ** 2) / (4 * t0)) - r
pprint(eq1)

print('\n', 'Solve f(r) for r. Take the positive root.', '\n')
eq2 = solve(eq1, r)
pprint(eq2[0])
print(eq2[0])

# Make a temporary variable z before defining L
eq3 = eq2[0].subs(L, z)

print('\n', 'L as a function of r0, r, t0, density of Ta (rhoTa),',
      '\n', 'sintered density of anode (Ds)', '\n')
L = (pi * r0 ** 2 * t0 * rhoTa) / (4 * r ** 2 * Ds) - t0
pprint(L)

print('\n', 'Substitute for L in equation for r to obtain an',
      '\n', 'f(r) = g(r, r0, t0, Ds)-r', '\n')
eq4 = eq3.subs(z, L)
eq5 = eq4 - r
pprint(eq5)
print(eq5)

print('\n', 'df(r)/dr', '\n')
eq6 = eq5.diff(r)
pprint(eq6)
print(eq6)
