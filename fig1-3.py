# data from:
# https://pydy.readthedocs.io/en/stable/examples/carvallo-whipple.html

from pydy.system import System
from bicycle import (
    KM,cons1,rf,rr,d1,d2,d3,l1,l2,l3,l4,mc,md,me,mf,g,
    q1,q2,q3,q4,q5,q7,u1,u2,u3,u4,u5,u6,u7,u8,
    ic11,ic22,ic31,ic33,id11,id22,
    ie11,ie22,ie31,ie33,if11,if22)

sys = System(KM)
sys.generate_ode_function(generator='cython')

# set constants
sys.constants = {
   rf: 0.35,
   rr: 0.3,
   d1: 0.9534570696121849,
   d3: 0.03207142672761929,
   d2: 0.2676445084476887,
   l1: 0.4707271515135145,
   l2: -0.47792881146460797,
   l4: -0.3699518200282974,
   l3: -0.00597083392418685,
   mc: 85.0,
   md: 2.0,
   me: 4.0,
   mf: 3.0,
   id11: 0.0603,
   id22: 0.12,
   if11: 0.1405,
   if22: 0.28,
   ic11: 7.178169776497895,
   ic22: 11.0,
   ic31: 3.8225535938357873,
   ic33: 4.821830223502103,
   ie11: 0.05841337700152972,
   ie22: 0.06,
   ie31: 0.009119225261946298,
   ie33: 0.007586622998470264,
   g: 9.81}

# fig1 ===========================================================
import sympy as sp
from scipy.optimize import newton

# set initial conditions
u1init = 4.6  # initial speed / m/s
u4init = 0.5  # initial roll rate / rad/s
q7init = 1e-8 # initial steering angle / rad

cons1 = cons1.subs(sys.constants)
cons1 = sp.lambdify((q5,q4,q7), cons1)
q5init = newton(cons1, 0, args=(0, q7init))

sys.initial_conditions = {
    q1: 0,
    q2: 0,
    q3: 0,
    q4: 0,
    q5: q5init,
    q7: q7init,
    u1: u1init,
    u2: 0,
    u3: 0,
    u4: u4init,
    u5: 0,
    u6: -u1init/sys.constants[rr],
    u7: 0,
    u8: -u1init/sys.constants[rf]}

# numerical integration
import numpy as np

sys.times = np.linspace(0, 5, 300) # seconds
y = sys.integrate().T # q1,...,q4,q7,q5,u4,u6,u7,u1,u2,u3,u5,u8
c = cons1(y[5],y[3],y[4]) # check constraint
print(np.max(np.abs(c)))

# plot graph
import matplotlib.pyplot as plt

label = [r'$q_1$  / m', r'$q_2$  / m', r'$q_3$  / rad',
         r'$q_4$  / rad', r'$q_7$  / rad', r'$q_5$  / rad']

plt.figure(figsize=(5,10))

for i,q in enumerate(y[:6]):
    plt.subplot(6,1,i+1)
    plt.plot(sys.times, q)
    plt.ylabel(label[i])

plt.xlabel(r'$t$ = time  / sec')
plt.tight_layout()
plt.savefig('fig1.eps')
plt.show()
plt.close()

# fig2 ===========================================================
u1init = 4.1  # initial speed / m/s

sys.initial_conditions[u1] = u1init
sys.initial_conditions[u6] = -u1init/sys.constants[rr]
sys.initial_conditions[u8] = -u1init/sys.constants[rf]

y = sys.integrate().T # q1,...,q4,q7,q5,u4,u6,u7,u1,u2,u3,u5,u8
c = cons1(y[5],y[3],y[4]) # check constraint
print(np.max(np.abs(c)))

plt.figure(figsize=(5,10))

for i,q in enumerate(y[:6]):
    plt.subplot(6,1,i+1)
    plt.plot(sys.times, q)
    plt.ylabel(label[i])

plt.xlabel(r'$t$ = time  / sec')
plt.tight_layout()
plt.savefig('fig2.eps')
plt.show()
plt.close()

# fig3 ===========================================================
u1init = 6  # initial speed / m/s

sys.initial_conditions[u1] = u1init
sys.initial_conditions[u6] = -u1init/sys.constants[rr]
sys.initial_conditions[u8] = -u1init/sys.constants[rf]

y = sys.integrate().T # q1,...,q4,q7,q5,u4,u6,u7,u1,u2,u3,u5,u8
c = cons1(y[5],y[3],y[4]) # check constraint
print(np.max(np.abs(c)))

plt.figure(figsize=(5,10))

for i,q in enumerate(y[:6]):
    plt.subplot(6,1,i+1)
    plt.plot(sys.times, q)
    plt.ylabel(label[i])

plt.xlabel(r'$t$ = time  / sec')
plt.tight_layout()
plt.savefig('fig3.eps')
plt.show()
