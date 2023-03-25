# reference:
#  J. K. Moore, "Human Control of a Bicycle"
#   PhD Thesis, University of California Davis

import sympy as sp
import sympy.physics.mechanics as mec

# coordinates
q1,q2,q3,q4 = mec.dynamicsymbols('q1 q2 q3 q4')
q5,q6,q7,q8 = mec.dynamicsymbols('q5 q6 q7 q8')

# bicycle orientation
N = mec.ReferenceFrame('N')
A = N.orientnew('A', 'Axis', (q3, N.z))# rear frame yaw
B = A.orientnew('B', 'Axis', (q4, A.x))# rear frame roll
C = B.orientnew('C', 'Axis', (q5, B.y))# rear frame pitch
D = C.orientnew('D', 'Axis', (q6, C.y))# rear whell spin
E = C.orientnew('E', 'Axis', (q7, C.z))# front frame steer
F = E.orientnew('F', 'Axis', (q8, E.y))# front wheel spin

# bicyle dimension
rf, rr = sp.symbols('rf rr') # wheel radius
d1, d2, d3 = sp.symbols('d1 d2 d3') # frame dimension
l1, l2, l3, l4 = sp.symbols('l1 l2 l3 l4') # position of frame centers
f_on = E.y.cross(A.z).cross(E.y).normalize() # unit vector from fo to fn

P = mec.Point('P')# point on ground touching rear wheel
dn = P.locatenew('dn', 0)# point on rear wheel touching ground
do = dn.locatenew('do', -rr*B.z)# rear wheel center
co = do.locatenew('co', l1 * C.x + l2 * C.z)# rear frame center
ce = do.locatenew('ce', d1 * C.x)# steer axis point
fo = ce.locatenew('fo', d2 * E.z + d3 * E.x)# front wheel center
eo = fo.locatenew('eo', l3 * E.x + l4 * E.z)# front frame center
fn = fo.locatenew('fn', rf*f_on)# point on front wheel touching ground

cons1 = fn.pos_from(dn).dot(A.z) # constraint among q5,q4,q7

# velocities
u1,u2,u3,u4 = mec.dynamicsymbols('u1 u2 u3 u4')
u5,u6,u7,u8 = mec.dynamicsymbols('u5 u6 u7 u8')

t = mec.dynamicsymbols._t # time

kde = [ # kinematical differential equaiton
    q1.diff(t) - u1,  # x translational
    q2.diff(t) - u2,  # y translational
    q3.diff(t) - u3,  # yaw
    q4.diff(t) - u4,  # roll
    q5.diff(t) - u5,  # pitch
    q7.diff(t) - u7]  # steer

# angular velocities
A.set_ang_vel(N, u3 * N.z)  # yaw rate
B.set_ang_vel(A, u4 * A.x)  # roll rate
C.set_ang_vel(B, u5 * B.y)  # pitch rate
D.set_ang_vel(C, u6 * C.y)  # rear wheel rate
E.set_ang_vel(C, u7 * C.z)  # steer rate
F.set_ang_vel(E, u8 * E.y)  # front wheel rate

# set velocities
P.set_vel(N, u1*N.x + u2*N.y)
do.v2pt_theory(P, N, B)
dn.v2pt_theory(do, N, D)# rear wheel contact velocity
co.v2pt_theory(do, N, C)
ce.v2pt_theory(do, N, C)
fo.v2pt_theory(ce, N, E)
eo.v2pt_theory(fo, N, E)
fn.v2pt_theory(fo, N, F)# front wheel contact velocity

# wheel contact does not slip
cons2 = [dn.vel(N).dot(a) for a in A][:2]# rear wheel
cons2+= [fn.vel(N).dot(a) for a in A]# front wheel

# mass and inertia
mc, md, me, mf = sp.symbols('mc md me mf') # mass
ic11, ic22, ic33, ic31 = sp.symbols('ic11 ic22 ic33 ic31') # rear frame
id11, id22 = sp.symbols('id11 id22') # rear wheel
ie11, ie22, ie33, ie31 = sp.symbols('ie11 ie22 ie33 ie31') # front frame
if11, if22 = sp.symbols('if11 if22') # front wheel

Ic = mec.inertia(C, ic11, ic22, ic33, 0, 0, ic31)
Id = mec.inertia(C, id11, id22, id11)
Ie = mec.inertia(E, ie11, ie22, ie33, 0, 0, ie31)
If = mec.inertia(E, if11, if22, if11)

RF = mec.RigidBody('Rear Frame', co, C, mc, (Ic, co))
RW = mec.RigidBody('Rear Wheel', do, D, md, (Id, do))
FF = mec.RigidBody('Front Frame', eo, E, me, (Ie, eo))
FW = mec.RigidBody('Front Wheel', fo, F, mf, (If, fo))

# gravity
g = sp.symbols('g') # acceleration of gravity
Fco = (co, mc*g*A.z)
Fdo = (do, md*g*A.z)
Feo = (eo, me*g*A.z)
Ffo = (fo, mf*g*A.z)

# torque
T4, T6, T7 = mec.dynamicsymbols('T4 T6 T7')
Tc = (C, T4*A.x - T6*B.y - T7*C.z)
Td = (D, T6*C.y)
Te = (E, T7*C.z)

# generate equation of motion
KM = mec.KanesMethod(
    N,
    [q1,q2,q3,q4,q7],# independent coordinates
    [u4,u6,u7],  # independent velocities
    kde, # kinematical differential equation
    [q5],  # dependent coordinates
    [cons1], # holonomic constraints
    [u1,u2,u3,u5,u8],  # dependent velocities
    cons2) # nonholonomic constraints

KM.kanes_equations(
    bodies=[RF, RW, FF, FW],
    loads=[Fco, Fdo, Feo, Ffo, Tc, Td, Te])
