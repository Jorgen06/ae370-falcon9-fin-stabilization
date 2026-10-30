"""Ten-state Euler-angle/rigid-body/fin-actuator model.

Preserves the governing equations and parameters from the submitted group notebook.
The vehicle and aerodynamic coefficients are assumptions in a Falcon 9-inspired model.
"""
import numpy as np

m, L, rb = 25000., 40., 1.85
a, ell = 2.6, 15.
rho, V, Af, CNd = 0.9, 200., 1.5, 2.
I0 = np.array([m*rb**2, m*(L**2/12 + rb**2/2),
               m*(L**2/12 + rb**2/2)])
D0 = np.array([1.e4, 3.e5, 3.e5])
Kf = 0.5*rho*V**2*Af*CNd
B0 = Kf*np.array([[a,a,a,a],[-ell,0,ell,0],[0,-ell,0,ell]])
G = np.array([[1,-1,0],[1,0,-1],[1,1,0],[1,0,1]])
wn = np.array([1.3,0.7,0.7])
Kp, Kd = I0*wn**2, 2*I0*wn
H = G / (Kf*np.array([4*a,2*ell,2*ell]))

def parameters(tau=.15, limit=15., force_scale=1., inertia_scale=1.,
               active=True, damping=True, inertia=None):
    return {'I': I0*inertia_scale if inertia is None else np.array(inertia),
            'D': D0.copy() if damping else np.zeros(3),
            'B': B0*force_scale, 'tau': tau,
            'limit': np.deg2rad(limit), 'active': active}

def rhs(u,t,P):
    phi,theta,psi,p,q,r = u[:6]
    s,c,ct = np.sin(phi),np.cos(phi),np.cos(theta)
    if abs(ct)<.15:
        raise ValueError('Motion is too close to the Euler pitch singularity.')
    angle_dot = np.array([p+(s*q+c*r)*np.tan(theta),
                          c*q-s*r,(s*q+c*r)/ct])
    w=u[3:6]
    moment=P['B']@u[6:]-P['D']*w
    rate_dot=(moment-np.cross(w,P['I']*w))/P['I']
    command=np.zeros(4)
    if P['active']:
        requested=-Kp*u[:3]-Kd*w
        command=P['limit']*np.tanh((H@requested)/P['limit'])
    fin_dot=(command-u[6:])/P['tau']
    return np.r_[angle_dot,rate_dot,fin_dot]

u0=np.r_[np.deg2rad([12.,-10.,8.,4.,-3.,2.]),np.zeros(4)]
P0=parameters()
