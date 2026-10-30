"""Linearized model and RK4 absolute-stability limit."""
import numpy as np
from .model import H, Kp, Kd

def jacobian(P):
    A=np.zeros((10,10))
    A[:3,3:6]=np.eye(3)
    A[3:6,3:6]=-np.diag(P['D']/P['I'])
    A[3:6,6:]=P['B']/P['I'][:,None]
    if P['active']:
        A[6:,:3]=-H*Kp/P['tau']
        A[6:,3:6]=-H*Kd/P['tau']
    A[6:,6:]=-np.eye(4)/P['tau']
    return A

def stability(z):
    return 1+z+z**2/2+z**3/6+z**4/24

def step_limit(P):
    lam=np.linalg.eigvals(jacobian(P))
    if np.max(lam.real)>=0:return np.nan
    lo,hi=0.,1.
    while np.max(np.abs(stability(hi*lam)))<=1:hi*=2
    for _ in range(60):
        mid=(lo+hi)/2
        if np.max(np.abs(stability(mid*lam)))<=1:lo=mid
        else:hi=mid
    return lo
