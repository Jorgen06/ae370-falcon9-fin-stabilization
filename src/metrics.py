"""Orientation-error and fin-effort measures from the original notebook."""
import numpy as np
from .integrators import rotation

def angle_error(u):
    matrices=np.array([rotation(angles) for angles in u[:,:3]])
    trace=np.trace(matrices,axis1=1,axis2=2)
    return np.arccos(np.clip((trace-1)/2,-1,1))

def measures(t,u,P):
    angle=angle_error(u)
    rate=np.linalg.norm(u[:,3:6],axis=1)
    outside=(angle>np.deg2rad(1))|(rate>np.deg2rad(1))
    bad=np.flatnonzero(outside)
    i=bad[-1]+1 if len(bad) else 0
    settling=t[i] if i<len(t) and t[-1]-t[i]>=2 else np.nan
    squared=np.sum(u[:,6:]**2,axis=1)
    effort=np.sum(.5*(squared[:-1]+squared[1:])*np.diff(t))
    return np.array([settling,np.rad2deg(np.max(angle)),effort,
                     np.rad2deg(np.max(np.abs(u[:,6:])))])
