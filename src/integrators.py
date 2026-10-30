"""From-scratch forward Euler and classical four-stage RK4 time integrators.

No third-party ODE time-marching routine is called.
"""
import numpy as np

def integrate(f,u0,tf,h,method='RK4'):
    n=int(round(tf/h))
    if not np.isclose(n*h,tf):
        raise ValueError('Choose h so it divides tf.')
    t=np.linspace(0,tf,n+1)
    u=np.empty((n+1,len(u0)))
    u[0]=u0
    for k in range(n):
        y1=f(u[k],t[k])
        if method=='Euler':
            u[k+1]=u[k]+h*y1
        elif method=='RK4':
            y2=f(u[k]+h*y1/2,t[k]+h/2)
            y3=f(u[k]+h*y2/2,t[k]+h/2)
            y4=f(u[k]+h*y3,t[k]+h)
            u[k+1]=u[k]+h*(y1+2*y2+2*y3+y4)/6
        else:
            raise ValueError('Choose Euler or RK4.')
    return t,u

def rotation(angles):
    phi,theta,psi=angles
    cp,sp=np.cos(phi),np.sin(phi)
    ct,st=np.cos(theta),np.sin(theta)
    cs,ss=np.cos(psi),np.sin(psi)
    return np.array([[cs*ct,cs*st*sp-ss*cp,cs*st*cp+ss*sp],
        [ss*ct,ss*st*sp+cs*cp,ss*st*cp-cs*sp],[-st,ct*sp,ct*cp]])
