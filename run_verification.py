"""Run independent checks, stability and convergence studies; save Figures 1-3."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from src.model import *
from src.integrators import integrate, rotation
from src.stability import jacobian, stability, step_limit
from scripts.plot_utils import save_figure

def main():
    columns=[]
    for beta in [0,np.pi/2,np.pi,3*np.pi/2]:
        position=np.array([ell,a*np.cos(beta),a*np.sin(beta)])
        force=Kf*np.array([0,-np.sin(beta),np.cos(beta)])
        columns.append(np.cross(position,force))
    assert np.allclose(np.array(columns).T,B0)
    assert np.linalg.matrix_rank(B0)==3

    hs=np.array([.12,.06,.03,.015,.0075])
    errors=np.zeros((3,len(hs)))
    for axis in range(3):
        P=parameters(tau=.3,active=False,damping=False)
        start=np.zeros(10)
        start[3+axis]=.01
        start[6:]=.01*G[:,axis]
        a0=(P['B']@start[6:])[axis]/P['I'][axis]
        for j,h in enumerate(hs):
            t,u=integrate(lambda u,t:rhs(u,t,P),start,1.2,h)
            decay=np.exp(-t/P['tau'])
            exact=np.zeros_like(u)
            exact[:,6:]=decay[:,None]*start[6:]
            exact[:,3+axis]=.01+a0*P['tau']*(1-decay)
            exact[:,axis]=.01*t+a0*P['tau']*(t-P['tau']*(1-decay))
            errors[axis,j]=np.max(np.abs(u-exact))
    orders=np.log2(errors[:,-2]/errors[:,-1])
    print('Exact-solution orders:',orders)
    assert np.all((orders>3.8)&(orders<4.3))
    for i,label in enumerate(['Roll','Pitch','Yaw']):
        plt.loglog(hs,errors[i],'o-',label=label)
    plt.loglog(hs,errors[0,0]*(hs/hs[0])**4,'k--',label='Slope 4')
    plt.xlabel('Time step (s)');plt.ylabel('Maximum scaled state error')
    plt.legend();plt.grid(True,which='both',alpha=.25);save_figure("figure-01.png")

    P=parameters(active=False,damping=False,inertia=[2.,3.,4.])
    start=np.r_[np.deg2rad([5.,-8.,6.]),[.13,.19,.23],np.zeros(4)]
    t,u=integrate(lambda u,t:rhs(u,t,P),start,4.,.01)
    energy=.5*np.sum(P['I']*u[:,3:6]**2,axis=1)
    momentum=np.array([rotation(row[:3])@(P['I']*row[3:6]) for row in u])
    energy_drift=np.max(np.abs(energy-energy[0]))/energy[0]
    momentum_drift=np.max(np.linalg.norm(momentum-momentum[0],axis=1))
    momentum_drift/=np.linalg.norm(momentum[0])
    print('Relative energy drift:',energy_drift)
    print('Relative inertial angular-momentum drift:',momentum_drift)
    assert energy_drift<1.e-7 and momentum_drift<1.e-7

    lam=np.linalg.eigvals(jacobian(P0))
    eps=1.e-7;basis=np.eye(10)
    Acheck=np.column_stack([(rhs(eps*e,0,P0)-rhs(-eps*e,0,P0))/(2*eps)
                            for e in basis])
    assert np.allclose(Acheck,jacobian(P0),rtol=1.e-7,atol=1.e-8)
    print('Nominal local step limit:',step_limit(P0))
    print('Fast-actuator local step limit:',step_limit(parameters(tau=.05)))
    assert np.max(lam.real)<0 and np.max(np.abs(stability(.01*lam)))<1
    delays=np.linspace(.03,3.,250)
    growth=np.array([np.max(np.linalg.eigvals(jacobian(parameters(tau=tau))).real)
                     for tau in delays])
    print('First locally unstable delay in scan:',delays[np.flatnonzero(growth>=0)[0]])
    fig,ax=plt.subplots(1,2,figsize=(11,4))
    ax[0].plot(delays,growth);ax[0].axhline(0,color='k',ls='--')
    ax[0].set(xlabel='Actuator time constant (s)',ylabel='Largest real eigenvalue (1/s)')
    xx,yy=np.meshgrid(np.linspace(-3.2,.3,250),np.linspace(-3.2,3.2,250))
    ax[1].contour(xx,yy,np.abs(stability(xx+1j*yy)),levels=[1],colors='k')
    for h in [.01,.5]:
        z=h*lam;ax[1].scatter(z.real,z.imag,label=f'h = {h} s')
    ax[1].set(xlabel='Real(h lambda)',ylabel='Imag(h lambda)');ax[1].legend()
    plt.tight_layout();save_figure("figure-02.png")

    tf,href=24.,.00125
    tr,ur=integrate(lambda u,t:rhs(u,t,P0),u0,tf,href)
    steps=np.array([.08,.04,.02,.01,.005])
    angle_errors={};work={}
    for method in ['Euler','RK4']:
        angle_errors[method]=[];work[method]=[]
        for h in steps:
            t,u=integrate(lambda u,t:rhs(u,t,P0),u0,tf,h,method)
            ref=ur[::int(round(h/href))]
            error=np.rad2deg(np.max(np.abs(u[:,:3]-ref[:,:3])))
            evaluations=int(round(tf/h))*(4 if method=='RK4' else 1)
            angle_errors[method].append(error);work[method].append(evaluations)
            print(f'{method:5s}  h={h:.3f}  error={error:.3e} deg  calls={evaluations}')
    full_order=np.log2(angle_errors['RK4'][1]/angle_errors['RK4'][2])
    print('Full-model observed order:',full_order)
    assert 3.7<full_order<4.5
    fig,ax=plt.subplots(1,2,figsize=(11,4))
    for method in ['Euler','RK4']:
        ax[0].loglog(steps,angle_errors[method],'o-',label=method)
        ax[1].loglog(work[method],angle_errors[method],'o-',label=method)
    ax[0].set(xlabel='Time step (s)',ylabel='Maximum angle difference (deg)')
    ax[1].set(xlabel='Derivative evaluations',ylabel='Maximum angle difference (deg)')
    for x in ax:x.legend();x.grid(True,which='both',alpha=.25)
    plt.tight_layout();save_figure("figure-03.png")
    for tau,limit in [(.05,25.),(1.2,5.)]:
        P=parameters(tau=tau,limit=limit)
        t,u=integrate(lambda u,t:rhs(u,t,P),u0,24.,.01)
        th,uh=integrate(lambda u,t:rhs(u,t,P),u0,24.,.005)
        error=np.rad2deg(np.max(np.abs(u[:,:3]-uh[::2,:3])))
        print(f'Corner tau={tau}, limit={limit}: {error:.3e} deg')
        assert error<.01

if __name__ == '__main__':
    main()
