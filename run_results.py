"""Run nominal, 49-case sweep, sensitivity studies; save Figures 4-6."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from src.model import *
from src.integrators import integrate, rotation
from src.metrics import angle_error, measures
from scripts.plot_utils import save_figure

def main():
    t,u=integrate(lambda u,t:rhs(u,t,P0),u0,24.,.01)
    Poff=parameters(active=False)
    _,uoff=integrate(lambda u,t:rhs(u,t,Poff),u0,24.,.01)
    nominal=measures(t,u,P0)
    print('Nominal [settling (s), peak error (deg), effort (rad^2 s), peak fin (deg)]:')
    print(nominal)
    print('Zero-command final angle error (deg):',np.rad2deg(angle_error(uoff)[-1]))
    assert abs(nominal[0]-7.27)<.011
    fig,ax=plt.subplots(2,2,figsize=(11,7))
    for i,name in enumerate(['Roll','Pitch','Yaw']):
        ax[0,0].plot(t,np.rad2deg(u[:,i]),label=name)
        ax[0,1].plot(t,np.rad2deg(u[:,i+3]),label=name)
    for j in range(4):ax[1,0].plot(t,np.rad2deg(u[:,j+6]),label=f'Fin {j+1}')
    for sign in [-1,1]:ax[1,0].axhline(sign*15,color='k',ls='--')
    ax[1,1].plot(t,np.rad2deg(angle_error(u)),label='Active fins')
    ax[1,1].plot(t,np.rad2deg(angle_error(uoff)),label='Zero commands')
    ax[1,1].axhline(1,color='k',ls=':')
    for x,ylabel in zip(ax.flat,['Angle (deg)','Body rate (deg/s)',
                                 'Fin angle (deg)','Total rotation error (deg)']):
        x.set(xlabel='Time (s)',ylabel=ylabel);x.legend();x.grid(alpha=.25)
    plt.tight_layout();save_figure("figure-04.png")

    fig=plt.figure(figsize=(7,5));ax=fig.add_subplot(111,projection='3d')
    for time,color in zip([0.,2.,5.,12.],['C0','C1','C2','C3']):
        R=rotation(u[int(round(time/.01)),:3])
        nose=R@np.array([1.,0.,0.]);side=R@np.array([0.,.24,0.])
        ax.plot([0,nose[0]],[0,nose[1]],[0,nose[2]],color=color,label=f't={time:g} s')
        ax.plot([0,side[0]],[0,side[1]],[0,side[2]],'--',color=color)
    ax.set(xlabel='Reference x',ylabel='Reference y',zlabel='Reference z',
           xlim=(-.15,1.05),ylim=(-.6,.6),zlim=(-.6,.6),
           title='Orientation about a fixed center of mass')
    ax.set_box_aspect([1,1,1]);ax.legend();save_figure("figure-05.png")

    limits=np.array([5.,7.5,10.,12.5,15.,20.,25.])
    taus=np.array([.05,.1,.15,.3,.5,.8,1.2])
    results=np.zeros((len(taus),len(limits),4))
    refinement=np.zeros((len(taus),len(limits),3))
    for i,tau in enumerate(taus):
        for j,limit in enumerate(limits):
            P=parameters(tau=tau,limit=limit)
            t,u=integrate(lambda u,t:rhs(u,t,P),u0,24.,.01)
            results[i,j]=measures(t,u,P)
            th,uh=integrate(lambda u,t:rhs(u,t,P),u0,24.,.005)
            fine=measures(th,uh,P)
            refinement[i,j]=[
                np.rad2deg(np.max(np.abs(u[:,:3]-uh[::2,:3]))),
                abs(results[i,j,0]-fine[0]),abs(results[i,j,2]-fine[2])]
            assert refinement[i,j,0]<.01
            assert np.isfinite(results[i,j,0])==np.isfinite(fine[0])
            assert np.max(np.abs(u[:,6:]))<=P['limit']+1.e-8
            assert np.max(np.abs(u[:,1]))<np.deg2rad(60)
    fig,ax=plt.subplots(1,2,figsize=(12,4.5))
    for x,index,title,label in zip(ax,[0,2],['Recovery time','Fin motion effort'],
        ['Settling time (s)','Squared-deflection effort (rad^2 s)']):
        im=x.imshow(np.ma.masked_invalid(results[:,:,index]),origin='lower',aspect='auto')
        x.set_xticks(range(len(limits)),[f'{v:g}' for v in limits])
        x.set_yticks(range(len(taus)),[f'{v:g}' for v in taus])
        x.set(xlabel='Maximum fin deflection (deg)',ylabel='Time constant (s)',title=title)
        fig.colorbar(im,ax=x,label=label)
    plt.tight_layout();save_figure("figure-06.png")
    for tau,limit in [(.15,5.),(.15,15.),(.15,25.),(.8,15.),(1.2,15.),(1.2,5.)]:
        i=np.flatnonzero(taus==tau)[0];j=np.flatnonzero(limits==limit)[0]
        print(f'tau={tau:.2f}, limit={limit:g}:',results[i,j,:3])
    print('Recovered cases:',np.count_nonzero(np.isfinite(results[:,:,0])),'of 49')
    print('Maximum refinement changes [angle (deg), settling (s), effort (rad^2 s)]:')
    print(np.max(refinement,axis=(0,1)))

    print('Parameter       factor   settling (s)   effort (rad^2 s)')
    for name,values in [('force_scale',[.5,1.,1.5]),('inertia_scale',[.75,1.,1.25])]:
        for value in values:
            P=parameters(**{name:value})
            t,u=integrate(lambda u,t:rhs(u,t,P),u0,24.,.01)
            output=measures(t,u,P)
            print(f'{name:14s} {value:5.2f}    {output[0]:8.2f}      {output[2]:.3f}')

if __name__ == '__main__':
    main()
