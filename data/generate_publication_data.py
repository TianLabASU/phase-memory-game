"""Seeded Python translation of the publication's MATLAB Figure 1 simulation.

Source: TianLab-ASU/PhaseSeparation_Growth_CircuitMemory,
Stochastic Simulation/Fig1_LLPS_StochasticSim_CellGrowth.m and dependencies.
The browser consumes the downsampled JSON produced by this script.
"""
from __future__ import annotations
import json, math, random
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

VT0=4e-19; VTF=1e-25; OMEGA=6e20; DROP_FACTOR=5e-4
GAMMA=1e-6; KBT=1.38e-23*305; MU0=4e-20
KM=0.003; K_M=4.; K_M0=.16; D_M=8.; K_P=10.; D_P=2.; GROWTH=1.3
CMIN=.6; CMAX=3.; NMAX=1.

def sa(lara): return CMIN+(CMAX-CMIN)*lara**3/(lara**3+KM**3)
def local_drop(): return 1/VTF/OMEGA*DROP_FACTOR

def split_tf(nt, vt, mu):
    if mu == 0 or nt <= 0: return nt, 0
    ndrop=np.arange(nt+1,dtype=float); ndil=nt-ndrop
    vd=ndrop*VTF; vl=vt-vd
    phi=np.maximum(ndil*VTF/vl,1e-300)
    f_dil=KBT/VTF*(phi*np.log(phi)-phi)
    area=4*np.pi*(3/(4*np.pi)*VTF*ndrop)**(2/3)
    energy=vl*f_dil-vd*mu/VTF+GAMMA*area
    idx=int(np.argmin(energy))
    return int(ndil[idx]),int(ndrop[idx])

def transcription(local, lara):
    s=sa(lara); return K_M*s*local**2/(s*local**2+1)+K_M0

def steady_state(mu):
    lara=.005
    def rhs(t,y):
        m,p,c=y; nt=max(0,int(math.ceil(p*VT0*OMEGA))); _,drop=split_tf(nt,VT0,mu)
        local=local_drop() if drop else p
        gr=GROWTH*(1-c)
        return [transcription(local,lara)-D_M*m-gr*m,K_P*m-D_P*p-gr*p,gr*c]
    sol=solve_ivp(rhs,[0,5000],[.1,.1,.01],method='BDF',rtol=1e-8,atol=1e-10)
    return np.ceil(sol.y[:2,-1]*VT0*OMEGA).astype(int)

def simulate(mu,seed,tmax=20.,step=.025):
    rng=random.Random(seed); x=steady_state(mu).astype(int)
    t=0.; cellnum=1/64; gr=GROWTH; doubling=math.log(2)/gr; since=0.; vt=VT0
    samples=[]; next_sample=0.; output=0.
    last=(0.,float(x[1]/vt/OMEGA),0.,0.,0.,0.)
    while t<tmax and x.sum()>0:
        vt=VT0*(since/doubling+1) if math.isfinite(doubling) else VT0
        ndil,ndrop=split_tf(int(x[1]),vt,mu)
        vd=ndrop*VTF; vl=vt-vd
        dilute=ndil/vl/OMEGA
        local=local_drop() if ndrop else dilute
        vm=transcription(local,.001)*OMEGA*vt
        rates=[vm,K_P*x[0],D_M*x[0],D_P*x[1]]; a0=sum(rates)
        if a0<=0: break
        tau=math.log(1/max(rng.random(),1e-15))/a0
        # Emit uniformly spaced browser data by holding the Gillespie state.
        avg=x[1]/vt/OMEGA; prod=transcription(local,.001)
        while next_sample<=min(t+tau,tmax):
            output += prod*(next_sample-last[0])
            samples.append(dict(time=round(next_sample,3),average_tf=avg,local_tf=local,
                                droplet_molecules=ndrop,cell_volume_um3=vt*1e18,
                                transcription_rate=prod,cumulative_output=output))
            last=(next_sample,avg,local,ndrop,vt,prod); next_sample+=step
        since+=tau; t+=tau
        if since>=doubling and math.isfinite(doubling):
            since=0.; x=np.ceil(x/2).astype(int); cellnum*=2; vt=VT0
            gr=max(0.,GROWTH*(1-cellnum/NMAX)); doubling=math.log(2)/gr if gr>0 else math.inf
        pick=rng.random()*a0; acc=0
        for i,r in enumerate(rates):
            acc+=r
            if pick<=acc:
                if i==0:x[0]+=1
                elif i==1:x[1]+=1
                elif i==2 and x[0]>0:x[0]-=1
                elif i==3 and x[1]>0:x[1]-=1
                break
    return samples

def main():
    data={
      "metadata":{
        "description":"Seeded trajectories generated from the publication MATLAB model",
        "source_repository":"https://github.com/TianLab-ASU/PhaseSeparation_Growth_CircuitMemory/tree/main/Stochastic%20Simulation",
        "source_script":"Fig1_LLPS_StochasticSim_CellGrowth.m",
        "induction_initial_mM":.005,"induction_after_dilution_mM":.001,
        "duration_h":20,"sample_interval_h":.025,"seeds":{"standard_sa":20251017,"drop_sa":20251018},
        "parameters":{"km":4,"km0":.16,"dm":8,"kp":10,"dp":2,"growthrate_per_h":1.3,"Vt_m3":4e-19,"V_TF_m3":1e-25,"Drop_Factor":.0005,"gamma_N_per_m":1e-6,"kBT_J":KBT,"mu_drop_J":MU0}
      },
      "standard_sa":simulate(0,20251017),"drop_sa":simulate(MU0,20251018)
    }
    out=Path(__file__).with_name('publication-trajectories.json')
    out.write_text(json.dumps(data,separators=(',',':')))
    browser=Path(__file__).with_name('publication-data.js')
    browser.write_text('window.PUBLICATION_SIM_DATA='+json.dumps(data,separators=(',',':'))+';\n')
    print(out, out.stat().st_size)
    for k in ('standard_sa','drop_sa'):
        a=data[k]; print(k,'n=',len(a),'TF start/min/end=',*[round(v,3) for v in (a[0]['average_tf'],min(x['average_tf'] for x in a),a[-1]['average_tf'])], 'drop max=',max(x['droplet_molecules'] for x in a))
if __name__=='__main__': main()
