"""Plot real QM ESP versus fitted RESP ESP after server calculations."""
import argparse
import numpy as np
from common import ROOT,read_json,read_esp,sha256


def plot(group):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    folder=ROOT/'runs/resp'/group
    report=read_json(folder/'review.json')
    template=report['template_charges']
    jobs=read_json(ROOT/'structures/manifest.json')['jobs']
    jobs=[j for j in jobs if (j['molecule'] in ['PMMA3','PMMA5'] if group=='PMMA' else j['molecule']==group)]
    output=folder/'figures';output.mkdir(exist_ok=True)
    for job in jobs:
        rows=read_json(ROOT/'structures'/(job['molecule']+'_atom_map.json'))
        q=np.array([template[r['equivalence_key']] for r in rows])
        actual=ROOT/'runs/qm'/job['id']
        metadata=read_json(actual/'result.json')
        if sha256(actual/'qm.esp')!=metadata['esp_sha256']:raise ValueError('ESP changed after extraction')
        atoms,points,observed=read_esp(actual/'qm.esp');predicted=np.zeros(len(points))
        for start in range(0,len(points),2048):
            d=np.linalg.norm(points[start:start+2048,None]-atoms[None],axis=2)
            predicted[start:start+2048]=(1/d)@q
        residual=predicted-observed
        ids=np.linspace(0,len(points)-1,min(3000,len(points)),dtype=int)
        fig,axes=plt.subplots(1,2,figsize=(9,3.8),constrained_layout=True)
        axes[0].scatter(observed[ids],predicted[ids],s=5,alpha=.4)
        lo=min(observed.min(),predicted.min());hi=max(observed.max(),predicted.max())
        axes[0].plot([lo,hi],[lo,hi],color='black',linewidth=1)
        axes[0].set(xlabel='QM ESP [atomic units]',ylabel='RESP point-charge ESP [atomic units]')
        axes[1].hist(residual,bins=60)
        axes[1].set(xlabel='RESP ESP - QM ESP [atomic units]',ylabel='Grid points')
        fig.suptitle(f"{job['id']} | {job['purpose']} | {len(points)} points | RMSE {np.sqrt(np.mean(residual**2)):.4g} a.u.")
        fig.savefig(output/(job['id']+'_esp.png'),dpi=180);plt.close(fig)
    print('Real ESP/RESP figures:',output)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('group',choices=['TPS','OTF','PMMA'])
    plot(p.parse_args().group)
