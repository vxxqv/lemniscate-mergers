"""Reproducible finite experiments; run from the project root.

python3 src/experiments.py --out computations/main-001
"""
import argparse, csv, json, time, platform
from pathlib import Path
import numpy as np
from scipy.ndimage import label
import mpmath as mp
from flint import arb, acb, ctx
from analytic import family, derivative, evaluate, find_certificate, mp_coeffs

FAMILIES = ['aligned4','opposed2','ring4','broken_ring4']
EPS = ['1/2','1/3','1/5','1/8','1/12','1/20','1/32','1/50']


def real(x):
    return float(x.mid())


def mp_eval(coeffs, w):
    return mp.polyval(list(reversed(coeffs)), w)


def write_csv(path, rows):
    with open(path, 'w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def asymptotics(out):
    ctx.prec=384; mp.mp.dps=95
    rows=[]; raw=[]
    for name in FAMILIES:
        for es in EPS:
            eps=arb(es); q,g,a=family(name,eps)
            ma=mp_coeffs(a,mp); md=mp_coeffs(derivative(a),mp)
            for gamma in [-1,1]:
                # Independent mp Newton computation; reference is NOT certificate.
                c=mp.findroot(lambda w:mp_eval(md,w),mp.mpc(gamma),solver='newton',tol=mp.mpf('1e-88'))
                actual=mp_eval(ma,c); f0=mp_eval(ma,mp.mpc(gamma))
                f1=mp_eval(md,mp.mpc(gamma))
                f2=mp_eval(mp_coeffs(derivative(a,2),mp),mp.mpc(gamma))
                vn=f0-f1*f1/(2*f2)
                cert=find_certificate(a,acb(gamma))
                # A 1e-80 reference box encloses the computed approximation only;
                # theorem localization supplies existence/uniqueness independently.
                check=False; err=None; radius=None; margin=None
                if cert:
                    cv=mp.mpc(cert['center'].real.mid().str(110,radius=False),
                              cert['center'].imag.mid().str(110,radius=False))
                    err=real(cert['error']); radius=real(cert['radius']); margin=real(cert['margin'])
                    check=bool(abs(actual-cv) <= mp.mpf(cert['error'].lower().str(110,radius=False)))
                    assert check
                    assert abs(c-gamma) < mp.mpf(cert['radius'].lower().str(110,radius=False))
                    raw.append(dict(family=name,epsilon=es,gamma=gamma,
                                    radius=cert['radius_exact'],rouche_margin=str(cert['margin']),
                                    value_center=str(cert['center']),value_radius=str(cert['error']),
                                    value_modulus_lower=str((abs(cert['center'])-cert['error']).lower()),
                                    value_modulus_upper=str((abs(cert['center'])+cert['error']).upper())))
                rows.append(dict(family=name,epsilon_exact=es,epsilon=real(eps),gamma=gamma,
                                 displacement=float(abs(c-gamma)),frozen_value_error=float(abs(actual-f0)),
                                 newton_value_error=float(abs(actual-vn)),certificate=bool(cert),
                                 certificate_radius=radius,certificate_error=err,rouche_margin=margin,
                                 reference_inside=check,reference_critical_height=float(abs(actual)),
                                 derivative_residual=float(abs(mp_eval(md,c)))))
    write_csv(out/'critical_values.csv',rows)
    (out/'ball_certificates.json').write_text(json.dumps(raw,indent=2))
    return rows


def partition_bounds(out):
    # All four unbroken configurations have four remote roots on |b|=sqrt(10).
    # For broken ring, the altered root lies farther out, so the shared
    # distance-only bracket below is still a safe conservative bracket.
    ctx.prec=384
    R=arb('5/2'); B=abs(acb(3,1))
    tau1=abs(acb(arb('-8/5'),arb('1/5')))
    tau2=abs(acb(arb('12/5'),arb('1/5')))
    # |w^3-3w+(2+i)/5| >= R^3-3R-|(2+i)/5| on |w|=R.
    barrier=R**3-3*R-abs(acb(arb('2/5'),arb('1/5')))
    rows=[]
    for name in FAMILIES:
        for es in EPS:
            e=arb(es); x=e*R/B
            dlo=(1-x)**4; dhi=(1+x)**4
            if name=='aligned4':
                delta=-4*(1-x).log()
            elif name=='opposed2':
                delta=-2*(1-x*x).log()
            elif name=='ring4':
                delta=-(1-x**4).log()
            else:
                # ring logarithm plus the exact ratio replacing one remote zero.
                # Difference series: |(1+d)^(-k)-1| <= k*d (d>0).
                delta=-(1-x**4).log()+arb('1/1000')*x/(1-x)
            low=(-delta).exp(); high=delta.exp()
            dv=bool(dhi*tau1 < dlo*tau2 and dhi*tau2 < dlo*barrier)
            mv=bool(high*tau1 < low*tau2 and high*tau2 < low*barrier)
            rows.append(dict(family=name,epsilon_exact=es,epsilon=real(e),
                             distance_certifies_hierarchy=dv,moment_certifies_hierarchy=mv,
                             distance_gap=real(dlo*tau2-dhi*tau1),
                             moment_gap=real(low*tau2-high*tau1),log_bound=real(delta),
                             boundary_lower=real(barrier)))
    write_csv(out/'hierarchy_bounds.csv',rows)
    return rows


def raster(out):
    # Test only the EXACT quadratic benchmark. Thresholds are regular. The
    # outcome is root-to-root connectivity in a sampled vertex graph, not
    # a theorem about a contouring library or exact topology certification.
    ratios=[0.125,0.25,0.5,1,2,4,8]
    etas=[0.0025,0.01,0.04]
    angles=[0,np.pi/16,np.pi/8,np.pi/4]
    shifts=[0,0.25,0.5,0.75]
    rows=[]
    for eta in etas:
      width=2*np.sqrt(eta)
      for ratio in ratios:
       h=ratio*width
       for theta in angles:
        for sx in shifts:
         for sy in shifts:
          # Cover the same fixed window for every shift, including coarse
          # grids. Otherwise translating the first vertex could crop a root.
          xs=(np.arange(np.floor(-1.8/h-sx),np.ceil(1.8/h-sx)+1)+sx)*h
          ys=(np.arange(np.floor(-1.8/h-sy),np.ceil(1.8/h-sy)+1)+sy)*h
          z=(xs[None,:]+1j*ys[:,None])*np.exp(-1j*theta)
          heights=np.abs(z*z-1)
          roots=np.array([-1,1])*np.exp(1j*theta)
          indices=[(int(np.argmin(abs(ys-r.imag))),int(np.argmin(abs(xs-r.real)))) for r in roots]
          for sign in [-1,1]:
           mask=heights <= 1+sign*eta
           for conn,structure in [(4,np.array([[0,1,0],[1,1,1],[0,1,0]])),(8,np.ones((3,3),int))]:
            lab,_=label(mask,structure=structure)
            ids=[int(lab[y,x]) for y,x in indices]
            seen=all(i>0 for i in ids)
            merged=seen and ids[0]==ids[1]
            status='roots_unresolved' if not seen else ('correct' if merged==(sign>0) else 'wrong')
            rows.append(dict(eta=eta,h_over_width=ratio,h=h,theta_radians=theta,
                             shift_x=sx,shift_y=sy,level_side=sign,adjacency=conn,
                             vertices=heights.size,status=status))
    write_csv(out/'raster_trials.csv',rows)
    return rows


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);args=p.parse_args()
    out=Path(args.out);out.mkdir(exist_ok=True,parents=True);start=time.perf_counter()
    cr=asymptotics(out); hb=partition_bounds(out); rr=raster(out)
    summary=dict(critical_cases=len(cr),certified_critical_cases=sum(x['certificate'] for x in cr),
                 reference_inside_all_valid_certificates=all(x['reference_inside'] for x in cr if x['certificate']),
                 hierarchy_cases=len(hb),distance_hierarchies=sum(x['distance_certifies_hierarchy'] for x in hb),
                 moment_hierarchies=sum(x['moment_certifies_hierarchy'] for x in hb),raster_trials=len(rr),
                 raster_outcomes={s:sum(x['status']==s for x in rr) for s in ['correct','wrong','roots_unresolved']},
                 elapsed_seconds=time.perf_counter()-start,python=platform.python_version())
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
