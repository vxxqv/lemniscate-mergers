"""Orthogonal checks for the mathematical implementation, not proof by sampling."""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from flint import arb,acb,ctx
import mpmath as mp
from analytic import family,evaluate,mul,stationary_certificate,find_certificate,mp_coeffs,derivative

ctx.prec=384;mp.mp.dps=90

# Exact quadratic at a nonstationary, nonreal translated center.
# F(z)=3+(z-2i)^2, c=2i, F(c)=3.
a=[acb(-1),acb(0,-4),acb(1)]
cert=stationary_certificate(a,acb(1),arb(3))
assert cert is not None and cert['error']==0 and cert['center']==3
assert stationary_certificate([acb(1),acb(0),acb(0),acb(1)],acb(0),arb('1/10')) is None

factor_checks=0
for name in ['aligned4','opposed2','ring4','broken_ring4']:
 for es in ['1/2','1/50']:
  q,g,a=family(name,arb(es));qc,gc,ac=[mp_coeffs(t,mp) for t in [q,g,a]]
  z=mp.mpc('0.12','0.25')
  assert abs(mp.polyval(ac[::-1],z)-mp.polyval(qc[::-1],z)*mp.polyval(gc[::-1],z))<mp.mpf('1e-85')
  factor_checks+=1

# Compare critical-height ordering under reflection with geometric root labels.
cp=[(-1-mp.sqrt(7))/3,(-1+mp.sqrt(7))/3]
pv=[abs((z+2)*z*(z-1)) for z in cp]
assert pv[1]<pv[0]  # first merger is middle/right, reversed after reflection

# Verify every included fine-grid trial, not only an aggregated headline.
import csv
root=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader(open(root/'computations/main-003/raster_trials.csv')))
fine=[x for x in rows if float(x['h_over_width'])<=.5]
assert len(rows)==5376 and len(fine)==2304
assert all(x['status']=='correct' for x in fine)

# Independent exact geometry identities at rational samples.
for numerator in [1,2,4]:
 eta=mp.mpf(numerator)/100
 for x in [mp.mpf(k)/10 for k in range(-10,11)]:
  value=abs((x+1j*mp.sqrt(eta))**2-1)**2-(1+eta)**2
  expected=x*x*(x*x+2*eta-2)
  assert abs(value-expected)<mp.mpf('1e-85') and value<=mp.mpf('1e-85')

summary={'exact_quadratic_stationary_value':'passed','degenerate_center_rejection':'passed',
         'independent_factorized_checks':factor_checks,'reflection_label_example':'passed',
         'fine_grid_trials_checked':len(fine),'tube_identity_samples':63,
         'non_claim':'These finite checks do not replace the universal proofs.'}
out=root/'records/orthogonal_checks.json';out.write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
