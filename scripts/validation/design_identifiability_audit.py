#!/usr/bin/env python3
"""Design audit for GSE178742 apiary experiment: prints design matrices and tests identifiability of the exposure effect
when colony enters as (a) a fixed effect, (b) a blocking/random effect, (c) pseudo-bulk colony units."""
import numpy as np,pandas as pd,sys
R=sys.argv[1] if len(sys.argv)>1 else '.'
m=pd.read_csv(f'{R}/data/metadata/GSE178742_sample_metadata_verified.tsv',sep='\t')
lines=[]
def P(s=''): lines.append(str(s)); print(s)
for ex in ['Apiary_Foragers','Apiary_Nurses']:
    d=m[m.experiment==ex].sort_values(['treatment','colony','geo_title']).reset_index(drop=True)
    trt=(d.treatment=='Exposed').astype(int).values; cols=sorted(d.colony.unique())
    P(f'===== {ex}: {len(d)} bees, colonies {cols}; exposed colonies {sorted(d[d.treatment=="Exposed"].colony.unique())}')
    # (a) fixed colony effects
    C=np.column_stack([(d.colony==c).astype(int).values for c in cols[1:]])  # treatment-coded, reference = first colony
    Xa=np.column_stack([np.ones(len(d)),trt,C]); ra=np.linalg.matrix_rank(Xa)
    P(f'(a) ~ exposure + colony (fixed): columns={Xa.shape[1]} rank={ra} -> {"NOT identifiable (exposure is a linear combination of colony dummies)" if ra<Xa.shape[1] else "identifiable"}')
    # show the exact linear dependence: exposure = sum of exposed-colony dummies (+ intercept if reference colony exposed)
    ref=cols[0]; refexp=int(d[d.colony==ref].treatment.iloc[0]=='Exposed')
    comb=refexp*np.ones(len(d))+sum(((d.colony==c)&(d.treatment=='Exposed')).astype(int).values - refexp*(d.colony==c).astype(int).values for c in cols[1:])
    P(f'    check: exposure column == {"intercept - " if refexp else ""}linear combination of colony dummies: {np.array_equal(comb,trt)}')
    # (b) model actually used in de_from_kallisto.R: ~ exposure, colony as blocking factor (duplicateCorrelation; random-intercept-like)
    Xb=np.column_stack([np.ones(len(d)),trt]); P(f'(b) ~ exposure; colony = blocking factor (NOT in design): columns=2 rank={np.linalg.matrix_rank(Xb)} -> identifiable')
    P('    Design matrix (b) [bee, colony, (Intercept), trtExposed]:')
    for i,r in d.iterrows(): P(f'    {r.geo_title:9s} {r.colony}  1  {trt[i]}')
    # (c) pseudo-bulk
    pb=d.groupby('colony').treatment.first(); Xc=np.column_stack([np.ones(len(pb)),(pb=='Exposed').astype(int)])
    P(f'(c) colony pseudo-bulk (PRIMARY): n={len(pb)} colonies, design rank={np.linalg.matrix_rank(Xc)}, residual df={len(pb)-2}')
    for c,t in pb.items(): P(f'    {c}  1  {int(t=="Exposed")}')
P('\nInterpretation: exposure varies only BETWEEN colonies. With colony as a fixed effect the exposure effect is not estimable.')
P('Model (b) estimates it only by treating colony as a random/blocking effect; the exposure contrast is still informed by 6 colonies,')
P('the consensus intra-colony correlation from duplicateCorrelation is estimated from 6 blocks (imprecise), and limma does not adjust')
P('the residual df to the number of colonies, so (b) can be anti-conservative. (b) is therefore a SENSITIVITY analysis only;')
P('(c) colony pseudo-bulk (3 v 3, residual df = 4) is the PRIMARY analysis. Field experiment: 1 pooled library per colony; exposure fully')
P('confounded with site (Wellington 2 v Toronto 2), so no site-independent exposure estimate exists.')
open(f'{R}/results/validation/design_identifiability_audit.txt','w').write('\n'.join(lines))
