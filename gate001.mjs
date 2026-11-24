// General Instrument — Gate 001. No dependencies. Node.js >=18.
function runGate() {
  const N=80, sigma=.12, rms=.60, peak=1.2, rounds=3, seeds=144;
  function rng(seed){let s=seed>>>0; return ()=>{s=(s+0x6D2B79F5)>>>0;let z=s;z=Math.imul(z^(z>>>15),z|1);z^=z+Math.imul(z^(z>>>7),z|61);return (((z^(z>>>14))>>>0)+.5)/4294967296}}
  function normal(r){return Math.sqrt(-2*Math.log(r()))*Math.cos(2*Math.PI*r())}
  const lp=(u,a)=>{let x=0;return u.map(v=>x=a*x+(1-a)*v)};
  const non=(u,c)=>u.map(v=>v+c*v*v*v);
  const models=[];
  for(const a of [.35,.65]) for(const b of [.82,.91]){
    models.push({kind:'linear',a,b,c:0});
    for(const c of [.25,.55]) for(const kind of ['early','middle','late-plant','late-sensor']) models.push({kind,a,b,c});
  }
  const group=m=>m.kind.startsWith('late')?'late (plant/sensor unresolved)':m.kind;
  function predict(m,u){const {a,b,c,kind}=m;
    if(kind==='early')return lp(lp(non(u,c),a),b);
    if(kind==='middle')return lp(non(lp(u,a),c),b);
    return non(lp(lp(u,a),b),c);
  }
  const power=u=>Math.sqrt(u.reduce((s,x)=>s+x*x,0)/u.length);
  const probes=[];
  function add(name,u){const scale=rms/power(u);u=u.map(x=>scale*x);if(Math.max(...u.map(Math.abs))>peak+1e-9)throw Error('peak '+name);probes.push({name,u});}
  add('constant',Array(N).fill(1));
  for(const period of [4,8,16,40,80]) add('square-'+period,Array.from({length:N},(_,t)=>t%period<period/2?1:-1));
  for(const bias of [0,.5,1]) add('bias-'+bias,Array.from({length:N},(_,t)=>bias+Math.sin(2*Math.PI*t/16)));
  add('burst',Array.from({length:N},(_,t)=>t%20<5?1:0));
  const means=probes.map(p=>models.map(m=>predict(m,p.u)));
  const d2=(a,b)=>a.reduce((s,x,i)=>s+(x-b[i])**2,0);
  const groups=[...new Set(models.map(group))];
  const gi=models.map(m=>groups.indexOf(group(m)));
  function masses(p){let q=groups.map(()=>0);p.forEach((v,i)=>q[gi[i]]+=v);return q;}
  const entropy=p=>-p.reduce((s,v)=>s+(v>0?v*Math.log(v):0),0);
  function posterior(p,y,mi){let l=p.map((v,i)=>Math.log(v)-d2(y,mi[i])/(2*sigma*sigma));let mx=Math.max(...l);let q=l.map(v=>Math.exp(v-mx)),z=q.reduce((a,b)=>a+b,0);return q.map(v=>v/z);}
  const prior=models.map(()=>1/models.length);
  // Standard Bayesian expected information gain, quadrature by seeded Monte Carlo.
  // Score group entropy, not arbitrary latent labels that are observational twins.
  const noise=Array.from({length:12},(_,j)=>{const r=rng(88000+j);return Array.from({length:N},()=>sigma*normal(r));});
  function eig(p,k){const h=entropy(masses(p));let expected=0;for(let i=0;i<p.length;i++){if(p[i]<1e-10)continue;let e=0;for(const n of noise){const y=means[k][i].map((x,t)=>x+n[t]);e+=entropy(masses(posterior(p,y,means[k])));}expected+=p[i]*e/noise.length;}return h-expected;}
  const initial=probes.map((p,k)=>eig(prior,k));
  const order=initial.map((_,k)=>k).sort((a,b)=>initial[b]-initial[a]);
  // Cheap alternative: posterior-weighted pair separation, saturating at noise scale.
  const sep=means.map(mi=>mi.map((a,i)=>mi.map((b,j)=>gi[i]===gi[j]?0:1-Math.exp(-d2(a,b)/(8*sigma*sigma)))));
  function pair(p,k){let s=0;for(let i=0;i<p.length;i++)for(let j=0;j<i;j++)s+=p[i]*p[j]*sep[k][i][j];return s;}
  function best(scores){return scores.indexOf(Math.max(...scores))}
  function decide(p){let q=masses(p);let k=best(q);return q[k]>=.95?groups[k]:'unresolved';}
  const arms=['constant','fixed-ranked','random','pair','EIG'];
  const totals=Object.fromEntries(arms.map(a=>[a,{correct:0,wrong:0,abstained:0,nll:0,lateTwinMaxGap:0,first:{},byRound:Array.from({length:rounds},()=>({correct:0,wrong:0,abstained:0}))}]));
  for(let s=0;s<seeds;s++){
    const truth=s%models.length, m=models[truth], target=group(m);
    for(const arm of arms){let p=prior.slice();const r=rng(10000+s), rr=rng(70000+s);
      for(let t=0;t<rounds;t++){
        let k=arm==='constant'?0:arm==='fixed-ranked'?order[t]:arm==='random'?Math.floor(rr()*probes.length):best(probes.map((_,k)=>arm==='pair'?pair(p,k):eig(p,k)));
        if(t===0){let f=totals[arm].first;f[probes[k].name]=(f[probes[k].name]||0)+1;}
        const y=means[k][truth].map(x=>x+sigma*normal(r));p=posterior(p,y,means[k]);
        const d=decide(p), key=d==='unresolved'?'abstained':d===target?'correct':'wrong';
        totals[arm].byRound[t][key]++;
        if(t===rounds-1){totals[arm][key]++;totals[arm].nll-=Math.log(Math.max(1e-300,masses(p)[groups.indexOf(target)]));let plant=0,sensor=0;p.forEach((v,i)=>{if(models[i].kind==='late-plant')plant+=v;if(models[i].kind==='late-sensor')sensor+=v});totals[arm].lateTwinMaxGap=Math.max(totals[arm].lateTwinMaxGap,Math.abs(plant-sensor));}
      }
    }
  }
  for(const a of arms)totals[a].nll/=seeds;
  const twinError=Math.max(...probes.map((p,k)=>Math.max(...models.map((m,i)=>m.kind==='late-plant'?d2(means[k][i],means[k][i+1]):0))));
  if(twinError!==0||arms.some(a=>totals[a].lateTwinMaxGap>1e-12))throw Error('Twin invariance failed');
  if(arms.some(a=>totals[a].correct+totals[a].wrong+totals[a].abstained!==seeds))throw Error('Bad ledger');
  return {protocol:{N,sigma,rms,peak,rounds,seeds,models:models.length,groups,decisionThreshold:.95},probes:probes.map((p,k)=>({name:p.name,priorEIG:initial[k]})),twinError,totals};
}
console.log(JSON.stringify(runGate(),null,2));
