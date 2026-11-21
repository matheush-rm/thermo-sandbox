'use strict';
// =====================================================================
//  ELETROMAGNETISMO
//  Elétrico (quase-estático, Euler implícito, unidades SI, célula = h³):
//    Σ_j h(ε_ij + σ_ij·dt)(V_i − V_j) = q_i        (Gauss + continuidade)
//    I_ij = h·σ_ij·(V_i − V_j)   q ← q − dt·ΣI    Joule: I·ΔV·dt
//  Ruptura dielétrica: |E| > E_rup (∝ densidade do gás) ioniza a célula → plasma condutor.
//  O canal cresce pela ponta (campo concentrado) = modelo de quebra dielétrica (DBM).
//  Magnético: Biot–Savart (correntes no plano → Bz) + dipolos (ímãs, ferro induzido → Bx,By,Bz),
//  convolução por FFT em grade 4× mais grossa. Força em ferromagnéticos soltos: f = |M|·∇|B|.
// =====================================================================
const EPS0=8.854e-12,KM=1e-7,MU0=4e-7*Math.PI,VOL=DX*DX*DX,SIGPL=1e4,DT_E=0.02,DT_F=2e-4,TMAX=30000,ITER=10,EST=2; // EST: solução elétrica a cada 2 passos (dt dobrado; implícito → estável)
const Vf=new Float64Array(N),aR=new Float64Array(N),aD=new Float64Array(N),eR=new Float64Array(N),eD=new Float64Array(N),
  sR=new Float64Array(N),sD=new Float64Array(N),dg=new Float64Array(N),sgE=new Float64Array(N),tmpG=new Float64Array(N);
const IxC=new Float32Array(N),IyC=new Float32Array(N),ionL=new Float32Array(N),fixV=new Uint8Array(N),dom=new Uint8Array(N);
for(let y=1;y<H-1;y++)for(let x=1;x<W-1;x++)dom[y*W+x]=1;
let hasSrc=false,ePl=0,eJoule=0,flash=0,nStrikes=0,strikeCD=0,bMax=0;

// ---------- coeficientes ----------
function elecCoef(dt){
  let live=storm;const vp=Vsrc/2;hasSrc=false;
  for(let i=0;i<N;i++){const sr=srcA[mat[i]];if(sr&&sr!==3&&dom[i]){hasSrc=true;break;}}
  // σ do plasma limitado a κ·ε0/dt (τ=ε/σ ≪ dt: o canal continua praticamente equipotencial, a transferência
  // de carga/energia é a mesma, e o sistema fica bem condicionado). ⚠️ assumido: com fonte ideal o arco fica
  // limitado a ~MΩ (faz o papel da impedância interna da fonte)
  const spl=Math.min(SIGPL,1e4*EPS0/dt);
  for(let i=0;i<N;i++){
    const m=mat[i];let s=sigA[m];
    if(alA[m]){const f=1+alA[m]*(T[i]-293.15);s/=f>.05?f:.05;}
    let io=ion[i];
    if(kind[m]===GAS&&T[i]>4000){let t=(T[i]-4000)/8000;t=t>1?1:t*t;if(t>io)io=t;} // ionização térmica
    ionL[i]=io; if(io>1e-3){s+=spl*io;live=true;}
    sgE[i]=s;
    const sr=srcA[m];
    if(sr&&dom[i]){fixV[i]=1;Vf[i]=sr===1?vp:sr===2?-vp:0;if(sr!==3)live=true;}else fixV[i]=0;
    if(!live&&(chg[i]>1e-13||chg[i]<-1e-13)) live=true;
  }
  if(!live) return false;
  const h=DX;
  for(let y=1;y<H-1;y++){const row=y*W;
    for(let x=1;x<W-1;x++){const i=row+x,ei=epsA[mat[i]]*EPS0,si=sgE[i];
      if(x<W-2){const j=i+1,ej=epsA[mat[j]]*EPS0,sj=sgE[j],e=2*h*ei*ej/(ei+ej),s=si+sj>0?2*h*si*sj/(si+sj):0;
        eR[i]=e;sR[i]=s;aR[i]=e+s*dt;}
      if(y<H-2){const j=i+W,ej=epsA[mat[j]]*EPS0,sj=sgE[j],e=2*h*ei*ej/(ei+ej),s=si+sj>0?2*h*si*sj/(si+sj):0;
        eD[i]=e;sD[i]=s;aD[i]=e+s*dt;}
      else{eD[i]=2*h*ei;sD[i]=2*h*si;aD[i]=eD[i]+sD[i]*dt;} // piso = terra (V=0)
    }}
  for(let y=1;y<H-1;y++){const row=y*W;for(let x=1;x<W-1;x++){const i=row+x;dg[i]=aR[i-1]+aR[i]+aD[i-W]+aD[i];}}
  return true;
}

// ---------- agregados condutores (correção grossa: condutor isolado acha seu potencial pela carga total) ----------
const LTH=[1e-8],NL=LTH.length; // condutor: σ·dt ≫ ε
const labA=[],ordA=[],cS=[],cE=[],cAn=[],nC=new Int32Array(NL);
for(let L=0;L<NL;L++){labA.push(new Int32Array(N).fill(-1));ordA.push(new Int32Array(N));cS.push(new Int32Array(N));cE.push(new Int32Array(N));cAn.push(new Uint8Array(N));}
let lGround=false,flashN=0;const cHi=new Uint8Array(N);
function label(){
  lGround=false;
  for(let L=0;L<NL;L++){
    const lab=labA[L],ord=ordA[L],S=cS[L],E=cE[L],An=cAn[L],th=LTH[L];
    lab.fill(-1);let nc=0,pos=0;
    for(let i=W;i<N-W;i++){
      if(!dom[i]||lab[i]>=0||sgE[i]<th) continue;
      const st=pos;ord[pos++]=i;lab[i]=nc;let an=fixV[i],hd=st,hi=0,tg=0,tp=0,tn=0; // tg/tp/tn: toca terra / polo + / polo −
      while(hd<pos){const a=ord[hd++];if(a>=(H-2)*W)an=1; // toca o piso = aterrado
        if(ion[a]>.6)hi=1;
        if(L===0){if(a>=(H-2)*W)tg=1;if(fixV[a]){const r=srcA[mat[a]];if(r===3)tg=1;else if(r===1)tp=1;else tn=1;}}
        let j=a-1;if(dom[j]&&lab[j]<0&&sgE[j]>=th){lab[j]=nc;ord[pos++]=j;an|=fixV[j];}
        j=a+1;if(dom[j]&&lab[j]<0&&sgE[j]>=th){lab[j]=nc;ord[pos++]=j;an|=fixV[j];}
        j=a-W;if(dom[j]&&lab[j]<0&&sgE[j]>=th){lab[j]=nc;ord[pos++]=j;an|=fixV[j];}
        j=a+W;if(dom[j]&&lab[j]<0&&sgE[j]>=th){lab[j]=nc;ord[pos++]=j;an|=fixV[j];}
      }
      S[nc]=st;E[nc]=pos;An[nc]=an;if(L===0)cHi[nc]=hi;if(hi&&((tg&&!tp&&!tn)||tg+tp+tn>=2))lGround=true;nc++; // canal completo: nuvem→terra ou entre potenciais distintos
    }
    nC[L]=nc;
  }
}
function aggCorrect(){
  for(let L=0;L<NL;L++){
    const lab=labA[L],ord=ordA[L],S=cS[L],E=cE[L],An=cAn[L],nc=nC[L];
    for(let c=0;c<nc;c++){
      const s=S[c],e=E[c];if(An[c]||e-s<2) continue;
      let num=0,den=0;
      for(let k=s;k<e;k++){const i=ord[k],v=Vf[i];num+=chg[i];
        let a=aR[i-1];if(a>0&&lab[i-1]!==c){num-=a*(v-Vf[i-1]);den+=a;}
        a=aR[i];if(a>0&&lab[i+1]!==c){num-=a*(v-Vf[i+1]);den+=a;}
        a=aD[i-W];if(a>0&&lab[i-W]!==c){num-=a*(v-Vf[i-W]);den+=a;}
        a=aD[i];if(a>0&&lab[i+W]!==c){num-=a*(v-Vf[i+W]);den+=a;}
      }
      if(den>0){const d=num/den;for(let k=s;k<e;k++)Vf[ord[k]]+=d;}
    }
  }
}
// ---------- SOR vermelho-preto ----------
function sor(n){
  const w=1.85;
  for(let it=0;it<n;it++){
    for(let c=0;c<2;c++){
      for(let y=1;y<H-1;y++){const row=y*W;
        for(let x=1+((y+c)&1);x<W-1;x+=2){const i=row+x;if(fixV[i])continue;
          const v=(chg[i]+aR[i-1]*Vf[i-1]+aR[i]*Vf[i+1]+aD[i-W]*Vf[i-W]+aD[i]*Vf[i+W])/dg[i];
          Vf[i]+=w*(v-Vf[i]);}}}
    if(it===(n>>1)||it===n-1) aggCorrect();
  }
}
// ---------- PCG deflacionado nas células condutoras (Jacobi + deflação por corpo condutor) ----------
// resolve com precisão a rede condutora (fios longos, contraste cobre/nicromo/eletrólito) mantendo o ar fixo;
// a deflação pelos indicadores de cada corpo resolve exatamente o potencial flutuante (carga total / capacitância).
const inS=new Uint8Array(N),Sx=new Int32Array(N),rr=new Float64Array(N),zz=new Float64Array(N),pq=new Float64Array(N),ww=new Float64Array(N),ext=new Float64Array(N);
let cEc=new Float64Array(1024),cMu=new Float64Array(1024),cAct=new Int32Array(1024);
function dpcg(maxIt,tol,onlyIon){ // onlyIon: só corpos com plasma (sub-passos rápidos)
  const lab=labA[0],ord=ordA[0],S=cS[0],E=cE[0],nc=nC[0];
  if(!nc) return 0;
  if(cEc.length<nc){cEc=new Float64Array(nc*2);cMu=new Float64Array(nc*2);}
  if(cAct.length<nc)cAct=new Int32Array(nc*2);
  let n=0,na=0; // cAct: corpos ativos (loops de deflação só sobre eles, não sobre todas as gotas da nuvem)
  for(let c=0;c<nc;c++){if(onlyIon&&!cHi[c])continue;const n0=n;for(let k=S[c];k<E[c];k++){const i=ord[k];if(!fixV[i]){inS[i]=1;Sx[n++]=i;}}if(n>n0)cAct[na++]=c;}
  if(!n) return 0;
  if(maxIt<0)maxIt=Math.min(1500,Math.max(-maxIt,(1.5e5/n)|0)); // orçamento de trabalho: conjuntos pequenos convergem por completo
  for(let a=0;a<na;a++){const c=cAct[a];cEc[c]=0;cMu[c]=0;}
  let bn=0;
  for(let t=0;t<n;t++){const i=Sx[t],c=lab[i];
    const e=dg[i]-(inS[i-1]?aR[i-1]:0)-(inS[i+1]?aR[i]:0)-(inS[i-W]?aD[i-W]:0)-(inS[i+W]?aD[i]:0);
    ext[i]=e;cEc[c]+=e;
    const v=Vf[i],q=chg[i];
    rr[i]=q+aR[i-1]*Vf[i-1]+aR[i]*Vf[i+1]+aD[i-W]*Vf[i-W]+aD[i]*Vf[i+W]-dg[i]*v;
    const b=q+(inS[i-1]?0:aR[i-1]*Vf[i-1])+(inS[i+1]?0:aR[i]*Vf[i+1])+(inS[i-W]?0:aD[i-W]*Vf[i-W])+(inS[i+W]?0:aD[i]*Vf[i+W]);
    bn+=b*b;cMu[c]+=rr[i];}
  // correção grossa (deflação): zera a carga residual líquida de cada corpo
  for(let a=0;a<na;a++){const c=cAct[a];cMu[c]=cEc[c]>0?cMu[c]/cEc[c]:0;}
  for(let t=0;t<n;t++){const i=Sx[t],m=cMu[lab[i]];Vf[i]+=m;rr[i]-=ext[i]*m;}
  const lim=tol*tol*bn;
  let rz=0,rn=0;for(let a=0;a<na;a++)cMu[cAct[a]]=0;
  for(let t=0;t<n;t++){const i=Sx[t],z=rr[i]/dg[i];zz[i]=z;rz+=rr[i]*z;rn+=rr[i]*rr[i];cMu[lab[i]]+=ext[i]*z;}
  for(let a=0;a<na;a++){const c=cAct[a];cMu[c]=cEc[c]>0?cMu[c]/cEc[c]:0;}
  for(let t=0;t<n;t++){const i=Sx[t];pq[i]=zz[i]-cMu[lab[i]];}
  let it=0;
  while(it<maxIt&&rn>lim&&rz>0){
    let pw=0;
    for(let t=0;t<n;t++){const i=Sx[t];
      const w=dg[i]*pq[i]-aR[i-1]*pq[i-1]-aR[i]*pq[i+1]-aD[i-W]*pq[i-W]-aD[i]*pq[i+W];ww[i]=w;pw+=pq[i]*w;}
    if(!(pw>0))break;
    const al=rz/pw;let rz2=0;rn=0;for(let a=0;a<na;a++)cMu[cAct[a]]=0;
    for(let t=0;t<n;t++){const i=Sx[t];Vf[i]+=al*pq[i];const r=rr[i]-al*ww[i];rr[i]=r;const z=r/dg[i];zz[i]=z;rz2+=r*z;rn+=r*r;cMu[lab[i]]+=ext[i]*z;}
    for(let a=0;a<na;a++){const c=cAct[a];cMu[c]=cEc[c]>0?cMu[c]/cEc[c]:0;}
    const be=rz2/rz;rz=rz2;
    for(let t=0;t<n;t++){const i=Sx[t];pq[i]=be*pq[i]+zz[i]-cMu[lab[i]];}
    it++;
  }
  for(let t=0;t<n;t++){const i=Sx[t];inS[i]=0;pq[i]=0;}
  return it;
}
// ---------- correntes, transporte de carga, Joule ----------
const jb=new Float64Array(N);
function joule(i,e){jb[i]+=e;}
// aplica o calor de Joule; sem fontes, limita a energia total à energia eletrostática armazenada ½Σq·V
// (impede que resíduo numérico do solver em canais de plasma crie energia do nada)
// energia do campo ½Σ hε·ΔV² (sempre ≥0, insensível ao resíduo nas ligações de plasma)
// ⚠️ assumido: a placa de 10 cm subestima a carga/energia de uma nuvem real (~5 C, 10⁸ V, ~10⁵ J por metro
// de canal); PLX escala o calor depositado no canal de plasma sem fonte para reproduzir ~10⁴ K
const PLX=3e5,EPLM=2e5; // EPLM: teto de energia do canal por sub-passo (J)
let Uprev=0;
function fieldU(){let u=0;
  for(let y=1;y<H-1;y++){const row=y*W;for(let x=1;x<W-1;x++){const i=row+x,v=Vf[i];let d;
    if(x<W-2){d=v-Vf[i+1];u+=eR[i]*d*d;}d=v-(y<H-2?Vf[i+W]:0);u+=eD[i]*d*d;}}
  return .5*u;}
function fluxJ(dt){
  const U=fieldU();
  const p0=ePl,e0=eJoule,E=flux(dt);let sc=1;
  if(!(E<1e15)){jb.fill(0);Uprev=U;return 0;} // overflow numérico: descarta o passo (o guarda de V reinicia o solver)
  if(!hasSrc){const av=Math.max(0,Uprev-U)+1e-9;if(E>av)sc=av/E;}
  Uprev=U;
  ePl=p0+(ePl-p0)*sc*(hasSrc?1:PLX);eJoule=e0+(eJoule-e0)*sc;
  // canal de raio (sem fonte): energia da coluna 3D (assumido), limitada por sub-passo
  let ep=0;if(!hasSrc)for(let i=0;i<N;i++)if(jb[i]&&ion[i]>.5&&kind[mat[i]]===GAS)ep+=jb[i]*sc*PLX;
  const px=ep>EPLM?PLX*EPLM/ep:PLX;
  for(let i=0;i<N;i++){const e=jb[i];if(e){jb[i]=0;const m=mat[i];
    const x=!hasSrc&&ion[i]>.5&&kind[m]===GAS?px:1;
    const t=T[i]+e*sc*x/(capOf(i,m)*1e3);T[i]=t>TMAX?TMAX:t;}}
  return sc;
}
function flux(dt){
  IxC.fill(0);IyC.fill(0);let E=0;
  for(let y=1;y<H-1;y++){const row=y*W;
    for(let x=1;x<W-1;x++){const i=row+x,fi=fixV[i];
      if(x<W-2){const s=sR[i];if(s>1e-16){const j=i+1,dv=Vf[i]-Vf[j],I=s*dv;
        if(I>1e-12||I<-1e-12){const dq=I*dt,fj=fixV[j];if(!fi)chg[i]-=dq;if(!fj)chg[j]+=dq;IxC[i]+=I*.5;IxC[j]+=I*.5;
          const e=I*dv*dt;if(e>1e-6){if(!fi)joule(i,e*.5);if(!fj)joule(j,e*.5);E+=e;if(ionL[i]>.1||ionL[j]>.1)ePl+=e;}}}}
      const s=sD[i];if(s>1e-16){const j=i+W,dv=Vf[i]-Vf[j],I=s*dv;
        if(I>1e-12||I<-1e-12){const dq=I*dt,g=!dom[j],fj=g||fixV[j];if(!fi)chg[i]-=dq;if(!fj)chg[j]+=dq;IyC[i]+=I*.5;if(!g)IyC[j]+=I*.5;
          const e=I*dv*dt;if(e>1e-6){if(!fi)joule(i,e*.5);if(!fj)joule(j,e*.5);E+=e;if(ionL[i]>.1||(!g&&ionL[j]>.1))ePl+=e;}}}
    }}
  eJoule+=E;return E;
}
// erro residual do solver dentro de condutores → redistribui uniformemente (relaxação instantânea de carga interna)
function smoothChg(){
  const L=NL-1,ord=ordA[L],S=cS[L],E=cE[L],An=cAn[L],nc=nC[L];
  for(let c=0;c<nc;c++){const s=S[c],e=E[c];if(e-s<2)continue;let sum=0,n=0;
    for(let k=s;k<e;k++){const i=ord[k];if(fixV[i])continue;const v=Vf[i];
      const g=eR[i-1]*(v-Vf[i-1])+eR[i]*(v-Vf[i+1])+eD[i-W]*(v-Vf[i-W])+eD[i]*(v-Vf[i+W]);
      tmpG[i]=g;sum+=chg[i]-g;n++;}
    if(!n)continue;const mr=An[c]?0:sum/n; // aterrado/ligado a fonte: carga = Gauss; isolado: conserva a carga total
    for(let k=s;k<e;k++){const i=ord[k];if(!fixV[i])chg[i]=tmpG[i]+mr;}
  }
}
// ---------- ruptura dielétrica ----------
function eMag(i){const v=Vf[i];
  const ex=(aR[i-1]>0&&aR[i]>0)?(Vf[i-1]-Vf[i+1])/(2*DX):aR[i]>0?(v-Vf[i+1])/DX:aR[i-1]>0?(Vf[i-1]-v)/DX:0;
  const ey=(aD[i-W]>0)?(Vf[i-W]-Vf[i+W])/(2*DX):(v-Vf[i+W])/DX;
  return Math.hypot(ex,ey);}
// ruptura estocástica (modelo de quebra dielétrica, DBM): avançam poucas células por sub-passo,
// sorteadas com peso ∝ (E/E_rup)^η → canal fino e ramificado que procura o caminho de menor rigidez
const VAPT=5000,QMAX=2e-6,PROP=.35,FLMAX=1500,BK=30,BETA=3,cdI=new Int32Array(N),cdW=new Float64Array(N);
// razão E/E_rup da célula i (0 se não rompe ou não é dielétrico); th: fração do campo de ruptura exigida
function brk(i,th){
  if(fixV[i]||ion[i]>.5)return 0;
  const m=mat[i],e0=ebA[m];if(!e0||sgE[i]>1)return 0;
  let eb=e0;if(kind[m]===GAS){let r=gd[i]/1.2;if(r<.05)r=.05;eb*=r;} // Paschen (pd alto): E_rup ∝ densidade
  const v=Vf[i];let d,dm=0;
  if(aR[i-1]>0){d=v-Vf[i-1];if(d<0)d=-d;if(d>dm)dm=d;}
  if(aR[i]>0){d=v-Vf[i+1];if(d<0)d=-d;if(d>dm)dm=d;}
  if(aD[i-W]>0){d=v-Vf[i-W];if(d<0)d=-d;if(d>dm)dm=d;}
  if(aD[i]>0){d=v-Vf[i+W];if(d<0)d=-d;if(d>dm)dm=d;}
  const r=dm/(eb*th*DX);return r>1?r:0;
}
const bkM=new Uint32Array(N);let bkS=0;
function breakdown(){
  if(flashN>=FLMAX||lGround)return 0; // canal completo (tocou a terra) ou líder esgotado
  // 1) ponta do líder: só vizinhos de canal ionizado (a ponta avança por streamers, ~1 MV/m no ar → PROP·E_rup)
  bkS++;let na=0,ta=0;
  for(let y=1;y<H-1;y++){const row=y*W;
    for(let x=1;x<W-1;x++){const c=row+x;if(!(ion[c]>.6))continue;
      for(let q=0;q<4;q++){const i=q===0?c-1:q===1?c+1:q===2?c-W:c+W;
        if(bkM[i]===bkS||i<W||i>=N-W)continue;bkM[i]=bkS;const xi=i%W;if(xi<1||xi>W-2)continue;
        const r=brk(i,PROP);if(r){const w=Math.pow(r,BETA);cdI[na]=i;cdW[na]=w;ta+=w;na++;}}}}
  // há líder ativo → cresce só a partir dele (canal único ramificado), poucas células sorteadas ∝ (E/E_rup)^η
  if(na){
    if(na<=BK){for(let k=0;k<na;k++)ion[cdI[k]]=1;flashN+=na;return na;}
    for(let k=1;k<na;k++)cdW[k]+=cdW[k-1];
    let n=0;
    for(let p=0;p<BK;p++){const u=rnd()*ta;let lo=0,hi=na-1;while(lo<hi){const md=(lo+hi)>>1;if(cdW[md]<u)lo=md+1;else hi=md;}
      const i=cdI[lo];if(ion[i]<1){ion[i]=1;n++;}}
    flashN+=n;return n;
  }
  // 2) iniciação: campo de ruptura completo em qualquer ponto; sorteia um
  let ni=0,ti=0;
  for(let y=1;y<H-1;y++){const row=y*W;
    for(let x=1;x<W-1;x++){const i=row+x;if(!ebA[mat[i]])continue;const r=brk(i,1);
      if(r){const w=Math.pow(r,BETA);cdI[ni]=i;cdW[ni]=w;ti+=w;ni++;}}}
  if(ni){let u=rnd()*ti,k=0;for(;k<ni-1;k++){u-=cdW[k];if(u<=0)break;}ion[cdI[k]]=1;return 1;}
  return 0;
}
let eCnt=0;
function elecStep(){
  if(++eCnt<EST)return;eCnt=0;const dt=DT_E*EST;
  if(!elecCoef(dt)){
    if(elecLive){Vf.fill(0);IxC.fill(0);IyC.fill(0);ionL.fill(0);elecLive=false;iMax=0;}
    return;
  }
  elecLive=true;
  label();sor(ITER);dpcg(-60,1e-9);sor(3);fluxJ(dt); // SOR de novo: dielétricos acompanham o novo potencial dos condutores
 smoothChg();
  // líder escalonado: sub-passos rápidos enquanto houver ruptura nova
  let nb=breakdown(),k=0;
  while(nb>0&&k<2){elecCoef(DT_F);label();sor(2);dpcg(-30,1e-9,true);sor(2);fluxJ(DT_F);smoothChg();nb=breakdown();k++;}
  // recombinação do plasma (τ≈50 ms) e resfriamento radiativo do canal
  const grow=nb>0; // líder ainda avançando: canal mantém-se ionizado até completar
  if(!grow)flashN=0;
  const f=Math.exp(-dt/.05);let vm=0,im=0;
  for(let i=0;i<N;i++){
    const q=chg[i];if(q>QMAX)chg[i]=QMAX;else if(q<-QMAX)chg[i]=-QMAX; // limite de corona: ε0·E_rup·área
    if(ion[i]>0&&!grow){const v=ion[i]*f;ion[i]=v<1e-3?0:v;}
    const t=T[i],mi=mat[i];if(t>1500&&kind[mi]===GAS)T[i]=t-(t-1500)*.3; // resfriamento radiativo/expansão do canal
    else if(t>VAPT&&kind[mi]!==GAS&&!fixV[i]&&mi!==PAREDE){convert(i,mi===AGUA||mi===SALMOURA?VAPOR:FUMACA);mxA[i]=0;myA[i]=0;} // vaporização explosiva (⚠️ vapor do material ≈ fumaça)
    if(dom[i]){const a=Vf[i]<0?-Vf[i]:Vf[i];if(a>vm)vm=a;const c=(IxC[i]<0?-IxC[i]:IxC[i])+(IyC[i]<0?-IyC[i]:IyC[i]);if(c>im)im=c;}
  }
  if(!(vm<1e10)){Vf.fill(0);chg.fill(0);vm=0;im=0;console.warn('solver elétrico reiniciado');} // proteção numérica
  vScale+=((vm>1?vm:1)-vScale)*(vm>vScale?.5:.05);iMax=im;
}

// =====================================================================
//  MAGNETISMO (placa fina de espessura h, grade grossa 100×100, FFT 256²)
// =====================================================================
const SH=2,CF=1<<SH,CW=W>>SH,CH=H>>SH,CN=CW*CH,HC=CF*DX,FB=7,FN=1<<FB,FS=FN*FN,SOFT2=(.5*HC)*(.5*HC);
const Bx=new Float32Array(CN),By=new Float32Array(CN),Bz=new Float32Array(CN),Bm=new Float32Array(CN),gBx=new Float32Array(CN),gBy=new Float32Array(CN);
const tmx=new Float32Array(CN),tmy=new Float32Array(CN),tmz=new Float32Array(CN),
  pmx=new Float32Array(CN),pmy=new Float32Array(CN),cxs=new Float32Array(CN),cys=new Float32Array(CN),magC=new Uint8Array(CN);
const Bfine=new Float32Array(N),poleT=new Float32Array(N).fill(-1);
const rv=new Uint16Array(FN),fcs=new Float64Array(FN/2),fsn=new Float64Array(FN/2);
for(let i=0;i<FN;i++){let r=0;for(let b=0;b<FB;b++)if(i&(1<<b))r|=1<<(FB-1-b);rv[i]=r;}
for(let k=0;k<FN/2;k++){fcs[k]=Math.cos(2*Math.PI*k/FN);fsn[k]=Math.sin(2*Math.PI*k/FN);}
const TW=[];for(let size=2;size<=FN;size<<=1){const h=size>>1,c=new Float64Array(h),s=new Float64Array(h);for(let k=0;k<h;k++){c[k]=Math.cos(2*Math.PI*k/size);s[k]=Math.sin(2*Math.PI*k/size);}TW.push([c,s]);}
// FFT radix-2; passe de colunas vetorizado (butterflies sobre linhas inteiras, acesso contíguo)
function rowF(re,im,o,inv){
  for(let i=0;i<FN;i++){const j=rv[i];if(j>i){const a=o+i,b=o+j;let t=re[a];re[a]=re[b];re[b]=t;t=im[a];im[a]=im[b];im[b]=t;}}
  for(let L=0,size=2;size<=FN;size<<=1,L++){const half=size>>1,wc=TW[L][0],ws=TW[L][1];
    for(let s0=o;s0<o+FN;s0+=size)for(let k=0;k<half;k++){const a=s0+k,b=a+half,wr=wc[k],wi=inv?ws[k]:-ws[k];
      const br=re[b],bi=im[b],xr=br*wr-bi*wi,xi=br*wi+bi*wr;re[b]=re[a]-xr;im[b]=im[a]-xi;re[a]+=xr;im[a]+=xi;}}
}
function colF(re,im,inv){
  for(let i=0;i<FN;i++){const j=rv[i];if(j>i){const a=i*FN,b=j*FN;for(let c=0;c<FN;c++){let t=re[a+c];re[a+c]=re[b+c];re[b+c]=t;t=im[a+c];im[a+c]=im[b+c];im[b+c]=t;}}}
  for(let L=0,size=2;size<=FN;size<<=1,L++){const half=size>>1,wc=TW[L][0],ws=TW[L][1];
    for(let s0=0;s0<FN;s0+=size)for(let k=0;k<half;k++){const a=(s0+k)*FN,b=a+half*FN,wr=wc[k],wi=inv?ws[k]:-ws[k];
      for(let c=0;c<FN;c++){const ia=a+c,ib=b+c,br=re[ib],bi=im[ib],xr=br*wr-bi*wi,xi=br*wi+bi*wr;re[ib]=re[ia]-xr;im[ib]=im[ia]-xi;re[ia]+=xr;im[ia]+=xi;}}}
}
function fft2(re,im,inv,nr){
  if(!inv){for(let r=0;r<nr;r++)rowF(re,im,r*FN,false);colF(re,im,false);}
  else{colF(re,im,true);for(let r=0;r<nr;r++)rowF(re,im,r*FN,true);const s=1/FS;for(let i=0;i<nr*FN;i++){re[i]*=s;im[i]*=s;}}
}
let Kcx,Kcy,Kxx,Kxy,Kyy,Kzz;
function mkK(f){const re=new Float64Array(FS),im=new Float64Array(FS);
  for(let dy=1-CH;dy<CH;dy++)for(let dx=1-CW;dx<CW;dx++){if(!dx&&!dy)continue;
    const rx=dx*HC,ry=dy*HC,r2=rx*rx+ry*ry+SOFT2,r3=r2*Math.sqrt(r2);
    re[((dy+FN)&(FN-1))*FN+((dx+FN)&(FN-1))]=f(rx,ry,r2,r3);}
  fft2(re,im,false,FN);return [re,im];}
function buildKernels(){ // r = campo − fonte; y para baixo, z para dentro da tela
  Kcx=mkK((x,y,r2,r3)=>KM*y/r3);Kcy=mkK((x,y,r2,r3)=>-KM*x/r3);                 // Biot–Savart: dB = μ0/4π · I dl × r / r³
  // dipolo de coluna 2D (corpo extrudado em z, ⚠️ assumido): B = μ0/2π·[2(λ·r̂)r̂ − λ]/r², λ = M·área = M·VOL/DX.
  // Com placa fina de 10 cm o campo de um ímã de 1 m cairia ~20× (0,03 T) e não atrairia nada.
  const K2=2*KM/DX;
  Kxx=mkK((x,y,r2)=>K2*(2*x*x/r2-1)/r2);Kxy=mkK((x,y,r2)=>K2*2*x*y/r2/r2);Kyy=mkK((x,y,r2)=>K2*(2*y*y/r2-1)/r2);
  Kzz=mkK((x,y,r2,r3)=>-KM/r3);
}
const Pr=new Float64Array(FS),Pi=new Float64Array(FS),Ar=new Float64Array(FS),Ai=new Float64Array(FS),Zr=new Float64Array(FS),Zi=new Float64Array(FS),
  O1r=new Float64Array(FS),O1i=new Float64Array(FS),O2r=new Float64Array(FS),O2i=new Float64Array(FS);
function magClear(){Bx.fill(0);By.fill(0);Bz.fill(0);Bm.fill(0);gBx.fill(0);gBy.fill(0);Bfine.fill(0);bMax=0;}
// convolução: P = momentos no plano (x+iy) → B no plano; A = correntes, Z = momentos z → Bz. Resultado em O1/O2.
function magConv(mx,my,ax,ay,mz,inP,inA,inZ){
  if(inP){Pr.fill(0);Pi.fill(0);
    for(let cy=0;cy<CH;cy++)for(let cx=0;cx<CW;cx++){const c=cy*CW+cx,k=cy*FN+cx;Pr[k]=mx[c];Pi[k]=my[c];}
    fft2(Pr,Pi,false,CH);}
  if(inA){Ar.fill(0);Ai.fill(0);for(let cy=0;cy<CH;cy++)for(let cx=0;cx<CW;cx++){const c=cy*CW+cx,k=cy*FN+cx;Ar[k]=ax[c];Ai[k]=ay[c];}fft2(Ar,Ai,false,CH);}
  if(inZ){Zr.fill(0);Zi.fill(0);for(let cy=0;cy<CH;cy++)for(let cx=0;cx<CW;cx++){const c=cy*CW+cx,k=cy*FN+cx;Zr[k]=mz[c];}fft2(Zr,Zi,false,CH);}
  const kxr=Kxx[0],kxi=Kxx[1],kyr=Kxy[0],kyi=Kxy[1],kzr=Kyy[0],kzi=Kyy[1],
    car=Kcx[0],cai=Kcx[1],cbr=Kcy[0],cbi=Kcy[1],zzr=Kzz[0],zzi=Kzz[1],MK=FN-1,oZ=inA||inZ;
  for(let k=0;k<FS;k++){
    const kn=(((FN-(k>>FB))&MK)<<FB)|((FN-(k&MK))&MK);
    if(inP){ // desempacota FFT de (mx + i·my)
      const pr=Pr[k],pi=Pi[k],qr=Pr[kn],qi=Pi[kn];
      const mxr=(pr+qr)*.5,mxi=(pi-qi)*.5,myr=(pi+qi)*.5,myi=(qr-pr)*.5;
      const bxr=kxr[k]*mxr-kxi[k]*mxi+kyr[k]*myr-kyi[k]*myi,bxi=kxr[k]*mxi+kxi[k]*mxr+kyr[k]*myi+kyi[k]*myr;
      const byr=kyr[k]*mxr-kyi[k]*mxi+kzr[k]*myr-kzi[k]*myi,byi=kyr[k]*mxi+kyi[k]*mxr+kzr[k]*myi+kzi[k]*myr;
      O1r[k]=bxr-byi;O1i[k]=bxi+byr;
    }
    if(oZ){let r=0,im=0;
      if(inA){const pr=Ar[k],pi=Ai[k],qr=Ar[kn],qi=Ai[kn];
        const cxr=(pr+qr)*.5,cxi=(pi-qi)*.5,cyr=(pi+qi)*.5,cyi=(qr-pr)*.5;
        r=car[k]*cxr-cai[k]*cxi+cbr[k]*cyr-cbi[k]*cyi;im=car[k]*cxi+cai[k]*cxr+cbr[k]*cyi+cbi[k]*cyr;}
      if(inZ){const zr=Zr[k],zi=Zi[k];r+=zzr[k]*zr-zzi[k]*zi;im+=zzr[k]*zi+zzi[k]*zr;}
      O2r[k]=r;O2i[k]=im;
    }
  }
  if(inP)fft2(O1r,O1i,true,CH);if(oZ)fft2(O2r,O2i,true,CH);
}
// 1) campo externo: ímãs permanentes + correntes; 2) magnetização induzida no ferro por esse campo (sem
// realimentação ferro→ferro: ⚠️ assumido; com χ alto a realimentação numa grade grossa diverge) somada por cima
function magStep(){
  tmx.fill(0);tmy.fill(0);tmz.fill(0);pmx.fill(0);pmy.fill(0);cxs.fill(0);cys.fill(0);magC.fill(0);
  let nP=0,nI=0,iT=0;
  for(let y=0;y<H;y++){const cr=(y>>SH)*CW;
    for(let x=0;x<W;x++){const i=y*W+x,m=mat[i],c=cr+(x>>SH);
      if(m===IMA){poleT[i]=-1;const a=mxA[i],b=myA[i];if(a||b){pmx[c]+=a*VOL;pmy[c]+=b*VOL;nP++;magC[c]=1;}}
      else if(mgA[m]===1&&T[i]<tcA[m])nI++;
      const ix=IxC[i],iy=IyC[i];if(ix!==0||iy!==0){cxs[c]+=ix*DX;cys[c]+=iy*DX;const a=(ix<0?-ix:ix)+(iy<0?-iy:iy);if(a>iT)iT=a;}
    }}
  const curOn=elecLive&&iT>.05;
  if(!nP&&!curOn){if(magLive){magClear();magLive=false;}return;}
  if(!Kcx) buildKernels();
  magLive=true;
  const SELF=MU0*.5/(VOL*CF*CF); // termo próprio da célula: B=μ0·M/2 (desmagnetização de quadrado 2D)
  magConv(pmx,pmy,cxs,cys,null,nP>0,curOn,false);
  for(let cy=0;cy<CH;cy++)for(let cx=0;cx<CW;cx++){const c=cy*CW+cx,k=cy*FN+cx;
    Bx[c]=nP?O1r[k]+SELF*pmx[c]:0;By[c]=nP?O1i[k]+SELF*pmy[c]:0;Bz[c]=curOn?O2r[k]:0;}
  if(nI){let any=false,anyZ=false;
    for(let y=0;y<H;y++){const cr=(y>>SH)*CW;
      for(let x=0;x<W;x++){const i=y*W+x,m=mat[i];if(mgA[m]!==1||!(T[i]<tcA[m]))continue;const c=cr+(x>>SH);
        const bx=Bx[c],by=By[c],bz=Bz[c],b=Math.sqrt(bx*bx+by*by+bz*bz);
        if(b>1e-9){let M=chiA[m]*b/MU0;if(M>msA[m])M=msA[m];const f=M*VOL/b;tmx[c]+=f*bx;tmy[c]+=f*by;tmz[c]+=f*bz;any=true;if(bz)anyZ=true;}}}
    if(any){magConv(tmx,tmy,null,null,tmz,true,false,anyZ);
      for(let cy=0;cy<CH;cy++)for(let cx=0;cx<CW;cx++){const c=cy*CW+cx,k=cy*FN+cx;
        Bx[c]+=O1r[k];By[c]+=O1i[k];if(anyZ)Bz[c]+=O2r[k];}}
  }
  let bm=0;
  for(let c=0;c<CN;c++){const bx=Bx[c],by=By[c],bz=Bz[c],b=Math.sqrt(bx*bx+by*by+bz*bz);Bm[c]=b;if(b>bm)bm=b;}
  bMax=bm;
  for(let cy=0;cy<CH;cy++)for(let cx=0;cx<CW;cx++){const c=cy*CW+cx;
    gBx[c]=cx>0&&cx<CW-1?(Bm[c+1]-Bm[c-1])/(2*HC):0;gBy[c]=cy>0&&cy<CH-1?(Bm[c+CW]-Bm[c-CW])/(2*HC):0;}
  // |B| fino (bilinear) e posição relativa dos polos nos ímãs
  for(let y=0;y<H;y++){let gy=(y+.5)/CF-.5;if(gy<0)gy=0;if(gy>CH-1.001)gy=CH-1.001;const y0=gy|0,fy=gy-y0;
    for(let x=0;x<W;x++){let gx=(x+.5)/CF-.5;if(gx<0)gx=0;if(gx>CW-1.001)gx=CW-1.001;const x0=gx|0,fx=gx-x0,c=y0*CW+x0;
      Bfine[y*W+x]=(Bm[c]*(1-fx)+Bm[c+1]*fx)*(1-fy)+(Bm[c+CW]*(1-fx)+Bm[c+CW+1]*fx)*fy;}}
  if(nP) for(let i=0;i<N;i++){if(mat[i]!==IMA)continue;const a=mxA[i],b=myA[i];if(!a&&!b)continue;
    const dx=Math.abs(a)>=Math.abs(b)?(a>0?1:-1):0,dy=dx?0:(b>0?1:-1),d=dx+dy*W;
    let n=0,j=i+d;while(n<80&&mat[j]===IMA){n++;j+=d;}let s=0;j=i-d;while(s<80&&mat[j]===IMA){s++;j-=d;}
    poleT[i]=(s+.5)/(n+s+1);}
}

// ---------- força externa em partícula solta: qE (carga) + |M|∇|B| (ferromagnético) ----------
function extMove(i,m,k,x,y){
  let ax=0,ay=0;
  if(magLive&&mgA[m]===1&&T[i]<tcA[m]){const c=(y>>SH)*CW+(x>>SH),b=Bm[c];
    if(b>1e-6){let M=chiA[m]*b/MU0;if(M>msA[m])M=msA[m];const r=rho0[m];ax+=M*gBx[c]/r;ay+=M*gBy[c]/r;}}
  const q=chg[i];
  if(elecLive&&(q>1e-11||q<-1e-11)){const f=q/(densOf(i)*VOL);ax+=f*(Vf[i-1]-Vf[i+1])/(2*DX);ay+=f*(Vf[i-W]-Vf[i+W])/(2*DX);}
  if(ax*ax+ay*ay<2.2) return false; // < 0,15 g: ignora
  const tx=ax,ty=ay+(gm>.02?G0*gm:0),atx=tx<0?-tx:tx,aty=ty<0?-ty:ty;
  let j,p;
  if(k!==MOV&&rnd()*(atx+aty)<atx){j=i+(tx>0?1:-1);p=atx;}else{j=i+(ty>0?W:-W);p=aty;}
  const mj=mat[j],kj=kind[mj];
  if(mv[j]!==stamp&&(kj===GAS||(kj===LIQ&&mj!==m))){
    let pr=p/(G0*1.5);if(pr>1)pr=1;if(k===LIQ)pr*=mobOf(i,m);
    if(rnd()<pr){swap(i,j);return true;}
  }
  return ty<0; // puxado para cima mais que a gravidade: fica preso; senão cai normalmente
}

// =====================================================================
//  CLIMA: tempestade (perfil térmico forçado, nuvem, graupel, chuva, eletrificação não-indutiva)
// =====================================================================
const RAIN=5e-4,HAIL=2.5e-4,QC=6e-10;
function weatherStep(){
  const CB=CBASE;
  for(let y=1;y<=CB+14;y++){const row=y*W,Tp=y<=CB?233.15+43*y/CB:276; // −40 °C no topo → +3 °C na base (lapse rate comprimido, assumido)
    for(let x=1;x<W-1;x++){const i=row+x,m=mat[i];
      if(m===VAPOR&&y<CB&&T[i]<360&&rnd()<.02){setCell(i,NUVEM,Tp);continue;} // recondensa
      if(m===AR){if(y>2&&y<CB&&rnd()<.004){const c=chg[i];setCell(i,NUVEM,Tp);chg[i]=c;}continue;}
      if(m===NUVEM){
        if(y>CB){if(rnd()<.06){mat[i]=AR;gd[i]=P[i]*molar[AR]/(RG*T[i]);}continue;} // evapora abaixo da base
        T[i]+=(Tp-T[i])*.08;
        if(rnd()<.35){const r=rnd(),j=r<.4?i-W:r<.6?i+W:r<.8?i-1:i+1,mj=mat[j]; // turbulência c/ corrente ascendente
          if((mj===NUVEM||mj===AR)&&dom[j]&&j>=2*W)swap(i,j);}
        const t=T[i];
        if(t>273.6&&rnd()<RAIN&&kind[mat[i+W]]===GAS)convert(i,AGUA);
        else if(t<258&&rnd()<HAIL)convert(i,GRANIZO);
        continue;
      }
      if(m===GRANIZO&&y<CB){ // colisão graupel–cristal: inversão de sinal a −15 °C (Takahashi)
        const j=nbOf(i,rnd());
        if(mat[j]===NUVEM&&T[j]<263&&rnd()<.6){const dq=T[i]<258?-QC:QC;chg[i]+=dq;
          // o cristal leve sobe com a corrente ascendente enquanto o graupel cai: separação macroscópica
          // (transporte comprimido no tempo: a carga do cristal vai para uma célula de nuvem 8–25 linhas acima — assumido)
          const yt=y-8-((rnd()*18)|0),xt=x+((rnd()*7)|0)-3,k2=yt*W+xt;
          if(yt>=2&&xt>0&&xt<W-1&&mat[k2]===NUVEM)chg[k2]-=dq;else chg[j]-=dq;}
      }
    }}
}

function emReset(){
  chg.fill(0);ion.fill(0);mxA.fill(0);myA.fill(0);Vf.fill(0);IxC.fill(0);IyC.fill(0);ionL.fill(0);poleT.fill(-1);
  Uprev=0;magClear();magLive=false;elecLive=false;vScale=1;iMax=0;nStrikes=0;flash=0;
  if(typeof setStorm==='function')setStorm(false);
}

// =====================================================================
//  OVERLAY: linhas de campo (B, E), ⊙/⊗ de Bz, partículas de corrente, brilho do plasma, clarão
// =====================================================================
const ov=document.getElementById('ov'),oc=ov.getContext('2d'),flashEl=document.getElementById('flash');
let ovW=1,ovH=1,dpr=1;
function fitOv(){dpr=Math.min(2,devicePixelRatio||1);const r=stage.getBoundingClientRect();ovW=Math.max(1,Math.round(r.width*dpr));ovH=Math.max(1,Math.round(r.height*dpr));ov.width=ovW;ov.height=ovH;}
addEventListener('resize',fitOv);
const OCC=CW,occ=new Int32Array(OCC*OCC),Exc=new Float32Array(CN),Eyc=new Float32Array(CN),Vc=new Float32Array(CN),smp=[0,0];
let lines=[],linesView=-1,lineMax=1;
function sampleXY(fx,fy,x,y){
  let gx=x-.5,gy=y-.5;if(gx<0)gx=0;if(gx>CW-1.001)gx=CW-1.001;if(gy<0)gy=0;if(gy>CH-1.001)gy=CH-1.001;
  const x0=gx|0,y0=gy|0,ax=gx-x0,ay=gy-y0,c=y0*CW+x0;
  smp[0]=(fx[c]*(1-ax)+fx[c+1]*ax)*(1-ay)+(fx[c+CW]*(1-ax)+fx[c+CW+1]*ax)*ay;
  smp[1]=(fy[c]*(1-ax)+fy[c+1]*ax)*(1-ay)+(fy[c+CW]*(1-ax)+fy[c+CW+1]*ax)*ay;
}
function trace(fx,fy,sx,sy,dir,id,fmin,out){
  let x=sx,y=sy;
  for(let n=0;n<320;n++){
    sampleXY(fx,fy,x,y);let m=Math.hypot(smp[0],smp[1]);if(m<fmin)break;
    const hx=x+dir*smp[0]/m*.25,hy=y+dir*smp[1]/m*.25;             // RK2 (ponto médio)
    sampleXY(fx,fy,hx,hy);m=Math.hypot(smp[0],smp[1]);if(m<fmin)break;
    x+=dir*smp[0]/m*.5;y+=dir*smp[1]/m*.5;
    if(x<0||y<0||x>=CW||y>=CH)break;
    const o=(y|0)*OCC+(x|0);if(occ[o]&&occ[o]!==id)break;occ[o]=id;
    out.push(x,y,m);
  }
}
function buildLines(fx,fy,fmin,seeds){
  occ.fill(0);lines=[];let id=1,mx=fmin;
  for(let s=0;s<seeds.length&&lines.length<320;s+=2){
    const sx=seeds[s],sy=seeds[s+1],o0=(sy|0)*OCC+(sx|0);if(o0<0||o0>=OCC*OCC||occ[o0])continue;
    const b=[],f=[];trace(fx,fy,sx,sy,-1,id,fmin,b);trace(fx,fy,sx,sy,1,id,fmin,f);
    const pts=[];for(let k=b.length-3;k>=0;k-=3)pts.push(b[k],b[k+1],b[k+2]);for(let k=0;k<f.length;k++)pts.push(f[k]);
    if(pts.length>=15){let sm=0;for(let k=2;k<pts.length;k+=3)sm+=Math.log(pts[k]);const g=Math.exp(sm/(pts.length/3));if(g>mx)mx=g;lines.push({p:pts,g});}
    id++;
  }
  lineMax=mx;
}
function seedsB(){
  const s=[];
  for(let c=0;c<CN;c++)if(magC[c]){const cx=c%CW,cy=(c/CW)|0;for(let k=0;k<4;k++)s.push(cx+.5+MDX[k]*.7,cy+.5+MDY[k]*.7);}
  for(let cy=1;cy<CH;cy+=2)for(let cx=1;cx<CW;cx+=2)s.push(cx+rnd(),cy+rnd());
  return s;
}
function seedsE(){
  const s=[];
  for(let cy=0;cy<CH;cy++)for(let cx=0;cx<CW;cx++){const i=(cy*CF+1)*W+cx*CF+1;
    if(fixV[i]||chg[i]>1e-9||chg[i]<-1e-9||ionL[i]>.3)s.push(cx+.5,cy+.5);}
  for(let k=s.length-2;k>0;k-=2){const j=((rnd()*(k/2+1))|0)*2;let t=s[k];s[k]=s[j];s[j]=t;t=s[k+1];s[k+1]=s[j+1];s[j+1]=t;}
  for(let cy=1;cy<CH;cy+=2)for(let cx=1;cx<CW;cx+=2)s.push(cx+rnd(),cy+rnd());
  return s;
}
function computeEc(){
  for(let cy=0;cy<CH;cy++)for(let cx=0;cx<CW;cx++){const i=(cy*CF+1)*W+cx*CF+1;Vc[cy*CW+cx]=(Vf[i]+Vf[i+1]+Vf[i+W]+Vf[i+W+1])*.25;}
  let em=0;
  for(let cy=0;cy<CH;cy++)for(let cx=0;cx<CW;cx++){const c=cy*CW+cx;
    const ex=cx>0&&cx<CW-1?-(Vc[c+1]-Vc[c-1])/(2*HC):0,ey=cy>0&&cy<CH-1?-(Vc[c+CW]-Vc[c-CW])/(2*HC):0;
    Exc[c]=ex;Eyc[c]=ey;const m=Math.hypot(ex,ey);if(m>em)em=m;}
  return em;
}
// partículas de corrente (sentido convencional, + → −)
const PN=700,ppx=new Float32Array(PN),ppy=new Float32Array(PN),ppl=new Int16Array(PN),curList=new Int32Array(8000);let curN=0;
function updCurList(ith){curN=0;for(let i=0;i<N;i++){if(!dom[i])continue;const a=(IxC[i]<0?-IxC[i]:IxC[i])+(IyC[i]<0?-IyC[i]:IyC[i]);
  if(a>ith){if(curN<curList.length)curList[curN++]=i;else curList[(rnd()*curN)|0]=i;}}}
function drawOverlay(fn){
  oc.setTransform(1,0,0,1,0,0);oc.clearRect(0,0,ovW,ovH);oc.globalCompositeOperation='source-over';
  const sx=ovW/W,sy=ovH/H,cs=CF*sx,ith=Math.max(1e-3,iMax*1e-4);
  // linhas de campo
  const wantB=view===4&&magLive,wantE=view===3&&elecLive;
  if(wantB||wantE){
    if(fn%10===0||linesView!==view){
      if(wantB){const fmin=Math.max(2e-6,bMax*2e-5);buildLines(Bx,By,fmin,seedsB());}
      else{const em=computeEc();buildLines(Exc,Eyc,Math.max(50,em*1e-3),seedsE());}
      linesView=view;
    }
    oc.lineWidth=1.2*dpr;oc.setLineDash([5*dpr,7*dpr]);oc.lineDashOffset=-(fn*.7*dpr)%(12*dpr);
    const lmin=Math.log(Math.max(1e-30,lineMax*1e-4)),lr=Math.log(lineMax)-lmin;
    for(const L of lines){let a=(Math.log(L.g)-lmin)/lr;a=a<0?0:a>1?1:a;a=.18+a*.7;
      oc.strokeStyle=wantB?'rgba(225,235,255,'+a.toFixed(2)+')':'rgba(255,236,190,'+a.toFixed(2)+')';
      oc.beginPath();const p=L.p;oc.moveTo(p[0]*cs,p[1]*cs);for(let k=3;k<p.length;k+=3)oc.lineTo(p[k]*cs,p[k+1]*cs);oc.stroke();}
    oc.setLineDash([]);
  }else{lines=[];linesView=-1;}
  // Bz: ⊙ saindo da tela, ⊗ entrando
  if(view===4&&magLive){const thr=Math.max(5e-6,bMax*2e-3);oc.lineWidth=1.1*dpr;
    for(let cy=1;cy<CH;cy+=3)for(let cx=1;cx<CW;cx+=3){const c=cy*CW+cx,bz=Bz[c],bp=Math.hypot(Bx[c],By[c]);
      if((bz<0?-bz:bz)<thr||(bz<0?-bz:bz)<bp*.7)continue;
      const r=(2.5+Math.min(4,Math.log10((bz<0?-bz:bz)/thr)*1.6))*dpr,X=(cx+.5)*cs,Y=(cy+.5)*cs;
      oc.strokeStyle=bz<0?'rgba(255,120,110,.85)':'rgba(120,200,255,.85)';oc.beginPath();oc.arc(X,Y,r,0,6.2832);oc.stroke();
      if(bz<0){oc.fillStyle=oc.strokeStyle;oc.beginPath();oc.arc(X,Y,r*.3,0,6.2832);oc.fill();}
      else{const q=r*.6;oc.beginPath();oc.moveTo(X-q,Y-q);oc.lineTo(X+q,Y+q);oc.moveTo(X+q,Y-q);oc.lineTo(X-q,Y+q);oc.stroke();}}}
  // partículas de corrente
  if(elecLive&&(view===0||view===3)){
    if(fn%8===0)updCurList(ith);
    const act=Math.min(PN,curN*2);oc.fillStyle='rgba(255,236,140,.9)';
    for(let p=0;p<act;p++){
      if(ppl[p]<=0){if(!curN)continue;const i=curList[(rnd()*curN)|0];ppx[p]=i%W+rnd();ppy[p]=(i/W|0)+rnd();ppl[p]=30+(rnd()*50|0);}
      const i=(ppy[p]|0)*W+(ppx[p]|0),ix=IxC[i],iy=IyC[i],a=Math.hypot(ix,iy);
      if(!(a>ith)){ppl[p]=0;continue;}
      let v=.25+.18*Math.log10(a/ith);if(v>1.1)v=1.1;
      ppx[p]+=ix/a*v;ppy[p]+=iy/a*v;ppl[p]--;
      oc.fillRect(ppx[p]*sx-1.2*dpr,ppy[p]*sy-1.2*dpr,2.4*dpr,2.4*dpr);
    }
  }else ppl.fill(0);
  // brilho do plasma
  if(plasmaN){oc.globalCompositeOperation='lighter';
    for(let k=0;k<plasmaN;k++){const i=plasmaL[k],x=i%W,y=i/W|0,io=Math.min(1,ionL[i]);
      oc.fillStyle='rgba(120,130,255,'+(.10*io).toFixed(3)+')';oc.fillRect((x-2)*sx,(y-2)*sy,5*sx,5*sy);
      oc.fillStyle='rgba(255,255,255,'+(.55*io).toFixed(3)+')';oc.fillRect(x*sx,y*sy,sx,sy);}
    oc.globalCompositeOperation='source-over';}
  // clarão (energia dissipada em plasma neste quadro)
  if(ePl>2e3){const f=Math.min(1,.25+Math.log10(ePl/2e3)*.25);if(f>flash)flash=f;if(strikeCD<=0&&f>.3){nStrikes++;strikeCD=15;}}
  ePl=0;eJoule=0;strikeCD--;flash*=.8;
  flashEl.style.opacity=flash>.01?(flash*.35).toFixed(3):'0';
}

// ================= inicialização =================
fitOv();clearAll();PRESETS.tank();updateLegend();
requestAnimationFrame(loop);
