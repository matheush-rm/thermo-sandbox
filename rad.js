'use strict';
// ================= radiação e física nuclear =================
// Decaimento radioativo emite α, β, γ (e nêutrons em fontes Am-Be / fissão espontânea).
//  α: alcance ~4 cm no ar → só atinge as células vizinhas (perigosa se inalada/em contato; papel/pele param).
//  β: alcance ~3 m no ar, poucos mm em sólidos/água: atenuação μβ ≈ ρ/0,4 g·cm⁻² (alumínio/plástico param).
//  γ: lei exponencial I = I₀·e^(−μx), μ = μm·ρ com μm ≈ 0,06 cm²/g (~1 MeV, quase igual p/ todo material;
//     chumbo ×1,4) → camada semirredutora: água ~11 cm, concreto ~5 cm, chumbo ~0,8 cm.
//  Geometria: raios de Monte Carlo isotrópicos a partir de cada fonte → fluência cai com a distância
//  (⚠️ assumido: grade 2D → cai com 1/r em vez de 1/r²).
// Nêutrons (partículas): fissão (U-235, Pu-239) libera ν≈2,4–2,9 nêutrons rápidos; espalhamento em núcleos leves
//  (água, grafite, concreto, tecido) modera → nêutron térmico, que fissiona o U-235 ~100× mais; boro/mercúrio absorvem.
//  Captura ativa o material (aço, concreto, água salgada ficam radioativos). k = nêutrons de uma geração / anterior:
//  k>1 → reação em cadeia. Massa crítica emerge da fuga pela superfície (forma, refletor, densidade: m_c ∝ 1/ρ²).
//  Cadeia RÁPIDA muito supercrítica (bomba) → detonação nuclear; cadeia térmica (reator) → só calor (fusão, explosão de vapor).
// Dose: campo em Sv/h ("se uma pessoa estivesse aqui"); vida acumula dose: >3 Sv começa a morrer, mutação ∝ dose,
//  campo alto esteriliza (não nasce).
// ⚠️ Assumido: seções de choque reduzidas p/ massa crítica caber na grade de 10 cm (real: 1 célula de Pu já é supercrítica);
//  1 passo ≈ 3 min de exposição; meias-vidas longas (Co, Cs, Ra, U, Pu) ≫ tempo de jogo → fontes constantes;
//  atividade induzida/produtos de fissão decaem com meia-vida comprimida única (~1500 passos).
const HPS=.05,RAYS=1500,NCAP=6000,NDOSE=5,EFIS=.35,FPK=.004,BUMAX=3e6,NUKE_F=6e3,PFIS_RA=6,RDEC=Math.pow(.5,1/1500),RNDEC=Math.pow(.5,1/900);
// por material: Γγ, Γβ (Sv/h a 1 m por célula), α local (Sv/h), nêutrons espontâneos/passo, calor de decaimento (K/passo)
const GAM=new Float32Array(NM),BET=new Float32Array(NM),ALF=new Float32Array(NM),NSP=new Float32Array(NM),DHEAT=new Float32Array(NM);
// nêutrons: prob. por célula percorrida de fissão/captura (rápido F, térmico T), espalhamento, moderação (térmaliza por choque), ativação
const NFF=new Float32Array(NM),NFT=new Float32Array(NM),NCF=new Float32Array(NM),NCT=new Float32Array(NM),NSC=new Float32Array(NM),
  NMOD=new Float32Array(NM),NACT=new Float32Array(NM),NI=new Uint8Array(NM),NU=new Float32Array(NM),MUG=new Float32Array(NM),MUB=new Float32Array(NM),
  FISS=new Float32Array(NM),CRIT=new Float32Array(NM);
const HYD=[AGUA,SALMOURA,OLEO,MEL,ALCOOL,MADEIRA,MADEIRAE,ISOPOR,ISOPORC,GELO,VIDA,VIDAM,NITRO,C4,HNO3]; // ricos em hidrogênio: moderam
for(let m=0;m<NM;m++){
  const k=kind[m],r=rho0[m]/1000;
  if(k!==GAS){MUG[m]=.06*r*10;MUB[m]=r*10/.4;NSC[m]=Math.min(.3,.03*r+.02);NCF[m]=.0005;NCT[m]=.004;NMOD[m]=.02;NACT[m]=.002;NI[m]=1;}
}
MUG[CHUMBO]*=1.4;MUG[CHUMBOF]*=1.4;
for(const m of HYD){NSC[m]=.3;NMOD[m]=.45;NCT[m]=.012;NCF[m]=.001;NACT[m]=.0005;}
NCT[SALMOURA]=.035;NACT[SALMOURA]=.02;               // Na-24
NSC[CONCRETO]=.25;NMOD[CONCRETO]=.15;NCT[CONCRETO]=.01;NACT[CONCRETO]=.006;
for(const m of [ACO,SUCATA,LIMALHA,IMA,ACOF,NICROMO]){NSC[m]=.25;NCT[m]=.02;NACT[m]=.03;} // Fe-59/Co-60 por ativação
NSC[COBRE]=.22;NCT[COBRE]=.03;NACT[COBRE]=.03;
NSC[GRAFITE]=.28;NMOD[GRAFITE]=.12;NCT[GRAFITE]=.0003;NCF[GRAFITE]=0;NACT[GRAFITE]=.0005;
NSC[CHUMBO]=.2;NMOD[CHUMBO]=.005;NCT[CHUMBO]=.001;NACT[CHUMBO]=.0005;NSC[CHUMBOF]=.2;NCT[CHUMBOF]=.001;
NCT[BORO]=.95;NCF[BORO]=.3;NSC[BORO]=.05;NACT[BORO]=0;   // B-10(n,α)
NCT[MERCURIO]=.3;NCT[VHG]=0;
NSC[AMBE]=.15;
// físseis
const fis=(m,ff,ft,cf,ct,sc,nu,crit)=>{NFF[m]=ff;NFT[m]=ft;NCF[m]=cf;NCT[m]=ct;NSC[m]=sc;NU[m]=nu;NMOD[m]=.005;NACT[m]=0;FISS[m]=1;CRIT[m]=crit;};
fis(U235,.04,.40,.004,.07,.18,2.43,260);
fis(PU239,.05,.45,.004,.16,.18,2.88,60);
fis(UNAT,.0018,.05,.03,.034,.18,2.43,0);
fis(CORIO,.014,.16,.004,.05,.12,2.45,0);
NCT[PFIS]=.12;NSC[PFIS]=.06;  // venenos de fissão (Xe-135, Sm-149)
// emissores
GAM[COBALTO]=1;BET[COBALTO]=.2;DHEAT[COBALTO]=.004;
GAM[CESIO]=.35;BET[CESIO]=1;DHEAT[CESIO]=.006;
GAM[RADIO]=.08;ALF[RADIO]=4;DHEAT[RADIO]=.01;
ALF[RADONIO]=2;
GAM[U235]=.0005;ALF[U235]=.05;NSP[U235]=.0004;
GAM[UNAT]=.0002;ALF[UNAT]=.01;NSP[UNAT]=.0002;
GAM[PU239]=.002;ALF[PU239]=.5;NSP[PU239]=.003;DHEAT[PU239]=.003;  // Pu-240: fissão espontânea; Pu é morno ao toque
GAM[CORIO]=3;BET[CORIO]=2;DHEAT[CORIO]=.35;NSP[CORIO]=.0005;
GAM[AMBE]=.02;NSP[AMBE]=.4;
const SRCM=new Uint8Array(NM);for(let m=0;m<NM;m++)SRCM[m]=(GAM[m]||BET[m]||ALF[m]||NSP[m]||DHEAT[m])?1:0;
SRCM[RADIO]=1;SRCM[CESIO]=1;

// estado por célula: atividade induzida/contaminação (Sv/h a 1 m), queima, compressão, dose acumulada (vida), campo
const rA=new Float32Array(N),bu=new Float32Array(N),cmp=new Float32Array(N),dAcc=new Float32Array(N),radF=new Float32Array(N),radAcc=new Float32Array(N);
// nêutrons
const nX=new Float32Array(NCAP),nY=new Float32Array(NCAP),nUX=new Float32Array(NCAP),nUY=new Float32Array(NCAP),nW=new Float32Array(NCAP),nFast=new Uint8Array(NCAP);
let cmpAny=false,nCnt=0,fisN=0,fisF=0,fisPrev=0,fisFp=0,fisLast=-1,radMax=0,radLive=false,rAny=false,nNukes=0,lastKt=0,fisShow=0,nukeCD=0,radSrcN=0;
const srcI=new Int32Array(N),srcG=new Float32Array(N),srcB=new Float32Array(N);

function radReset(){nCnt=0;nNukes=0;lastKt=0;radF.fill(0);radAcc.fill(0);fisN=fisF=fisPrev=fisFp=fisShow=0;radMax=0;radLive=false;nukeCD=0;}
function radClearCell(i){rA[i]=0;bu[i]=0;cmp[i]=0;dAcc[i]=0;}

function nSpawn(x,y,w,fast){
  if(nCnt>=NCAP){ // roleta russa global: mata metade, dobra o peso das sobreviventes (sem viés na média)
    let k=0;for(let n=0;n<nCnt;n++)if(rnd()<.5){nX[k]=nX[n];nY[k]=nY[n];nUX[k]=nUX[n];nUY[k]=nUY[n];nW[k]=Math.min(1e12,nW[n]*2);nFast[k]=nFast[n];k++;}
    nCnt=k; if(rnd()<.5) return; w*=2;
  }
  const a=rnd()*6.2832;nX[nCnt]=x;nY[nCnt]=y;nUX[nCnt]=Math.cos(a);nUY[nCnt]=Math.sin(a);nW[nCnt]=w;nFast[nCnt]=fast;nCnt++;
}
function fission(j,w,fast){
  const m=mat[j];fisN+=w;if(fast)fisF+=w;fisLast=j;
  let dT=EFIS*w/CvA[m];if(dT>800)dT=800;T[j]+=dT;if(T[j]>5000)T[j]=5000;
  bu[j]+=w;rA[j]+=w*FPK;if(rA[j]>200)rA[j]=200;rAny=true;
  radAcc[j]+=w*.5; // γ imediato
  const nu=NU[m],n=nu|0,x=j%W,y=(j/W)|0,c=n+(rnd()<nu-n?1:0);
  for(let k=0;k<c;k++) nSpawn(x+rnd(),y+rnd(),w,1);
  if(bu[j]>BUMAX){mat[j]=PFIS;lat[j]=0;life[j]=0;rA[j]=Math.max(rA[j],PFIS_RA);} // combustível queimado → produtos de fissão
}
function neutronStep(){
  for(let n=nCnt-1;n>=0;n--){
    let x=nX[n],y=nY[n],ux=nUX[n],uy=nUY[n],w=nW[n],f=nFast[n],dead=false;
    const sub=f?12:2; // rápido ≫ térmico; geração de nêutrons ≪ tempo mecânico
    for(let s=0;s<sub;s++){
      x+=ux;y+=uy;
      if(x<1||y<1||x>=W-1||y>=H-1){dead=true;break;}
      const j=(y|0)*W+(x|0),m=mat[j];
      if(m===PAREDE){dead=true;break;}
      radAcc[j]+=NDOSE*w;
      if(!NI[m]) continue;
      const c=1+cmp[j];let r=rnd();
      const pf=(f?NFF[m]:NFT[m])*c; if(r<pf){fission(j,w,f);dead=true;break;} r-=pf;
      const pc=(f?NCF[m]:NCT[m])*c; if(r<pc){rA[j]+=w*NACT[m];if(rA[j]>200)rA[j]=200;if(NACT[m])rAny=true;radAcc[j]+=w*2;dead=true;break;} r-=pc;
      if(r<NSC[m]*c){const a=rnd()*6.2832;ux=Math.cos(a);uy=Math.sin(a);if(f&&rnd()<NMOD[m])f=0;}
    }
    if(dead){const l=--nCnt;nX[n]=nX[l];nY[n]=nY[l];nUX[n]=nUX[l];nUY[n]=nUY[l];nW[n]=nW[l];nFast[n]=nFast[l];}
    else{nX[n]=x;nY[n]=y;nUX[n]=ux;nUY[n]=uy;nFast[n]=f;}
  }
}
// raio γ (t=0) ou β (t=1) a partir da célula i, ângulo a, peso w
function castRay(i,a,w,t){
  const ux=Math.cos(a),uy=Math.sin(a),s=Math.max(ux<0?-ux:ux,uy<0?-uy:uy),dx=ux/s,dy=uy/s,L=1/s,wmin=w*.003;
  let x=i%W+.5,y=((i/W)|0)+.5;
  for(let q=0;q<300;q++){
    x+=dx;y+=dy;if(x<1||y<1||x>=W-1||y>=H-1)return;
    const j=(y|0)*W+(x|0),m=mat[j];if(m===PAREDE)return;
    radAcc[j]+=w*L;
    let mu;
    if(kind[m]===GAS) mu=t?gd[j]*.025:0; else mu=t?MUB[m]:MUG[m];
    if(mu>0){w*=Math.exp(-mu*L);if(w<wmin)return;}
  }
}
function radStep(){
  // 1) fontes: varredura; efeitos locais (α, calor, radônio, césio solúvel, decaimento da atividade induzida)
  let ns=0,S=0;const dec=rAny;rAny=false;cmpAny=false;
  for(let y=1;y<H-1;y++){const row=y*W;
    for(let x=1;x<W-1;x++){const i=row+x,m=mat[i];
      let g=0,b=0;
      if(dec){let a=rA[i];if(a>0){a*=RDEC;if(a<1e-5)a=0;rA[i]=a;if(a>0){rAny=true;g=a*.5;b=a*.5;if(a>.05){let dT=a*.004/(CvA[m]||1);T[i]+=dT>5?5:dT;}}}}
      if(SRCM[m]){
        g+=GAM[m];b+=BET[m];
        const al=ALF[m];if(al){radAcc[i]+=al;radAcc[i-1]+=al;radAcc[i+1]+=al;radAcc[i-W]+=al;radAcc[i+W]+=al;}
        if(DHEAT[m])T[i]+=DHEAT[m];
        if(NSP[m]&&rnd()<NSP[m]*(1+cmp[i]))nSpawn(x+rnd(),y+rnd(),1,1);
        if(m===AMBE&&(cmp[i-1]>.3||cmp[i+1]>.3||cmp[i-W]>.3||cmp[i+W]>.3)&&rnd()<.5)for(let k=0;k<20;k++)nSpawn(x+rnd(),y+rnd(),1,1); // iniciador esmagado: rajada
        if(m===RADIO){if(rnd()<.004){const j=nbOf(i,rnd());if(mat[j]===AR)setCell(j,RADONIO,T[j]);}} // emanação de radônio
        else if(m===RADONIO){if(rnd()<1-RNDEC){mat[i]=AR;gd[i]=P[i]*molar[AR]/(RG*T[i]);}} // → Po-218 (sólido, desprezado)
        else if(m===CESIO){const j=nbOf(i,rnd()),mj=mat[j];
          if((mj===AGUA||mj===SALMOURA)&&rnd()<.01){const t=T[i];setCell(i,mj,t);rA[i]=1.2;rA[j]+=.6;rAny=true;}} // CsCl dissolve: água contaminada
      }
      if(cmp[i]>0){cmpAny=true;cmp[i]*=.996;if(cmp[i]<.01)cmp[i]=0;else if(cmp[i]>.3){pvx[i]=0;pvy[i]=0;}} // inércia do tamper segura o núcleo comprimido
      if(g>0||b>0){srcI[ns]=i;srcG[ns]=g;srcB[ns]=b;ns++;S+=g+b;}
    }}
  radSrcN=ns;
  // 2) raios γ/β: orçamento fixo, peso igual; fonte de Γ (Sv/h a 1 m) emite 2π·10·Γ por passo (campo 2D = Γ·10/r)
  if(S>0){const k=RAYS/S,wr=62.83/k;
    for(let s=0;s<ns;s++){const i=srcI[s];
      let e=srcG[s]*k,n=e|0;if(rnd()<e-n)n++;for(let q=0;q<n;q++)castRay(i,rnd()*6.2832,wr,0);
      e=srcB[s]*k;n=e|0;if(rnd()<e-n)n++;for(let q=0;q<n;q++)castRay(i,rnd()*6.2832,wr,1);}}
  // 3) nêutrons
  fisPrev=fisPrev*.7+fisN*.3;fisN=0;fisF=0;
  // cadeia rápida: vários ciclos de transporte por passo (geração ~10 ns ≪ desmontagem mecânica ~µs)
  const cyc=cmpAny||(fisPrev>4&&fisFp>.6*fisPrev)?6:1;
  for(let c=0;c<cyc&&nCnt;c++) neutronStep();
  fisFp=fisFp*.7+fisF*.3;
  fisShow=fisShow*.9+fisN*.1;
  if(nukeCD>0) nukeCD--;
  else if(fisN>NUKE_F&&fisF>.6*fisN) nuke(fisLast);
  // 4) campo suavizado (média móvel)
  if(ns||nCnt||radLive){let mx=0;
    for(let i=0;i<N;i++){let v=radF[i]+(radAcc[i]-radF[i])*.12;if(v<1e-9)v=0;radF[i]=v;radAcc[i]=0;if(v>mx)mx=v;}
    radMax=mx;radLive=mx>1e-9;}
}

// ================= vida × radiação =================
// chamado por lifeHazard (sim.js). true = célula morreu
function radLife(i){
  const r=radF[i];let d=dAcc[i]*.9995;
  if(r>0){d+=r*HPS;
    if(mat[i]===VIDA&&rnd()<r*HPS*.08){mat[i]=VIDAM;shade[i]=(rnd()*8)|0;} // mutação ∝ dose (modelo linear sem limiar)
  }
  dAcc[i]=d;
  if(d>3&&rnd()<(d-3)*.02){lifeKill(i);return true;} // síndrome aguda da radiação (DL50 ≈ 4–5 Sv)
  return false;
}
const RAD_STERILE=1; // Sv/h: acima disso não nasce vida

// ================= detonação nuclear =================
function nuke(c){
  if(c<0) return;
  nukeCD=120;
  const cx=c%W,cy=(c/W)|0;
  // massa físsil e compressão perto do centro → eficiência (fração fissionada): ~1–2% canhão, ~10–20% implosão
  let cells=0,cs=0,crit=0;
  for(let y=Math.max(1,cy-25);y<=Math.min(H-2,cy+25);y++)for(let x=Math.max(1,cx-25);x<=Math.min(W-2,cx+25);x++){
    const i=y*W+x,m=mat[i];if(!CRIT[m])continue;cells++;const k=1+cmp[i];cs+=k*k;crit+=CRIT[m];}
  if(!cells){cells=1;cs=1;crit=150;}
  // sr = massas críticas (compressão conta ρ²); ⚠️ massa efetiva 0,25 kg/célula; 1 kg fissionado ≈ 17 kt
  const sr=cs/(crit/cells),eff=Math.min(.3,Math.max(.004,.008*sr*sr)),kt=Math.max(.05,cells*.25*17*eff);
  lastKt=kt;nNukes++;
  const Wkg=kt*1e6,Rfb=Math.min(60,6+9*Math.cbrt(kt)),Rfb2=Rfb*Rfb; // ⚠️ bola de fogo limitada (real: ~66 m·kt^0,4)
  // bola de fogo: vaporiza tudo (menos as paredes da câmara); parte vira precipitação radioativa
  for(let y=Math.max(1,cy-Rfb|0);y<=Math.min(H-2,cy+Rfb|0);y++)for(let x=Math.max(1,cx-Rfb|0);x<=Math.min(W-2,cx+Rfb|0);x++){
    const dx=x-cx,dy=y-cy,d2=dx*dx+dy*dy;if(d2>Rfb2)continue;
    const i=y*W+x,m=mat[i];if(m===PAREDE)continue;
    const f=1-Math.sqrt(d2)/Rfb,t=3000+6000*f;
    if(FISS[m]||m===PFIS||(kind[m]!==GAS&&rnd()<.12+.3*f)){setCell(i,PFIS,t);rA[i]=PFIS_RA*(2+8*f);}
    else{setCell(i,FOGO,t);life[i]=60+(rnd()*150|0);rA[i]=.3*f;}
    rAny=true;
  }
  // onda de choque (Mills via explosivos) e sobrepressão no campo de ar
  detonate(c,Wkg);
  if(airOn){const R=Rfb/AC;for(let y=0;y<AH;y++)for(let x=0;x<AW;x++){const dx=x-cx/AC,dy=y-cy/AC;if(dx*dx+dy*dy<R*R)aInj[y*AW+x]+=AP_MAX*1.5;}}
  // pulso térmico e radiação imediata (γ + nêutrons) em linha de visada; nêutrons ativam o que atingem
  const Rth=Math.min(260,40+60*Math.cbrt(kt));
  for(let a=0;a<1440;a++){const an=a/1440*6.2832,ux=Math.cos(an),uy=Math.sin(an);
    let x=cx+.5,y=cy+.5,wg=2e5,wn=1e4,heat=true;
    for(let q=0;q<Rth;q++){x+=ux;y+=uy;if(x<1||y<1||x>=W-1||y>=H-1)break;
      const j=(y|0)*W+(x|0),m=mat[j];if(m===PAREDE)break;
      const r=q+1,amp=1/r;
      radAcc[j]+=(wg+wn*NDOSE)*amp*.02;
      if(kind[m]!==GAS){
        if(heat){heat=false;const h=Math.min(2500,8e4*amp*Math.cbrt(kt)/(CvA[m]||1));T[j]+=h;if(ignA[m]>0&&T[j]>ignA[m])ignite(j);}
        rA[j]+=Math.min(50,wn*amp*NACT[m]*.05);wg*=Math.exp(-MUG[m]);wn*=Math.exp(-3*(NSC[m]+NCT[m]));
      }
    }}
  rAny=true;radLive=true;
  for(let n=0;n<nCnt;n++)nW[n]=1;nCnt=Math.min(nCnt,200); // cadeia termina: núcleo desmontado
  flash=8;
  toast('☢ Detonação nuclear: ~'+(kt>=1?kt.toFixed(kt<10?1:0)+' kt':(kt*1000).toFixed(0)+' t')+' de TNT');
}

// ================= vista: nêutrons na sobreposição =================
function drawRad(oc,sx,sy,dpr){
  if(!nCnt) return;
  const st=Math.max(1,(nCnt/1500)|0),s=1.6*dpr;
  oc.fillStyle='rgba(120,255,240,.85)';
  for(let n=0;n<nCnt;n+=st){if(!nFast[n])continue;oc.fillRect(nX[n]*sx-s/2,nY[n]*sy-s/2,s,s);}
  oc.fillStyle='rgba(255,255,255,.75)';
  for(let n=0;n<nCnt;n+=st){if(nFast[n])continue;oc.fillRect(nX[n]*sx-s/2,nY[n]*sy-s/2,s,s);}
}

// ================= cenários =================
function discM(cx,cy,R,m,t){for(let y=-R;y<=R;y++)for(let x=-R;x<=R;x++)if(x*x+y*y<=R*R+R*.8)rectM(cx+x,cy+y,cx+x,cy+y,m,t);}
function soup(x0,y0,x1,y1,p){for(let y=y0;y<=y1;y++)for(let x=x0;x<=x1;x++)if(rnd()<p)rectM(x,y,x,y,VIDA);}
function townM(gy){ // cidade-alvo: prédios com janelas, árvores, lago e colônias de vida
  rectM(1,gy,198,198,AREIA); rectM(60,gy,84,gy+6,AGUA);
  boxM(10,gy-40,34,gy-1,2,CONCRETO); rectM(10,gy-28,11,gy-22,VIDRO); rectM(33,gy-28,34,gy-22,VIDRO);
  boxM(160,gy-52,186,gy-1,2,CONCRETO); rectM(160,gy-36,161,gy-30,VIDRO);
  rectM(46,gy-22,47,gy-1,MADEIRA); rectM(41,gy-30,52,gy-23,MADEIRA);
  rectM(140,gy-18,141,gy-1,MADEIRA); rectM(136,gy-25,145,gy-19,MADEIRA);
  soup(14,gy-30,30,gy-4,.35); soup(164,gy-44,182,gy-4,.35); soup(110,gy-20,132,gy-2,.35);
}
Object.assign(PRESETS,{
  nukeImp(){ // implosão (tipo Fat Man/Trinity): núcleo de Pu subcrítico + tamper de U natural + lente de C-4; choque comprime ρ → supercrítico
    clearAll(); townM(176);
    rectM(99,122,101,175,ACO);                                     // torre de aço
    discM(100,108,12,ACO); discM(100,108,10,C4); discM(100,108,6,UNAT); discM(100,108,3,PU239); rectM(100,108,100,108,AMBE); // iniciador Am-Be
    heaterSet(700); rectM(99,95,101,95,RESIST); rectM(100,96,100,97,C4); // detonador: resistência aquece o C-4 (>247 °C)
  },
  nukeArm(){ // implosão com disparo manual: bunker com bateria; feche o vão do fio com Cobre (ou Faísca/Aquecer no C-4) → nicromo esquenta → C-4 → implosão
    clearAll(); setVolt(60); const gy=176;                       // cidade só à direita: nada condutor (vida, água) perto do circuito
    rectM(1,gy,198,198,AREIA); rectM(120,gy,150,gy+6,AGUA);
    boxM(160,gy-52,186,gy-1,2,CONCRETO); rectM(160,gy-36,161,gy-30,VIDRO);
    rectM(140,gy-18,141,gy-1,MADEIRA); rectM(136,gy-25,145,gy-19,MADEIRA);
    soup(164,gy-44,182,gy-4,.35); soup(112,gy-20,132,gy-2,.35);
    rectM(99,122,101,175,ACO);                                     // torre
    discM(100,108,12,ACO); discM(100,108,10,C4); discM(100,108,6,UNAT); discM(100,108,3,PU239); rectM(100,108,100,108,AMBE);
    rectM(100,94,100,97,C4); rectM(97,92,103,93,NICROMO);          // detonador de fio-ponte (nicromo encostado no C-4)
    boxM(18,146,62,176,2,CONCRETO);                                // bunker de controle
    rectM(30,158,34,162,FONTEP); rectM(30,166,34,170,FONTEN);      // bateria
    rectM(25,167,29,168,COBRE); rectM(24,79,25,168,COBRE); rectM(24,79,105,80,COBRE); rectM(104,79,105,93,COBRE); // fio −
    rectM(35,159,38,160,COBRE); rectM(42,159,44,160,COBRE); rectM(43,85,44,160,COBRE); rectM(43,85,96,86,COBRE); rectM(95,85,96,93,COBRE); // fio + (vão em x 39–41 = chave aberta)
    rectM(39,162,41,162,LATA);                                     // marca sob a chave
    toast('☢ Feche a chave no bunker: pinte Cobre no vão do fio (acima da marca)',false);
  },
  nukeGun(){ // tipo canhão (Little Boy): projétil de U-235 disparado contra o alvo; juntos passam da massa crítica
    clearAll(); townM(178);
    rectM(87,40,89,172,ACO); rectM(111,40,113,172,ACO); rectM(84,150,116,175,ACO); rectM(87,150,113,172,UNAT); rectM(100,170,100,170,AMBE); // cano + tamper/refletor de U natural + iniciador
    rectM(90,150,110,169,AR); rectM(90,162,110,169,U235);           // alvo (subcrítico) dentro do refletor
    rectM(90,60,110,66,U235); rectM(90,67,110,67,VIDRO);           // projétil (subcrítico) sobre placa de vidro
    rectM(96,52,104,53,C4); heaterSet(700); rectM(99,51,101,51,RESIST); // carga propelente + detonador
  },
  reactor(){ // reator de piscina: varetas de U-235 moderadas pela água; barras de boro controlam (apague-as para subir a potência)
    clearAll();
    rectM(1,188,198,198,CONCRETO);
    boxM(30,80,170,187,4,CONCRETO); rectM(34,84,166,183,AGUA);
    for(let x=54,k=0;x<150;x+=8,k++){rectM(x,138,x,183,U235); if(k%2===0)rectM(x+4,112,x+4,183,BORO);}
    rectM(100,184,100,184,AMBE);                                    // fonte de partida
    rectM(176,150,196,187,CONCRETO); soup(4,150,26,186,.35);        // vida perto da piscina
  },
  radiation(){ // blindagem: fonte de Co-60 atrás de chumbo, concreto, água e nada; porão com rádio emanando radônio; césio na água
    clearAll();
    rectM(1,184,198,198,CONCRETO);
    rectM(14,12,16,130,COBALTO);
    rectM(1,40,60,41,PAREDE); rectM(1,72,60,73,PAREDE); rectM(1,104,60,105,PAREDE);
    rectM(40,14,43,39,CHUMBO); rectM(40,42,43,71,CONCRETO); boxM(39,74,45,103,1,VIDRO); rectM(40,75,44,102,AGUA);
    soup(70,14,190,38,.3); soup(70,44,190,70,.3); soup(70,76,190,102,.3); soup(70,108,190,130,.3);
    rectM(70,40,198,41,PAREDE); rectM(70,72,198,73,PAREDE); rectM(70,104,198,105,PAREDE);
    boxM(10,146,80,183,3,CONCRETO); rectM(14,177,30,179,RADIO); soup(40,160,76,179,.35); // porão fechado
    rectM(110,160,111,183,VIDRO); rectM(189,160,190,183,VIDRO); rectM(110,182,190,183,VIDRO); rectM(112,164,188,181,AGUA); rectM(115,150,121,159,CESIO); // cápsula de césio caindo na piscina
  }
});

// ================= UI =================
function fmtSv(d){return d>=1?d.toPrecision(3)+' Sv':d>=1e-3?(d*1e3).toPrecision(3)+' mSv':(d*1e6).toPrecision(3)+' µSv';}
function radReadout(i,m){let h='';
  if(radLive||radF[i]>1e-7) h+=kv('dose',fmtSv(radF[i])+'/h'+(radF[i]>RAD_STERILE?' · esteriliza':''));
  if((m===VIDA||m===VIDAM)&&dAcc[i]>1e-4) h+=kv('dose acum.',fmtSv(dAcc[i])+(dAcc[i]>3?' · letal':''));
  if(rA[i]>.01) h+=kv('atividade',rA[i].toFixed(2)+' (contaminado)');
  if(cmp[i]>.01) h+=kv('compressão',(1+cmp[i]).toFixed(2)+'× ρ');
  if(bu[i]>0&&FISS[m]) h+=kv('queima',(bu[i]/BUMAX*100).toPrecision(2)+'%');
  return h;}
function radStats(){if(!radLive&&!nCnt&&!nNukes)return '';
  return '<span>dose máx <b>'+fmtSv(radMax)+'/h</b></span>'+(nCnt?'<span>nêutrons <b>'+nCnt+'</b></span><span>fissões <b>'+fisShow.toFixed(0)+'</b>/passo</span>':'')+(nNukes?'<span>☢ <b>'+nNukes+'</b> · '+lastKt.toFixed(1)+' kt</span>':'');}
