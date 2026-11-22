'use strict';
// ================= ar dinâmico =================
// Euler compressível linearizado numa grade grossa (AC×AC células), escalonada (velocidade nas faces):
//   ∂p'/∂t = −γ·P₀·∇·u + fontes      ∂u/∂t = −∇p'/ρ − (u·∇)u + g·(ρ−ρ₀)/ρ
// p' = sobrepressão dinâmica (Pa) somada à pressão hidrostática/selada de pressureStep.
// Fontes: aquecimento a volume constante (p·ΔT/T), vaporização (expansão do vapor) e explosões.
// Acoplamentos: compressão/expansão adiabática T·(p₂/p₁)^((γ−1)/γ); arrasto em gases, líquidos e pós.
// ⚠️ Assumido: velocidade do som reduzida (C_AIR células grossas/passo) para caber na escala de tempo do jogo.
const AC=2,AW=W/AC,AH=H/AC,AN=AW*AH,C_AIR=.3,KP=1.4*PATM,KV=C_AIR*C_AIR/KP,RHO_AIR=PATM*.029/(RG*293.15);
const AP_MAX=20*PATM,AP_MIN=-.9*PATM,THERM=.15,VAP_SRC=.25*PATM/(AC*AC),BUOY=.003,BLAST_KG=.6*PATM,PGK=.016; // PGK: empurrão −∇p/ρ em grãos/líquidos (células finas/passo por Pa/célula grossa·kg/m³)
const GX_AIR=.4/1.4,AF=new Uint16Array(N); for(let i=0;i<N;i++) AF[i]=(((i/W)|0)>>1)*AW+((i%W)>>1);
const ap=new Float32Array(AN),avx=new Float32Array(AN),avy=new Float32Array(AN),tvx=new Float32Array(AN),tvy=new Float32Array(AN),
  aInj=new Float32Array(AN),aQ=new Float32Array(AN),aRho=new Float32Array(AN),aKv=new Float32Array(AN),aCmp=new Float32Array(AN),
  ucx=new Float32Array(AN),ucy=new Float32Array(AN),agx=new Float32Array(AN),agy=new Float32Array(AN),apS=new Float32Array(AN);
const aBlk=new Uint8Array(AN),aGas=new Uint8Array(AN),aDen=new Uint8Array(AN),aAct=new Uint8Array(AN);
let airOn=true,blastOn=true,apMax=0,avMax=0;
// arrasto por material: gás segue o ar; pós/líquidos/blocos com acoplamento ~1/√ρ e limiar
const cplA=new Float32Array(NM),thrA=new Float32Array(NM),pgA=new Float32Array(NM);
for(let m=0;m<NM;m++){const k=kind[m],r=Math.max(1,rho0[m]);pgA[m]=k===GAS||k===FIX?0:PGK/r;
  if(k===GAS){cplA[m]=1;thrA[m]=0;}
  else if(k===POW){cplA[m]=Math.min(1,50/Math.sqrt(r));thrA[m]=.05;}
  else if(k===LIQ){cplA[m]=Math.min(1,40/Math.sqrt(r));thrA[m]=.08;}
  else if(k===MOV){cplA[m]=Math.min(1,20/Math.sqrt(r));thrA[m]=.15;}}

function airReset(){ap.fill(0);avx.fill(0);avy.fill(0);aInj.fill(0);apS.fill(0);apMax=0;avMax=0;}
// sobrepressão vista por uma célula fina; perto de parede usa a célula grossa do lado do gás (d = direção para fora)
function apAt(j,d){let c=AF[j];if(aBlk[c]){const n=j+d;if(n>=0&&n<N&&!aBlk[AF[n]])c=AF[n];}return ap[c];}
function bil(a,fx,fy){
  if(fx<0)fx=0;else if(fx>AW-1.001)fx=AW-1.001; if(fy<0)fy=0;else if(fy>AH-1.001)fy=AH-1.001;
  const x0=fx|0,y0=fy|0,tx=fx-x0,ty=fy-y0,c=y0*AW+x0;
  return (a[c]*(1-tx)+a[c+1]*tx)*(1-ty)+(a[c+AW]*(1-tx)+a[c+AW+1]*tx)*ty;
}
function airStep(){
  if(!airOn) return;
  const AV_MAX=Math.max(.6,1.05*windF/343*C_AIR); // teto da grade: ≥ explosões (.6) e ≥ força máx. do vento
  // 1) agrega células finas: sólidos bloqueiam; fonte térmica Σ ΔT/T dos gases (T−tOld: calor desde o último passo, sem transporte)
  aBlk.fill(0);aGas.fill(0);aDen.fill(0);aRho.fill(0);aQ.fill(0);
  for(let i=0;i<N;i++){
    const m=mat[i],k=kind[m],c=AF[i];
    if(k===FIX){aBlk[c]=1;continue;}
    if(k===GAS){aGas[c]++;aRho[c]+=gd[i];const t0=tOld[i];if(t0>1){let r=(T[i]-t0)/t0;if(r>3)r=3;else if(r<-.5)r=-.5;aQ[c]+=r;}}
    else aDen[c]++;
  }
  const qs=PATM*THERM/(AC*AC);
  for(let c=0;c<AN;c++){
    const g=aGas[c],rho=g?aRho[c]/g:RHO_AIR; aRho[c]=rho;
    let kv=RHO_AIR/rho; kv=kv>4?4:kv<.25?.25:kv; aKv[c]=KV*kv;
    if(aBlk[c]){aInj[c]=0;continue;}
    let p=ap[c]+aInj[c]+qs*aQ[c]; aInj[c]=0;
    ap[c]=p>AP_MAX?AP_MAX:p<AP_MIN?AP_MIN:p;
  }
  // 2) momento: gradiente de pressão + empuxo; arrasto de líquidos/pós na célula
  const gb=BUOY*gm;
  for(let cy=0;cy<AH;cy++){const row=cy*AW;
    for(let cx=0;cx<AW;cx++){const c=row+cx;
      if(aBlk[c]){avx[c]=0;avy[c]=0;continue;}
      const dr=1-.05*aDen[c];
      if(cx<AW-1&&!aBlk[c+1]){let u=avx[c]*dr*.998-(ap[c+1]-ap[c])*.5*(aKv[c]+aKv[c+1]);
        if(gAny&&aGas[c]&&aGas[c+1]){const rf=.5*(aRho[c]+aRho[c+1]);let b=(rf-RHO_AIR)/rf;b=b>1?1:b<-4?-4:b;u+=BUOY*.5*(gCx[c]+gCx[c+1])*b;} // empuxo na direção de g local (grav.js)
        avx[c]=u>AV_MAX?AV_MAX:u<-AV_MAX?-AV_MAX:u;} else avx[c]=0;
      if(cy<AH-1&&!aBlk[c+AW]){const n=c+AW;let v=avy[c]*dr*.998-(ap[n]-ap[c])*.5*(aKv[c]+aKv[n]);
        if((gb>0||gAny)&&aGas[c]&&aGas[n]){const rf=.5*(aRho[c]+aRho[n]);let b=(rf-RHO_AIR)/rf;b=b>1?1:b<-4?-4:b;v+=(gb+BUOY*.5*(gCy[c]+gCy[n]))*b;}
        avy[c]=v>AV_MAX?AV_MAX:v<-AV_MAX?-AV_MAX:v;} else avy[c]=0;
    }}
  // 3) advecção semi-lagrangiana da velocidade (termo (u·∇)u: vórtices e jatos se transportam)
  for(let cy=0;cy<AH;cy++){const row=cy*AW;
    for(let cx=0;cx<AW;cx++){const c=row+cx;
      let u=avx[c];
      if(u!==0){const v=.25*(avy[c]+(cx<AW-1?avy[c+1]:0)+(cy>0?avy[c-AW]+(cx<AW-1?avy[c-AW+1]:0):0));
        tvx[c]=(u<1e-3&&u>-1e-3&&v<1e-3&&v>-1e-3)?u:bil(avx,cx-u,cy-v);} else tvx[c]=0;
      const v=avy[c];
      if(v!==0){const u2=.25*(avx[c]+(cx>0?avx[c-1]:0)+(cy<AH-1?avx[c+AW]+(cx>0?avx[c+AW-1]:0):0));
        tvy[c]=(v<1e-3&&v>-1e-3&&u2<1e-3&&u2>-1e-3)?v:bil(avy,cx-u2,cy-v);} else tvy[c]=0;
    }}
  for(let c=0;c<AN;c++){avx[c]=aBlk[c]?0:tvx[c];avy[c]=aBlk[c]?0:tvy[c];}
  // 4) continuidade: p' −= γP₀·∇·u (compressão isentrópica); viscosidade artificial leve
  let pm=0,vm=0;
  for(let cy=0;cy<AH;cy++){const row=cy*AW;
    for(let cx=0;cx<AW;cx++){const c=row+cx;
      if(aBlk[c]){aCmp[c]=0;continue;}
      const div=avx[c]-(cx>0?avx[c-1]:0)+avy[c]-(cy>0?avy[c-AW]:0);
      const p0=ap[c];let p=p0*.998-KP*div;
      p=p>AP_MAX?AP_MAX:p<AP_MIN?AP_MIN:p; apS[c]=p;
      aQ[c]=p0; // pressão antes da divergência (já com fontes isocóricas)
    }}
  const open=!lid;
  for(let cy=0;cy<AH;cy++){const row=cy*AW;
    for(let cx=0;cx<AW;cx++){const c=row+cx;
      if(aBlk[c]) continue;
      if(open&&cy===0){ap[c]=0;aCmp[c]=Math.pow(PATM/(PATM+aQ[c]),GX_AIR)-1;continue;} // topo aberto: reservatório atmosférico
      let s=0,n=0;
      if(cx>0&&!aBlk[c-1]){s+=apS[c-1];n++;} if(cx<AW-1&&!aBlk[c+1]){s+=apS[c+1];n++;}
      if(cy>0&&!aBlk[c-AW]){s+=apS[c-AW];n++;} if(cy<AH-1&&!aBlk[c+AW]){s+=apS[c+AW];n++;}
      const p=n?apS[c]+.04*(s/n-apS[c]):apS[c]; ap[c]=p;
      // fator adiabático EXATO a partir da pressão final armazenada: o produto telescopa (ciclo fechado volta à mesma T);
      // aproximação linear 1+γ'r deixava +γ'r² por ciclo → bomba de calor numérica sob vento/ciclone
      aCmp[c]=Math.pow((PATM+p)/(PATM+aQ[c]),GX_AIR)-1;
      const a=p<0?-p:p; if(a>pm)pm=a;
    }}
  // células bloqueadas herdam a média dos vizinhos livres (leituras perto de paredes)
  for(let c=0;c<AN;c++){if(!aBlk[c])continue;const cx=c%AW;let s=0,n=0;
    if(cx>0&&!aBlk[c-1]){s+=ap[c-1];n++;} if(cx<AW-1&&!aBlk[c+1]){s+=ap[c+1];n++;}
    if(c>=AW&&!aBlk[c-AW]){s+=ap[c-AW];n++;} if(c<AN-AW&&!aBlk[c+AW]){s+=ap[c+AW];n++;}
    ap[c]=n?s/n:0;}
  // 5) velocidade no centro (células finas/passo) e máscara de atividade
  for(let c=0;c<AN;c++){
    if(aBlk[c]){ucx[c]=0;ucy[c]=0;aAct[c]=0;continue;}
    const cx=c%AW,u=(avx[c]+(cx>0?avx[c-1]:0))*.5*AC,v=(avy[c]+(c>=AW?avy[c-AW]:0))*.5*AC;
    ucx[c]=u;ucy[c]=v;const s=(u<0?-u:u)+(v<0?-v:v);if(s>vm)vm=s;
    const cy=(c/AW)|0,pc=ap[c],pl=cx>0&&!aBlk[c-1]?ap[c-1]:pc,pr=cx<AW-1&&!aBlk[c+1]?ap[c+1]:pc,
      pu=cy>0&&!aBlk[c-AW]?ap[c-AW]:pc,pd=cy<AH-1&&!aBlk[c+AW]?ap[c+AW]:pc,gx=(pl-pr)*.5,gy=(pu-pd)*.5;
    agx[c]=gx;agy[c]=gy; // −∇p (Pa/célula grossa)
    const r=aCmp[c];aAct[c]=(s>.04||r>3e-4||r<-3e-4||gx>3e3||gx<-3e3||gy>3e3||gy<-3e3)?1:0;
  }
  apMax=pm;avMax=vm;
  // 6) acoplamento com as partículas: adiabática nos gases + transporte/arrasto
  stamp=stamp>=255?1:stamp+1;
  for(let y=1;y<H-1;y++){
    const row=y*W,ltr=rnd()<.5;
    for(let xx=1;xx<W-1;xx++){
      const x=ltr?xx:W-1-xx,i=row+x;let c=AF[i];
      // superfície de líquido/pó: célula grossa sem gás → sente o vento da célula de cima (camada limite)
      if(!aGas[c]&&c>=AW&&aGas[c-AW]&&kind[mat[i-W]]===GAS) c-=AW;
      if(!aAct[c]) continue;
      const m=mat[i],k=kind[m];
      if(k===FIX) continue;
      if(k===GAS){const r=aCmp[c];if(r>1e-7||r<-1e-7){ // aCmp = (p₂/p₁)^((γ−1)/γ)−1 do ar; outros γ: reescala o expoente
        let t=T[i]*(gexp[m]===GX_AIR?1+r:Math.pow(1+r,gexp[m]/GX_AIR));T[i]=t>6000?6000:t<2?2:t;}}
      if(mv[i]===stamp) continue;
      const f=cplA[m]; if(!f) continue;
      let u=ucx[c],v=ucy[c];
      if(k!==GAS){ // arrasto dv/dt = k·(u_ar − v) + empurrão −∇p/ρ; gravidade e voo em flyStep
        const pg=pgA[m];u+=agx[c]*pg;v+=agy[c]*pg;
        if((u<0?-u:u)+(v<0?-v:v)<thrA[m]) continue;
        const a=f*.6;pvx[i]+=(u-pvx[i])*a;pvy[i]+=(v-pvy[i])*a;flyOn=true;
        // grão/gota apoiado: arrasto também rola/escorrega direto (prob. ∝ excesso sobre o limiar), subindo na diagonal se preciso
        const au=u<0?-u:u; if(au<=thrA[m]||mv[i]===stamp) continue;
        if(rnd()<(au-thrA[m])*f){const d=u>0?1:-1,j=i+d;
          if(kind[mat[j]]===GAS&&mv[j]!==stamp) swap(i,j);
          else if(kind[mat[j-W]]===GAS&&mv[j-W]!==stamp&&rnd()<.6) swap(i,j-W);}
        continue;}
      const ax=u<0?-u:u,ay=v<0?-v:v,s=ax+ay; if(s<.02) continue;
      // um sorteio: escolhe o eixo ∝ |componente| e reaproveita a fração como 2º número uniforme
      const rs=rnd()*s;let p,d,r2;
      if(rs<ax){p=ax;d=u>0?1:-1;r2=rs/ax;}else{p=ay;d=v>0?W:-W;r2=(rs-ax)/ay;}
      if(k!==GAS){p=(p-thrA[m])*f;if(p<=0)continue;}
      if(r2>=p) continue;
      const j=i+d,mj=mat[j];
      if(kind[mj]!==GAS||mv[j]===stamp) continue;
      if(k===GAS&&kind[mj]===GAS){ // troca é simétrica: só a célula "diferente" é advectada (senão o ar de trás a devolve)
        if(m===AR){if(mj!==AR)continue;let a=T[i]-293.15,b=T[j]-293.15;a=a<0?-a:a;b=b<0?-b:b;if(a<b+.5)continue;}
        else if(mj===m)continue;}
      swap(i,j);
    }
  }
}

// ================= ferramentas de ar =================
// Vento: arrastar empurra o ar na direção do gesto. Ciclone: vórtice de Rankine (anti-horário; Shift = horário).
// Vento = jato sustentado (ventilador): enquanto o botão está pressionado, relaxa o ar do pincel para
// a velocidade-alvo na direção do último arrasto (|alvo| cresce com a rapidez do gesto, até a força máx.)
let windDX=0,windDY=0;
function airWindDir(dx,dy){windDX=windDX*.6+dx*.4;windDY=windDY*.6+dy*.4;}
function airWind(x,y){
  const d=Math.hypot(windDX,windDY); if(!airOn||d<.05) return;
  const WM=windF/343*C_AIR,sp=WM*Math.min(1,.35+d/3),tu=windDX/d*sp,tv=windDY/d*sp;
  const R=Math.max(1.5,(brush+2)/AC),X=(x+.5)/AC,Y=(y+.5)/AC;
  const x0=Math.max(0,X-R|0),x1=Math.min(AW-1,X+R+1|0),y0=Math.max(0,Y-R|0),y1=Math.min(AH-1,Y+R+1|0);
  for(let cy=y0;cy<=y1;cy++)for(let cx=x0;cx<=x1;cx++){const c=cy*AW+cx;if(aBlk[c])continue;
    const r=Math.hypot(cx+.5-X,cy+.5-Y);if(r>R)continue;const w=.5*(1-r/R*.5);
    avx[c]+=(tu-avx[c])*w;avy[c]+=(tv-avy[c])*w;}
}
function airCyclone(x,y,cw){
  if(!airOn) return;
  const R=Math.max(2,(brush+2)/AC),RR=R*1.6,X=(x+.5)/AC,Y=(y+.5)/AC,S=.83*windF/343*C_AIR*(cw?-1:1);
  const x0=Math.max(0,X-RR|0),x1=Math.min(AW-1,X+RR+1|0),y0=Math.max(0,Y-RR|0),y1=Math.min(AH-1,Y+RR+1|0);
  const vt=r=>r<R?S*r/R:S*R/r*(1-(r-R)/(RR-R)); // Rankine, cauda zerando em RR
  for(let cy=y0;cy<=y1;cy++)for(let cx=x0;cx<=x1;cx++){const c=cy*AW+cx;if(aBlk[c])continue;
    let rx=cx+1-X,ry=cy+.5-Y,r=Math.hypot(rx,ry); // face x em (cx+1, cy+.5)
    if(r>.05&&r<RR){const t=vt(r)*ry/r;avx[c]+=(t-avx[c])*.25;}
    rx=cx+.5-X;ry=cy+1-Y;r=Math.hypot(rx,ry);   // face y em (cx+.5, cy+1)
    if(r>.05&&r<RR){const t=-vt(r)*rx/r;avy[c]+=(t-avy[c])*.25;}
  }
}

// ================= explosivos =================
// Carga conectada detona inteira no mesmo passo (velocidade de detonação ~8 km/s ≫ 1 célula/passo).
// Pico de sobrepressão da onda de choque: Mills (1987), ΔP[kPa] = 1772/Z³ − 114/Z² + 108/Z, Z = R/W^(1/3) (m/kg^(1/3)),
// W em kg de TNT equivalente; sólido rompe se a pressão refletida (2·ΔP) supera a resistência.
// A frente se propaga por busca em largura (contorna obstáculos, para em paredes que resistem).
const exW=new Float32Array(NM),exTd=new Float32Array(NM),exSh=new Float32Array(NM);
for(let m=0;m<NM;m++){const e=D[m].ex;if(e){exW[m]=e[0];exTd[m]=e[1];exSh[m]=e[2];}}
function mills(Z){if(Z<.1)Z=.1;const iz=1/Z;return (1772*iz*iz*iz-114*iz*iz+108*iz)*1e3;}
const ZCUT=5.2; // 2·ΔP(Z) < ~60 kPa: abaixo do material mais fraco
const exQ=new Int32Array(N),exO=new Int32Array(N),exSeen=new Uint16Array(N),detL=new Int32Array(8192),detWq=new Float32Array(8192);
let exStamp=0,detN=0,nBlasts=0;
function detonate(i,w){if(detN<detL.length){detL[detN]=i;detWq[detN]=w;detN++;}}
function explodeStep(){
  if(!detN) return;
  exStamp=exStamp>=65535?1:exStamp+1;const S=exStamp;
  let qh=0,qt=0,Wt=0;
  for(let n=0;n<detN;n++){const i=detL[n];if(exSeen[i]===S)continue;exSeen[i]=S;exQ[qt++]=i;Wt+=detWq[n];}
  detN=0;
  while(qh<qt){const i=exQ[qh++];
    for(let n=0;n<4;n++){const j=n===0?i-1:n===1?i+1:n===2?i-W:i+W;
      if(exSeen[j]!==S&&exW[mat[j]]>0){exSeen[j]=S;exQ[qt++]=j;}}}
  const ns=qt,inj=airOn&&blastOn;
  for(let a=0;a<ns;a++){const i=exQ[a],w=exW[mat[i]]||.4;Wt+=exW[mat[i]];
    setCell(i,FOGO,2800+rnd()*800);life[i]=25+(rnd()*40|0);exO[i]=i;
    if(inj) aInj[AF[i]]+=BLAST_KG*w;}
  nBlasts++;
  if(!blastOn||Wt<=0) return;
  const cw=Math.cbrt(Wt),Rmax=ZCUT*cw/DX,Rmax2=Rmax*Rmax,Rf=3*Math.cbrt(ns)+2,Rf2=Rf*Rf;
  qh=0;
  while(qh<qt){
    const i=exQ[qh++],o=exO[i],ox=o%W,oy=(o/W)|0;
    for(let n=0;n<4;n++){
      const j=n===0?i-1:n===1?i+1:n===2?i-W:i+W;
      if(exSeen[j]===S) continue; exSeen[j]=S;
      const x=j%W,y=(j/W)|0; if(x<1||x>W-2||y<1||y>H-2) continue;
      const dx=x-ox,dy=y-oy,d2=dx*dx+dy*dy; if(d2>Rmax2) continue;
      const m=mat[j];
      if(exW[m]>0){detonate(j,0);continue;} // detonação por simpatia no passo seguinte
      if(kind[m]===FIX){
        const s=strA[m]; if(!s) continue;
        if(2*mills((Math.sqrt(d2)*DX+.05)/cw)<s) continue;
        breakFixed(j,m);
      }
      exO[j]=o; exQ[qt++]=j;
      if(FISS[mat[j]]){const c=Math.min(1.6,2*mills((Math.sqrt(d2)*DX+.05)/cw)/2e9);if(c>cmp[j])cmp[j]=c;} // implosão: choque comprime o núcleo (ρ↑ → massa crítica ∝ 1/ρ²)
      {const mm=mat[j],k=kind[mm]; // impulso nos detritos: Δv ∝ ΔP/ρ, radial
        if(k!==GAS&&k!==FIX){const d=Math.sqrt(d2)+.5;let v=mills((d*DX+.05)/cw)/(Math.max(rho0[mm],20)*FLYK);
          if(v>.08){if(v>FLY_MAX)v=FLY_MAX;pvx[j]+=v*dx/d;pvy[j]+=v*dy/d-v*.25;flyOn=true;}}}
      if(d2<Rf2){const f=1-Math.sqrt(d2)/Rf,mm=mat[j],k=kind[mm];
        if(k===GAS){setCell(j,FOGO,1300+1800*f);life[j]=8+(rnd()*30|0);}
        else if(ignA[mm]>0) ignite(j);
        else T[j]+=150*f;}
    }
  }
}

// ================= detritos balísticos =================
// grãos/líquidos com momento próprio: voam em linha reta por células de gás, gravidade + arrasto, ricocheteiam/perdem energia ao colidir
const FLYK=160,FLY_MAX=4;let flyOn=false;
function flyStep(){
  if(!flyOn) return;
  stamp=stamp>=255?1:stamp+1;
  const g=.04*gm;let any=false;
  for(let y=H-2;y>=1;y--){const row=y*W;
    for(let x=1;x<W-1;x++){let i=row+x;
      let vx=pvx[i],vy=pvy[i]; if(vx===0&&vy===0) continue;
      const m=mat[i],k=kind[m];
      if(k===GAS||k===FIX){pvx[i]=0;pvy[i]=0;continue;}
      if(mv[i]===stamp) continue;
      vx+=gfx[i];vy+=g+gfy[i];vx*=.97;vy*=.97; // + gravidade local (grav.js)
      const s=Math.abs(vx)>Math.abs(vy)?Math.abs(vx):Math.abs(vy);
      if(s<.08){pvx[i]=0;pvy[i]=0;continue;}
      any=true;pvx[i]=vx;pvy[i]=vy;
      let n=s|0; if(rnd()<s-n) n++;
      const sx=vx/s,sy=vy/s;let fx=x+.5,fy=y+.5;
      for(let q=0;q<n;q++){
        fx+=sx;fy+=sy;const j=((fy|0)*W)+(fx|0);
        if(j===i) continue;
        if(fx<1||fx>=W-1||fy<1||fy>=H-1||kind[mat[j]]!==GAS||mv[j]===stamp){
          // arraste de superfície: bloqueado na horizontal → sobe na diagonal se houver gás (ondas, spray, monte a jusante)
          const ju=j-W;if(fx>=1&&fx<W-1&&fy>=2&&(vx<0?-vx:vx)>(vy<0?-vy:vy)&&kind[mat[ju]]===GAS&&mv[ju]!==stamp){swap(i,ju);i=ju;fy-=1;pvx[i]*=.85;continue;}
          pvx[i]*=.45;pvy[i]=pvy[i]>0&&(fy|0)!==y?-pvy[i]*.45:pvy[i]*.3;break;} // ricocheteia no chão
        swap(i,j);i=j;
      }
    }}
  flyOn=any;
}

// ================= vista do fluxo de ar (sobreposição) =================
const TRN=1400,trX=new Float32Array(TRN),trY=new Float32Array(TRN),trL=new Int16Array(TRN);
function drawAir(oc,sx,sy,dpr,fn){
  if(!airOn) return;
  const k=Math.max(.5,steps)*1.6;oc.lineCap='round';
  const B=[[],[],[]];
  for(let p=0;p<TRN;p++){
    if(trL[p]<=0){trX[p]=1+rnd()*(W-2);trY[p]=1+rnd()*(H-2);trL[p]=40+(rnd()*60|0);}
    const c=AF[(trY[p]|0)*W+(trX[p]|0)];
    if(aBlk[c]){trL[p]=0;continue;}
    const u=ucx[c],v=ucy[c],s=Math.hypot(u,v);trL[p]--;
    if(s<.01){if(rnd()<.05)trL[p]=0;continue;}
    const x0=trX[p],y0=trY[p];trX[p]+=u*k;trY[p]+=v*k;
    if(trX[p]<1||trX[p]>W-1||trY[p]<1||trY[p]>H-1){trL[p]=0;continue;}
    const b=s<.08?0:s<.35?1:2,e=Math.min(4,2+s*6);
    B[b].push(x0*sx,y0*sy,(x0+u*k*e)*sx,(y0+v*k*e)*sy);
  }
  const col=['rgba(150,200,255,.35)','rgba(190,235,255,.6)','rgba(255,255,255,.85)'];
  oc.lineWidth=1.3*dpr;
  for(let b=0;b<3;b++){const L=B[b];if(!L.length)continue;oc.strokeStyle=col[b];oc.beginPath();
    for(let q=0;q<L.length;q+=4){oc.moveTo(L[q],L[q+1]);oc.lineTo(L[q+2],L[q+3]);}oc.stroke();}
  oc.lineCap='butt';
}
