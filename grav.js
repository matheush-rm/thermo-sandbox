'use strict';
// ================= gravidade local (Newton) e buracos negros (Schwarzschild) =================
// Fontes pontuais de massa M (ferramenta Gravidade: segurar aumenta M exponencialmente; Shift diminui).
// Campo real 3D: g = Σ G·M·r̂/r² (amolecido em ε = ½ célula) + gravidade uniforme do painel (gm·g₀, para baixo).
// Efeitos: queda/escorregamento/espalhamento na direção local de g (autômato generalizado em 8 direções),
// voo balístico (detritos) sob g local → órbitas/estilingue; pressão hidrostática ∫ρ g·dl (atmosfera e oceano
// em volta da massa, gás comprimido e aquecido adiabaticamente); empuxo do ar dinâmico na direção de g;
// forças de maré ΔF ≈ 2GM/r³ rasgam sólidos (espaguetificação) quando ρ·(2GM/r³)·L² > resistência.
// Singularidade: quando o raio de Schwarzschild r_s = 2GM/c² alcança ½ célula (M ≳ 3,4·10²⁵ kg) a fonte vira buraco
// negro: absorve tudo que cruza o horizonte (massa somada a M), suga o ar, aquece o disco de acreção (r < 3 r_s, ISCO),
// mostra sombra/anel de fótons (2,6 r_s) e lente gravitacional; dois buracos negros que se tocam fundem
// perdendo ~5% da massa em ondas gravitacionais.
// A autogravidade da matéria da câmara é desprezível (1 L de aço = 7,85 kg → 5·10⁻⁸ m/s² a 10 cm) e não é calculada.
// ⚠️ Assumido: hidrostática limitada a 100 g (gás ideal com escala de altura RT/Mg ≪ 1 célula divergiria);
// empuxo do ar limitado a 30 g; aquecimento do disco limitado (+≤400 K/passo, ≤ 5500 K) — a energia real
// (~6% de mc²) vaporizaria tudo instantaneamente; maré avaliada no comprimento de uma célula (L = 10 cm);
// sombra/lente desenhadas como vista de um observador distante (o plano da câmara é uma fatia 2D).
const GN=6.674e-11,CLUZ=2.998e8,HBAR=1.0546e-34,KB=1.380649e-23,GSOFT2=(.5*DX)*(.5*DX),GH_MAX=100*G0,GAIR_MAX=30,
  GM0=1e9,GMMAX=2e27,GRATE=1.07,GMAXS=8,BH_RS=.5,GW_LOSS=.05,BH_SINK=.08*PATM;
const gS=[]; // fontes {x,y (células, contínuo), M (kg), bh, abs (kg absorvidos)}
const gTx=new Float32Array(N),gTy=new Float32Array(N),gTa=new Float32Array(N),gFr=new Float32Array(N),gO=new Uint8Array(N),gLoc=new Uint8Array(N),
  gfx=new Float32Array(N),gfy=new Float32Array(N),gHx=new Float32Array(N),gHy=new Float32Array(N),gTid=new Float32Array(N),gPhi=new Float32Array(N),
  gCx=new Float32Array(AN),gCy=new Float32Array(AN);
const GDIR=[1,1+W,W,W-1,-1,-1-W,-W,1-W]; // octantes (y para baixo): 0 →, 2 ↓, 4 ←, 6 ↑
let gAny=false,gBH=false,gDirty=false,gLastGm=-1,gTidMax=0;
const gWaves=[]; // ondas gravitacionais (fusões) para a sobreposição
let TID_MIN=1e30; // menor maré (s⁻²) que rasga algum sólido: σ/(ρ·L²)
for(let m=0;m<NM;m++)if(kind[m]===FIX&&strA[m]>0&&crushA[m]&&rho0[m]>0)TID_MIN=Math.min(TID_MIN,strA[m]/(rho0[m]*DX*DX));

const rsOf=M=>2*GN*M/(CLUZ*CLUZ); // m
function gravReset(){gS.length=0;gWaves.length=0;gDirty=true;gravField();}
function gravGet(){return gS.map(s=>({x:s.x,y:s.y,M:s.M,abs:s.abs||0}));}
function gravSet(a){gS.length=0;if(Array.isArray(a))for(const s of a)if(s&&isFinite(s.M)&&s.M>0)gS.push({x:+s.x,y:+s.y,M:+s.M,abs:+s.abs||0,bh:false});gDirty=true;gravField();}

function gravField(){
  gDirty=false;gLastGm=gm;gAny=gS.length>0;gBH=false;
  for(const s of gS){s.bh=rsOf(s.M)/DX>=BH_RS;if(s.bh)gBH=true;}
  if(!gAny){gLoc.fill(0);gfx.fill(0);gfy.fill(0);gHx.fill(0);gHy.fill(0);gTid.fill(0);gPhi.fill(0);gCx.fill(0);gCy.fill(0);gTidMax=0;return;}
  const gu=gm*G0,thr=.02*G0*Math.max(gm,.1),ns=gS.length,SX=new Float64Array(ns),SY=new Float64Array(ns),SM=new Float64Array(ns);
  for(let s=0;s<ns;s++){SX[s]=gS[s].x;SY[s]=gS[s].y;SM[s]=GN*gS[s].M;}
  let tm=0;
  for(let y=0;y<H;y++)for(let x=0;x<W;x++){const i=y*W+x;
    let ax=0,ay=0,td=0,ph=0;
    for(let s=0;s<ns;s++){const dx=(SX[s]-x-.5)*DX,dy=(SY[s]-y-.5)*DX,r2=dx*dx+dy*dy+GSOFT2,r=Math.sqrt(r2),a=SM[s]/r2;
      ax+=a*dx/r;ay+=a*dy/r;td+=2*SM[s]/(r2*r);ph-=SM[s]/r;}
    const tx=ax,ty=ay+gu,ta=Math.hypot(tx,ty),as=Math.hypot(ax,ay);
    gTx[i]=tx;gTy[i]=ty;gTa[i]=ta;gTid[i]=td;gPhi[i]=ph;if(td>tm)tm=td;
    gLoc[i]=as>thr?1:0;
    let ang=Math.atan2(ty,tx)/(Math.PI/4);if(ang<0)ang+=8;const o=Math.floor(ang);gO[i]=o&7;gFr[i]=ang-o;
    let fx=.04*ax/G0,fy=.04*ay/G0;fx=fx>FLY_MAX?FLY_MAX:fx<-FLY_MAX?-FLY_MAX:fx;fy=fy>FLY_MAX?FLY_MAX:fy<-FLY_MAX?-FLY_MAX:fy;gfx[i]=fx;gfy[i]=fy;
    const ch=as>GH_MAX?GH_MAX/as:1;gHx[i]=.5*DX*ax*ch;gHy[i]=.5*DX*ay*ch;
    if(!(x&1)&&!(y&1)){const c=(y>>1)*AW+(x>>1);let cx=ax/G0,cy=ay/G0;const cm=Math.hypot(cx,cy);if(cm>GAIR_MAX){cx*=GAIR_MAX/cm;cy*=GAIR_MAX/cm;}gCx[c]=cx;gCy[c]=cy;}
  }
  gTidMax=tm;
}

// ---------- movimento sob g local (substitui moveStep nas células com gLoc) ----------
function gravMove(i,x,y){
  const m=mat[i],k=kind[m];if(k===FIX) return true;
  const ga=gTa[i];if(ga<1e-4) return true; // imponderável
  const gf=ga>=G0?1:ga/G0,o0=gO[i],o=(rnd()<gFr[i]?o0+1:o0)&7,j=i+GDIR[o];
  if(j<0||j>=N) return true;
  const mj=mat[j],kj=kind[mj];
  if(k===GAS){
    if(kj===GAS&&mv[j]!==stamp&&!(mj===m&&T[j]-T[i]<2)){const ri=gd[i],rj=gd[j];if(ri>rj&&rnd()<(ri-rj)/ri*40*gf){swap(i,j);return true;}}
    if(m!==AR){const jj=i+GDIR[(o+(rnd()<.5?2:6))&7],mm=mat[jj];if(jj>=0&&jj<N&&mv[jj]!==stamp&&kind[mm]===GAS&&mm!==m&&rnd()<.35)swap(i,jj);}
    return true;
  }
  if(y===0) return true;
  if(pvx[i]!==0||pvy[i]!==0){flyOn=true;return true;} // em voo: flyStep cuida (gravidade local inclusa)
  const ri=densOf(i);
  if(mv[j]!==stamp&&j>=W){
    if(kj===GAS){
      const a=.04*ga/G0;
      if(a>=.04){pvx[i]+=gTx[i]/ga*a;pvy[i]+=gTy[i]/ga*a;flyOn=true;return true;} // queda livre balística
      if(ri>gd[j]&&rnd()<gf*(k===LIQ?mobOf(i,m):1)){swap(i,j);return true;}
    }else if(kj===LIQ||((kj===POW||kj===MOV)&&k===LIQ)){
      const rj=densOf(j);
      if(ri>rj){let p=(ri-rj)/ri*40;if(p>1)p=1;if(k===LIQ)p*=mobOf(i,m);if(kj===LIQ)p*=mobOf(j,mj);if(rnd()<p*gf){swap(i,j);return true;}}
    }
  }
  if(k===MOV) return true;
  // escorrega nas diagonais da direção de g (pilhas cônicas / esféricas em volta da massa)
  const sd=rnd()<.5?1:-1;
  for(let s=0;s<2;s++){const e=s?-sd:sd,od=(o+e+8)&7,jd=i+GDIR[od];if(jd<W||jd>=N||mv[jd]===stamp)continue;
    const md=mat[jd],kd=kind[md];if(kd!==GAS&&kd!==LIQ)continue;
    if((od&1)&&kj===FIX&&kind[mat[i+GDIR[(o+2*e+8)&7]]]===FIX)continue; // não vaza por quinas
    if(ri<=densOf(jd))continue;
    let p=k===POW?.65:mobOf(i,m);if(kd===LIQ)p*=mobOf(jd,md);
    if(rnd()<p*gf){swap(i,jd);return true;}}
  // líquidos espalham perpendicularmente a g
  if(k===LIQ&&rnd()<mobOf(i,m)*gf){
    const D2=GDIR[(o+(rnd()<.5?2:6))&7],dg=GDIR[o],dis=dispA[m];let last=-1;
    for(let s=1;s<=dis;s++){const jj=i+D2*s;if(jj<W||jj>=N)break;const mm=mat[jj],kk=kind[mm];
      if(mv[jj]===stamp)break;
      if(kk===GAS){last=jj;const jb=jj+dg;if(jb>=0&&jb<N&&kind[mat[jb]]===GAS)break;continue;}
      if(kk===LIQ&&mm!==m&&densOf(jj)<ri-1){last=jj;continue;}
      break;}
    if(last>=0)swap(i,last);
  }
  return true;
}

// ---------- passo: campo, fusões, horizonte, disco de acreção, maré ----------
function gravStep(){
  if(gm!==gLastGm&&gAny) gDirty=true;
  if(gDirty) gravField();
  if(!gAny) return;
  // fusões: horizontes que se tocam (ou fontes na mesma célula)
  for(let a=0;a<gS.length;a++)for(let b=a+1;b<gS.length;b++){const A=gS[a],B=gS[b];
    const d=Math.hypot(A.x-B.x,A.y-B.y),ra=rsOf(A.M)/DX,rb=rsOf(B.M)/DX;
    if(d<Math.max(1,ra+rb)){const M=A.M+B.M,both=A.bh&&B.bh;
      A.x=(A.x*A.M+B.x*B.M)/M;A.y=(A.y*A.M+B.y*B.M)/M;A.M=both?M*(1-GW_LOSS):M;A.abs=(A.abs||0)+(B.abs||0);
      if(both)gWaves.push({x:A.x,y:A.y,r:ra+rb,a:1,E:M*GW_LOSS*CLUZ*CLUZ});
      gS.splice(b,1);b--;gDirty=true;}}
  if(gDirty) gravField();
  // maré: sólidos fixos rasgam quando ρ·(2GM/r³)·L² supera a resistência
  if(gTidMax>TID_MIN){
    const L2=DX*DX;
    for(let y=1;y<H-1;y++){const row=y*W;for(let x=1;x<W-1;x++){const i=row+x,m=mat[i];
      if(kind[m]!==FIX||!strA[m]||!crushA[m])continue;
      if(rho0[m]*gTid[i]*L2>strA[m]&&rnd()<.2)breakFixed(i,m);}}
  }
  if(!gBH) return;
  for(const s of gS){if(!s.bh)continue;
    const R=rsOf(s.M)/DX,Ra=Math.max(3,3*R),Rb=Math.ceil(Ra)+1,R2=Math.max(R,.75)**2,Ra2=Ra*Ra;
    const x0=Math.max(1,(s.x-Rb)|0),x1=Math.min(W-2,(s.x+Rb)|0),y0=Math.max(1,(s.y-Rb)|0),y1=Math.min(H-2,(s.y+Rb)|0);
    let dm=0;
    for(let y=y0;y<=y1;y++)for(let x=x0;x<=x1;x++){const i=y*W+x,m=mat[i];if(m===PAREDE)continue;
      const dx=x+.5-s.x,dy=y+.5-s.y,d2=dx*dx+dy*dy;if(d2>Ra2)continue;
      const k=kind[m];
      if(d2<=R2){ // cruzou o horizonte: massa-energia vai para M (inclui energia térmica E = c·m·T, desprezível)
        const ms=(k===GAS?gd[i]:densOf(i))*DX*DX*DX;dm+=ms;
        if(k===GAS&&airOn)aInj[AF[i]]-=BH_SINK*(gd[i]/RHO_AIR)/(AC*AC); // gás removido: sumidouro de pressão (ar é sugado)
        setCell(i,AR,293.15);
        continue;}
      // disco de acreção (dentro da ISCO, 3 r_s): atrito viscoso/compressão aquecem a matéria
      const d=Math.sqrt(d2),f=Ra/d-1;let t=T[i]+(f*120>400?400:f*120);T[i]=t>5500?5500:t;
    }
    s.M+=dm;s.abs=(s.abs||0)+dm;
  }
}

// ---------- ferramenta ----------
function gravPaint(cx,cy,shrink){
  const X=cx+.5,Y=cy+.5,R=Math.max(2,brush);let best=-1,bd=1e9;
  for(let s=0;s<gS.length;s++){const d=Math.hypot(gS[s].x-X,gS[s].y-Y);if(d<bd){bd=d;best=s;}}
  if(best>=0&&bd<=Math.max(R,rsOf(gS[best].M)/DX+1)){const s=gS[best];
    if(shrink){s.M/=GRATE;if(s.M<GM0*.5)gS.splice(best,1);} // ⚠️ ferramenta: na física real um buraco negro só perde massa por Hawking
    else s.M=Math.min(GMMAX,s.M*GRATE);
    gDirty=true;return;}
  if(shrink||gS.length>=GMAXS||cx<1||cx>W-2||cy<1||cy>H-2) return;
  gS.push({x:X,y:Y,M:GM0,abs:0,bh:false});gDirty=true;
}

// ---------- inspetor ----------
function sciM(v,u){if(!(v>0))return '0 '+u;const e=Math.floor(Math.log10(v)),mnt=v/Math.pow(10,e);
  const SUP='⁰¹²³⁴⁵⁶⁷⁸⁹';const ex=(e<0?'⁻':'')+String(Math.abs(e)).split('').map(c=>SUP[+c]).join('');
  return (e>=-2&&e<4?(+v.toPrecision(3)).toString().replace('.',','):mnt.toFixed(2).replace('.',',')+'·10'+ex)+' '+u;}
function gravReadout(i){
  if(!gAny) return '';
  const a=gTa[i],x=i%W+.5,y=(i/W|0)+.5;let h=kv('g local',sciM(a,'m/s²')+' ('+sciM(a/G0,'g')+')');
  const phi=gPhi[i];if(phi<0)h+=kv('v escape',sciM(Math.sqrt(-2*phi),'m/s'));
  let s0=null,bd=1e9;for(const s of gS){const d=Math.hypot(s.x-x,s.y-y);if(d<bd){bd=d;s0=s;}}
  if(s0){const rs=rsOf(s0.M),r=bd*DX;
    h+=kv('fonte',sciM(s0.M,'kg')+' a '+r.toFixed(1)+' m');
    h+=kv('r_s',sciM(rs,'m')+(s0.bh?' · buraco negro':''));
    if(s0.bh){h+=kv('tempo local',r>rs?'×'+Math.sqrt(1-rs/r).toFixed(4)+' (dilatação)':'dentro do horizonte');
      h+=kv('T Hawking',sciM(HBAR*CLUZ*CLUZ*CLUZ/(8*Math.PI*GN*s0.M*KB),'K'));
      if(s0.abs>0)h+=kv('absorvido',sciM(s0.abs,'kg'));}
    else{const df=1+phi/(CLUZ*CLUZ);if(df<.999999)h+=kv('tempo local','×'+df.toFixed(6));}
    if(gTid[i]>1e-3)h+=kv('maré',sciM(gTid[i],'s⁻²'));}
  return h;
}

// ---------- renderização: vista de gravidade e lente do buraco negro ----------
const GVIEW=7,GST2=[[0,[10,8,20]],[.2,[48,20,110]],[.4,[120,30,160]],[.6,[220,60,120]],[.8,[255,150,50]],[1,[255,250,220]]];
const glut=new Uint32Array(256);for(let k=0;k<256;k++){const c=lerpStops(GST2,k/255);glut[k]=pack(c[0]|0,c[1]|0,c[2]|0);}
const G_LO=-2,G_HI=14,lensB=new Uint32Array(N),BLK=pack(0,0,0);
function gravPost(px){
  if(view===GVIEW){const gu=gm*G0,iv=255/(G_HI-G_LO);
    for(let i=0;i<N;i++){const a=gAny?gTa[i]:gu;let k=a>0?((Math.log10(a)-G_LO)*iv)|0:0;k=k<0?0:k>255?255:k;
      const m=mat[i],kk=kind[m];px[i]=kk===GAS?glut[k]:blend(dim(px[i],.55),glut[k],kk===FIX?.25:.4);}}
  if(!gBH) return;
  for(const s of gS){if(!s.bh)continue;
    const R=rsOf(s.M)/DX,Rsh=2.6*R,TE2=Rsh*(Rsh-R),RL=Math.max(6,14*R),RL2=RL*RL;
    const x0=Math.max(0,(s.x-RL)|0),x1=Math.min(W-1,(s.x+RL)|0),y0=Math.max(0,(s.y-RL)|0),y1=Math.min(H-1,(s.y+RL)|0);
    for(let y=y0;y<=y1;y++){const r0=y*W;for(let x=x0;x<=x1;x++)lensB[r0+x]=px[r0+x];}
    for(let y=y0;y<=y1;y++)for(let x=x0;x<=x1;x++){
      const dx=x+.5-s.x,dy=y+.5-s.y,t2=dx*dx+dy*dy;if(t2>RL2)continue;const i=y*W+x,t=Math.sqrt(t2);
      if(t<Rsh){px[i]=BLK;continue;} // sombra (r < 2,6 r_s: fótons capturados)
      let f=(t-.6*RL)/(.4*RL);f=f<0?1:f>1?0:1-f*f*(3-2*f); // desvanece a lente na borda
      const b=t-TE2/t*f,sx=s.x+dx/t*b,sy=s.y+dy/t*b; // lente β = θ − θ_E²/θ, θ_E² = 2,6·1,6 r_s²: borda da sombra mostra o horizonte (anel comprimido)
      const xi=sx|0,yi=sy|0;let c=(xi>=0&&xi<W&&yi>=0&&yi<H)?lensB[yi*W+xi]:BLK;
      const pr=(t-Rsh)/Math.max(.6,.25*R);if(pr<1)c=blend(c,pack(255,214,150),.85*(1-pr)); // anel de fótons
      px[i]=c;}
  }
}
// sobreposição: malha do poço de potencial, setas de g, marcadores, disco de acreção e ondas gravitacionais
function drawGrav(oc,sx,sy,dpr,fn){
  if(!gAny&&!gWaves.length) return;
  const showM=view===GVIEW||tool===GRAVT||tool2===GRAVT;
  if(view===GVIEW&&gAny){
    const q=gS.map(s=>Math.max(0,Math.min(1,Math.log10(s.M/GM0)/(Math.log10(GMMAX/GM0))))),R0=40,ST=8;
    const warp=(x,y)=>{let X=x,Y=y;for(let s=0;s<gS.length;s++){const S=gS[s],dx=S.x-x,dy=S.y-y,r=Math.hypot(dx,dy);
      if(r<1e-6)continue;const w=.85*q[s]*R0/(r+R0);X+=dx*w;Y+=dy*w;}return [X,Y];};
    oc.lineWidth=1*dpr;oc.strokeStyle='rgba(190,170,255,.32)';oc.beginPath();
    for(let g=0;g<=W;g+=ST){for(let t=0;t<=H;t+=2){const [X,Y]=warp(g,t);t?oc.lineTo(X*sx,Y*sy):oc.moveTo(X*sx,Y*sy);}
      for(let t=0;t<=W;t+=2){const [X,Y]=warp(t,g);t?oc.lineTo(X*sx,Y*sy):oc.moveTo(X*sx,Y*sy);}}
    oc.stroke();
    // setas: direção de g, comprimento ∝ log|g|
    oc.strokeStyle='rgba(255,235,200,.7)';oc.lineWidth=1.2*dpr;oc.beginPath();
    for(let y=6;y<H;y+=12)for(let x=6;x<W;x+=12){const i=y*W+x,a=gTa[i];if(a<.05)continue;
      let L=(Math.log10(a)-G_LO)/(G_HI-G_LO)*9+1.5;L=L>10?10:L;const ux=gTx[i]/a,uy=gTy[i]/a,X=(x+.5)*sx,Y=(y+.5)*sy,ex=X+ux*L*sx,ey=Y+uy*L*sy;
      oc.moveTo(X,Y);oc.lineTo(ex,ey);oc.moveTo(ex,ey);oc.lineTo(ex-(ux*3-uy*2)*dpr,ey-(uy*3+ux*2)*dpr);oc.moveTo(ex,ey);oc.lineTo(ex-(ux*3+uy*2)*dpr,ey-(uy*3-ux*2)*dpr);}
    oc.stroke();
  }
  for(const s of gS){const X=s.x*sx,Y=s.y*sy,R=rsOf(s.M)/DX;
    if(s.bh){ // brilho do disco de acreção (dentro da ISCO)
      oc.globalCompositeOperation='lighter';const Ra=Math.max(3,3*R)*1.6*sx,gr=oc.createRadialGradient(X,Y,2.6*R*sx,X,Y,Ra);
      gr.addColorStop(0,'rgba(255,170,80,.45)');gr.addColorStop(1,'rgba(255,90,20,0)');oc.fillStyle=gr;
      oc.beginPath();oc.arc(X,Y,Ra,0,6.2832);oc.fill();oc.globalCompositeOperation='source-over';
      oc.strokeStyle='rgba(255,220,170,.9)';oc.lineWidth=1.2*dpr;oc.beginPath();oc.arc(X,Y,Math.max(2.6*R*sx,2*dpr),0,6.2832);oc.stroke();}
    if(showM){const pr=(4+2*Math.sin(fn*.15))*dpr;oc.strokeStyle=s.bh?'rgba(255,200,140,.9)':'rgba(200,170,255,.9)';oc.lineWidth=1.2*dpr;
      oc.beginPath();oc.arc(X,Y,pr,0,6.2832);oc.moveTo(X-pr*1.8,Y);oc.lineTo(X-pr*.6,Y);oc.moveTo(X+pr*.6,Y);oc.lineTo(X+pr*1.8,Y);
      oc.moveTo(X,Y-pr*1.8);oc.lineTo(X,Y-pr*.6);oc.moveTo(X,Y+pr*.6);oc.lineTo(X,Y+pr*1.8);oc.stroke();
      const lbl=sciM(s.M,'kg')+(s.bh?' · BN r_s '+sciM(rsOf(s.M),'m'):'');
      oc.font=(10.5*dpr)+'px ui-monospace,monospace';const tw=oc.measureText(lbl).width;
      oc.fillStyle='rgba(0,0,0,.6)';oc.fillRect(X+8*dpr,Y-20*dpr,tw+8*dpr,15*dpr);oc.fillStyle='#fff';oc.fillText(lbl,X+12*dpr,Y-9*dpr);}
  }
  // ondas gravitacionais da fusão (anel que se expande à velocidade da luz, desenhado lento)
  for(let w=gWaves.length-1;w>=0;w--){const g=gWaves[w];g.r+=2.5;g.a*=.965;if(g.a<.03){gWaves.splice(w,1);continue;}
    oc.strokeStyle='rgba(200,180,255,'+g.a.toFixed(3)+')';oc.lineWidth=2*dpr;
    for(let k=0;k<3;k++){const r=g.r-k*7;if(r<=0)continue;oc.beginPath();oc.arc(g.x*sx,g.y*sy,r*sx,0,6.2832);oc.stroke();}}
}
function gravLegend(){lgGrad.style.background='linear-gradient(90deg,'+GST2.map(s=>css(s[1])+' '+s[0]*100+'%').join(',')+')';
  lgLo.textContent='10⁻² m/s²';lgHi.textContent='10¹⁴ m/s²';}
document.getElementById('bGravZ').onclick=()=>{pushUndo();gravReset();};
