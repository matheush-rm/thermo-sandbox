'use strict';
// ================= constantes físicas =================
const W=200,H=200,N=W*H,DX=0.1,PATM=101325,RG=8.314,G0=9.81,STEAM_MASS=0.05,HS=0.07,FIRE_Q=0.5;
const GAS=1,LIQ=2,POW=3,MOV=4,FIX=5;
const AR=0,PAREDE=1,AGUA=2,SALMOURA=3,OLEO=4,MERCURIO=5,MEL=6,ALCOOL=7,LAVA=8,N2L=9,ALUF=10,ACOF=11,
  ACO=12,CONCRETO=13,VIDRO=14,LATA=15,BASALTO=16,RESIST=17,MADEIRA=18,MADEIRAE=19,ISOPOR=20,ISOPORC=21,
  GELO=22,HGS=23,AREIA=24,SAL=25,CACOS=26,ENTULHO=27,SUCATA=28,CINZA=29,VAPOR=30,VALC=31,N2G=32,VHG=33,
  FUMACA=34,FOGO=35,HELIO=36,COBRE=37,COBREF=38,NICROMO=39,IMA=40,LIMALHA=41,FERROFL=42,FONTEP=43,FONTEN=44,
  TERRA=45,NUVEM=46,GRANIZO=47,VIDA=48,C4=49,NITRO=50,HNO3=51,
  U235=52,PU239=53,UNAT=54,CORIO=55,PFIS=56,RADIO=57,RADONIO=58,COBALTO=59,CESIO=60,CHUMBO=61,CHUMBOF=62,GRAFITE=63,BORO=64,AMBE=65,VIDAM=66,NM=67;
const ERASE=100,HEAT=101,COOL=102,CHGP=103,CHGN=104,FLIP=105,SPARK=106,WIND=107,CYCL=108;
const GRAVT=120; // ferramenta Gravidade (grav.js)
const CONT={[GRAVT]:1,[HEAT]:1,[COOL]:1,[CHGP]:1,[CHGN]:1,[SPARK]:1,[CYCL]:1}; // ferramentas aplicadas continuamente

// hi/lo: [tipo(1 const,2 ebulição C-C,3 Clapeyron linear), T0(K), Lm(J/mol) ou dT/dP, destino, latente(K eq.), Tcrítica]
const D=[];
D[AR]       ={n:'Ar',kind:GAS,M:.029,gam:1.4,cpg:1005,k:.026,col:[14,17,24]};
D[PAREDE]   ={n:'Parede',kind:FIX,rho:7850,k:.05,C:3,col:[44,48,58]};
D[AGUA]     ={n:'Água',kind:LIQ,rho:1000,k:.6,C:4.18,beta:2.1e-4,kap:4.6e-10,mob:1,disp:5,col:[38,104,214],var:8,hi:[2,373.15,40660,VAPOR,540,647],lo:[3,273.15,-7.4e-8,GELO,80]};
D[SALMOURA] ={n:'Água salgada',kind:LIQ,rho:1100,k:.58,C:3.9,beta:3e-4,kap:4e-10,mob:.95,disp:5,col:[44,140,165],var:8,hi:[2,381.5,40660,VAPOR,560,650],lo:[1,252,0,GELO,80]};
D[OLEO]     ={n:'Óleo',kind:LIQ,rho:920,k:.17,C:1.7,beta:7e-4,kap:6e-10,mob:.45,disp:3,col:[196,150,30],var:10,ign:600,fl:60,fsp:.25};
D[MERCURIO] ={n:'Mercúrio',kind:LIQ,rho:13534,k:8.3,C:1.9,beta:1.8e-4,kap:4e-11,mob:1,disp:3,col:[176,182,194],var:18,hi:[2,629.9,59110,VHG,2100,1750],lo:[1,234.3,0,HGS,81]};
D[MEL]      ={n:'Mel',kind:LIQ,rho:1420,k:.5,C:3.2,beta:4.8e-4,kap:3e-10,mob:.08,disp:2,col:[214,138,16],var:10};
D[ALCOOL]   ={n:'Álcool',kind:LIQ,rho:789,k:.17,C:1.92,beta:1.1e-3,kap:1.1e-9,mob:1,disp:6,col:[160,180,228],var:6,hi:[2,351.4,38560,VALC,345,514],ign:636,fl:24,fsp:.6};
D[LAVA]     ={n:'Lava',kind:LIQ,rho:2600,k:1.5,C:3.0,beta:1e-4,kap:1e-10,mob:.08,disp:2,col:[255,84,10],var:40,lo:[1,1000,0,BASALTO,130],T0:1450};
D[N2L]      ={n:'N₂ líquido',kind:LIQ,rho:808,k:.14,C:1.65,beta:4.6e-3,kap:2e-9,mob:1,disp:6,col:[196,226,255],var:6,hi:[2,77.36,5570,N2G,97,126],T0:77};
D[ALUF]     ={n:'Alumínio fundido',kind:LIQ,rho:2375,k:90,C:2.6,beta:1.2e-4,kap:2e-11,mob:.8,disp:2,col:[232,196,180],var:10,lo:[1,933,0,SUCATA,336],T0:1000};
D[ACOF]     ={n:'Aço fundido',kind:LIQ,rho:6900,k:30,C:5.5,beta:1e-4,kap:1e-11,mob:.7,disp:2,col:[255,190,110],var:10,lo:[1,1783,0,ACO,300],T0:1900};
D[ACO]      ={n:'Aço',kind:FIX,rho:7850,k:50,C:3.6,str:250e6,crush:SUCATA,col:[118,128,142],var:6,hi:[1,1783,0,ACOF,300]};
D[CONCRETO] ={n:'Concreto',kind:FIX,rho:2400,k:1.7,C:2.0,str:30e6,crush:ENTULHO,col:[148,146,140],var:22};
D[VIDRO]    ={n:'Vidro',kind:FIX,rho:2500,k:1.0,C:2.1,str:2e6,crush:CACOS,col:[140,200,218],var:6,shock:120};
D[LATA]     ={n:'Lata (alumínio)',kind:FIX,rho:2700,k:237,C:2.42,str:110e3,crush:SUCATA,col:[196,202,214],var:8,hi:[1,933,0,ALUF,336]};
D[BASALTO]  ={n:'Basalto',kind:FIX,rho:2900,k:2,C:2.4,str:150e6,crush:ENTULHO,col:[64,58,58],var:16,hi:[1,1473,0,LAVA,130]};
D[RESIST]   ={n:'Resistência',kind:FIX,rho:8000,k:20,C:3.5,col:[150,62,40],var:4};
D[MADEIRA]  ={n:'Madeira',kind:MOV,rho:600,k:.12,C:1.0,col:[150,98,52],var:22,ign:573,fl:140,fsp:.04};
D[MADEIRAE] ={n:'Madeira encharcada',kind:MOV,rho:1100,k:.35,C:2.6,col:[88,62,40],var:14,ign:680,fl:140,fsp:.015};
D[ISOPOR]   ={n:'Isopor',kind:MOV,rho:20,k:.033,C:.026,str:100e3,col:[240,240,234],var:12,hi:[1,513,0,FUMACA,40],ign:700,fl:14,fsp:.3};
D[ISOPORC]  ={n:'Isopor comprimido',kind:MOV,rho:160,k:.04,C:.2,col:[200,200,186],var:10,hi:[1,513,0,FUMACA,40]};
D[GELO]     ={n:'Gelo',kind:MOV,rho:917,k:2.2,C:1.93,col:[186,226,250],var:12,hi:[3,273.15,-7.4e-8,AGUA,159],T0:258};
D[HGS]      ={n:'Mercúrio sólido',kind:MOV,rho:14190,k:30,C:1.9,col:[140,146,158],var:10,hi:[1,234.3,0,MERCURIO,81],T0:200};
D[AREIA]    ={n:'Areia',kind:POW,rho:1600,k:.3,C:1.3,col:[218,186,118],var:26,hi:[1,1973,0,VIDRO,200]};
D[SAL]      ={n:'Sal',kind:POW,rho:2160,k:6,C:1.9,col:[236,236,242],var:12};
D[CACOS]    ={n:'Cacos de vidro',kind:POW,rho:2500,k:1,C:2.1,col:[150,205,220],var:40};
D[ENTULHO]  ={n:'Entulho',kind:POW,rho:2000,k:1,C:1.8,col:[118,114,108],var:28};
D[SUCATA]   ={n:'Sucata metálica',kind:POW,rho:2700,k:237,C:2.42,col:[160,166,178],var:30,hi:[1,933,0,ALUF,336]};
D[CINZA]    ={n:'Cinza',kind:POW,rho:600,k:.1,C:.5,col:[86,86,86],var:18};
D[VAPOR]    ={n:"Vapor d'água",kind:GAS,M:.018,gam:1.33,cpg:2010,k:.025,col:[205,215,228],a:.5,lo:[2,373.15,40660,AGUA,1124,647],T0:390};
D[VALC]     ={n:'Vapor de álcool',kind:GAS,M:.046,gam:1.13,cpg:1420,k:.015,col:[175,185,225],a:.4,lo:[2,351.4,38560,ALCOOL,600,514],ign:636,fl:12,fsp:.7,T0:360};
D[N2G]      ={n:'N₂ gasoso',kind:GAS,M:.028,gam:1.4,cpg:1040,k:.026,col:[150,175,215],a:.28,lo:[2,77.36,5570,N2L,191,126]};
D[VHG]      ={n:'Vapor de mercúrio',kind:GAS,M:.2006,gam:1.67,cpg:104,k:.008,col:[150,175,150],a:.4,lo:[2,629.9,59110,MERCURIO,2800,1750],T0:650};
D[FUMACA]   ={n:'Fumaça',kind:GAS,M:.029,gam:1.4,cpg:1005,k:.026,col:[80,80,84],a:.6};
D[FOGO]     ={n:'Fogo',kind:GAS,M:.029,gam:1.4,cpg:1005,k:.06,col:[255,130,30],a:1,T0:1300};
D[HELIO]    ={n:'Hélio',kind:GAS,M:.004,gam:1.67,cpg:5193,k:.15,col:[235,170,215],a:.3};
D[COBRE]    ={n:'Cobre',kind:FIX,rho:8960,k:400,C:3.45,col:[196,112,62],var:10,hi:[1,1358,0,COBREF,530]};
D[COBREF]   ={n:'Cobre fundido',kind:LIQ,rho:8000,k:160,C:4.2,beta:1e-4,kap:2e-11,mob:.8,disp:2,col:[255,150,90],var:10,lo:[1,1358,0,COBRE,530],T0:1450};
D[NICROMO]  ={n:'Nicromo',kind:FIX,rho:8400,k:11.3,C:3.8,col:[126,116,104],var:6,hi:[1,1673,0,ACOF,300]};
D[IMA]      ={n:'Ímã (NdFeB)',kind:FIX,rho:7500,k:9,C:3.3,col:[78,80,90],var:4,hi:[1,1583,0,ACOF,300]};
D[LIMALHA]  ={n:'Limalha de ferro',kind:POW,rho:3900,k:5,C:1.8,col:[74,76,84],var:22,hi:[1,1811,0,ACOF,250]};
D[FERROFL]  ={n:'Ferrofluido',kind:LIQ,rho:1250,k:.2,C:2.0,beta:7e-4,kap:6e-10,mob:.5,disp:3,col:[26,24,30],var:6};
D[FONTEP]   ={n:'Polo + (fonte)',kind:FIX,rho:8000,k:50,C:3.5,col:[214,52,48],var:4};
D[FONTEN]   ={n:'Polo − (fonte)',kind:FIX,rho:8000,k:50,C:3.5,col:[48,92,222],var:4};
D[TERRA]    ={n:'Terra (aterramento)',kind:FIX,rho:2000,k:1,C:2,col:[72,98,58],var:10};
D[NUVEM]    ={n:'Nuvem',kind:GAS,M:.029,gam:1.4,cpg:1005,k:.026,col:[196,202,214],a:.8,var:14,hi:[2,373.15,40660,VAPOR,30,647]};
D[VIDA]     ={n:'Vida (Conway)',kind:FIX,rho:1050,k:.5,C:3.6,str:400e3,crush:CINZA,col:[70,214,110],var:24,ign:450,fl:20,fsp:.5};
// ex: [kg de TNT equivalente por célula (1 L), T de detonação (K), sensibilidade a choque (Pa, 0 = insensível)]
D[C4]       ={n:'C-4',kind:FIX,rho:1600,k:.25,C:1.6,col:[222,212,176],var:10,ex:[2.1,520,0]};
D[NITRO]    ={n:'Nitroglicerina',kind:LIQ,rho:1600,k:.29,C:1.9,beta:8e-4,kap:5e-10,mob:.55,disp:3,col:[226,232,168],var:8,ex:[2.4,490,1.5e5]};
D[HNO3]     ={n:'Ácido nítrico',kind:LIQ,rho:1510,k:.3,C:2.6,beta:1.2e-3,kap:7e-10,mob:1,disp:5,col:[232,196,92],var:10,hi:[1,356,0,FUMACA,200]};
D[GRANIZO]  ={n:'Granizo (graupel)',kind:POW,rho:700,k:2,C:1.4,col:[226,236,246],var:14,hi:[1,273.15,0,AGUA,80],T0:258};
// radioativos / nucleares (emissão e nêutrons em rad.js)
D[U235]     ={n:'Urânio-235 (enriquecido)',kind:MOV,rho:19050,k:27,C:2.2,col:[96,116,86],var:10,hi:[1,1405,0,CORIO,330]};
D[PU239]    ={n:'Plutônio-239',kind:MOV,rho:19800,k:6.7,C:2.6,col:[150,138,128],var:10,hi:[1,913,0,CORIO,90]};
D[UNAT]     ={n:'Urânio natural',kind:MOV,rho:19050,k:27,C:2.2,col:[72,78,72],var:10,hi:[1,1405,0,CORIO,330]};
D[CORIO]    ={n:'Cório (combustível fundido)',kind:LIQ,rho:8000,k:3,C:4,beta:1e-4,kap:1e-10,mob:.12,disp:2,col:[255,118,40],var:36,lo:[1,1300,0,PFIS,200],hi:[1,3500,0,FUMACA,300],T0:2000};
D[PFIS]     ={n:'Produtos de fissão (fallout)',kind:POW,rho:1200,k:.3,C:1.2,col:[118,128,96],var:30};
D[RADIO]    ={n:'Rádio-226',kind:POW,rho:5000,k:18.6,C:.6,col:[196,226,206],var:12};
D[RADONIO]  ={n:'Radônio-222',kind:GAS,M:.222,gam:1.67,cpg:94,k:.0036,col:[140,220,170],a:.35};
D[COBALTO]  ={n:'Cobalto-60',kind:FIX,rho:8900,k:100,C:3.7,str:200e6,crush:SUCATA,col:[112,124,178],var:6};
D[CESIO]    ={n:'Cloreto de césio-137',kind:POW,rho:3990,k:1,C:1.9,col:[224,232,255],var:14};
D[CHUMBO]   ={n:'Chumbo',kind:FIX,rho:11340,k:35,C:1.45,str:12e6,crush:SUCATA,col:[92,96,110],var:6,hi:[1,600.6,0,CHUMBOF,180]};
D[CHUMBOF]  ={n:'Chumbo fundido',kind:LIQ,rho:10660,k:16,C:1.5,beta:1.2e-4,kap:2e-11,mob:.9,disp:3,col:[176,178,190],var:8,lo:[1,600.6,0,CHUMBO,180],T0:700};
D[GRAFITE]  ={n:'Grafite (moderador)',kind:FIX,rho:1700,k:120,C:1.2,str:20e6,crush:CINZA,col:[50,50,56],var:8,ign:1020,fl:200,fsp:.02};
D[BORO]     ={n:'Carbeto de boro (absorvedor)',kind:FIX,rho:2520,k:30,C:2.4,str:300e6,crush:ENTULHO,col:[36,42,36],var:6};
D[AMBE]     ={n:'Fonte de nêutrons (Am-Be)',kind:FIX,rho:3000,k:5,C:2,col:[196,118,220],var:6};
D[VIDAM]    ={n:'Vida mutante (B36/S23)',kind:FIX,rho:1050,k:.5,C:3.6,str:400e3,crush:CINZA,col:[206,96,226],var:24,ign:450,fl:20,fsp:.5};

// elétrico: [σ condutividade S/m, εr permissividade relativa, rigidez dielétrica MV/m (0 = não rompe)]
const EL={
 [AR]:[3e-15,1,3],[PAREDE]:[0,1,0],[AGUA]:[5e-3,80,0],[SALMOURA]:[5,75,0],[OLEO]:[1e-12,2.2,15],[MERCURIO]:[1.04e6,1,0],
 [MEL]:[2e-2,20,0],[ALCOOL]:[1.35e-7,24,10],[LAVA]:[10,10,0],[N2L]:[1e-16,1.43,20],[ALUF]:[4e6,1,0],[ACOF]:[7e5,1,0],
 [ACO]:[6e6,1,0],[CONCRETO]:[1e-6,5,3],[VIDRO]:[1e-12,6,10],[LATA]:[3.5e7,1,0],[BASALTO]:[1e-6,8,15],[RESIST]:[9e5,1,0],
 [MADEIRA]:[1e-12,3,5],[MADEIRAE]:[1e-4,20,0],[ISOPOR]:[1e-16,1.03,2],[ISOPORC]:[1e-15,1.2,4],[GELO]:[1e-8,3.2,10],[HGS]:[4e6,1,0],
 [AREIA]:[1e-9,4,3],[SAL]:[1e-12,6,10],[CACOS]:[1e-12,6,10],[ENTULHO]:[1e-7,5,3],[SUCATA]:[1e5,1,0],[CINZA]:[1e-8,3,3],
 [VAPOR]:[1e-14,1,2.7],[VALC]:[1e-14,1,3],[N2G]:[1e-15,1,3.4],[VHG]:[1e-14,1,1],[FUMACA]:[1e-13,1,2.5],[FOGO]:[1e-4,1,.6],[HELIO]:[1e-15,1,.45],
 [COBRE]:[5.96e7,1,0],[COBREF]:[5e6,1,0],[NICROMO]:[9.1e5,1,0],[IMA]:[6.7e5,1,0],[LIMALHA]:[1e3,1,0],[FERROFL]:[1e-8,5,10],
 [FONTEP]:[6e7,1,0],[FONTEN]:[6e7,1,0],[TERRA]:[6e7,1,0],[NUVEM]:[1e-13,1,1],[GRANIZO]:[1e-9,3.2,10],[VIDA]:[.5,50,0],
 [C4]:[1e-12,3,15],[NITRO]:[1e-9,19,10],[HNO3]:[50,50,0],
 [U235]:[3.6e6,1,0],[PU239]:[6.7e5,1,0],[UNAT]:[3.6e6,1,0],[CORIO]:[1e5,1,0],[PFIS]:[1e-6,5,3],[RADIO]:[1e-6,6,5],[RADONIO]:[1e-14,1,2.5],
 [COBALTO]:[1.6e7,1,0],[CESIO]:[1e-8,6,10],[CHUMBO]:[4.8e6,1,0],[CHUMBOF]:[1e6,1,0],[GRAFITE]:[1e5,1,0],[BORO]:[1e2,8,0],[AMBE]:[1e-6,5,5],[VIDAM]:[.5,50,0]
};
// coef. térmico da resistividade (1/K)
const ALPHA={[COBRE]:.0039,[LATA]:.0039,[ACO]:.005,[NICROMO]:.0004,[MERCURIO]:.0009,[IMA]:.001,[ALUF]:.001,[COBREF]:.001,[ACOF]:.001};
// magnético: [1 = ferromagneto mole / 2 = ímã permanente, Ms saturação A/m, χ efetivo (com desmagnetização), T Curie K]
const MG={[ACO]:[1,1.6e6,2.9,1043],[LIMALHA]:[1,8e5,5,1043],[FERROFL]:[1,4.5e4,1.2,850],[IMA]:[2,9.5e5,0,583]};
const SRC={[FONTEP]:1,[FONTEN]:2,[TERRA]:3}; // potencial fixo: +V/2, −V/2, 0

// compila em arrays tipados
const kind=new Uint8Array(NM),rho0=new Float32Array(NM),betaA=new Float32Array(NM),bkA=new Float32Array(NM),
  kc=new Float32Array(NM),CvA=new Float32Array(NM),molar=new Float32Array(NM),gexp=new Float32Array(NM),cpgA=new Float32Array(NM),
  mobA=new Float32Array(NM),dispA=new Uint8Array(NM),strA=new Float32Array(NM),crushA=new Uint8Array(NM),
  ignA=new Float32Array(NM),flA=new Float32Array(NM),fspA=new Float32Array(NM),T0A=new Float32Array(NM),shockA=new Float32Array(NM),
  hiTy=new Uint8Array(NM),hiT=new Float32Array(NM),hiLm=new Float64Array(NM),hiTo=new Uint8Array(NM),hiL=new Float32Array(NM),hiCap=new Float32Array(NM),
  loTy=new Uint8Array(NM),loT=new Float32Array(NM),loLm=new Float64Array(NM),loTo=new Uint8Array(NM),loL=new Float32Array(NM),loCap=new Float32Array(NM),
  actA=new Uint8Array(NM),fillA=new Float32Array(NM),
  sigA=new Float64Array(NM),epsA=new Float32Array(NM),ebA=new Float32Array(NM),alA=new Float32Array(NM),srcA=new Uint8Array(NM),
  mgA=new Uint8Array(NM),msA=new Float32Array(NM),chiA=new Float32Array(NM),tcA=new Float32Array(NM),mgMov=new Uint8Array(NM);
for(let m=0;m<NM;m++){const d=D[m],e=EL[m],g=MG[m];
  sigA[m]=e[0]; epsA[m]=e[1]; ebA[m]=e[2]*1e6; alA[m]=ALPHA[m]||0; srcA[m]=SRC[m]||0;
  if(g){mgA[m]=g[0];msA[m]=g[1];chiA[m]=g[2];tcA[m]=g[3];mgMov[m]=g[0]===1&&d.kind!==FIX?1:0;}
  kind[m]=d.kind; rho0[m]=d.rho||0; betaA[m]=d.beta||1e-5; bkA[m]=d.kap?d.beta/d.kap:0;
  kc[m]=d.k*HS; CvA[m]=d.C||0; molar[m]=d.M||0; gexp[m]=d.gam?(d.gam-1)/d.gam:0; cpgA[m]=d.cpg||0;
  mobA[m]=d.mob||1; dispA[m]=d.disp||1; strA[m]=d.str||0; crushA[m]=d.crush||0;
  ignA[m]=d.ign||0; flA[m]=d.fl||0; fspA[m]=d.fsp||0; T0A[m]=d.T0||293.15; shockA[m]=d.shock||0;
  if(d.hi){hiTy[m]=d.hi[0];hiT[m]=d.hi[1];hiLm[m]=d.hi[2];hiTo[m]=d.hi[3];hiL[m]=d.hi[4];hiCap[m]=d.hi[5]||9999;}
  if(d.lo){loTy[m]=d.lo[0];loT[m]=d.lo[1];loLm[m]=d.lo[2];loTo[m]=d.lo[3];loL[m]=d.lo[4];loCap[m]=d.lo[5]||9999;}
  actA[m]=(d.hi||d.lo||d.str||d.ign||d.shock||d.ex||m===FOGO||m===FUMACA||m===RESIST||m===SAL||m===MADEIRA)?1:0;
  fillA[m]=d.kind===POW?.55:d.kind===LIQ?.85:d.kind===GAS?.6:1;
}

// ================= estado (estrutura de arrays) =================
const mat=new Uint8Array(N),shade=new Uint8Array(N),life=new Uint8Array(N),mv=new Uint8Array(N),vis=new Uint8Array(N);
const T=new Float32Array(N),tOld=new Float32Array(N),pp=new Float32Array(N),gd=new Float32Array(N),lat=new Float32Array(N),st=new Float32Array(N),rhoA=new Float32Array(N);
let P=new Float32Array(N),Pn=new Float32Array(N);
const pvx=new Float32Array(N),pvy=new Float32Array(N); // velocidade balística de grãos/líquidos arremessados (células/passo)
// carga (C por célula), ionização (0..1), magnetização permanente (A/m)
const chg=new Float64Array(N),ion=new Float32Array(N),mxA=new Float32Array(N),myA=new Float32Array(N),flipMk=new Uint16Array(N);
let Vsrc=100,qAmt=1e-7,magDir=0,storm=false,windF=400,elecLive=false,magLive=false,vScale=1,iMax=0,strokeId=1;
const CBASE=100; // base da nuvem (linha) no modo tempestade
function si(v,u,d){const a=Math.abs(v);if(a===0)return '0 '+u;
  const P=[[1e9,'G'],[1e6,'M'],[1e3,'k'],[1,''],[1e-3,'m'],[1e-6,'µ'],[1e-9,'n'],[1e-12,'p'],[1e-15,'f']];
  for(const [f,p] of P) if(a>=f*.9995) return (v/f).toFixed(d===undefined?(a/f<10?2:a/f<100?1:0):d)+' '+p+u;
  return v.toExponential(1)+' '+u;}
const fmtV=v=>si(v,'V'),fmtQ=v=>si(v,'C'),fmtI=v=>si(v,'A'),fmtB=v=>si(v,'T'),fmtE=v=>si(v,'V/m');
const fmtS=s=>s===0?'0':s.toExponential(1).replace('e+','·10^').replace('e-','·10^-');
const q=new Int32Array(N);
let stamp=1,gm=1,heaterT=773.15,power=4,lid=false,paused=false,steps=2,view=0;

let seed=0x2545F491;
function rnd(){seed^=seed<<13;seed^=seed>>>17;seed^=seed<<5;return (seed>>>0)/4294967296;}

function densOf(i){
  const m=mat[i];
  if(kind[m]===GAS) return gd[i];
  if(m===AGUA){ // Thiesen: máximo em ~4 °C
    let c=T[i]-273.15; if(c<-10)c=-10; else if(c>150)c=150;
    const a=c-3.983035; return 999.975*(1-a*a*(c+301.797)/(522528.9*(c+69.34881)));
  }
  return rho0[m]*(1-betaA[m]*(T[i]-293.15));
}
function mobOf(i,m){
  if(m===MEL){const c=T[i]-273.15;const v=0.8/(10*Math.exp(-0.075*(c-20)));return v>1?1:v<.003?.003:v;}
  if(m===LAVA){const v=8/(100*Math.exp(0.02*(1473-T[i])));return v>1?1:v<.002?.002:v;}
  return mobA[m];
}
function tb(T0,Lm,cap,p){ // Clausius–Clapeyron
  if(p<200)p=200; const d=1/T0-RG*Math.log(p/PATM)/Lm;
  return d<=1/cap?cap:1/d;
}

function setCell(i,m,t){
  mat[i]=m;T[i]=t;tOld[i]=t;lat[i]=0;life[i]=0;shade[i]=(rnd()*8)|0;st[i]=0;chg[i]=0;ion[i]=0;mxA[i]=0;myA[i]=0;pvx[i]=0;pvy[i]=0;
  const p=P[i]>0?P[i]:PATM; pp[i]=p;
  gd[i]=kind[m]===GAS?p*molar[m]/(RG*t):0;
  if(m===FOGO){life[i]=60+(rnd()*60|0);lat[i]=0;}
  else if(m===FUMACA) life[i]=150;
  radClearCell(i);
}
function convert(i,to){
  const from=mat[i];
  if(from===SALMOURA&&to===VAPOR&&rnd()<.25) to=SAL;
  mat[i]=to;lat[i]=0;life[i]=0;
  if(to===PFIS&&rA[i]<PFIS_RA){rA[i]=PFIS_RA;rAny=true;} // cório solidificado: detrito radioativo
  const p=P[i]>0?P[i]:PATM; pp[i]=p; tOld[i]=T[i];
  if(kind[to]===GAS){
    if(kind[from]!==GAS) aInj[AF[i]]+=VAP_SRC; // vapor novo se expande (ebulição, explosão de vapor)
    const eq=p*molar[to]/(RG*T[i]);
    gd[i]=(kind[from]===LIQ&&vis[i]!==1)?Math.max(eq,rho0[from]*STEAM_MASS):eq;
    if(to===FUMACA) life[i]=120+(rnd()*100|0);
  }
}
function swap(i,j){
  let t=mat[i];mat[i]=mat[j];mat[j]=t;
  t=shade[i];shade[i]=shade[j];shade[j]=t;
  t=life[i];life[i]=life[j];life[j]=t;
  let f=T[i];T[i]=T[j];T[j]=f;
  f=tOld[i];tOld[i]=tOld[j];tOld[j]=f;
  f=pp[i];pp[i]=pp[j];pp[j]=f;
  f=gd[i];gd[i]=gd[j];gd[j]=f;
  f=lat[i];lat[i]=lat[j];lat[j]=f;
  f=chg[i];chg[i]=chg[j];chg[j]=f; // carga viaja com a matéria (chuva, granizo)
  f=ion[i];ion[i]=ion[j];ion[j]=f;
  f=pvx[i];pvx[i]=pvx[j];pvx[j]=f; f=pvy[i];pvy[i]=pvy[j];pvy[j]=f; // momento balístico (detritos)
  f=rA[i];rA[i]=rA[j];rA[j]=f; f=bu[i];bu[i]=bu[j];bu[j]=f; f=cmp[i];cmp[i]=cmp[j];cmp[j]=f; // contaminação viaja com a matéria
  mv[i]=stamp;mv[j]=stamp;
}
function copyFrom(i,j){mat[i]=mat[j];T[i]=T[j];tOld[i]=tOld[j];pp[i]=pp[j];gd[i]=gd[j];lat[i]=0;life[i]=0;chg[i]=0;}

// ================= movimento (empuxo, gravidade, viscosidade) =================
function moveStep(){
  stamp=stamp>=255?1:stamp+1;
  const gOn=gm>0.02,gf=gm<1?gm:1;
  const ext=elecLive||magLive;
  for(let y=H-2;y>=0;y--){
    const row=y*W,ltr=rnd()<.5;
    for(let xx=1;xx<W-1;xx++){
      const x=ltr?xx:W-1-xx,i=row+x;
      if(mv[i]===stamp) continue;
      if(gLoc[i]&&gravMove(i,x,y)) continue; // gravidade local (grav.js)
      const m=mat[i],k=kind[m];
      if(k===FIX) continue;
      const j=i+W,mj=mat[j],kj=kind[mj];
      if(m===AR){ // ar: só afunda em gás mais leve (convecção / gases leves sobem)
        if(!gOn||kj!==GAS||mv[j]===stamp) continue;
        if(mj===AR&&T[j]-T[i]<2) continue;
        const ri=gd[i],rj=gd[j];
        if(ri>rj&&rnd()<(ri-rj)/ri*40*gf) swap(i,j);
        continue;
      }
      if(y===0) continue;
      // forças elétrica (qE) e magnética (M·∇|B|)
      if(ext&&k!==GAS&&(mgMov[m]||chg[i]>1e-11||chg[i]<-1e-11)&&extMove(i,m,k,x,y)) continue;
      if(m===GRANIZO&&storm&&y<CBASE&&rnd()<.8) continue; // corrente ascendente sustenta o graupel
      const ri=densOf(i);
      if(gOn&&mv[j]!==stamp&&(kj===GAS||kj===LIQ||((kj===POW||kj===MOV)&&(k===LIQ||k===GAS)))){
        const rj=densOf(j);
        if(ri>rj){
          let p=(ri-rj)/ri*40; if(p>1)p=1;
          if(k===LIQ)p*=mobOf(i,m); if(kj===LIQ)p*=mobOf(j,mj);
          if(rnd()<p*gf){swap(i,j);continue;}
        }
      }
      if(k===MOV) continue;
      const d=rnd()<.5?1:-1;
      if(k===GAS){
        const jj=i+d,mm=mat[jj];
        if(mv[jj]!==stamp&&kind[mm]===GAS&&mm!==m&&rnd()<.35) swap(i,jj);
        continue;
      }
      if(gOn){
        let done=false;
        for(let s=0;s<2&&!done;s++){
          const dd=s?-d:d,jd=j+dd,md=mat[jd],kd=kind[md];
          if(mv[jd]===stamp||(kd!==GAS&&kd!==LIQ)) continue;
          if(kind[mj]===FIX&&kind[mat[i+dd]]===FIX) continue; // não vaza por quinas
          if(ri<=densOf(jd)) continue;
          let p=k===POW?.65:mobOf(i,m); if(kd===LIQ)p*=mobOf(jd,md);
          if(rnd()<p){swap(i,jd);done=true;}
        }
        if(done) continue;
      }
      if(k===LIQ&&rnd()<mobOf(i,m)){
        const dis=dispA[m]; let last=-1;
        for(let s=1;s<=dis;s++){
          const jj=i+d*s,mm=mat[jj],kk=kind[mm];
          if(mv[jj]===stamp) break;
          if(kk===GAS){last=jj; if(kind[mat[jj+W]]===GAS) break; continue;}
          if(kk===LIQ&&mm!==m&&densOf(jj)<ri-1){last=jj;continue;}
          break;
        }
        if(last>=0) swap(i,last);
      }
    }
  }
}

// ================= pressão =================
// BFS a partir da atmosfera: P_j = P_i + ρ̄·g·Δh. Regiões não alcançadas = seladas, com equação de estado própria.
function pressureStep(){
  const hg=0.5*G0*gm*DX;
  for(let i=0;i<N;i++){
    if(kind[mat[i]]===FIX){vis[i]=3;Pn[i]=0;rhoA[i]=0;}
    else{vis[i]=0;rhoA[i]=densOf(i);}
  }
  let qh=0,qt=0;
  if(!lid) for(let x=1;x<W-1;x++){ if(vis[x]===0){vis[x]=1;Pn[x]=PATM;q[qt++]=x;} }
  qh=bfs(qh,qt,1,hg); qt=qh;
  // regiões seladas
  for(let s=0;s<N;s++){
    if(vis[s]!==0) continue;
    const start=qt; vis[s]=2; Pn[s]=0; q[qt++]=s;
    qt=bfs(start,qt,2,hg);
    let gasN=0,snT=0,sumHg=0,cnt=0,sumPrev=0,liqN=0,sumDT=0;
    for(let a=start;a<qt;a++){
      const i=q[a],m=mat[i],k=kind[m],h=Pn[i];
      if(P[i]>0){sumPrev+=P[i]-h;cnt++;}
      if(k===GAS){gasN++;snT+=gd[i]/molar[m]*T[i];sumHg+=h;}
      else if(k===LIQ){liqN++;sumDT+=bkA[m]*(T[i]-tOld[i]);}
    }
    let base;
    if(gasN>0) base=(RG*snT-sumHg)/gasN;                 // gás ideal, conserva mols
    else base=(cnt?sumPrev/cnt:PATM)+(liqN?sumDT/liqN:0); // líquido selado: pressão térmica β/κ·ΔT
    for(let a=start;a<qt;a++){const i=q[a];let p=base+Pn[i];if(p<1000)p=1000;Pn[i]=p;}
  }
  // pós: adiabático em gases abertos, densidade dos gases, aquecimento por compressão em líquidos
  for(let i=0;i<N;i++){
    const v=vis[i]; if(v===3) continue;
    const m=mat[i],k=kind[m];let p=Pn[i];
    if(p<1000){p=1000;Pn[i]=p;} // coluna sustentada por força externa (ímã) acima de bolsão aberto: hidrostática daria P<0
    if(k===GAS){
      if(v===1){const r=p/pp[i]; if(r>1.00002||r<0.99998) T[i]*=Math.pow(r,gexp[m]);}
      gd[i]=p*molar[m]/(RG*T[i]);
    }else if(k===LIQ&&v===1){
      T[i]+=T[i]*betaA[m]*(p-pp[i])/(CvA[m]*1e6);
    }
    pp[i]=p; tOld[i]=T[i];
  }
  const t=P;P=Pn;Pn=t;
}
function bfs(qh,qt,mark,hg){
  while(qh<qt){
    const i=q[qh++],pi=Pn[i],ri=rhoA[i],x=i%W;
    let j=i+W;
    if(j<N&&vis[j]===0){vis[j]=mark;Pn[j]=pi+(hg+gHy[i])*(ri+rhoA[j]);q[qt++]=j;}
    j=i-W;
    if(j>=0&&vis[j]===0){vis[j]=mark;Pn[j]=pi-(hg+gHy[i])*(ri+rhoA[j]);q[qt++]=j;}
    if(x>0){j=i-1;if(vis[j]===0){vis[j]=mark;Pn[j]=pi-gHx[i]*(ri+rhoA[j]);q[qt++]=j;}}
    if(x<W-1){j=i+1;if(vis[j]===0){vis[j]=mark;Pn[j]=pi+gHx[i]*(ri+rhoA[j]);q[qt++]=j;}}
  }
  return qt;
}

// ================= calor (condução conservativa) =================
function capOf(i,m){if(kind[m]===GAS){const c=gd[i]*cpgA[m]*1e-6;return c<1e-6?1e-6:c;}return CvA[m];}
function xfer(i,j,mi,a,Ci){
  const dT=T[j]-T[i];
  if(dT<.02&&dT>-.02) return;
  const mj=mat[j];
  if(mi===mj&&(mi===AR||mi===PAREDE)&&dT<.3&&dT>-.3) return;
  const b=kc[mj],G=2*a*b/(a+b);
  const Cj=capOf(j,mj);
  let E=G*dT; const Em=dT*Ci*Cj/(Ci+Cj)*.24;
  if(dT>0){if(E>Em)E=Em;}else if(E<Em)E=Em;
  T[i]+=E/Ci; T[j]-=E/Cj;
}
function heatStep(){
  const odd=stamp&1;
  for(let yy=0;yy<H;yy++){
    const y=odd?yy:H-1-yy,row=y*W;
    for(let x=0;x<W;x++){
      const i=row+x,mi=mat[i],a=kc[mi],Ci=capOf(i,mi);
      if(x<W-1) xfer(i,i+1,mi,a,Ci);
      if(y<H-1) xfer(i,i+W,mi,a,Ci);
    }
  }
  for(let x=0;x<W;x++){T[x+(H-1)*W]=293.15; if(lid||x===0||x===W-1) T[x]=293.15;}
  for(let y=0;y<H;y++){T[y*W]=293.15;T[y*W+W-1]=293.15;}
}

// ================= estado: fases, fogo, esmagamento =================
function hasOxy(i){const a=mat[i-1],b=mat[i+1],c=mat[i-W],d=mat[i+W];return a===AR||b===AR||c===AR||d===AR||a===FOGO||b===FOGO||c===FOGO||d===FOGO;}
function nbOf(i,r){return r<.25?i-1:r<.5?i+1:r<.75?i-W:i+W;}
function ignite(i){
  const src=mat[i];
  if(kind[src]!==GAS){ // combustível sólido/líquido queima na superfície
    if(life[i]===0){const L=flA[src]*(.6+.8*rnd());life[i]=L>250?250:(L|0)||1;}
    if(T[i]<900)T[i]=900; return;
  }
  mat[i]=FOGO; lat[i]=src;
  const L=flA[src]*(.6+.8*rnd()); life[i]=L>250?250:L|0;
  if(T[i]<1200)T[i]=1200;
  const p=P[i]>0?P[i]:PATM; pp[i]=p; gd[i]=p*molar[FOGO]/(RG*T[i]); // tOld mantido: calor da ignição expande o ar (airStep)
}
function burn(i,m){ // combustível em brasa: emite chama no ar adjacente, aquece, propaga
  let l=life[i];
  if(T[i]<850)T[i]=850+rnd()*250;
  let oxy=false;
  const nb=[i-1,i+1,i-W,i+W];
  for(let n=0;n<4;n++){
    const j=nb[n],mj=mat[j];
    if(mj===AR){oxy=true; if(rnd()<.35){setCell(j,FOGO,1150+rnd()*200);life[j]=8+(rnd()*18|0);lat[j]=0;}}
    else if(mj===FOGO) oxy=true;
    else if(ignA[mj]>0&&life[j]===0&&rnd()<fspA[mj]) ignite(j);
    else if(mj===AGUA||mj===SALMOURA){ l=l>10?l-10:0; }
    else if(kind[mj]!==GAS&&mj!==PAREDE){const dT=FIRE_Q/CvA[mj];T[j]+=dT>15?15:dT;}
  }
  if(l>0&&rnd()<(oxy?.6:.08)) l--;
  if(l<=0){
    life[i]=0;
    if(!oxy&&kind[m]===MOV){T[i]-=200;return;} // abafado: apaga
    if((m===MADEIRA||m===MADEIRAE)&&rnd()<.2){mat[i]=CINZA;lat[i]=0;}
    else {mat[i]=FOGO;lat[i]=m;life[i]=10+(rnd()*20|0);gd[i]=P[i]*molar[FOGO]/(RG*T[i]);}
  } else life[i]=l;
}
function breakFixed(i,m){
  const to=crushA[m]; mat[i]=to; lat[i]=0; st[i]=0; T[i]+=3;
  const p=P[i]>0?P[i]:PATM; P[i]=p; pp[i]=p; tOld[i]=T[i];
}
function fire(i){
  let l=life[i];
  if(T[i]<1150) T[i]=1150+rnd()*250;
  const nb=[i-1,i+1,i-W,i+W];
  for(let n=0;n<4;n++){
    const j=nb[n],mj=mat[j];
    if(ignA[mj]>0){ if(life[j]===0&&rnd()<fspA[mj]*(T[j]>ignA[mj]*.8?2:1)) ignite(j); }
    else if(mj===AGUA||mj===SALMOURA) l=l>8?l-8:0;
    else if(kind[mj]!==GAS&&mj!==PAREDE){ const dT=FIRE_Q/CvA[mj]; T[j]+=dT>15?15:dT; }
  }
  if(rnd()<.7&&l>0) l--;
  if(l<=0){
    const src=lat[i]; lat[i]=0; life[i]=0;
    if((src===MADEIRA||src===MADEIRAE)&&rnd()<.15){mat[i]=CINZA;}
    else if(rnd()<.45){mat[i]=FUMACA;life[i]=120+(rnd()*100|0);}
    else mat[i]=AR;
    if(kind[mat[i]]===GAS) gd[i]=P[i]*molar[mat[i]]/(RG*T[i]);
  } else life[i]=l;
}
// ================= vida (Conway B3/S23) acoplada ao mundo =================
// nasce só em ar habitável; morre por calor, frio, fogo (queima), raio/plasma, choque elétrico, afogamento e pressão (str)
const LIFE_TMIN=268.15,LIFE_TMAX=333.15,LIFE_IMAX=2e-2,lifeN=new Uint8Array(N);
let lifeEvery=3,lifeTick=0;
function lifeKill(i){mat[i]=CINZA;lat[i]=0;life[i]=0;st[i]=0;const p=P[i]>0?P[i]:PATM;P[i]=p;pp[i]=p;tOld[i]=T[i];}
function lifeHazard(i){
  const t=T[i];
  if(t>ignA[VIDA]&&hasOxy(i)){ignite(i);return true;}
  if(t>LIFE_TMAX||t<LIFE_TMIN||ionL[i]>.05||(IxC[i]<0?-IxC[i]:IxC[i])+(IyC[i]<0?-IyC[i]:IyC[i])>LIFE_IMAX){lifeKill(i);return true;}
  if(rnd()<.05){let l=0;if(kind[mat[i-1]]===LIQ)l++;if(kind[mat[i+1]]===LIQ)l++;if(kind[mat[i-W]]===LIQ)l++;if(kind[mat[i+W]]===LIQ)l++;
    if(l>=3){lifeKill(i);return true;}} // afogada
  return radLife(i); // dose de radiação: morte e mutação (rad.js)
}
// mutante (rad.js) segue HighLife B36/S23; vizinhos contam as duas espécies (nibble alto = nº de mutantes)
function lifeStep(){
  if(++lifeTick<lifeEvery) return; lifeTick=0;
  lifeN.fill(0); let any=false;
  for(let y=1;y<H-1;y++){const row=y*W;
    for(let x=1;x<W-1;x++){const i=row+x,m=mat[i]; if((m!==VIDA&&m!==VIDAM)||life[i]>0) continue; any=true; // em chamas não conta
      const a=m===VIDAM?17:1;
      lifeN[i-W-1]+=a;lifeN[i-W]+=a;lifeN[i-W+1]+=a;lifeN[i-1]+=a;lifeN[i+1]+=a;lifeN[i+W-1]+=a;lifeN[i+W]+=a;lifeN[i+W+1]+=a;}}
  if(!any) return;
  for(let y=1;y<H-1;y++){const row=y*W;
    for(let x=1;x<W-1;x++){const i=row+x,v=lifeN[i],n=v&15,nm=v>>4,m=mat[i];
      if(m===VIDA||m===VIDAM){ if(life[i]===0&&(n<2||n>3)) setCell(i,AR,T[i]); }
      else if((n===3||(n===6&&nm>=3))&&m===AR&&T[i]>LIFE_TMIN&&T[i]<LIFE_TMAX&&ionL[i]<.01&&radF[i]<RAD_STERILE)
        setCell(i,nm>=2||n===6?VIDAM:VIDA,T[i]);
    }}
}
const acidR=new Uint8Array(NM); // 1 metal, 2 orgânico, 3 álcool, 4 alcalino (concreto), 5 água (diluição)
for(const m of [ACO,LATA,COBRE,SUCATA,LIMALHA,IMA,MERCURIO])acidR[m]=1;
for(const m of [MADEIRA,MADEIRAE,ISOPOR,ISOPORC,OLEO,MEL,VIDA,VIDAM])acidR[m]=2;
acidR[ALCOOL]=3;acidR[CONCRETO]=4;acidR[AGUA]=5;acidR[SALMOURA]=5;
function stateStep(){
  for(let y=1;y<H-1;y++){
    const row=y*W;
    for(let x=1;x<W-1;x++){
      const i=row+x,m=mat[i];
      if(!actA[m]) continue;
      if(m===FOGO){fire(i);continue;}
      if(m===RESIST){T[i]=heaterT;continue;}
      if(exW[m]>0){ // explosivo: calor, chama, faísca/plasma ou choque (nitro) detonam
        let tr=T[i]>exTd[m]||ion[i]>.3||ionL[i]>.05;
        if(!tr){const a=mat[i-1],b=mat[i+1],c=mat[i-W],d=mat[i+W];tr=a===FOGO||b===FOGO||c===FOGO||d===FOGO;}
        if(!tr&&exSh[m]>0){const pa=ap[AF[i]];tr=pa>exSh[m]||pa<-exSh[m];}
        if(tr) detonate(i,0);
        continue;
      }
      if(m===IMA&&T[i]>tcA[m]&&(mxA[i]||myA[i])){mxA[i]=0;myA[i]=0;} // acima de Curie: desmagnetiza (irreversível)
      if(m===FUMACA){ if(rnd()<.5){ if(life[i]>0) life[i]--; else {mat[i]=AR;gd[i]=P[i]*molar[AR]/(RG*T[i]);} } continue; }
      const k=kind[m],p=P[i]+ap[AF[i]]; // + sobrepressão dinâmica (ar)
      if(m===SAL){
        const j=nbOf(i,rnd()),mj=mat[j];
        if(mj===AGUA&&rnd()<.05){mat[j]=SALMOURA;mat[i]=SALMOURA;lat[i]=lat[j]=0;T[i]-=1;T[j]-=1;continue;}
        if(mj===GELO&&rnd()<.02){mat[j]=SALMOURA;mat[i]=SALMOURA;lat[i]=lat[j]=0;T[i]-=4;T[j]-=4;continue;}
      }
      else if(m===MADEIRA){ // encharca sob pressão
        if(p-PATM>150e3&&rnd()<.004){const j=nbOf(i,rnd());if(mat[j]===AGUA||mat[j]===SALMOURA){mat[i]=MADEIRAE;continue;}}
      }
      else if(m===HNO3){ // oxidante forte: dissolve metais (NO₂), inflama orgânicos, reage de forma explosiva com álcool
        const j=nbOf(i,rnd()),mj=mat[j],r=acidR[mj];
        if(r===1&&rnd()<.04){setCell(j,FUMACA,T[j]+40);life[j]=220;T[i]+=25;if(rnd()<.35){mat[i]=AGUA;lat[i]=0;}continue;}
        if(r===2&&rnd()<.05){T[j]+=120;T[i]+=40;if(ignA[mj]>0&&hasOxy(j))ignite(j);if(rnd()<.3){mat[i]=AGUA;lat[i]=0;}continue;}
        if(r===3&&rnd()<.1){detonate(i,.4);detonate(j,.4);continue;}
        if(r===4&&rnd()<.006){breakFixed(j,mj);T[i]+=5;if(rnd()<.5){mat[i]=AGUA;lat[i]=0;}continue;}
        if(r===5&&rnd()<.02){T[i]+=3;T[j]+=3;}
      }
      else if(m===ISOPOR){ // espuma de células fechadas colapsa
        if(p-PATM>strA[m]&&rnd()<.06){
          T[i]+=4;
          if(rnd()<.15) mat[i]=ISOPORC;
          else{const j=nbOf(i,rnd()); if(kind[mat[j]]===LIQ) copyFrom(i,j); else mat[i]=ISOPORC;}
          continue;
        }
      }
      if((m===VIDA||m===VIDAM)&&life[i]===0&&lifeHazard(i)) continue;
      if(ignA[m]>0){ if(life[i]>0&&kind[m]!==GAS){burn(i,m);continue;} if(T[i]>ignA[m]&&hasOxy(i)){ignite(i);continue;} }
      if(k===FIX&&strA[m]>0){
        let mx=-1,mn=1e15,c=0;
        let j=i-1; if(vis[j]!==3){const q1=P[j]+apAt(j,-1);if(q1>mx)mx=q1;if(q1<mn)mn=q1;c++;}
        j=i+1; if(vis[j]!==3){const q1=P[j]+apAt(j,1);if(q1>mx)mx=q1;if(q1<mn)mn=q1;c++;}
        j=i-W; if(vis[j]!==3){const q1=P[j]+apAt(j,-W);if(q1>mx)mx=q1;if(q1<mn)mn=q1;c++;}
        j=i+W; if(vis[j]!==3){const q1=P[j]+apAt(j,W);if(q1>mx)mx=q1;if(q1<mn)mn=q1;c++;}
        const s=c>=2?(mx-mn)/strA[m]:0; st[i]=s;
        if(s>1&&rnd()<.2){breakFixed(i,m);continue;}
        if(shockA[m]>0){ // choque térmico
          const sh=shockA[m],t=T[i];
          if((mat[i-1]===m&&Math.abs(T[i-1]-t)>sh)||(mat[i+1]===m&&Math.abs(T[i+1]-t)>sh)||
             (mat[i-W]===m&&Math.abs(T[i-W]-t)>sh)||(mat[i+W]===m&&Math.abs(T[i+W]-t)>sh)){
            if(rnd()<.03){breakFixed(i,m);continue;}
          }
        }
      }
      // mudanças de fase com calor latente (lat>0 rumo a hi, lat<0 rumo a lo)
      const ht=hiTy[m];
      if(ht){
        const th=ht===1?hiT[m]:ht===2?tb(hiT[m],hiLm[m],hiCap[m],p):hiT[m]+hiLm[m]*(p-PATM);
        const t=T[i];
        if(t>th){const l=lat[i]+t-th;T[i]=th;if(l>=hiL[m]){T[i]=th;convert(i,hiTo[m]);continue;}lat[i]=l;}
        else if(lat[i]>0){const dd=Math.min(lat[i],th-t);lat[i]-=dd;T[i]+=dd;}
      }
      const lt=loTy[m];
      if(lt){
        const th=lt===1?loT[m]:lt===2?tb(loT[m],loLm[m],loCap[m],p):loT[m]+loLm[m]*(p-PATM);
        const t=T[i];
        if(t<th){const l=lat[i]-(th-t);T[i]=th;if(-l>=loL[m]){convert(i,loTo[m]);continue;}lat[i]=l;}
        else if(lat[i]<0){const dd=Math.min(-lat[i],t-th);lat[i]+=dd;T[i]-=dd;}
      }
    }
  }
  if(!lid){ // topo aberto: reservatório atmosférico
    const g0=PATM*molar[AR]/(RG*293.15);
    for(let x=1;x<W-1;x++){ if(mat[x]!==AR) mat[x]=AR; T[x]=293.15; gd[x]=g0; pp[x]=PATM; lat[x]=0; chg[x]=0; ion[x]=0; }
  }
}

// ================= ferramentas =================
let brush=4,tool=AGUA,drawing=false,erasing=false,shiftKey=false,curX=-1,curY=-1,lastX=-1,lastY=-1,hoverI=-1;
let brushSq=false; // pincel quadrado
function paintAt(cx,cy,tl){
  if(tl===GRAVT){gravPaint(cx,cy,shiftKey);return;}
  if(tl===CYCL){airCyclone(cx,cy,shiftKey);return;} if(tl===WIND){airWind(cx,cy);return;}
  const r=brush,r2=r*r+r*.8;
  for(let dy=-r;dy<=r;dy++)for(let dx=-r;dx<=r;dx++){
    if(!brushSq&&dx*dx+dy*dy>r2) continue;
    const x=cx+dx,y=cy+dy;
    if(x<1||x>W-2||y<1||y>H-2) continue;
    paintCell(y*W+x,tl,false);
  }
}
// aplica ferramenta/material numa célula; solid: ignora preenchimento aleatório (formas, balde)
function paintCell(i,tl,solid){
    const cm=mat[i];
    if(cm===PAREDE||tl===WIND||tl===CYCL||tl===GRAVT) return;
    if(tl===ERASE){ if(cm!==AR) setCell(i,AR,293.15); else chg[i]=0; return; }
    if(tl===CHGP||tl===CHGN){ if(!srcA[cm]) chg[i]+=tl===CHGP?qAmt:-qAmt; return; }
    if(tl===SPARK){ ion[i]=1; return; }
    if(tl===FLIP){ // inverte carga, polo do ímã e polaridade da fonte (uma vez por traço)
      if(flipMk[i]===strokeId) return; flipMk[i]=strokeId;
      chg[i]=-chg[i];
      if(cm===IMA){mxA[i]=-mxA[i];myA[i]=-myA[i];}
      else if(cm===FONTEP) mat[i]=FONTEN; else if(cm===FONTEN) mat[i]=FONTEP;
      return;
    }
    if(tl===HEAT||tl===COOL){
      const C=capOf(i,cm); let dT=power/C; if(dT>80)dT=80;
      if(tl===HEAT){T[i]+=dT;if(T[i]>4000)T[i]=4000;} else {T[i]-=dT;if(T[i]<2)T[i]=2;}
      return;
    }
    const k=kind[tl],ck=kind[cm];
    if(cm===tl) return;
    if(k===FIX){ if(!shiftKey&&ck===FIX) return; if((tl===VIDA||tl===VIDAM)&&!solid&&rnd()>.3) return; }
    else { if(!shiftKey&&ck!==GAS) return; if(!solid&&rnd()>fillA[tl]) return; }
    setCell(i,tl,tl===RESIST?heaterT:T0A[tl]);
    if(tl===IMA) magnetize(i,magDir);
}
// formas: retângulo (cheio/contorno), elipse, linha; balde = preenche região conexa de mesmo material
function paintShape(md,x0,y0,x1,y1,tl){
  if(md==='line'){paintLine(x0,y0,x1,y1,tl);return;}
  const xa=Math.max(1,Math.min(x0,x1)),xb=Math.min(W-2,Math.max(x0,x1)),ya=Math.max(1,Math.min(y0,y1)),yb=Math.min(H-2,Math.max(y0,y1));
  const cx=(x0+x1)/2,cy=(y0+y1)/2,rx=Math.abs(x1-x0)/2+.5,ry=Math.abs(y1-y0)/2+.5,th=brush+1;
  for(let y=ya;y<=yb;y++)for(let x=xa;x<=xb;x++){
    if(md==='rectO'&&x-Math.min(x0,x1)>=th&&Math.max(x0,x1)-x>=th&&y-Math.min(y0,y1)>=th&&Math.max(y0,y1)-y>=th)continue;
    if(md==='ell'){const u=(x-cx)/rx,v=(y-cy)/ry;if(u*u+v*v>1)continue;}
    paintCell(y*W+x,tl,true);
  }
}
function bucket(cx,cy,tl){
  if(cx<1||cx>W-2||cy<1||cy>H-2)return;const s0=cy*W+cx,m=mat[s0];if(m===PAREDE)return;
  const seen=new Uint8Array(N),qq=new Int32Array(N);let h=0,t=0;qq[t++]=s0;seen[s0]=1;
  while(h<t){const c=qq[h++],x=c%W,y=c/W|0;
    for(const n of [c-1,c+1,c-W,c+W]){const nx=n%W,ny=n/W|0;if(nx<1||nx>W-2||ny<1||ny>H-2||seen[n]||mat[n]!==m)continue;seen[n]=1;qq[t++]=n;}
    paintCell(c,tl,true);}
}
const MDX=[0,1,0,-1],MDY=[-1,0,1,0]; // direção do polo N: ↑ → ↓ ←
function magnetize(i,d){mxA[i]=MDX[d]*msA[IMA];myA[i]=MDY[d]*msA[IMA];}
function paintLine(x0,y0,x1,y1,tl){
  if(tl===WIND){airWindDir(x1-x0,y1-y0);airWind(x1,y1);return;}
  const n=Math.max(Math.abs(x1-x0),Math.abs(y1-y0),1);
  const stepN=Math.max(1,Math.floor(brush/2));
  for(let s=0;s<=n;s+=stepN) paintAt(Math.round(x0+(x1-x0)*s/n),Math.round(y0+(y1-y0)*s/n),tl);
  paintAt(x1,y1,tl);
}

// ================= cenários =================
function rectM(x0,y0,x1,y1,m,t){for(let y=y0;y<=y1;y++)for(let x=x0;x<=x1;x++){if(x<1||x>W-2||y<1||y>H-2)continue;setCell(y*W+x,m,t===undefined?(m===RESIST?heaterT:T0A[m]):t);}}
function boxM(x0,y0,x1,y1,th,m){rectM(x0,y0,x1,y0+th-1,m);rectM(x0,y1-th+1,x1,y1,m);rectM(x0,y0,x0+th-1,y1,m);rectM(x1-th+1,y0,x1,y1,m);}
function magRect(x0,y0,x1,y1,d){rectM(x0,y0,x1,y1,IMA);for(let y=y0;y<=y1;y++)for(let x=x0;x<=x1;x++)magnetize(y*W+x,d);}
function clearAll(){
  P.fill(PATM); emReset(); airReset(); gravReset(); radReset(); rA.fill(0); bu.fill(0); cmp.fill(0); dAcc.fill(0);
  for(let i=0;i<N;i++) setCell(i,AR,293.15);
  for(let x=0;x<W;x++){setCell(x+(H-1)*W,PAREDE,293.15); if(lid) setCell(x,PAREDE,293.15);}
  for(let y=0;y<H;y++){setCell(y*W,PAREDE,293.15);setCell(y*W+W-1,PAREDE,293.15);}
}
function updLidBtn(){const bl=document.getElementById('bLid');bl.textContent='Tampa '+(lid?'fechada':'aberta');bl.classList.toggle('on',lid);}
function setLid(v){
  lid=v; updLidBtn();
  for(let x=1;x<W-1;x++) setCell(x,v?PAREDE:AR,293.15);
}
const PRESETS={
  tank(){
    clearAll();
    rectM(1,182,198,198,MERCURIO); rectM(1,45,198,181,AGUA); rectM(1,36,198,44,OLEO);
    rectM(18,24,52,28,MADEIRA); rectM(150,20,162,31,GELO);
    rectM(22,156,46,172,ISOPOR);
    boxM(78,138,100,174,1,LATA); rectM(79,139,99,173,AR);
    boxM(128,150,162,180,1,VIDRO); rectM(129,151,161,179,AR);
    rectM(170,120,180,130,MADEIRA);
  },
  cooker(){
    clearAll(); heaterSet(500);
    boxM(30,70,112,160,3,ACO); rectM(33,112,109,157,AGUA);
    rectM(34,161,108,163,RESIST);
    boxM(132,104,176,160,2,VIDRO); rectM(134,128,174,158,AGUA);
    rectM(136,161,172,163,RESIST);
    rectM(1,180,198,198,AGUA);
  },
  fire(){
    clearAll();
    rectM(1,160,118,198,AGUA); rectM(20,150,60,159,GELO); rectM(30,18,50,34,LAVA);
    rectM(70,40,100,50,N2L);
    rectM(120,140,122,198,VIDRO); rectM(196,140,198,198,VIDRO);
    rectM(123,170,195,198,ALCOOL); rectM(135,160,180,166,MADEIRA); rectM(150,158,158,160,MADEIRA,720);
  },
  bell(){
    clearAll();
    rectM(1,40,198,198,AGUA);
    // sino invertido (aço) preso por haste: ar aprisionado comprime
    rectM(60,40,140,43,ACO); rectM(60,44,63,150,ACO); rectM(137,44,140,150,ACO); rectM(64,44,136,150,AR);
    rectM(98,1,102,39,ACO);
    rectM(80,186,120,190,ISOPOR); rectM(150,175,160,185,HELIO);
  },
  storm(){ // tempestade: nuvem se forma, graupel carrega, raio procura o caminho de maior campo
    clearAll(); setStorm(true);
    rectM(1,186,198,198,AREIA); rectM(110,180,170,198,AGUA);
    rectM(28,150,31,185,MADEIRA); rectM(18,136,42,152,MADEIRA);        // árvore
    rectM(150,148,175,179,CONCRETO); rectM(161,128,163,185,ACO);        // prédio + para-raios até o solo
    rectM(161,186,163,198,ACO); rectM(156,196,168,198,TERRA);
    rectM(70,176,90,185,LATA);
  },
  circuit(){ // fonte, fios de cobre, nicromo aquecendo, eletrólito e limalha ao redor do fio
    clearAll(); setVolt(150);
    rectM(20,40,24,52,FONTEP); rectM(20,56,24,68,FONTEN);
    rectM(25,44,180,46,COBRE); rectM(25,62,90,64,COBRE);
    rectM(178,47,180,150,COBRE);
    rectM(110,148,177,150,COBRE); rectM(90,148,109,150,NICROMO);      // nicromo: maior resistência → esquenta
    rectM(60,148,89,150,COBRE); rectM(58,65,60,150,COBRE); rectM(61,62,90,64,AR);
    rectM(58,62,60,64,COBRE);
    boxM(120,110,170,140,1,VIDRO); rectM(121,111,169,139,SALMOURA);   // eletrólito em paralelo
    rectM(130,47,132,112,COBRE); rectM(158,112,160,147,COBRE);
    rectM(60,90,120,110,LIMALHA);
  },
  magnets(){ // ímãs, limalha mostrando o campo, ferrofluido e barra de aço
    clearAll();
    magRect(40,60,70,72,1); magRect(130,60,160,72,3);                // polos N frente a frente: repelem campo
    magRect(92,166,108,181,0);
    boxM(76,184,124,198,1,VIDRO); rectM(77,187,123,197,FERROFL);
    rectM(30,120,70,124,ACO);
    rectM(20,20,180,34,LIMALHA);
  },
  blast(){ // explosivos: bunker de concreto com janelas de vidro, C-4 com detonador, nitro, ácido × álcool, areia e isopor
    clearAll();
    rectM(1,170,198,198,AREIA);
    boxM(20,110,80,169,4,CONCRETO); rectM(20,128,23,140,VIDRO); rectM(77,128,80,140,VIDRO);
    rectM(47,160,53,165,C4); rectM(50,154,50,159,MADEIRA); rectM(50,153,50,153,MADEIRA,720); // pavio de madeira aceso
    boxM(110,130,140,169,2,VIDRO); rectM(112,150,138,167,NITRO);                        // frasco de nitro (sensível ao choque)
    rectM(150,120,195,123,CONCRETO); rectM(160,100,185,119,ISOPOR);
    boxM(150,150,196,169,2,ACO); rectM(152,160,172,167,HNO3); rectM(174,160,194,167,ALCOOL); rectM(173,152,173,167,VIDRO); // divisória de vidro
    rectM(90,20,110,40,MADEIRA);
  },
  life(){ // jogo da vida: canhão de planadores rumo a aquecedor, sopa entre eletrodos, lago (afogamento) e fogueira
    clearAll(); heaterSet(400); setVolt(30);
    const gun=['........................O...........','......................O.O...........','............OO......OO............OO',
      '...........O...O....OO............OO','OO........O.....O...OO..............','OO........O...O.OO....O.O...........',
      '..........O.....O.......O...........','...........O...O....................','............OO......................'];
    gun.forEach((r,y)=>{for(let x=0;x<r.length;x++)if(r[x]==='O')rectM(4+x,6+y,4+x,6+y,VIDA);}); // Gosper → planadores para SE
    rectM(150,118,152,160,RESIST);                                     // parede quente: planadores morrem no ar aquecido
    rectM(1,184,198,198,AGUA);                                         // lago
    rectM(10,120,12,170,FONTEP); rectM(70,120,72,170,FONTEN);           // eletrodos: colônia que os conectar leva choque
    for(let y=122;y<=168;y++)for(let x=13;x<=69;x++)if(rnd()<.35)rectM(x,y,x,y,VIDA); // sopa aleatória
    rectM(100,150,120,166,VIDA); rectM(104,174,116,183,MADEIRA); rectM(108,172,112,173,MADEIRA,720); // colônia sobre fogueira
  },
  tesla(){ // eletrodo de alta tensão: descargas buscam a ponta mais próxima/condutora
    clearAll(); setVolt(1e7);
    rectM(80,30,120,40,FONTEP); rectM(95,41,105,46,COBRE); rectM(100,47,100,62,COBRE); // agulha: ponta concentra o campo
    rectM(30,150,33,198,ACO); rectM(100,175,130,198,CONCRETO); rectM(150,140,170,198,MADEIRA);
    rectM(1,190,29,198,AGUA); rectM(171,192,198,198,AREIA);
  },
  teslacoil(){ // bobina de Tesla: toroide em alta tensão; espiras isoladas do secundário assumem o gradiente de potencial
    // ⚠️ assumido: simulação quase-estática (sem RF/ressonância) → toroide alimentado por fonte DC equivalente ao pico
    clearAll(); setVolt(1e7);
    rectM(1,188,198,198,CONCRETO); rectM(90,194,110,198,TERRA);         // piso + haste de aterramento
    rectM(76,176,124,187,CONCRETO);                                    // base
    for(let k=0;k<4;k++){const r=16+k*3;rectM(100-r,173-k,100-r+1,174-k,COBRE);rectM(100+r-1,173-k,100+r,174-k,COBRE);} // primário cônico
    rectM(96,80,104,175,VIDRO);                                       // forma do secundário (tubo isolante)
    for(let y=82;y<=174;y+=2){rectM(95,y,95,y,COBRE);rectM(105,y,105,y,COBRE);} // espiras (seção) isoladas entre si
    rectM(105,175,106,175,COBRE); rectM(106,176,107,193,COBRE);        // base do secundário → terra
    rectM(99,74,101,79,COBRE);                                        // ligação topo do secundário
    // toroide (seção: dois lóbulos + disco central). Núcleo = fonte HV, casca de alumínio: a corrente de carga
    // entra pela superfície toda (fonte pequena/linha concentra a corrente → Joule vaporiza o metal)
    const tor=(x,y,d)=>{const a=(x-81)/(15-d),b=(x-119)/(15-d),c=(y-70)/(8-d);
      return a*a+c*c<=1||b*b+c*c<=1||(x>=81&&x<=119&&y>=66+d&&y<=74-d);};
    for(let y=60;y<=80;y++)for(let x=66;x<=134;x++)if(tor(x,y,0))rectM(x,y,x,y,tor(x,y,2)?FONTEP:LATA);
    rectM(135,70,146,70,COBRE);                                       // ponta de descarga (breakout point)
    rectM(166,110,168,187,ACO); rectM(163,108,171,110,ACO);            // haste aterrada (alvo)
    boxM(14,140,48,187,1,ACO); rectM(26,176,36,187,MADEIRA);           // gaiola de Faraday com madeira dentro
    rectM(132,180,142,187,MADEIRA);                                    // madeira no chão: pode inflamar
  }
};

// ================= renderização =================
const cv=document.getElementById('cv'),ctx=cv.getContext('2d',{alpha:false});
const img=ctx.createImageData(W,H),u32=new Uint32Array(img.data.buffer);
const pack=(r,g,b)=>(255<<24)|(b<<16)|(g<<8)|r;
const BG=D[AR].col;
const pal=new Uint32Array(NM*8);
for(let m=0;m<NM;m++){const d=D[m],v=d.var||0,a=d.a===undefined?1:d.a;
  for(let s=0;s<8;s++){const o=(s/7-.5)*v;
    let r=d.col[0]+o,g=d.col[1]+o,b=d.col[2]+o;
    if(d.kind===GAS){r=BG[0]+(r-BG[0])*a;g=BG[1]+(g-BG[1])*a;b=BG[2]+(b-BG[2])*a;}
    pal[m*8+s]=pack(Math.max(0,Math.min(255,r))|0,Math.max(0,Math.min(255,g))|0,Math.max(0,Math.min(255,b))|0);}}
function lerpStops(stops,x){
  if(x<=stops[0][0]) return stops[0][1];
  for(let s=1;s<stops.length;s++){ if(x<=stops[s][0]){const a=stops[s-1],b=stops[s],t=(x-a[0])/(b[0]-a[0]);
    return [a[1][0]+(b[1][0]-a[1][0])*t,a[1][1]+(b[1][1]-a[1][1])*t,a[1][2]+(b[1][2]-a[1][2])*t];}}
  return stops[stops.length-1][1];
}
const PST=[[0,[28,52,255]],[.25,[0,165,255]],[.5,[40,215,130]],[.75,[255,215,40]],[1,[255,36,28]]];
const plut=new Uint32Array(256); for(let i=0;i<256;i++){const c=lerpStops(PST,i/255);plut[i]=pack(c[0]|0,c[1]|0,c[2]|0);}
const TST=[[0,[30,0,70]],[77,[110,40,255]],[200,[40,70,255]],[273,[30,170,255]],[293,[22,34,46]],[330,[120,40,96]],[373,[230,40,40]],[600,[255,130,20]],[1200,[255,225,80]],[2500,[255,255,255]]];
const tlut=new Uint32Array(4001); for(let t=0;t<=4000;t++){const c=lerpStops(TST,t);tlut[t]=pack(c[0]|0,c[1]|0,c[2]|0);}
const GST=[[700,[90,0,0]],[900,[190,30,0]],[1200,[255,100,10]],[1600,[255,175,70]],[2200,[255,230,160]],[3000,[255,250,240]]];
const glR=new Uint8Array(4001),glG=new Uint8Array(4001),glB=new Uint8Array(4001),glA=new Float32Array(4001);
for(let t=0;t<=4000;t++){const c=lerpStops(GST,t);glR[t]=c[0];glG[t]=c[1];glB[t]=c[2];let a=(t-700)/700;glA[t]=a<0?0:a>.95?.95:a;}
function blend(c1,c2,a){
  const r1=c1&255,g1=(c1>>>8)&255,b1=(c1>>>16)&255,r2=c2&255,g2=(c2>>>8)&255,b2=(c2>>>16)&255;
  return pack((r1+(r2-r1)*a)|0,(g1+(g2-g1)*a)|0,(b1+(b2-b1)*a)|0);
}
function dim(c,f){return pack(((c&255)*f)|0,(((c>>>8)&255)*f)|0,(((c>>>16)&255)*f)|0);}
const BGC=pack(BG[0],BG[1],BG[2]),RED=pack(255,40,30);
let pRange=20000,cntP=0,maxT=0,maxG=0;
const PLC=pack(228,230,255),VPOS=pack(255,78,40),VNEG=pack(40,124,255),YEL=pack(255,226,90),QPOS=pack(255,186,160),QNEG=pack(160,200,255),
  NCOL=pack(228,48,44),SCOL=pack(46,98,232),MCOL=pack(96,96,108);
const MST=[[0,[0,0,4]],[.25,[80,18,123]],[.5,[183,55,121]],[.75,[251,135,97]],[1,[252,253,191]]];
const mlut=new Uint32Array(256); for(let i=0;i<256;i++){const c=lerpStops(MST,i/255);mlut[i]=pack(c[0]|0,c[1]|0,c[2]|0);}
const plasmaL=new Int32Array(3000); let plasmaN=0,airRange=5000;
const AIRP=pack(255,70,50),AIRN=pack(60,130,255),AIRV=pack(40,90,110);
// dose: log10 1e-7..1e4 Sv/h (fundo natural ~1e-7 · letal em horas >1 Sv/h)
const RST=[[0,[14,17,24]],[.2,[20,60,90]],[.4,[30,150,110]],[.55,[200,220,60]],[.7,[255,150,30]],[.85,[240,40,40]],[1,[255,220,255]]];
const rlut=new Uint32Array(256); for(let i=0;i<256;i++){const c=lerpStops(RST,i/255);rlut[i]=pack(c[0]|0,c[1]|0,c[2]|0);}
const CHER=pack(90,150,255),RAGL=pack(150,255,170);
function poleC(i,c){const t=poleT[i]; if(t<0) return c; // polo N vermelho, S azul
  return t>.5?blend(MCOL,NCOL,Math.min(1,(t-.5)*2.4)):blend(MCOL,SCOL,Math.min(1,(.5-t)*2.4));}
function render(){
  const v=view,pr=pRange,iar=1/airRange; let mg=0,cnt=0,mt=0;
  const ivs=1/vScale,ith=Math.max(1e-3,iMax*1e-4),qth=1e-9; plasmaN=0;
  for(let i=0;i<N;i++){
    const m=mat[i],k=kind[m];
    let c=pal[(m<<3)|(shade[i]&7)];
    const t=T[i];
    if(m!==AR&&m!==PAREDE){cnt++; if(t>mt)mt=t;}
    if(k===LIQ){const g=P[i]-PATM+ap[AF[i]]; if(g>mg)mg=g;}
    if(v===0){
      if(m===IMA) c=poleC(i,c);
      if(m===FOGO){const f=(life[i]&7)/7; c=pack(255,(110+f*120)|0,(20+f*40)|0);}
      else if(t>700&&k!==GAS){const ti=t>4000?4000:t|0,a=glA[ti];c=blend(c,pack(glR[ti],glG[ti],glB[ti]),a);}
      else if(m===AR&&t>330){let a=(t-330)/900;if(a>.5)a=.5;c=blend(c,pack(150,60,20),a);}
      if(radLive){const d=radF[i]; // Cherenkov na água sob fluxo intenso; radioluminescência do rádio
        if(d>2&&(m===AGUA||m===SALMOURA)){let a=Math.log10(d/2)*.25+.15;if(a>.8)a=.8;c=blend(c,CHER,a);}}
      if(m===RADIO) c=blend(c,RAGL,.25+.15*Math.random());
    }else if(v===6){ // taxa de dose (Sv/h, escala log)
      const d=radF[i];let r=d>1e-7?(Math.log10(d)+7)/11:0;if(r>1)r=1;
      const rc=rlut[(r*255)|0];
      c=k===GAS?rc:blend(dim(c,.45),rc,r>0?.6:0);
      if(GAM[m]||BET[m]||ALF[m]||rA[i]>.05) c=blend(c,RAGL,.35);
    }else if(v===3){ // potencial elétrico (divergente), corrente e carga
      const vv=Vf[i]*ivs; let s=vv<0?-vv:vv; s=s>1?1:Math.sqrt(s);
      const vc=blend(BGC,vv>=0?VPOS:VNEG,s);
      c=k===GAS?vc:blend(dim(c,.5),vc,.55);
      const ia=(IxC[i]<0?-IxC[i]:IxC[i])+(IyC[i]<0?-IyC[i]:IyC[i]);
      if(ia>ith){let a=Math.log10(ia/ith)*.2+.3;if(a>.9)a=.9;c=blend(c,YEL,a);}
      const qq=chg[i]; if(qq>qth||qq<-qth){let a=Math.log10((qq<0?-qq:qq)/qth)*.25+.45;if(a>1)a=1;c=blend(c,qq>0?QPOS:QNEG,a);}
    }else if(v===4){ // |B| (escala log) com polos
      const b=Bfine[i]; let tb=b>1e-5?(Math.log10(b)+5)*.2:0; if(tb>1)tb=1;
      const mc=mlut[(tb*255)|0];
      if(m===IMA) c=blend(poleC(i,c),mc,.12);
      else if(mgA[m]) c=blend(c,mc,.3);
      else c=k===GAS?mc:blend(dim(c,.4),mc,.6);
    }else if(v===1){
      if(k===FIX){
        c=dim(c,.55); const s=st[i]; if(s>.02) c=blend(c,RED,s>1?.9:s*.9);
      }else{
        const g=P[i]-PATM+ap[AF[i]]; let ix=g<=0?0:(g/pr*255)|0; if(ix>255)ix=255;
        const pc=plut[ix];
        if(k===LIQ) c=blend(pc,c,.16);
        else if(k===GAS) c=m===FOGO?blend(BGC,pc,.6):blend(BGC,pc,m===AR?.14:.4);
        else c=blend(dim(c,.7),pc,.3);
      }
    }else if(v===5){ // fluxo de ar: sobrepressão dinâmica (vermelho +, azul −) e velocidade (brilho)
      const ca=AF[i];
      if(k===FIX) c=dim(c,.45);
      else{const pa=ap[ca]*iar;let a=pa<0?-pa:pa;a=a>1?1:Math.sqrt(a);
        const u=ucx[ca],w=ucy[ca];let s=(u<0?-u:u)+(w<0?-w:w);s=s>1?1:s;
        const base=k===GAS?blend(BGC,AIRV,s*.55):dim(c,.6);
        c=a>.02?blend(base,pa>0?AIRP:AIRN,a*.85):base;}
    }else{
      let ti=t|0; if(ti<0)ti=0; else if(ti>4000)ti=4000;
      c=k===GAS?blend(BGC,tlut[ti],m===AR?.7:.85):(k===LIQ?tlut[ti]:blend(tlut[ti],c,.18));
    }
    const io=ionL[i]; // canal de plasma (raio/arco) visível em todas as vistas
    if(io>.03){c=blend(c,PLC,io>1?1:io); if(io>.2&&plasmaN<plasmaL.length) plasmaL[plasmaN++]=i;}
    u32[i]=c;
  }
  gravPost(u32); // vista de gravidade + lente do buraco negro (grav.js)
  ctx.putImageData(img,0,0);
  cntP=cnt; maxT=mt; maxG=mg;
  const target=mg>2000?mg:2000;
  pRange+= (target>pRange?(target-pRange)*.3:(target-pRange)*.04);
  const at=apMax>500?apMax:500; airRange+=(at>airRange?(at-airRange)*.3:(at-airRange)*.03);
}

// ================= UI =================
const IC={
 [AGUA]:'<svg viewBox="0 0 32 32"><path d="M16 3C13 8 7 14 7 20a9 9 0 0 0 18 0c0-6-6-12-9-17z" fill="#2f7be0"/><path d="M11.5 20a4.5 4.5 0 0 0 3.5 4.6" stroke="#c4e3ff" stroke-width="2" fill="none" stroke-linecap="round"/></svg>',
 [SALMOURA]:'<svg viewBox="0 0 32 32"><path d="M16 3C13 8 7 14 7 20a9 9 0 0 0 18 0c0-6-6-12-9-17z" fill="#2a8fa8"/><rect x="11.5" y="15" width="4" height="4" fill="#fff" transform="rotate(20 13.5 17)"/><rect x="17" y="19" width="3.4" height="3.4" fill="#eef" transform="rotate(-18 18.7 20.7)"/><rect x="12.5" y="22.5" width="2.8" height="2.8" fill="#fff"/></svg>',
 [OLEO]:'<svg viewBox="0 0 32 32"><path d="M16 3C13 8 7 14 7 20a9 9 0 0 0 18 0c0-6-6-12-9-17z" fill="#c8901a"/><path d="M10 19c3-2 9-2 12 0" stroke="#7fe8c8" stroke-width="1.4" fill="none"/><path d="M10.5 21.5c3-1.6 8-1.6 11 0" stroke="#c88af5" stroke-width="1.4" fill="none"/><circle cx="12.5" cy="15" r="1.6" fill="#ffe7a6"/></svg>',
 [MERCURIO]:'<svg viewBox="0 0 32 32"><defs><radialGradient id="gHg" cx=".35" cy=".3" r=".8"><stop offset="0" stop-color="#fff"/><stop offset=".35" stop-color="#c3c9d3"/><stop offset="1" stop-color="#4b525e"/></radialGradient></defs><ellipse cx="16" cy="19" rx="11" ry="8.5" fill="url(#gHg)"/><ellipse cx="25" cy="9" rx="3" ry="2.4" fill="url(#gHg)"/></svg>',
 [MEL]:'<svg viewBox="0 0 32 32"><path d="M16 3l9 5v10l-9 5-9-5V8z" fill="#e3a01a" stroke="#8a5806" stroke-width="1.5"/><path d="M16 9l4 2.3v4.6L16 18.2 12 15.9v-4.6z" fill="none" stroke="#8a5806" stroke-width="1.2"/><path d="M14.5 22.5v4a1.8 1.8 0 0 0 3.6 0v-5" fill="#e3a01a" stroke="#8a5806" stroke-width="1.2"/></svg>',
 [ALCOOL]:'<svg viewBox="0 0 32 32"><path d="M13 4v9l-7 12.5A2 2 0 0 0 7.8 28h16.4a2 2 0 0 0 1.8-2.5L19 13V4" fill="none" stroke="#cfd8ee" stroke-width="2"/><path d="M9.5 20h13l3 5.5a1.4 1.4 0 0 1-1.2 1.8H7.7a1.4 1.4 0 0 1-1.2-1.8z" fill="#8fa6ea"/><path d="M11 4h10" stroke="#cfd8ee" stroke-width="2" stroke-linecap="round"/><circle cx="14" cy="23.5" r="1" fill="#e6ecff"/><circle cx="18" cy="22" r=".8" fill="#e6ecff"/></svg>',
 [LAVA]:'<svg viewBox="0 0 32 32"><defs><linearGradient id="gLv" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="#ffd23a"/><stop offset=".5" stop-color="#ff6a10"/><stop offset="1" stop-color="#c41d00"/></linearGradient></defs><path d="M16 3C13 8 7 14 7 20a9 9 0 0 0 18 0c0-6-6-12-9-17z" fill="url(#gLv)"/><path d="M10 21l3 1 1 3-3 0zM19 17l3 1-1 3-2-1zM16 25l3-.5.5 2-3 .3z" fill="#3a1408" opacity=".85"/></svg>',
 [N2L]:'<svg viewBox="0 0 32 32"><rect x="9" y="11" width="14" height="17" rx="3" fill="#1f3550" stroke="#bcd8ff" stroke-width="1.5"/><rect x="11" y="17" width="10" height="9" rx="1.5" fill="#bfe2ff"/><path d="M12 9c-2-2 1-3 0-5M16 9c-2-2 1-3 0-5M20 9c-2-2 1-3 0-5" stroke="#dbeeff" stroke-width="1.4" fill="none" stroke-linecap="round"/><text x="16" y="24.5" font-size="5.5" text-anchor="middle" fill="#1f3550" font-family="sans-serif" font-weight="700">N₂</text></svg>',
 [ACO]:'<svg viewBox="0 0 32 32"><path d="M5 5h22v5h-8v12h8v5H5v-5h8V10H5z" fill="#8794a6"/><path d="M5 5h22v1.6H5zM5 22h8v1.4H5zM19 22h8v1.4h-8z" fill="#c9d2de"/></svg>',
 [CONCRETO]:'<svg viewBox="0 0 32 32"><rect x="4" y="7" width="24" height="18" rx="1.5" fill="#9a9890"/><circle cx="9" cy="12" r="1.6" fill="#6e6c66"/><circle cx="16" cy="18" r="2" fill="#77756e"/><circle cx="22" cy="11" r="1.3" fill="#bdbab2"/><circle cx="11" cy="20" r="1.1" fill="#c3c0b8"/><circle cx="24" cy="20" r="1.7" fill="#6e6c66"/><circle cx="18" cy="12" r=".9" fill="#5f5d58"/></svg>',
 [VIDRO]:'<svg viewBox="0 0 32 32"><rect x="6" y="4" width="20" height="24" rx="1.5" fill="rgba(140,205,225,.28)" stroke="#a6e0ef" stroke-width="1.6"/><path d="M10 14l6-6M10 20l11-11M17 24l5-5" stroke="#eaffff" stroke-width="1.4" stroke-linecap="round" opacity=".85"/></svg>',
 [LATA]:'<svg viewBox="0 0 32 32"><defs><linearGradient id="gCan" x1="0" x2="1"><stop offset="0" stop-color="#7f8796"/><stop offset=".45" stop-color="#eef1f5"/><stop offset="1" stop-color="#8d95a3"/></linearGradient></defs><rect x="9" y="6" width="14" height="22" rx="3" fill="url(#gCan)"/><rect x="9" y="12" width="14" height="9" fill="#c8323c"/><ellipse cx="16" cy="6.5" rx="7" ry="2" fill="#d6dbe3" stroke="#7f8796"/><path d="M14 13.5l4 6" stroke="#fff" stroke-width="1.2" opacity=".7"/></svg>',
 [BASALTO]:'<svg viewBox="0 0 32 32"><path d="M4 22l4-10 7-5 9 3 4 9-6 7H9z" fill="#4c4646"/><path d="M8 12l7-5 2 8-6 4zM17 15l7-5 4 9z" fill="#625a5a"/><path d="M11 19l6-4 5 9" stroke="#2c2828" fill="none"/></svg>',
 [RESIST]:'<svg viewBox="0 0 32 32"><rect x="3" y="23" width="26" height="4" rx="1" fill="#4a4f5c"/><path d="M4 18l3-8 3 8 3-8 3 8 3-8 3 8 3-8 3 8" fill="none" stroke="#ff6a3d" stroke-width="2.2" stroke-linejoin="round"/><path d="M4 18l3-8 3 8 3-8 3 8 3-8 3 8 3-8 3 8" fill="none" stroke="#ffd08a" stroke-width=".8" stroke-linejoin="round"/></svg>',
 [MADEIRA]:'<svg viewBox="0 0 32 32"><circle cx="16" cy="16" r="12" fill="#a8703c"/><circle cx="16" cy="16" r="9" fill="none" stroke="#7a4c24" stroke-width="1.2"/><circle cx="16" cy="16" r="6" fill="none" stroke="#7a4c24" stroke-width="1.1"/><circle cx="16" cy="16" r="3" fill="none" stroke="#7a4c24"/><path d="M16 16l7-9" stroke="#5a3616" stroke-width="1.2"/></svg>',
 [ISOPOR]:'<svg viewBox="0 0 32 32"><g fill="#f4f4ef" stroke="#a9aba3" stroke-width=".9"><circle cx="10" cy="11" r="4"/><circle cx="18" cy="9" r="4"/><circle cx="23" cy="16" r="4"/><circle cx="14" cy="17" r="4.2"/><circle cx="9" cy="22" r="4"/><circle cx="18" cy="24" r="4"/><circle cx="25" cy="24" r="3"/></g></svg>',
 [GELO]:'<svg viewBox="0 0 32 32"><path d="M16 4l11 6-11 6-11-6z" fill="#e8f7ff"/><path d="M5 10l11 6v12L5 22z" fill="#a9d8f5"/><path d="M27 10l-11 6v12l11-6z" fill="#86c1ea"/><path d="M9 13l3 1.7v4" stroke="#fff" stroke-width="1.2" fill="none" stroke-linecap="round"/></svg>',
 [AREIA]:'<svg viewBox="0 0 32 32"><path d="M3 27C8 14 12 10 16 10s8 4 13 17z" fill="#dcb874"/><g fill="#a9853f"><circle cx="11" cy="21" r=".9"/><circle cx="16" cy="16" r=".9"/><circle cx="20" cy="22" r=".9"/><circle cx="14" cy="24" r=".8"/><circle cx="23" cy="25" r=".8"/></g><g fill="#f2dca6"><circle cx="13" cy="18" r=".8"/><circle cx="19" cy="19" r=".7"/><circle cx="8" cy="25" r=".8"/></g></svg>',
 [SAL]:'<svg viewBox="0 0 32 32"><path d="M15 8l8 4-8 4-8-4z" fill="#ffffff"/><path d="M7 12l8 4v9l-8-4z" fill="#d8dbe6"/><path d="M23 12l-8 4v9l8-4z" fill="#b9bdcc"/><path d="M24 5v4M22 7h4" stroke="#fff" stroke-width="1.2" stroke-linecap="round"/><path d="M25 24l2.5 1.2-2.5 1.2-2.5-1.2z" fill="#e8ebf3"/></svg>',
 [VAPOR]:'<svg viewBox="0 0 32 32"><path d="M10 27c-3-4 3-6 0-10s3-6 0-10M16 27c-3-4 3-6 0-10s3-6 0-10M22 27c-3-4 3-6 0-10s3-6 0-10" stroke="#dfe8f2" stroke-width="2" fill="none" stroke-linecap="round"/></svg>',
 [HELIO]:'<svg viewBox="0 0 32 32"><ellipse cx="16" cy="12" rx="8" ry="9.5" fill="#ec9bd0"/><path d="M15 21.3h2l-1 1.6z" fill="#ec9bd0"/><path d="M16 23c-2 2 2 3 0 6" stroke="#c9ccd6" fill="none"/><ellipse cx="13" cy="8.5" rx="2" ry="3" fill="#fbd6ee"/><text x="16" y="15" font-size="6.5" text-anchor="middle" fill="#7a2a5f" font-family="sans-serif" font-weight="700">He</text></svg>',
 [FOGO]:'<svg viewBox="0 0 32 32"><path d="M16 3c1 6 9 9 9 17a9 9 0 0 1-18 0c0-5 3-7 4-10 1 3 2 4 3 4-1-4 0-8 2-11z" fill="#ff6a1a"/><path d="M16 14c.6 3 5 5 5 9a5 5 0 0 1-10 0c0-2.5 1.5-3.5 2-5 .8 1.4 1.4 2 2 2-.5-2 0-4 1-6z" fill="#ffd34a"/></svg>',
 [ERASE]:'<svg viewBox="0 0 32 32"><path d="M18 5l9 9-12 12H9l-5-5z" fill="#f08fa8"/><path d="M11.5 11.5l9 9" stroke="#fff" stroke-width="3"/><path d="M8 27h18" stroke="#7d8597" stroke-width="1.5" stroke-linecap="round"/></svg>',
 [HEAT]:'<svg viewBox="0 0 32 32"><rect x="13" y="4" width="6" height="17" rx="3" fill="none" stroke="#e6e9ef" stroke-width="1.5"/><circle cx="16" cy="24" r="5" fill="#ff4a3a"/><rect x="14.6" y="10" width="2.8" height="12" fill="#ff4a3a"/><path d="M24 8c-1.5 1.5 1.5 3 0 4.5M27 6c-1.5 1.5 1.5 3 0 4.5" stroke="#ff9a6a" stroke-width="1.4" fill="none" stroke-linecap="round"/></svg>',
 [COOL]:'<svg viewBox="0 0 32 32"><g stroke="#7fd6ff" stroke-width="2" stroke-linecap="round"><path d="M16 4v24M5.6 10l20.8 12M5.6 22l20.8-12"/><path d="M13 6.5l3 2.5 3-2.5M13 25.5l3-2.5 3 2.5M6.5 14l3.7-.8-1-3.6M25.5 18l-3.7.8 1 3.6M6.5 18l3.7.8-1 3.6M25.5 14l-3.7-.8 1-3.6"/></g></svg>',
 [COBRE]:'<svg viewBox="0 0 32 32"><g fill="none" stroke="#d27a42" stroke-width="2.6"><ellipse cx="10" cy="16" rx="3" ry="8"/><ellipse cx="16" cy="16" rx="3" ry="8"/><ellipse cx="22" cy="16" rx="3" ry="8"/></g><path d="M3 16h4M25 16h4" stroke="#f0a070" stroke-width="2.2" stroke-linecap="round"/></svg>',
 [NICROMO]:'<svg viewBox="0 0 32 32"><path d="M3 20h4l2-8 3 12 3-12 3 12 3-12 3 12 2-4h3" fill="none" stroke="#ff7a3a" stroke-width="2.2" stroke-linejoin="round"/><path d="M3 20h4l2-8 3 12 3-12 3 12 3-12 3 12 2-4h3" fill="none" stroke="#ffd9a0" stroke-width=".8" stroke-linejoin="round"/></svg>',
 [IMA]:'<svg viewBox="0 0 32 32"><rect x="4" y="10" width="12" height="12" fill="#e03530"/><rect x="16" y="10" width="12" height="12" fill="#2e62e8"/><text x="10" y="19.5" font-size="8" text-anchor="middle" fill="#fff" font-family="sans-serif" font-weight="700">N</text><text x="22" y="19.5" font-size="8" text-anchor="middle" fill="#fff" font-family="sans-serif" font-weight="700">S</text></svg>',
 [LIMALHA]:'<svg viewBox="0 0 32 32"><g stroke="#9aa0ac" stroke-width="1.6" stroke-linecap="round"><path d="M6 8l3 2M13 6l3 1M21 7l3 2M5 15l3 1M12 13l3 2M20 14l3 1M26 16l2 2M7 22l3 1M14 21l3 2M22 22l3 1M10 27l3 1M18 27l3 1"/></g></svg>',
 [FERROFL]:'<svg viewBox="0 0 32 32"><path d="M4 26c3 0 4-3 6-3 1-6 2-12 3-15 1 3 2 9 3 15 1-5 2-10 3-12 1 2 2 7 3 12 2 0 3 3 6 3z" fill="#1c1a22" stroke="#6a6680"/></svg>',
 [FONTEP]:'<svg viewBox="0 0 32 32"><circle cx="16" cy="16" r="11" fill="#d6342f"/><path d="M16 10v12M10 16h12" stroke="#fff" stroke-width="3" stroke-linecap="round"/></svg>',
 [FONTEN]:'<svg viewBox="0 0 32 32"><circle cx="16" cy="16" r="11" fill="#305cde"/><path d="M10 16h12" stroke="#fff" stroke-width="3" stroke-linecap="round"/></svg>',
 [TERRA]:'<svg viewBox="0 0 32 32"><path d="M16 4v12M6 16h20M10 21h12M14 26h4" stroke="#8fc070" stroke-width="2.4" stroke-linecap="round"/></svg>',
 [NUVEM]:'<svg viewBox="0 0 32 32"><path d="M8 24a5 5 0 0 1 0-10 7 7 0 0 1 13-3 5.5 5.5 0 0 1 3 13z" fill="#b9c0cc"/><path d="M14 24l-2 5M19 24l-2 5" stroke="#5aa0ff" stroke-width="1.6" stroke-linecap="round"/></svg>',
 [GRANIZO]:'<svg viewBox="0 0 32 32"><g fill="#eef6ff" stroke="#9cc0dd"><circle cx="10" cy="10" r="4"/><circle cx="21" cy="13" r="3.2"/><circle cx="13" cy="21" r="3.6"/><circle cx="23" cy="24" r="2.6"/></g></svg>',
 [VIDA]:'<svg viewBox="0 0 32 32"><g fill="#46d66e"><rect x="13" y="5" width="6" height="6" rx="1"/><rect x="20" y="12" width="6" height="6" rx="1"/><rect x="6" y="19" width="6" height="6" rx="1"/><rect x="13" y="19" width="6" height="6" rx="1"/><rect x="20" y="19" width="6" height="6" rx="1"/></g></svg>',
 [C4]:'<svg viewBox="0 0 32 32"><rect x="5" y="9" width="22" height="15" rx="2" fill="#ddd3b0" stroke="#8a8166" stroke-width="1.2"/><rect x="5" y="14" width="22" height="3" fill="#3c3a34"/><path d="M20 9V5h4" stroke="#d8402c" stroke-width="1.6" fill="none"/><circle cx="24" cy="5" r="1.6" fill="#ffd24a"/></svg>',
 [NITRO]:'<svg viewBox="0 0 32 32"><path d="M12 3h8v7l6 13a3 3 0 0 1-2.7 4.3H8.7A3 3 0 0 1 6 23l6-13z" fill="none" stroke="#c9d1e0" stroke-width="1.6"/><path d="M9 18h14l2.3 5a1.6 1.6 0 0 1-1.5 2.3H8.2A1.6 1.6 0 0 1 6.7 23z" fill="#e2e8a8"/><path d="M15 12l-2 4h3l-2 4" stroke="#d8402c" stroke-width="1.5" fill="none" stroke-linejoin="round"/></svg>',
 [HNO3]:'<svg viewBox="0 0 32 32"><path d="M13 3h6v8l7 14a2 2 0 0 1-1.8 3H7.8A2 2 0 0 1 6 25l7-14z" fill="none" stroke="#c9d1e0" stroke-width="1.6"/><path d="M9.5 18h13l3 6.5a1.4 1.4 0 0 1-1.3 2H7.8a1.4 1.4 0 0 1-1.3-2z" fill="#e8c45c"/><path d="M18 9c3-2 6-1 7 1M17 6c2-2 5-2 7 0" stroke="#c0602a" stroke-width="1.5" fill="none" stroke-linecap="round"/></svg>',
 [WIND]:'<svg viewBox="0 0 32 32"><path d="M4 12h16a4 4 0 1 0-4-4M4 18h21a4 4 0 1 1-4 4M4 24h9" stroke="#9fd8ff" stroke-width="2.4" fill="none" stroke-linecap="round"/></svg>',
 [CYCL]:'<svg viewBox="0 0 32 32"><path d="M26 9c-3-4-12-5-17 0M27 15c-4-3-15-3-20 1M24 21c-3-2-11-2-14 1M20 26c-2-1-5-1-7 1" stroke="#9fd8ff" stroke-width="2.2" fill="none" stroke-linecap="round"/></svg>',
 [CHGP]:'<svg viewBox="0 0 32 32"><circle cx="16" cy="16" r="11" fill="none" stroke="#ff7a60" stroke-width="2"/><path d="M16 10v12M10 16h12" stroke="#ffb09a" stroke-width="3" stroke-linecap="round"/></svg>',
 [CHGN]:'<svg viewBox="0 0 32 32"><circle cx="16" cy="16" r="11" fill="none" stroke="#5a9cff" stroke-width="2"/><path d="M10 16h12" stroke="#a8c8ff" stroke-width="3" stroke-linecap="round"/></svg>',
 [FLIP]:'<svg viewBox="0 0 32 32"><path d="M7 12a9 9 0 0 1 16-3l2-2v7h-7l2.5-2.5A6 6 0 0 0 10 13z" fill="#e8665c"/><path d="M25 20a9 9 0 0 1-16 3l-2 2v-7h7l-2.5 2.5A6 6 0 0 0 22 19z" fill="#4f86ec"/></svg>',
 [SPARK]:'<svg viewBox="0 0 32 32"><path d="M18 2L7 18h7l-3 12 13-17h-7l4-11z" fill="#e6e8ff" stroke="#8f9cff" stroke-width="1.2" stroke-linejoin="round"/></svg>'
};
const fmtP=p=>p>=1e6?(p/1e6).toFixed(2)+' MPa':(p/1e3).toFixed(p<1e4?1:0)+' kPa';
const DESC={
 [AGUA]:'ρ 1000 kg/m³ · máx. densidade a 4 °C · ferve 100 °C @1 atm',
 [SALMOURA]:'ρ 1100 · ferve ~108 °C · congela −21 °C',
 [OLEO]:'ρ 920 · viscoso · inflama ~327 °C',
 [MERCURIO]:'ρ 13 534 · quase tudo flutua · ferve 357 °C, congela −39 °C',
 [MEL]:'ρ 1420 · viscosidade cai muito com calor',
 [ALCOOL]:'ρ 789 · ferve 78 °C · inflamável',
 [LAVA]:'ρ 2600 · ~1180 °C · solidifica em basalto <727 °C',
 [N2L]:'ρ 808 · −196 °C · ferve ao contato',
 [ACO]:'σ 250 MPa · funde 1510 °C · ótimo condutor',
 [CONCRETO]:'σ 30 MPa · vira entulho',
 [VIDRO]:'σ 2 MPa · quebra por pressão ou choque térmico',
 [LATA]:'σ 110 kPa · parede fina, amassa · funde 660 °C',
 [BASALTO]:'σ 150 MPa · funde ~1200 °C',
 [RESIST]:'elemento aquecedor na temperatura ajustada',
 [MADEIRA]:'ρ 600 · flutua · encharca sob >150 kPa · queima 300 °C',
 [ISOPOR]:'ρ 20 · colapsa >100 kPa manométrica · derrete 240 °C',
 [GELO]:'ρ 917 · flutua · fusão cai com pressão',
 [AREIA]:'ρ 1600 · vira vidro a 1700 °C',
 [SAL]:'ρ 2160 · dissolve na água · derrete gelo (endotérmico)',
 [VAPOR]:'gás M 18 g/mol · condensa abaixo de Teb(P)',
 [HELIO]:'gás M 4 g/mol · sobe rápido',
 [FOGO]:'combustão · gás quente sobe',
 [ERASE]:'remove (botão direito)',[HEAT]:'injeta calor: ΔT = Q/(ρc)',[COOL]:'retira calor',
 [COBRE]:'σ 5,96·10⁷ S/m · ótimo condutor · funde 1085 °C',
 [NICROMO]:'σ 9,1·10⁵ S/m · resistivo → aquece por efeito Joule',
 [IMA]:'NdFeB · Br ≈ 1,2 T · perde magnetização acima de 310 °C (Curie) · direção N no seletor "Polo N"',
 [LIMALHA]:'ferro em pó · atraído por ∇|B| · mostra o campo',
 [FERROFL]:'líquido magnético · Ms 45 kA/m · sobe até o ímã',
 [FONTEP]:'terminal positivo: potencial fixo +V/2 (ajuste "Tensão")',
 [FONTEN]:'terminal negativo: potencial fixo −V/2',
 [TERRA]:'aterramento: potencial 0 (o piso da câmara também é terra)',
 [NUVEM]:'gotículas/cristais suspensos · com graupel separa cargas',
 [GRANIZO]:'graupel · colisões com cristais de gelo transferem carga (− abaixo de −15 °C, + acima)',
 [VIDA]:'jogo da vida B3/S23 · nasce só no ar entre −5 e 60 °C · morre por calor, frio, fogo, raio, choque (>20 mA), afogamento e esmagamento',
 [CHGP]:'deposita carga positiva (C/célula ajustável)',[CHGN]:'deposita carga negativa',
 [FLIP]:'inverte: sinal da carga, polo do ímã, polaridade da fonte',
 [SPARK]:'ioniza o meio (gatilho de arco)',
 [C4]:'explosivo plástico (RDX) · 2,1 kg TNT/célula · insensível a choque; detona com chama, faísca/raio ou >250 °C · gruda (sólido fixo)',
 [NITRO]:'explosivo líquido · 2,4 kg TNT/célula · sensível a choque (onda >150 kPa), chama, faísca ou >215 °C',
 [HNO3]:'oxidante forte · dissolve metais soltando NO₂, inflama orgânicos, reação explosiva com álcool · ferve a 83 °C',
 [WIND]:'arraste para dar a direção; segurar parado mantém o jato (ventilador) · força no painel Ar',
 [CYCL]:'segure para criar um vórtice (anti-horário; Shift = horário)'
};
const GROUPS=[
 ['Líquidos',[AGUA,SALMOURA,OLEO,MERCURIO,MEL,ALCOOL,LAVA,N2L]],
 ['Sólidos fixos',[ACO,CONCRETO,VIDRO,LATA,BASALTO,RESIST]],
 ['Soltos',[MADEIRA,ISOPOR,GELO,AREIA,SAL]],
 ['Gases',[VAPOR,HELIO,FOGO]],
 ['Elétrica',[COBRE,NICROMO,FONTEP,FONTEN,TERRA]],
 ['Magnetismo',[IMA,LIMALHA,FERROFL]],
 ['Clima',[NUVEM,GRANIZO]],
 ['Vida',[VIDA,VIDAM]],
 ['Radiação',[U235,PU239,UNAT,AMBE,RADIO,COBALTO,CESIO,GRAFITE,BORO,CHUMBO,PFIS,CORIO,RADONIO]],
 ['Explosivos',[C4,NITRO,HNO3]],
 ['Ferramentas',[ERASE,HEAT,COOL,CHGP,CHGN,FLIP,SPARK,WIND,CYCL,GRAVT]]
];
const SUB={[ERASE]:'',[HEAT]:'+Q',[COOL]:'−Q',[RESIST]:'aquecedor',[FOGO]:'~1100 °C',[VAPOR]:'gás',[HELIO]:'gás',
  [COBRE]:'6·10⁷ S/m',[NICROMO]:'9·10⁵ S/m',[FONTEP]:'+V/2',[FONTEN]:'−V/2',[TERRA]:'0 V',[IMA]:'1,2 T',[LIMALHA]:'ferromag.',
  [FERROFL]:'ferromag.',[NUVEM]:'gás',[GRANIZO]:'ρ 700',[VIDA]:'B3/S23',[CHGP]:'+q',[CHGN]:'−q',[FLIP]:'±',[SPARK]:'plasma',[C4]:'2,1 kg TNT',[NITRO]:'2,4 kg TNT',[HNO3]:'oxidante',[WIND]:'ar',[CYCL]:'vórtice'};
const NAMES={[ERASE]:'Borracha',[HEAT]:'Aquecer',[COOL]:'Resfriar',[CHGP]:'Carga +',[CHGN]:'Carga −',[FLIP]:'Inverter polo',[SPARK]:'Faísca',[WIND]:'Vento',[CYCL]:'Ciclone'};
NAMES[GRAVT]='Gravidade';SUB[GRAVT]='massa';
DESC[GRAVT]='segure para aumentar a massa local (×1,07/quadro, g = GM/r²); Shift diminui · r_s = 2GM/c² ≥ 5 cm → buraco negro';
IC[GRAVT]='<svg viewBox="0 0 32 32"><ellipse cx="16" cy="16" rx="13" ry="5" stroke="#c8b4ff" stroke-width="2" fill="none"/><path d="M5 10c4 3 7 4 11 4s7-1 11-4M5 22c4-3 7-4 11-4s7 1 11 4" stroke="#8f7ad8" stroke-width="1.4" fill="none"/><circle cx="16" cy="16" r="4.5" fill="#0b0b10" stroke="#ffcf8a" stroke-width="1.6"/></svg>';
// ---- radiação (rad.js) ----
const trefoil=(bg,fg,x)=>'<svg viewBox="0 0 32 32"><circle cx="16" cy="16" r="14" fill="'+bg+'"/><g fill="'+fg+'"><circle cx="16" cy="16" r="2.4"/>'+[0,120,240].map(a=>'<path d="M16 16L'+(16+11*Math.cos((a-30)*Math.PI/180)).toFixed(1)+' '+(16+11*Math.sin((a-30)*Math.PI/180)).toFixed(1)+'A11 11 0 0 1 '+(16+11*Math.cos((a+30)*Math.PI/180)).toFixed(1)+' '+(16+11*Math.sin((a+30)*Math.PI/180)).toFixed(1)+'z" transform="rotate(-90 16 16)" opacity=".9"/>').join('')+'</g><circle cx="16" cy="16" r="4" fill="'+bg+'"/><circle cx="16" cy="16" r="2.2" fill="'+fg+'"/>'+(x||'')+'</svg>';
const blk=(c,e,t)=>'<svg viewBox="0 0 32 32"><rect x="5" y="7" width="22" height="18" rx="2" fill="'+c+'" stroke="'+e+'" stroke-width="1.5"/><text x="16" y="20.5" text-anchor="middle" font-size="9" font-family="monospace" font-weight="700" fill="'+e+'">'+t+'</text></svg>';
Object.assign(IC,{[U235]:trefoil('#3d5a2a','#b7f27a'),[PU239]:trefoil('#5a3a2a','#ffb27a'),[UNAT]:trefoil('#3a3f3a','#9aa89a'),
 [AMBE]:trefoil('#2a2f4a','#9ab0ff'),[RADIO]:trefoil('#1a3a28','#8affb0'),[COBALTO]:trefoil('#2a3448','#8fb4ff'),[CESIO]:trefoil('#4a4020','#f2dc7a'),
 [CORIO]:trefoil('#5a2a10','#ffa040'),[PFIS]:trefoil('#3a3530','#c8b890'),[RADONIO]:trefoil('#202830','#a0ffd0'),
 [GRAFITE]:blk('#2b2d31','#8a8d95','C'),[BORO]:blk('#3a3530','#d8c08a','B'),[CHUMBO]:blk('#4a4e58','#c8ccd8','Pb'),[CHUMBOF]:blk('#6a6e78','#e8ecf8','Pb'),
 [VIDAM]:'<svg viewBox="0 0 32 32"><g fill="#ce60e2"><rect x="12" y="5" width="7" height="7"/><rect x="19" y="12" width="7" height="7"/><rect x="5" y="19" width="7" height="7"/><rect x="12" y="19" width="7" height="7"/><rect x="19" y="19" width="7" height="7"/></g></svg>'});
Object.assign(SUB,{[U235]:'físsil',[PU239]:'físsil',[UNAT]:'0,7% ²³⁵U',[AMBE]:'fonte n',[RADIO]:'α γ',[COBALTO]:'γ forte',[CESIO]:'β γ solúvel',
 [CORIO]:'lava nuclear',[PFIS]:'precipitação',[RADONIO]:'gás α',[GRAFITE]:'moderador',[BORO]:'absorve n',[CHUMBO]:'blinda γ',[VIDAM]:'B36/S23'});
Object.assign(DESC,{
 [U235]:'urânio altamente enriquecido · ρ 19 050 · físsil: nêutron → fissão (~200 MeV, 2,4 nêutrons) · massa crítica nua ~R6 células; menos com refletor/compressão · funde 1132 °C → cório',
 [PU239]:'plutônio · ρ 19 800 · mais físsil que U-235 (2,9 n/fissão) e emite nêutrons espontâneos · massa crítica ~R4,5 · núcleo de bomba de implosão',
 [UNAT]:'urânio natural (0,7% U-235) · não sustenta cadeia rápida; com moderador (grafite/água) quase crítico · ótimo refletor/tamper denso',
 [AMBE]:'fonte de nêutrons amerício-berílio (α + Be → n) · iniciador: esmagada por implosão solta rajada de nêutrons',
 [RADIO]:'rádio-226 · emissor α (parado por papel/pele) + γ fraco · brilha · emana radônio',
 [COBALTO]:'cobalto-60 · γ forte 1,2 MeV · atravessa metros de ar; blindar com chumbo/concreto/água',
 [CESIO]:'cloreto de césio-137 · β + γ 0,66 MeV · pó solúvel: contamina água e se espalha (Goiânia 1987)',
 [CORIO]:'combustível derretido (~2500 °C) · altamente radioativo · ainda fissiona com moderador · solidifica em resíduo',
 [PFIS]:'produtos de fissão / precipitação radioativa · β γ · decai com o tempo',
 [RADONIO]:'radônio-222 · gás nobre pesado α · acumula em porões fechados · meia-vida curta',
 [GRAFITE]:'carbono puro · modera nêutrons (rápido → térmico) sem absorver · queima a ~750 °C',
 [BORO]:'boro-10 · absorve nêutrons térmicos (barra de controle)',
 [CHUMBO]:'ρ 11 340 · melhor blindagem γ por espessura · funde 327 °C',
 [VIDAM]:'vida mutante (radiação) · regra HighLife B36/S23 (nasce com 3 ou 6 vizinhos)'});
// ================= paleta de materiais (abas + busca) e trilho de ferramentas =================
const MODES=[
 ['free','Pincel redondo','arrastar para desenhar (B)','<circle cx="12" cy="12" r="7"/>'],
 ['sq','Pincel quadrado','arrastar para desenhar','<rect x="5" y="5" width="14" height="14" rx="1"/>'],
 ['line','Linha','arrastar do início ao fim (L)','<path d="M5 19L19 5"/><circle cx="5" cy="19" r="1.5" fill="#c9d1e0"/><circle cx="19" cy="5" r="1.5" fill="#c9d1e0"/>'],
 ['rect','Retângulo cheio','arrastar a diagonal (R)','<rect x="4" y="6" width="16" height="12" rx="1" fill="#c9d1e0" fill-opacity=".35"/>'],
 ['rectO','Retângulo contorno','espessura = pincel','<rect x="4" y="6" width="16" height="12" rx="1"/><rect x="7.5" y="9.5" width="9" height="5" stroke-dasharray="2 2" stroke-width="1"/>'],
 ['ell','Elipse cheia','arrastar a diagonal (E)','<ellipse cx="12" cy="12" rx="8" ry="6" fill="#c9d1e0" fill-opacity=".35"/>'],
 ['fill','Balde','preenche a região conexa de mesmo material (G)','<path d="M6 10l6-6 7 7-6 6a2 2 0 0 1-3 0l-4-4a2 2 0 0 1 0-3z"/><path d="M5 11h13"/><path d="M20 15s-1.5 2-1.5 3a1.5 1.5 0 0 0 3 0c0-1-1.5-3-1.5-3z" fill="#c9d1e0"/>']];
let drawMode='free',shapePrev=null;
const palEl=document.getElementById('palette'),tabsEl=document.getElementById('tabs'),railEl=document.getElementById('rail');
const searchEl=document.getElementById('matSearch'),selInfo=document.getElementById('selInfo'),tipEl=document.getElementById('tip');
const MAT_GROUPS=GROUPS.filter(g=>g[0]!=='Ferramentas'),TOOL_LIST=GROUPS.find(g=>g[0]==='Ferramentas')[1];
const nameOf=id=>NAMES[id]||D[id].n;
const subOf=id=>SUB[id]!==undefined?SUB[id]:(D[id].kind===FIX?'σ '+fmtP(D[id].str):'ρ '+D[id].rho);
const norm=s=>s.normalize('NFD').replace(/[̀-ͯ]/g,'').toLowerCase();
let matTab=0,lastMat=tool;
function selectTool(id){
  tool=id; if(id<100) lastMat=id;
  document.querySelectorAll('[data-id]').forEach(e=>e.classList.toggle('on',+e.dataset.id===id));
  selInfo.innerHTML='<div class="ic">'+IC[id]+'</div><div><div class="t">'+nameOf(id)+' <span style="color:var(--mut);font-weight:400;font-family:var(--mono);font-size:10px">'+subOf(id)+'</span></div><div class="d">'+(DESC[id]||'')+'</div></div>';
  const cur=document.getElementById('railCur'); if(cur) cur.innerHTML=IC[lastMat];
}
function tile(id){
  const b=document.createElement('button');b.className='it'+(id===tool?' on':'');b.dataset.id=id;b.title=DESC[id]||'';
  b.innerHTML=IC[id]+'<span class="nm">'+nameOf(id)+'</span><span class="sb">'+subOf(id)+'</span>';
  b.onclick=()=>selectTool(id);return b;
}
function renderPalette(){
  const q=norm(searchEl.value.trim());palEl.innerHTML='';let n=0;
  tabsEl.querySelectorAll('.tab').forEach(t=>t.classList.toggle('on',!q&&+t.dataset.k===matTab));
  for(let g=0;g<MAT_GROUPS.length;g++){
    const [gn,list]=MAT_GROUPS[g];
    const ids=q?list.filter(id=>norm(nameOf(id)+' '+(DESC[id]||'')+' '+gn).includes(q)):(matTab<0||g===matTab?list:[]);
    if(!ids.length)continue;
    if(q||matTab<0){const h=document.createElement('div');h.className='gh';h.textContent=gn;palEl.appendChild(h);}
    for(const id of ids){palEl.appendChild(tile(id));n++;}
  }
  if(!n)palEl.innerHTML='<div class="empty">Nenhum material encontrado</div>';
}
[['Todos'],...MAT_GROUPS].forEach(([gn],k)=>{k--;const t=document.createElement('button');t.className='tab';t.textContent=gn;t.dataset.k=k;
  t.onclick=()=>{matTab=k;searchEl.value='';renderPalette();};tabsEl.appendChild(t);});
document.getElementById('matCount').textContent=MAT_GROUPS.reduce((a,g)=>a+g[1].length,0);
searchEl.oninput=renderPalette;
searchEl.onkeydown=e=>{if(e.key==='Escape'){searchEl.value='';renderPalette();searchEl.blur();}
  if(e.key==='Enter'){const f=palEl.querySelector('.it');if(f){f.click();searchEl.blur();}}};
// trilho: material atual + ferramentas
{const cur=document.createElement('button');cur.className='tl';cur.id='railCur';cur.dataset.tip='Material atual|voltar a desenhar com ele';
  cur.onclick=()=>selectTool(lastMat);railEl.appendChild(cur);railEl.appendChild(document.createElement('hr'));
  for(const id of TOOL_LIST){const b=document.createElement('button');b.className='tl';b.dataset.id=id;b.innerHTML=IC[id];
    b.dataset.tip=nameOf(id)+'|'+(DESC[id]||'');b.onclick=()=>selectTool(id);railEl.appendChild(b);}
  railEl.appendChild(document.createElement('hr'));
  for(const [md,nm,ds,ic] of MODES){const b=document.createElement('button');b.className='tl md'+(md===drawMode?' on':'');b.dataset.md=md;
    b.innerHTML='<svg viewBox="0 0 24 24" fill="none" stroke="#c9d1e0" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">'+ic+'</svg>';
    b.dataset.tip=nm+'|'+ds;b.onclick=()=>setMode(md);railEl.appendChild(b);}}
function setMode(md){drawMode=md;brushSq=md==='sq';railEl.querySelectorAll('.md').forEach(e=>e.classList.toggle('on',e.dataset.md===md));}
// tooltip
document.addEventListener('pointerover',e=>{const t=e.target.closest('[data-tip]');if(!t){tipEl.style.opacity=0;return;}
  const [a,b]=t.dataset.tip.split('|');tipEl.innerHTML='<b>'+a+'</b>'+(b?'<span>'+b+'</span>':'');
  const r=t.getBoundingClientRect();tipEl.style.opacity=1;
  let x=r.right+10,y=r.top+r.height/2-tipEl.offsetHeight/2;
  if(x+tipEl.offsetWidth>innerWidth-8){x=r.left;y=r.bottom+8;}
  tipEl.style.left=x+'px';tipEl.style.top=Math.max(8,y)+'px';});
renderPalette();selectTool(tool);

// ================= controles =================
function bindRange(id,out,fn,fmt){const el=document.getElementById(id),o=document.getElementById(out);
  const up=()=>{fn(+el.value);o.textContent=fmt(+el.value);el.style.setProperty('--p',((el.value-el.min)/(el.max-el.min)*100)+'%');};
  el.oninput=up;up();return el;}
const brushEl0=bindRange('sBrush','oBrush',v=>brush=v,v=>v*2+1+' px');
bindRange('sGrav','oGrav',v=>gm=v,v=>v.toFixed(1)+' g');
const heaterEl=bindRange('sHeater','oHeater',v=>heaterT=v+273.15,v=>v+' °C');
function heaterSet(c){heaterEl.value=c;heaterEl.oninput();}
bindRange('sPow','oPow',v=>power=v,v=>v+' MJ/m³');
const SPEEDS=[.125,.25,.5,1,2,3,4];let stepAcc=0;
bindRange('sSpeed','oSpeed',v=>steps=SPEEDS[v],v=>SPEEDS[v]<1?'1/'+(1/SPEEDS[v])+' (lento)':'×'+SPEEDS[v]);
const voltEl=bindRange('sVolt','oVolt',v=>Vsrc=Math.pow(10,v/10),v=>fmtV(Math.pow(10,v/10)));
function setVolt(V){voltEl.value=Math.round(10*Math.log10(V));voltEl.oninput();}
bindRange('sChg','oChg',v=>qAmt=Math.pow(10,v/10),v=>fmtQ(Math.pow(10,v/10)));
const windMs=v=>50*Math.pow(200,v/100); // slider 0–100 → 50 m/s … 10 km/s (log)
bindRange('sWind','oWind',v=>windF=windMs(v),v=>{const s=windMs(v);return '≤ '+(s>=1000?(s/1000).toFixed(1)+' km/s':Math.round(s)+' m/s');}); // .06·v células grossas/passo ÷ som (.3) × 343 m/s (sem consts do air.js: ainda não carregado) // força máx. das ferramentas Vento/Ciclone
const poleEl=document.getElementById('sPole');
poleEl.onclick=e=>{const b=e.target.closest('button');if(!b)return;magDir=+b.dataset.d;poleEl.querySelectorAll('button').forEach(x=>x.classList.toggle('on',x===b));};
const bStorm=document.getElementById('bStorm');
function setStorm(v){storm=v;bStorm.textContent='Tempestade '+(v?'ligada':'desligada');bStorm.classList.toggle('on',v);}
bStorm.onclick=()=>setStorm(!storm);setStorm(storm);
const bAir=document.getElementById('bAir'),bBlast=document.getElementById('bBlast');
function setAir(v){airOn=v;if(!v)airReset();bAir.textContent='Ar '+(v?'ligado':'desligado');bAir.classList.toggle('on',v);updateLegend();}
function setBlast(v){blastOn=v;bBlast.textContent='Choque '+(v?'ligado':'desligado');bBlast.classList.toggle('on',v);}
bAir.onclick=()=>setAir(!airOn);bBlast.onclick=()=>setBlast(!blastOn);
document.getElementById('bAirZ').onclick=()=>airReset();
const viewsEl=document.getElementById('views');
function setView(v){view=v;viewsEl.querySelectorAll('button').forEach(b=>b.classList.toggle('on',+b.dataset.v===v));updateLegend();}
viewsEl.onclick=e=>{const b=e.target.closest('button');if(b)setView(+b.dataset.v);};
const bPause=document.getElementById('bPause');
const IPLAY='<svg viewBox="0 0 24 24"><path d="M7 5l12 7-12 7z" fill="currentColor"/></svg>',IPAUSE='<svg viewBox="0 0 24 24"><path d="M8 5v14M16 5v14"/></svg>';
function updPause(){bPause.innerHTML=paused?IPLAY+'Continuar':IPAUSE+'Pausar';bPause.classList.toggle('primary',!paused);}
function togglePause(){paused=!paused;updPause();}
bPause.onclick=togglePause;updPause();
document.getElementById('bLid').onclick=()=>{pushUndo();setLid(!lid);};
document.getElementById('bClear').onclick=()=>{pushUndo();clearAll();};
document.getElementById('sPreset').onchange=e=>{const p=e.target.value;if(p){pushUndo();if(lid)setLid(false);PRESETS[p]();zoomTo(1,0,0);}e.target.value='';};
const helpDlg=document.getElementById('help');
document.getElementById('bHelp').onclick=()=>helpDlg.showModal();
document.getElementById('bHelpX').onclick=()=>helpDlg.close();
helpDlg.onclick=e=>{if(e.target===helpDlg)helpDlg.close();};
addEventListener('keydown',e=>{
  if(e.key==='Shift'){shiftKey=true;constrainShape();}
  if(e.target.tagName==='INPUT'||e.target.tagName==='SELECT'||helpDlg.open)return;
  const mod=e.ctrlKey||e.metaKey,k=e.key.toLowerCase();
  if(mod){
    if(k==='z'&&!e.shiftKey){e.preventDefault();undo();}
    else if(k==='y'||(k==='z'&&e.shiftKey)){e.preventDefault();redo();}
    else if(k==='s'){e.preventDefault();saveScene();}
    else if(k==='o'){e.preventDefault();fileIn.click();}
    return;
  }
  if(e.key==='.'){stepOnce();return;}
  if(k==='p'){screenshot();return;}
  const MK={b:'free',l:'line',r:'rect',e:'ell',g:'fill'};if(MK[k]){setMode(MK[k]);return;}
  if(e.code==='Space'){e.preventDefault();togglePause();}
  if(e.key>='1'&&e.key<='8')setView(+e.key-1);
  if(e.key==='/'){e.preventDefault();searchEl.focus();searchEl.select();}
  if(e.key==='['||e.key===']'){brushEl0.value=+brushEl0.value+(e.key===']'?1:-1);brushEl0.oninput();}
  if(e.key==='+'||e.key==='=')zoomBy(1.5);
  if(e.key==='-')zoomBy(1/1.5);
  if(e.key==='0'||(e.key==='Escape'&&Z>1.001))zoomTo(1,0,0);
  if(e.key==='Escape'){pinI=-1;if(shapePrev){shapePrev=null;drawing=false;}}
  if(e.key==='f'||e.key==='F')focusAt(hoverI>=0?hoverI:pinI);
});
addEventListener('keyup',e=>{if(e.key==='Shift'){shiftKey=false;constrainShape();}});

// legenda
const lgGrad=document.getElementById('lgGrad'),lgLo=document.getElementById('lgLo'),lgHi=document.getElementById('lgHi'),lgView=document.getElementById('lgView');
const VIEWN=['material','pressão manom.','temperatura','potencial','|B|','fluxo de ar','taxa de dose'];VIEWN[7]='gravidade'; // GVIEW (grav.js)
const css=c=>'rgb('+(c[0]|0)+','+(c[1]|0)+','+(c[2]|0)+')';
function updateLegend(){
  lgView.textContent=VIEWN[view];
  if(view===GVIEW){gravLegend();return;}
  if(view===1){lgGrad.style.background='linear-gradient(90deg,'+PST.map(s=>css(s[1])+' '+s[0]*100+'%').join(',')+')';lgLo.textContent='0 kPa';lgHi.textContent='+'+fmtP(pRange);}
  else if(view===2){const lo=0,hi=1500;lgGrad.style.background='linear-gradient(90deg,'+TST.filter(s=>s[0]<=hi).map(s=>css(s[1])+' '+(s[0]-lo)/(hi-lo)*100+'%').join(',')+')';lgLo.textContent='−273 °C';lgHi.textContent='1227 °C';}
  else if(view===3){lgGrad.style.background='linear-gradient(90deg,rgb(40,124,255),rgb(14,17,24) 50%,rgb(255,78,40))';lgLo.textContent='−'+fmtV(vScale);lgHi.textContent='+'+fmtV(vScale)+' · amarelo = corrente';}
  else if(view===4){lgGrad.style.background='linear-gradient(90deg,'+MST.map(s=>css(s[1])+' '+s[0]*100+'%').join(',')+')';lgLo.textContent='10 µT';lgHi.textContent='1 T · ⊙ sai ⊗ entra';}
  else if(view===6){lgGrad.style.background='linear-gradient(90deg,'+RST.map(s=>css(s[1])+' '+s[0]*100+'%').join(',')+')';lgLo.textContent='0,1 µSv/h';lgHi.textContent='10 kSv/h · pontos = nêutrons';}
  else if(view===5){lgGrad.style.background='linear-gradient(90deg,rgb(60,130,255),rgb(14,17,24) 50%,rgb(255,70,50))';lgLo.textContent='−'+fmtP(airRange);lgHi.textContent=(airOn?'+'+fmtP(airRange)+' · linhas = vento':'ar desativado');}
  else{lgGrad.style.background='linear-gradient(90deg,#000,'+GST.map(s=>css(s[1])).join(',')+')';lgLo.textContent='brilho térmico';lgHi.textContent='>2000 °C';}
}

// ================= zoom / vista =================
// Z = ampliação; ZX,ZY = canto superior esquerdo visível (células)
let Z=1,ZX=0,ZY=0,zT=null,pinI=-1,panning=null;
const ZMAX=16,stage=document.getElementById('stage'),brushEl=document.getElementById('brush'),zLvl=document.getElementById('zLvl');
function clampView(){if(Z<1)Z=1;if(Z>ZMAX)Z=ZMAX;const v=W/Z;ZX=Math.max(0,Math.min(W-v,ZX));ZY=Math.max(0,Math.min(H-v,ZY));}
function applyZoom(){clampView();
  cv.style.transform=Z===1?'':'translate('+(-ZX/W*Z*100)+'%,'+(-ZY/H*Z*100)+'%) scale('+Z+')';
  zLvl.textContent=Math.round(Z*100)+'%';updRuler();}
function zoomAround(nz,fx,fy){ // fx,fy: fração [0,1] da tela que fica fixa
  const cx=ZX+fx*W/Z,cy=ZY+fy*H/Z;Z=Math.max(1,Math.min(ZMAX,nz));ZX=cx-fx*W/Z;ZY=cy-fy*H/Z;zT=null;applyZoom();}
function zoomBy(f){zoomAround(Z*f,.5,.5);}
function zoomTo(z,x,y){zT={z:Math.max(1,Math.min(ZMAX,z)),x,y};}
function zoomAnim(){if(!zT)return;const a=.2;Z+=(zT.z-Z)*a;ZX+=(zT.x-ZX)*a;ZY+=(zT.y-ZY)*a;
  if(Math.abs(zT.z-Z)<.002&&Math.abs(zT.x-ZX)<.02&&Math.abs(zT.y-ZY)<.02){Z=zT.z;ZX=zT.x;ZY=zT.y;zT=null;}applyZoom();}
function focusAt(i){ // enquadra o corpo conexo de mesmo material
  if(i<0)return;const m=mat[i];if(m===AR||m===PAREDE){zoomTo(1,0,0);return;}
  const seen=new Uint8Array(N),q=new Int32Array(N);let h=0,t=0,x0=W,x1=0,y0=H,y1=0;q[t++]=i;seen[i]=1;
  while(h<t){const c=q[h++],x=c%W,y=c/W|0;if(x<x0)x0=x;if(x>x1)x1=x;if(y<y0)y0=y;if(y>y1)y1=y;
    for(const n of [c-1,c+1,c-W,c+W])if(n>=0&&n<N&&!seen[n]&&mat[n]===m){seen[n]=1;q[t++]=n;}}
  const sz=Math.max(x1-x0+1,y1-y0+1)*1.35+8,z=Math.max(1,Math.min(ZMAX,W/sz)),v=W/z;
  zoomTo(z,(x0+x1+1)/2-v/2,(y0+y1+1)/2-v/2);pinI=i;
}
stage.addEventListener('wheel',e=>{e.preventDefault();const r=stage.getBoundingClientRect();
  zoomAround(Z*Math.exp(-e.deltaY*(e.deltaMode?0.05:0.0015)),(e.clientX-r.left)/r.width,(e.clientY-r.top)/r.height);updBrush(e);},{passive:false});
document.getElementById('zIn').onclick=()=>zoomBy(1.5);
document.getElementById('zOut').onclick=()=>zoomBy(1/1.5);
document.getElementById('zFit').onclick=()=>zoomTo(1,0,0);
document.getElementById('zFocus').onclick=()=>focusAt(pinI>=0?pinI:lastPaintI);

// régua de profundidade (acompanha o zoom)
const ruler=document.getElementById('ruler');
function updRuler(){
  const top=ZY*DX,span=H/Z*DX,steps=[.1,.2,.5,1,2,5];let st=5;for(const s of steps)if(span/s<=8){st=s;break;}
  let h='';for(let d=Math.ceil(top/st)*st;d<=top+span+1e-9;d+=st){const p=(d-top)/span*100;
    const tr=p<3?'none':p>97?'translateY(-100%)':'';h+='<div style="top:'+p+'%'+(tr?';transform:'+tr:'')+'">'+(st<1?d.toFixed(1):d.toFixed(0))+' m</div>';}
  ruler.innerHTML=h;
}

// ================= ponteiro =================
let lastPaintI=-1;
const isFree=()=>drawMode==='free'||drawMode==='sq';
function cellOf(e){const r=cv.getBoundingClientRect();return [Math.floor((e.clientX-r.left)/r.width*W),Math.floor((e.clientY-r.top)/r.height*H),r];}
function updBrush(e){const [x,y,r]=cellOf(e),sr=stage.getBoundingClientRect();const s=(brush*2+1)/W*r.width;
  brushEl.style.display=panning||stage.classList.contains('pan')?'none':'block';brushEl.style.borderRadius=drawMode==='sq'?'2px':'50%';brushEl.style.width=brushEl.style.height=(isFree()||drawMode==='rectO'?s:Math.max(8,r.width/W))+'px';
  brushEl.style.left=(e.clientX-sr.left)+'px';brushEl.style.top=(e.clientY-sr.top)+'px';
  curX=x;curY=y;hoverI=(x>=0&&x<W&&y>=0&&y<H)?y*W+x:-1;}
stage.addEventListener('contextmenu',e=>e.preventDefault());
stage.addEventListener('pointerdown',e=>{stage.setPointerCapture(e.pointerId);updBrush(e);
  if(e.button===1||e.altKey){e.preventDefault();panning={x:e.clientX,y:e.clientY,zx:ZX,zy:ZY};zT=null;stage.classList.add('pan','drag');brushEl.style.display='none';return;}
  if(e.ctrlKey||e.metaKey){pinI=hoverI;return;}
  if(tool===WIND){windDX=0;windDY=0;}
  shiftKey=e.shiftKey;erasing=e.button===2;lastX=curX;lastY=curY;lastPaintI=hoverI;
  strokeId=strokeId>=65535?1:strokeId+1;pushUndo();
  const tl=erasing?ERASE:tool;
  if(drawMode==='fill'){bucket(curX,curY,tl);return;}
  drawing=true;
  if(!isFree()){shapePrev={md:drawMode,x0:curX,y0:curY,x1:curX,y1:curY,cx:curX,cy:curY,ow:shiftKey};return;}
  if(!CONT[tl]) paintAt(curX,curY,tl);});
stage.addEventListener('pointermove',e=>{
  if(panning){const r=stage.getBoundingClientRect();ZX=panning.zx-(e.clientX-panning.x)/r.width*W/Z;ZY=panning.zy-(e.clientY-panning.y)/r.height*H/Z;applyZoom();return;}
  stage.classList.toggle('pan',e.altKey);
  updBrush(e);shiftKey=e.shiftKey;
  if(shapePrev){shapePrev.cx=curX;shapePrev.cy=curY;constrainShape();return;}
  if(drawing){const tl=erasing?ERASE:tool;if(!CONT[tl])paintLine(lastX,lastY,curX,curY,tl);lastX=curX;lastY=curY;lastPaintI=hoverI;}});
// Shift durante arrasto: linha só horiz./vert.; retângulo/elipse quadrado/círculo. Sobrescrever em formas = Shift no clique inicial
function constrainShape(){const sp=shapePrev;if(!sp)return;let dx=sp.cx-sp.x0,dy=sp.cy-sp.y0;
  if(shiftKey){if(sp.md==='line'){if(Math.abs(dx)>=Math.abs(dy))dy=0;else dx=0;}
    else{const d=Math.max(Math.abs(dx),Math.abs(dy));dx=(dx<0?-1:1)*d;dy=(dy<0?-1:1)*d;}}
  sp.x1=sp.x0+dx;sp.y1=sp.y0+dy;}
const endDraw=()=>{if(shapePrev){const sp=shapePrev;shapePrev=null;const sk=shiftKey;shiftKey=sp.ow;paintShape(sp.md,sp.x0,sp.y0,sp.x1,sp.y1,erasing?ERASE:tool);shiftKey=sk;}drawing=false;if(panning){panning=null;stage.classList.remove('drag','pan');}};
stage.addEventListener('pointerup',endDraw);stage.addEventListener('pointercancel',endDraw);
stage.addEventListener('pointerleave',()=>{brushEl.style.display='none';hoverI=-1;});
stage.addEventListener('auxclick',e=>{if(e.button===1)e.preventDefault();});

// ================= inspetor =================
const readout=document.getElementById('readout'),statsEl=document.getElementById('stats'),inspMode=document.getElementById('inspMode');
const kv=(k,v)=>'<span class="k">'+k+'</span><span class="v">'+v+'</span>';
function updReadout(){
  const i=hoverI>=0?hoverI:pinI;
  inspMode.textContent=hoverI>=0?'cursor':pinI>=0?'fixado':'—';
  if(i<0){readout.innerHTML='<div class="hint">Passe o mouse sobre a câmara. <b>Ctrl+clique</b> fixa uma célula; <b>F</b> foca o objeto.</div>';return;}
  const m=mat[i],k=kind[m],x=i%W,y=i/W|0;
  let h='<div class="nm"><i style="background:'+css(D[m].col)+'"></i>'+D[m].n+'</div><div class="kv">';
  h+=kv('T',(T[i]-273.15).toFixed(1)+' °C');
  if(k!==FIX){
    h+=kv('P',fmtP(P[i]));
    h+=kv('P manom.',fmtP(P[i]-PATM).replace('-','−'));
    if(airOn){const c=AF[i],pa=ap[c];h+=kv('P din. (ar)',(pa<0?'−':'+')+fmtP(pa<0?-pa:pa));
      const s=Math.hypot(ucx[c],ucy[c])/AC;if(s>1e-3)h+=kv('vento','Mach '+(s/C_AIR).toFixed(2)+' ≈ '+(s/C_AIR*343).toFixed(0)+' m/s');}
    h+=kv('ρ',densOf(i).toFixed(k===GAS?3:0)+' kg/m³');
    h+=kv('região',vis[i]===1?'aberta (atm.)':'selada');
    if(hiTy[m]===2) h+=kv('ebulição',(tb(hiT[m],hiLm[m],hiCap[m],P[i])-273.15).toFixed(1)+' °C');
    if(loTy[m]===2) h+=kv('condensa <',(tb(loT[m],loLm[m],loCap[m],P[i])-273.15).toFixed(1)+' °C');
    if(hiTy[m]===3) h+=kv('fusão',(hiT[m]+hiLm[m]*(P[i]-PATM)-273.15).toFixed(3)+' °C');
    if(loTy[m]===3) h+=kv('congela <',(loT[m]+loLm[m]*(P[i]-PATM)-273.15).toFixed(3)+' °C');
  }else if(strA[m]>0){
    h+=kv('tensão',(st[i]*100).toFixed(1)+'% da resist.');
  }
  h+=kv('σ',fmtS(elecLive?sgE[i]:sigA[m])+' S/m');
  if(elecLive){
    h+=kv('V',fmtV(Vf[i]))+kv('|E|',fmtE(eMag(i)));
    if(chg[i]>1e-13||chg[i]<-1e-13) h+=kv('q',fmtQ(chg[i]));
    const I=Math.hypot(IxC[i],IyC[i]); if(I>1e-6) h+=kv('I',fmtI(I));
    if(ionL[i]>.02) h+=kv('ionização',(ionL[i]*100).toFixed(0)+'%');
  }
  if(magLive) h+=kv('|B|',fmtB(Bfine[i]));
  if(m===IMA) h+=kv('M',mxA[i]||myA[i]?si(Math.hypot(mxA[i],myA[i]),'A/m')+' '+(myA[i]<0?'↑':myA[i]>0?'↓':mxA[i]>0?'→':'←'):'desmagnetizado');
  h+=gravReadout(i);
  h+=radReadout(i,m);
  h+=kv('posição',(x*DX).toFixed(1)+' m · prof. '+(y*DX).toFixed(1)+' m');
  readout.innerHTML=h+'</div>';
}


// ================= histórico (desfazer/refazer) =================
const snapArrs=()=>({mat,shade,life,T,tOld,pp,gd,lat,st,P,chg,ion,mxA,myA,Vf,ionL,rA,bu,dAcc}); // P é trocado com Pn → pegar referência atual
const UNDO_MAX=20,undoS=[],redoS=[];
function snapshot(){const o={lid,storm,grav:gravGet()},a=snapArrs();for(const k in a)o[k]=a[k].slice();return o;}
function restore(sn){const a=snapArrs();for(const k in a)a[k].set(sn[k]);
  lid=sn.lid;updLidBtn();setStorm(sn.storm);gravSet(sn.grav);IxC.fill(0);IyC.fill(0);airReset();radReset();updHist();}
function pushUndo(){undoS.push(snapshot());if(undoS.length>UNDO_MAX)undoS.shift();redoS.length=0;updHist();}
function undo(){if(!undoS.length)return;redoS.push(snapshot());restore(undoS.pop());}
function redo(){if(!redoS.length)return;undoS.push(snapshot());restore(redoS.pop());}
const bUndo=document.getElementById('bUndo'),bRedo=document.getElementById('bRedo');
function updHist(){bUndo.disabled=!undoS.length;bRedo.disabled=!redoS.length;}
bUndo.onclick=undo;bRedo.onclick=redo;updHist();

// ================= salvar / carregar / captura =================
const SLIDERS=['sBrush','sGrav','sHeater','sPow','sSpeed','sVolt','sChg','sWind'];
function download(blob,name){const a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(a.href),2000);}
const stampName=()=>new Date().toISOString().slice(0,19).replace(/[:T]/g,'-');
async function saveScene(){
  const a=snapArrs(),meta=[],bufs=[];
  for(const k in a){meta.push([k,a[k].constructor.name,a[k].length]);bufs.push(a[k]);}
  const set={};for(const id of SLIDERS)set[id]=+document.getElementById(id).value;
  const hdr=new TextEncoder().encode(JSON.stringify({app:'camara',v:1,W,H,lid,storm,magDir,view,set,meta,air:airOn,blast:blastOn,grav:gravGet()}));
  const raw=new Blob([new Uint32Array([hdr.length]),hdr,...bufs]);
  const gz=await new Response(raw.stream().pipeThrough(new CompressionStream('gzip'))).blob();
  download(gz,'camara-'+stampName()+'.camara');
  toast('Cenário salvo');
}
async function loadScene(file){
  try{
    const buf=await new Response(file.stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
    const hl=new Uint32Array(buf.slice(0,4))[0],h=JSON.parse(new TextDecoder().decode(buf.slice(4,4+hl)));
    if(h.app!=='camara'||h.W!==W||h.H!==H)throw new Error('arquivo incompatível');
    const a=snapArrs();let off=4+hl;const TY={Uint8Array,Float32Array,Float64Array};
    pushUndo();
    for(const [k,ty,len] of h.meta){const C=TY[ty],nb=len*C.BYTES_PER_ELEMENT;if(!a[k]||len!==N)throw new Error('campo '+k);
      a[k].set(new C(buf.slice(off,off+nb)));off+=nb;}
    for(let i=0;i<N;i++)if(mat[i]>=NM)mat[i]=AR;
    for(const id in h.set){const el=document.getElementById(id);if(el){el.value=h.set[id];el.oninput();}}
    lid=!!h.lid;updLidBtn();setStorm(!!h.storm);setView(h.view|0);setAir(h.air!==false);setBlast(h.blast!==false);airReset();radReset();gravSet(h.grav);
    magDir=h.magDir|0;poleEl.querySelectorAll('button').forEach(x=>x.classList.toggle('on',+x.dataset.d===magDir));
    zoomTo(1,0,0);toast('Cenário carregado: '+file.name);
  }catch(err){toast('Falha ao carregar: '+err.message,true);}
}
function screenshot(){
  const S=1200,c=document.createElement('canvas');c.width=c.height=S;const x=c.getContext('2d');
  x.imageSmoothingEnabled=false;x.drawImage(cv,ZX,ZY,W/Z,H/Z,0,0,S,S);x.imageSmoothingEnabled=true;x.drawImage(ov,0,0,S,S);
  c.toBlob(b=>{download(b,'camara-'+stampName()+'.png');toast('Captura salva');},'image/png');
}
const fileIn=document.getElementById('fileIn');
fileIn.onchange=()=>{if(fileIn.files[0])loadScene(fileIn.files[0]);fileIn.value='';};
document.getElementById('bSave').onclick=saveScene;
document.getElementById('bLoad').onclick=()=>fileIn.click();
document.getElementById('bShot').onclick=screenshot;
stage.addEventListener('dragover',e=>{e.preventDefault();stage.classList.add('drop');});
stage.addEventListener('dragleave',()=>stage.classList.remove('drop'));
stage.addEventListener('drop',e=>{e.preventDefault();stage.classList.remove('drop');const f=e.dataTransfer.files[0];if(f)loadScene(f);});
const toastEl=document.getElementById('toast');let toastT=0;
function toast(msg,err){toastEl.textContent=msg;toastEl.className='toast show'+(err?' err':'');clearTimeout(toastT);toastT=setTimeout(()=>toastEl.className='toast',2400);}
function stepOnce(){if(!paused)togglePause();simStep();magStep();}
document.getElementById('bStep').onclick=stepOnce;
applyZoom();
// loop
function simStep(){gravStep();moveStep();flyStep();airStep();pressureStep();heatStep();stateStep();explodeStep();lifeStep();radStep();if(storm)weatherStep();elecStep();}
let last=performance.now(),fps=60,frameN=0,simMs=0;
function loop(now){
  const dt=now-last;last=now;fps+=(1000/Math.max(dt,1)-fps)*.05;
  const t0=performance.now();zoomAnim();
  if(!paused){
    if(drawing&&curX>=0&&isFree()){const tl=erasing?ERASE:tool;paintAt(curX,curY,tl);}
    stepAcc+=steps;while(stepAcc>=1){stepAcc--;simStep();}
  }else if(drawing&&curX>=0&&CONT[tool]&&!erasing&&isFree()){paintAt(curX,curY,tool);}
  if((frameN&3)===0) magStep();
  const t1=performance.now();simMs+=(t1-t0-simMs)*.1;
  render(); drawOverlay(frameN);
  if((frameN++&7)===0){
    updReadout();
    if(view===1) lgHi.textContent='+'+fmtP(pRange);
    statsEl.innerHTML='<span><b>'+fps.toFixed(0)+'</b> fps</span><span>sim <b>'+simMs.toFixed(1)+'</b> ms</span><span><b>'+cntP+'</b> partículas</span><span>Pmáx <b>+'+fmtP(maxG)+'</b></span><span>Tmáx <b>'+(cntP?(maxT-273.15).toFixed(0):'—')+'</b> °C</span>'+(elecLive?'<span>Vmáx <b>'+fmtV(vScale)+'</b></span><span>Imáx <b>'+fmtI(iMax)+'</b></span>':'')+(nStrikes?'<span>descargas <b>'+nStrikes+'</b></span>':'')+(magLive?'<span>Bmáx <b>'+fmtB(bMax)+'</b></span>':'')+radStats();
    if(view===3||view===4) updateLegend();
  }
  requestAnimationFrame(loop);
}
// inicialização em em.js (carregado depois)
