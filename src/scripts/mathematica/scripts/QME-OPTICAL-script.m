(* ::Package:: *)

Off[General::spell1];
Off[Flatten::normal];
$HistoryLength=0; (*saves memory and enables deallocation*)
meter=10^10angstrom;second = 10^15 fs;joule=5.0341182*10^(22)inverseCm;
kilogram = joule second^2/meter^2;

angstrom=1;fs=1;inverseCm=1;kelvin=1;

\[HBar]=6.626068*10^(-34)/2/\[Pi] meter^2 kilogram/second;
c=299792458.0 meter/second;
kB = 1.3806503*10^(-23) meter^2 kilogram second^(-2) kelvin^(-1);

c \[HBar] 2 \[Pi]/(meter/angstrom/100)

(*creates list of all pairs => used to build 2-excitonic hamiltonian*)
ListOfPairs:=Module[{i,j,tmp},
tmp={};
For[i=1,i<=MoleculeNumber,i++,
For[j=1,j<=MoleculeNumber,j++,
If[i<j,tmp=Union[tmp,{{i,j}}]]
];
];
Return[tmp];
];


\[CapitalLambda]\[CapitalLambda][s_]=1/\[Tau]c[s];

Print["{kB T, \[Lambda]\[Lambda], 2\[Pi]\[HBar]c/\[Tau]c, J, \[CapitalDelta]e} = ",{kB T, \[Lambda]\[Lambda], 2\[Pi] \[HBar] c/\[Tau]c, J, \[CapitalDelta]e}];
Print["Memory: ",MemoryInUse[]];

HsTmp = N[DiagonalMatrix[Join[Table[\[Epsilon][i],{i,MoleculeNumber}],Table[\[Epsilon][ListOfPairs[[i]][[1]]]+\[Epsilon][ListOfPairs[[i]][[2]]],{i,MoleculeNumber (MoleculeNumber-1)/2}]]]+Table[Table[If[i<=MoleculeNumber && j<=MoleculeNumber&&i!=j,J[i,j],0]+If[i>MoleculeNumber && j>MoleculeNumber&&i!=j,J2[i-MoleculeNumber,j-MoleculeNumber],0],{i,1,MoleculeNumber+MoleculeNumber (MoleculeNumber-1)/2}],{j,1,MoleculeNumber+MoleculeNumber (MoleculeNumber-1)/2}]];
Hs=Table[Table[If[i>0&&j>0,HsTmp[[i]][[j]],0],{i,0,Length[HsTmp]}],{j,0,Length[HsTmp]}];
HsTmp=.;
Print["Hs zaveden"];
U[i_]:=DiagonalMatrix[Table[If[j==i,1,0],{j,1,Length[Hs]}]];
Print["Projektory U[i] zavedeny"];


(*CrossCorel[r_,s_]:=
Module[{},
If[Max[r,s]>Length[ListOfPairs]||Min[r,s]<1,Return[0]];
Return[Intersection[Part[ListOfPairs,r],Part[ListOfPairs,s]]];
];*)
CrossCorel[r_,s_]:=
Module[{},
If[Max[r,s]>Length[ListOfPairs]||Min[r,s]<1,Return[{}]];
Return[Intersection[Part[ListOfPairs,r],Part[ListOfPairs,s]]];
];
CorWithoutOC=Table[Table[Union[{},CrossCorel[i-MoleculeNumber-1,j-MoleculeNumber-1],If[i==j&&j>1&&j<=MoleculeNumber+1&&i>1&&i<=MoleculeNumber+1,{j-1},{},{}]*KroneckerDelta[i,j]],{i,1,Length[Hs]}],{j,1,Length[Hs]}];
Cortmp[m_,n_]=Table[Table[Intersection[CorWithoutOC[[i]][[i]],CorWithoutOC[[j]][[j]]],{i,1,Length[Hs]}],{j,1,Length[Hs]}][[m]][[n]];

\[Beta][s_]:=1/kB/T[s];
\[Nu][n_,\[Beta]_]:=2\[Pi]/\[HBar]/\[Beta] *n;
suma\[Nu][\[Beta]\[Beta]_,\[CapitalLambda]\[CapitalLambda]\[CapitalLambda]_,t_]=Sum[\[Nu][n,\[Beta]\[Beta]]E^(-\[Nu][n,\[Beta]\[Beta]] t)/(\[Nu][n,\[Beta]\[Beta]]^2-\[CapitalLambda]\[CapitalLambda]\[CapitalLambda]^2),{n,1,MezCC}];
CCC[t_,\[CapitalLambda]\[CapitalLambda]_,\[Lambda]\[Lambda]_,\[Beta]_]:=-I \[Lambda]\[Lambda] \[CapitalLambda]\[CapitalLambda] E^(-\[CapitalLambda]\[CapitalLambda] Abs[t]) Sign[t] \[HBar];
(*\[HBar] zde predpokladam, v literature neni, jina definice?*)
CC[t_,\[CapitalLambda]\[CapitalLambda]_,\[Lambda]\[Lambda]_,\[Beta]_]:=(\[Lambda]\[Lambda] \[CapitalLambda]\[CapitalLambda] Coth[\[Beta] \[HBar] \[CapitalLambda]\[CapitalLambda]/2]E^(-\[CapitalLambda]\[CapitalLambda] Abs[t])+4\[CapitalLambda]\[CapitalLambda] \[Lambda]\[Lambda]/\[HBar]/\[Beta] suma\[Nu][\[Beta],\[CapitalLambda]\[CapitalLambda],Abs[t]])\[HBar];
(*\[HBar] zde predpokladam, v literature neni, jina definice?*)


(*Zavadeni korelacnich funkci*)
If[POUZITjPROcor,

Module[{s,NN,dt,dom, \[Omega]max,TRange,\[Omega]list,timelist,pokusOmega,pokusTime,Fft,iFft,FftShift},
Fft[x_]:=Fourier[x 10^12,FourierParameters->{1,-1}]/10^12;
iFft[x_]:=InverseFourier[x 10^12,FourierParameters->{1,-1}]/10^12;
FftShift[x_]:=RotateRight[x,(Length[x]-Mod[Length[x],2])/2];

NN=1024*8*2*4;
dt=1/1000/2;

TRange=dt NN ;
dom=2 \[Pi]/TRange;
\[Omega]max=NN/2*dom;
Print["\[Omega]max=",N[\[Omega]max]];
Print["TRange=",N[TRange \[HBar]]];

\[Omega]list=(Table[i-1,{i,1,NN}]-NN/2)dom;
timelist=(Table[dt i,{i,0,NN-1}]-NN/2 dt);

For[s=1,s<=MoleculeNumber,s=s+1,

pokusOmega[s]=2*\[Pi]*Table[2\[Pi] (\[Omega])^2((1+nTermo[-\[Omega],s])SpectralDensity[-\[Omega],s]UnitStep[-\[Omega]]+UnitStep[\[Omega]]nTermo[\[Omega],s]SpectralDensity[\[Omega],s])/.{\[Omega]->(i/NN*\[Omega]max-\[Omega]max/2)*2},{i,0,NN-1}]/.{Indeterminate->0};
(*pokusOmega[s]=2*\[Pi]*Table[2\[Pi] (\[Omega])^2((nTermo[-\[Omega],s])SpectralDensity[-\[Omega],s]UnitStep[-\[Omega]]+UnitStep[\[Omega]](1+nTermo[\[Omega],s])SpectralDensity[\[Omega],s])/.{\[Omega]->(i/NN*\[Omega]max-\[Omega]max/2)*2},{i,0,NN-1}]/.{Indeterminate->0};*)
pokusTime[s]=FftShift[iFft[FftShift[pokusOmega[s]]]]/dt/2/\[Pi];

NovaCor[s]=Interpolation[Table[{\[HBar] timelist[[i]], pokusTime[s][[i]]/\[HBar]},{i,1,NN}]];
];
];

Cor[m_,n_,t_]:=Sum[(NovaCor[Part[Cortmp[m,n],i]][t]),{i,1,Length[Cortmp[m,n]]}];
Print["Korelacni funkce zavedeny"];


,
Cor[m_,n_,t_]:=Sum[(CC[t,\[CapitalLambda]\[CapitalLambda][Part[Cortmp[m,n],i]],\[Lambda]\[Lambda][Part[Cortmp[m,n],i]],\[Beta][Part[Cortmp[m,n],i]]]+CCC[t,\[CapitalLambda]\[CapitalLambda][Part[Cortmp[m,n],i]],\[Lambda]\[Lambda][Part[Cortmp[m,n],i]],\[Beta][Part[Cortmp[m,n],i]]]),{i,1,Length[Cortmp[m,n]]}];
Print["Korelacni funkce zavedeny"];
,
Print["Neni korelacni fce"];
];


Commutator[A_,B_]:=A.B-B.A;
ToExcitonBasis[X_]:=Module[{tmp,vect,i,j,\[CapitalOmega],Hs1,Hs2,tmp1,tmp2,vect1,vect2,order1,order2,eig1,eig2,check},
Hs1=Table[Table[Hs[[j]][[i]],{i,2,MoleculeNumber+1}],{j,2,MoleculeNumber+1}];
Hs2=Table[Table[Hs[[j]][[i]],{i,MoleculeNumber+1+1,Length[Hs]}],{j,MoleculeNumber+1+1,Length[Hs]}]; 
eig1=Eigensystem[Hs1];
eig2=Eigensystem[Hs2];
vect1=N[eig1[[2]]];
vect2=N[eig2[[2]]];
tmp1=Table[(vect1)[[Length[vect1]-k+1]]/Norm[(vect1)[[Length[vect1]-k+1]]]*If[Abs[Max[(vect1)[[Length[vect1]-k+1]]]]>Abs[Min[(vect1)[[Length[vect1]-k+1]]]],1,-1,1],{k,1,Length[vect1]}]; tmp2=Table[(vect2)[[Length[vect2]-k+1]]/Norm[(vect2)[[Length[vect2]-k+1]]]*If[Abs[Max[(vect2)[[Length[vect2]-k+1]]]]>Abs[Min[(vect2)[[Length[vect2]-k+1]]]],1,-1,1],{k,1,Length[vect2]}]; 

check=Table[(tmp1.Hs1.Inverse[tmp1])[[i]][[i]],{i,1,Length[Hs1]}];
While[!(Norm[check-Sort[check]]<10^-6),
i=RandomInteger[{1,Length[tmp1]}];
j=RandomInteger[{1,Length[tmp1]}];
If[i==j,Continue[]];
tmp1=Table[tmp1[[If[k==i,j,If[k==j,i,k]]]],{k,1,Length[tmp1]}];
check=Table[(tmp1.Hs1.Inverse[tmp1])[[i]][[i]],{i,1,Length[Hs1]}];
];
check=Table[(tmp2.Hs2.Inverse[tmp2])[[i]][[i]],{i,1,Length[Hs2]}];
While[!(Norm[check-Sort[check]]<10^-6),
i=RandomInteger[{1,Length[tmp2]}];
j=RandomInteger[{1,Length[tmp2]}];
If[i==j,Continue[]];
tmp2=Table[tmp2[[If[k==i,j,If[k==j,i,k]]]],{k,1,Length[tmp2]}];
check=Table[(tmp2.Hs2.Inverse[tmp2])[[i]][[i]],{i,1,Length[Hs2]}];
];
(*Print[tmp1];

order1=Table[{eig1[[1]][[i]],tmp1[[i]]},{i,1,Length[tmp1]}];
Print[order1];
order1=Sort[order1,#1[[1]]<#2[[1]]&];
Print[order1];
tmp1=Table[order1[[i]][[2]],{i,1,Length[order1]}];
Print[tmp1];

order2=Table[{eig2[[1]][[i]],tmp2[[i]]},{i,1,Length[tmp2]}];
order2=Sort[order2,#1[[1]]<#2[[1]]&];
tmp2=Table[order2[[i]][[2]],{i,1,Length[order2]}];*)
tmp=Table[Table[If[i==1&&j==1,1,0,0]+
If[j>1&&i>1&&(i<=MoleculeNumber+1)&&(j<=MoleculeNumber+1),tmp1[[j-1]][[i-1]],0,0]+If[(i>MoleculeNumber+1)&&(j>MoleculeNumber+1),tmp2[[j-1-MoleculeNumber]][[i-1-MoleculeNumber]],0,0],{i,1,Length[Hs]}],{j,1,Length[Hs]}];
FromExciton2x2=Inverse[tmp];
ToExciton2x2=tmp;
Return[ToExciton2x2.X.FromExciton2x2];
];
ToExcitonBasis[Hs];
(*FromExcitonBasis[X_]:=Module[{tmp,vect,i,j}, (*prevadi do excitonove baze*)
vect=N[Eigensystem[Hs][[2]]];     (*Protoze podobnostni transformace je nejednoznacna    *)
vect=SortBy[N[Eigensystem[Hs][[2]]],Dot[Abs[#],Join[Table[1,{j,1,MoleculeNumber}],Table[10^6,{j,1,Length[Hs]-MoleculeNumber}]]]&]; (*sorted by 1/2 excitons*)
tmp=Transpose[Table[(vect)[[k]]/Norm[(vect)[[k]]],{k,1,Length[vect]}]]; (*normuji vl. vektory*)
Return[tmp.X.Inverse[tmp]];
];*)
PartMatr[A_,i_,j_]:=(IdentityMatrix[Length[A]][[i]]).A.(IdentityMatrix[Length[A]][[j]]);
Print["To/FromExcitonBasis zavedeny"];

\[Rho]\[Rho][k_]=Table[Table[If[k==(jjj-1)+Length[Hs]*(jj-1) + 1,1,0],{jjj,1,Length[Hs]}],{jj,1,Length[Hs]}];
exc[i_]:=If[i==1,0,If[i-1<=MoleculeNumber,1,2]];
To4x4FormalismVector[matice_,exc1_,exc2_]:=Module[{pom,eps,ep,ret},  
pom=Table[Table[If[(exc[l]==exc1&&exc[k]==exc2)||exc1==-1||exc1==-2,1,0],{k,1,Length[Hs]}],{l,1,Length[Hs]}];
eps=Table[Table[ep,{k,1,Length[Hs]}],{l,1,Length[Hs]}];
ret=DeleteCases[Flatten[pom (matice+eps)],0];
ep=0;
Return[ret];
];
 (*Prevadi lin. operaci do 4x4 formalismu (=Flatten)*) 
To4x4Formalism[\[Rho]_,operace_,exc1_,exc2_]:=Module[{tabindex,ep,len},                    
tabindex=DeleteCases[Table[k Norm[To4x4FormalismVector[\[Rho]\[Rho][k],exc1,exc2]],{k,1,Length[Hs]^2}],0];
len=Length[To4x4FormalismVector[Hs,exc1,exc2]];
Return[                                                                      
Table[Table[Flatten[((operace)/.{\[Rho]->Evaluate[\[Rho]\[Rho][tabindex[[k]]]] })][[tabindex[[l]]]],{k,1,Length[tabindex]}],{l,1,Length[tabindex]}]];
];
\[Rho]\[Rho]=.;
Print["To4x4Formalism zavedeno"];


time1=AbsoluteTime[];
(*JordanDecomposition slow - preventing doing it for each superoperator element*)
vect=N[Eigensystem[Hs][[2]]];    
vect=SortBy[N[Eigensystem[Hs][[2]]],Dot[Abs[#],Join[Table[1,{j,1,MoleculeNumber}],Table[10^6,{j,1,Length[Hs]-MoleculeNumber}]]]&]; (*sorted by 1/2 excitons*)
tmp=Transpose[Table[(vect)[[k]]/Norm[(vect)[[k]]],{k,1,Length[vect]}]]; 


Module[{i},
For[i=1,i<=Length[Hs],i++,
Ue[i]=Evaluate[Chop[ToExciton2x2.U[i].FromExciton2x2]];
]];
Us[t_]=Chop[tmp.Evaluate[E^(-Inverse[tmp].Hs.tmp*I/\[HBar] t)DiagonalMatrix[Table[1,{j,1,Length[Hs]}]]].Inverse[tmp],10^-11];
Use[t_]=Chop[Expand[ToExciton2x2.Us[t].FromExciton2x2]];
Print["Us[t] zavedeno"];


ToExciton=To4x4Formalism[S,Inverse[tmp].S.tmp];
FromExciton=Inverse[ToExciton];
vect=.;tmp=.;
time2=AbsoluteTime[];
Print["To/FromExciton,Us[t] zavedeno, ",time2-time1," s"];


(*Pro minimalizaci poctu chyb jednoduse vynecham integraci pres cas u lambd, cimz ziskam primo superoperator M*)
time1=AbsoluteTime[];
exc11[i_]:=If[i==1,0,If[i-1<=MoleculeNumber,1,1]]; 
\[CapitalLambda][i_]=Table[Table[\[Lambda][i][j][k][t],{j,1,Length[Hs]}],{k,1,Length[Hs]}]Table[Table[If[exc11[i]==exc11[j]&&exc11[j]==exc11[k],1,0],{j,1,Length[Hs]}],{k,1,Length[Hs]}];
(*Projektory na jednotlive komponenty matic, pro sestaveni rovnic.*)
Projector[i_]=Table[If[k==i,1,0],{k,1,Length[Hs]}];

epsilon=0 10.0^-2;
Equations\[CapitalLambda][s_]:=Module[{i,ret,tmp,\[CapitalLambda]t,rs},
ret={};
rs=Sum[Cor[s,n,t+epsilon]Use[t].Ue[n].Use[-t],{n,1,Length[Hs]}];
For[i=1,i<=Length[Hs]^2,i++,
If[exc11[Mod[i-1,Length[Hs]] +1]==exc11[(i-1-Mod[i-1,Length[Hs]])/Length[Hs] +1]&&exc11[Mod[i-1,Length[Hs]] +1]==exc11[s],
ret=Join[ret,
{Projector[Mod[i-1,Length[Hs]] +1].\[CapitalLambda][s].Projector[(i-1-Mod[i-1,Length[Hs]])/Length[Hs] +1]->Projector[Mod[i-1,Length[Hs]] +1].rs.Projector[(i-1-Mod[i-1,Length[Hs]])/Length[Hs] +1]}]
];
];
Return[ret];
];
Module[{i},
Equations={};
For[i=1,i<=Length[Hs],i++,
Equations=Join[Equations,Equations\[CapitalLambda][i]]]];
Print["\[Lambda]ambdas zavedeno"];

Hm:=-Sum[I/\[HBar] Ue[i].\[CapitalLambda][i],{i,1,Length[Hs]}];
exxc1=0;
exxc2=1;
MMM[t_]=Evaluate[Simplify[((To4x4Formalism[\[Rho][t],-I/\[HBar] (Hm.\[Rho][t]-\[Rho][t].Transpose[Conjugate[Hm]])+1/\[HBar]^2Sum[Ue[j].\[Rho][t].Transpose[Conjugate[\[CapitalLambda][j]]]+\[CapitalLambda][j].\[Rho][t].Ue[j],{j,1,Length[Hs]}],exxc1,exxc2]/.Equations).To4x4Formalism[S,Use[t].S.Use[-t],exxc1,exxc2]),t>0]];
MMM[t_]=If[SECULAR,Table[Table[If[i==j,1,0],{i,1,Length[MMM[t]]}],{j,1,Length[MMM[t]]}]*MMM[t],MMM[t],MMM[t]];


Print["MMM[t] zavedeno"];
LLL=Chop[Evaluate[Expand[To4x4Formalism[\[Rho][t],-I/\[HBar] (ToExciton2x2.Hs.FromExciton2x2.\[Rho][t]-\[Rho][t].Transpose[Conjugate[ToExciton2x2.Hs.FromExciton2x2]]),exxc1,exxc2]]]];
(*oprava na e^I\[CapitalOmega]*)
LLL=LLL-IdentityMatrix[Length[LLL]]I/\[HBar]*\[CapitalOmega];
MMM[t_]=Chop[Expand[MMM[t]E^(-I/\[HBar] \[CapitalOmega] t)]];

Print["LLL zavedeno"];


Fft[x_]:=Fourier[x 10^12,FourierParameters->{1,-1}]/10^12;
iFft[x_]:=InverseFourier[x 10^12,FourierParameters->{1,-1}]/10^12;
FftShift[x_]:=RotateRight[x,(Length[x]-Mod[Length[x],2])/2];
FftMat[ListOfA_]:=Module[{k,i,j,ret},
For[i=1,i<=Length[ListOfA[[1]]],i++,
For[j=1,j<=Length[ListOfA[[1]]],j++,
ret[i,j]=Fft[Module[{k},Table[ListOfA[[k]][[i]][[j]],{k,1,Length[ListOfA]}]]];
];
];
Table[Table[Table[ret[j,i][[k]],{i,1,Length[ListOfA[[1]]]}],{j,1,Length[ListOfA[[1]]]}],{k,1,Length[ListOfA]}]
];
iFftMat[ListOfA_]:=Module[{k,i,j,ret},
For[i=1,i<=Length[ListOfA[[1]]],i++,
For[j=1,j<=Length[ListOfA[[1]]],j++,
ret[i,j]=iFft[Module[{k},Table[ListOfA[[k]][[i]][[j]],{k,1,Length[ListOfA]}]]];
];
];
Table[Table[Table[ret[j,i][[k]],{i,1,Length[ListOfA[[1]]]}],{j,1,Length[ListOfA[[1]]]}],{k,1,Length[ListOfA]}]
];

Gg[gg_,n_]:=Inverse[(I Part[\[Omega]list,n]+gg)IdentityMatrix[Length[LLL]]-Table[Table[Part[AaA[r][s],n],{s,1,Length[LLL]}],{r,1,Length[LLL]}]]/dt


Tmax=dt NN;
dom=2 \[Pi]/Tmax;

\[Omega]list=(Table[i-1,{i,1,NN}]-NN/2)dom;
timelist=FftShift[Table[dt i,{i,0,NN-1}]-NN/2 dt];
gg=1/(Tmax/8);


Print["Fft initiated"];
MemoryInUse[]
(*casove zavisle MMM jako list, vcetne E^(gg t)*)
(*listOfMMM=Table[If[0<=timelist[[i]]<20\[Tau]c[1],Chop[MMM[timelist[[i]]],10^-12]If[timelist[[i]]>=0,1,0],0IdentityMatrix[Length[MMM[t]]^2]] E^(-gg timelist[[i]]),{i,1,NN}];*)
time1=AbsoluteTime[];
Module[{r,s,len,vyraz},
len=Length[MMM[t]];
For[r=1,r<=len,r++,For[s=1,s<=len,s++,
vyraz[t_]=MMM[t][[r]][[s]];
listOfM[r][s]=Table[If[0<timelist[[i]],Chop[vyraz[timelist[[i]]],10^-12],0] E^(-gg timelist[[i]]),{i,1,NN}];
];];
];
time2=AbsoluteTime[];
Print["MMM[t]\[ExponentialE]^(-gg t) zavedeno, ",time2-time1," s, "];
MemoryInUse[]

time1=AbsoluteTime[];
(*listOfFourieredMMM=FftShift[FftMat[listOfMMM]];*)
Module[{r,s,len,vyraz},
len=Length[MMM[t]];
For[r=1,r<=len,r++,For[s=1,s<=len,s++,
listOfFourieredM[r][s]=FftShift[Fft[listOfM[r][s]]];
];];
];
(*listOfMMM=.;*)
time2=AbsoluteTime[];
Print["MMM[\[Omega]] zavedeno, ",time2-time1," s, ",ByteCount[listOfFourieredMMM]];
(*LLLlist=Table[LLL,{l,1,NN}];*)

(*Initial conditions for laser pulse*)
(*Module[{\[Mu],i,j,\[Rho]G,E0},
E0=0.4;
\[Rho]G=Table[Table[If[i==1&&j==1,1,0],{i,1,Length[Hs]}],{j,1,Length[Hs]}];
\[Mu]=Table[Table[If[(exc[i]==1&&exc[j]==0)||(exc[j]==1&&exc[i]==0),1,0],{i,1,Length[Hs]}],{j,1,Length[Hs]}];
\[Rho]0=\[Rho]G+E0 I Commutator[-\[Mu],\[Rho]G]-E0^2/2 Commutator[-\[Mu],Commutator[-\[Mu],\[Rho]G]];
\[Rho]0=Chop[ToExcitonBasis[\[Rho]0]];
Print[MatrixForm[\[Rho]0]];
\[Rho]0=To4x4FormalismVector[\[Rho]0,exxc1,exxc2];
];
\[Rho]0={1,0,0}*)


AaA[r_][s_]:=listOfFourieredM[r][s] dt+LLL[[r]][[s]];

lenred=Length[To4x4FormalismVector[Hs,exxc1,exxc2]];
time1=AbsoluteTime[];
Module[{i},For[i=1,i<=NN,i++,
MemoryWaste[i]=(Gg[gg, i].\[Rho]0);
];];
Print["*, ",AbsoluteTime[]-time1," s"];
time1=AbsoluteTime[];
For[s=1,s<=Length[\[Rho]0],s++,
(*If[Length[Intersection[{s},Table[(lenred^(1/2)+1)(s-1)+1,{s,1,lenred^(1/2)}]]]==0,Continue[]];*)
res[s]=iFft[FftShift[Table[(MemoryWaste[i])[[s]],{i,1,NN}]]];
res[s]=res[s] E^(timelist gg);
(*oprava na e^I\[CapitalOmega]*)
ress[s]=res[s]*Table[E^(-I/\[HBar]  \[CapitalOmega] i dt),{i,1,NN}];
(*Print["@@@Oq",s," ",Re[Drop[res[s],-NN/2]],Im[Drop[res[s],-NN/2]]];*)
If[SECULAR,
Print["@@@OQ",s," ",Re[Drop[res[s],-NN/2]],Im[Drop[res[s],-NN/2]]],
Print["@@@Oq",s," ",Re[Drop[res[s],-NN/2]],Im[Drop[res[s],-NN/2]]],
Print["@@@Oq",s," ",Re[Drop[res[s],-NN/2]],Im[Drop[res[s],-NN/2]]]
];

Print[s," ",MemoryInUse[]," ",time2-time1," s"];
]; s=.;
time2=AbsoluteTime[];
Print["Fourier, ",time2-time1," s"];

Quit[];


(*ListPlot[ Re[res[1]],PlotRange->All,Joined->True]*)
(*ListPlot[ Im[res[1]],PlotRange->All,Joined->True]*)
(*ListPlot[ Re[res[2]],PlotRange->All,Joined->True]*)
(*ListPlot[ Im[res[2]],PlotRange->All,Joined->True]*)
(*ListPlot[ Re[res[3]],PlotRange->All,Joined->True]*)
(*ListPlot[ Im[res[3]],PlotRange->All,Joined->True]*)
(**)
(*ListPlot[ Re[ress[1]],PlotRange->All,Joined->True]*)
(*ListPlot[ Im[ress[1]],PlotRange->All,Joined->True]*)
(*ListPlot[ Re[ress[2]],PlotRange->All,Joined->True]*)
(*ListPlot[ Im[ress[2]],PlotRange->All,Joined->True]*)
(*ListPlot[ Re[ress[3]],PlotRange->All,Joined->True]*)
(*ListPlot[ Im[ress[3]],PlotRange->All,Joined->True]*)
(**)
(*ListPlot[ Re[Drop[res[1],-3 NN/4]],PlotRange->All,Joined->True]*)
(*ListPlot[ Re[Drop[res[2],-3 NN/4]],PlotRange->All,Joined->True]*)
(*ListPlot[ Re[Drop[res[3],-3 NN/4]],PlotRange->All,Joined->True]*)



