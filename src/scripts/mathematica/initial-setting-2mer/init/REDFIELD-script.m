(* ::Package:: *)

Off[General::spell1];
Off[Flatten::normal];
Off[$MaxPrecision::"preclck"];
Off[$MinPrecision::"preclck"];
Off[N::"preclg"];
Off[SetPrecision::"preclg"];
$HistoryLength=2; (*saves memory and enables deallocation*)
meter=10^10angstrom;second = 10^15 fs;joule=5.0341182*10^(22)inverseCm;
kilogram = joule second^2/meter^2;

angstrom=1;fs=1;inverseCm=1;kelvin=1;

\[HBar]=6.626068*10^(-34)/2/\[Pi] meter^2 kilogram/second;
c=299792458.0 meter/second;
kB = 1.3806503*10^(-23) meter^2 kilogram second^(-2) kelvin^(-1);

c \[HBar] 2 \[Pi]/(meter/angstrom);

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


HsTmp = N[DiagonalMatrix[Join[Table[\[Epsilon][i],{i,MoleculeNumber}],Table[\[Epsilon][ListOfPairs[[i]][[1]]]+\[Epsilon][ListOfPairs[[i]][[2]]],{i,MoleculeNumber (MoleculeNumber-1)/2}]]]+Table[Table[If[i<=MoleculeNumber && j<=MoleculeNumber&&i!=j,J[i,j],0]+If[i>MoleculeNumber && j>MoleculeNumber&&i!=j,J2[i-MoleculeNumber,j-MoleculeNumber],0],{i,1,MoleculeNumber+MoleculeNumber (MoleculeNumber-1)/2}],{j,1,MoleculeNumber+MoleculeNumber (MoleculeNumber-1)/2}]];
Hs=Table[Table[If[i>0&&j>0,HsTmp[[i]][[j]],0],{i,0,Length[HsTmp]}],{j,0,Length[HsTmp]}];
HsTmp=.;
Print["Hs zaveden"];
U[i_]:=DiagonalMatrix[Table[If[j==i,1,0],{j,1,Length[Hs]}]];
Print["Projektory U[i] zavedeny"];


(*
CrossCorel[r_,s_]:=
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
tmpCor[m_,n_]=Table[Table[Intersection[CorWithoutOC[[i]][[i]],CorWithoutOC[[j]][[j]]],{i,1,Length[Hs]}],{j,1,Length[Hs]}][[m]][[n]];

\[Beta][s_]:=1/kB/T[s];
\[Nu][n_,\[Beta]_]:=2\[Pi]/\[HBar]/\[Beta] *n;
suma\[Nu][\[Beta]\[Beta]_,\[CapitalLambda]\[CapitalLambda]\[CapitalLambda]_,t_]=Sum[\[Nu][n,\[Beta]\[Beta]]E^(-\[Nu][n,\[Beta]\[Beta]] t)/(\[Nu][n,\[Beta]\[Beta]]^2-\[CapitalLambda]\[CapitalLambda]\[CapitalLambda]^2),{n,1,5}];
CCC[t_,\[CapitalLambda]\[CapitalLambda]_,\[Lambda]\[Lambda]_,\[Beta]_]:=-I \[Lambda]\[Lambda] \[CapitalLambda]\[CapitalLambda] E^(-\[CapitalLambda]\[CapitalLambda] Abs[t]) Abs[t]/t \[HBar];
(*\[HBar] zde predpokladam, v literature neni, jina definice?*)
CC[t_,\[CapitalLambda]\[CapitalLambda]_,\[Lambda]\[Lambda]_,\[Beta]_]:=(\[Lambda]\[Lambda] \[CapitalLambda]\[CapitalLambda] Coth[\[Beta] \[HBar] \[CapitalLambda]\[CapitalLambda]/2]E^(-\[CapitalLambda]\[CapitalLambda] Abs[t])+4\[CapitalLambda]\[CapitalLambda] \[Lambda]\[Lambda]/\[HBar]/\[Beta] suma\[Nu][\[Beta],\[CapitalLambda]\[CapitalLambda],Abs[t]])\[HBar];
(*\[HBar] zde predpokladam, v literature neni, jina definice?*)
(*tmpCor[m_,n_]:=If[m==1&&n==1,0,
If[Min[m,n]-1>MoleculeNumber,CrossCorel[m-1-MoleculeNumber,n-1-MoleculeNumber],DeleteCases[{(m-1) KroneckerDelta[m-1,n-1]},0]]];*)
Cor[m_,n_,t_]:=Sum[(CC[t,\[CapitalLambda]\[CapitalLambda][Part[tmpCor[m, n],i]],\[Lambda]\[Lambda][Part[tmpCor[m, n],i]],\[Beta][Part[tmpCor[m, n],i]]]+CCC[t,\[CapitalLambda]\[CapitalLambda][Part[tmpCor[m, n],i]],\[Lambda]\[Lambda][Part[tmpCor[m, n],i]],\[Beta][Part[tmpCor[m, n],i]]]),{i,1,Length[tmpCor[m, n]]}];
tmp=.
Print["Korelacni funkce zavedeny"];

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
(*ToExcitonBasis[X_]:=Module[{tmp,vect,i,j},
vect=N[Eigensystem[Hs][[2]]];     (*Protoze podobnostni transformace je nejednoznacna    *)
(*vect=SortBy[N[Eigensystem[Hs][[2]]],Dot[Abs[#],Join[Table[1,{j,1,MoleculeNumber}],Table[10^6,{j,1,Length[Hs]-MoleculeNumber}]]]&]; (*sorted by 1/2 excitons*)*)
tmp=Transpose[Table[(-vect)[[Length[vect]-k+1]]/Norm[(vect)[[Length[vect]-k+1]]],{k,1,Length[vect]}]]; (*normuji vl. vektory*)
ToExciton2x2=Inverse[tmp];
FromExciton2x2=tmp;
Return[Inverse[tmp].X.tmp];
];
ToExcitonBasis[Hs];*)
(*FromExcitonBasis[X_]:=Module[{tmp,vect,i,j}, (*prevadi do excitonove baze*)
vect=N[Eigensystem[Hs][[2]]];     (*Protoze podobnostni transformace je nejednoznacna    *)
vect=SortBy[N[Eigensystem[Hs][[2]]],Dot[Abs[#],Join[Table[1,{j,1,MoleculeNumber}],Table[10^6,{j,1,Length[Hs]-MoleculeNumber}]]]&]; (*sorted by 1/2 excitons*)
tmp=Transpose[Table[(vect)[[k]]/Norm[(vect)[[k]]],{k,1,Length[vect]}]]; (*normuji vl. vektory*)
Return[tmp.X.Inverse[tmp]];
];*)
PartMatr[A_,i_,j_]:=(IdentityMatrix[Length[A]][[i]]).A.(IdentityMatrix[Length[A]][[j]]);
Print["To/FromExcitonBasis zavedeny"];

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
Return[                                                                      Table[Table[Flatten[((operace)/.{\[Rho]->Evaluate[\[Rho]\[Rho][tabindex[[k]]]] })][[tabindex[[l]]]],{k,1,Length[tabindex]}],{l,1,Length[tabindex]}]];
];

Module[{i},
For[i=1,i<=Length[Hs],i++,
Ue[i]=Evaluate[Chop[ToExciton2x2.U[i].FromExciton2x2]];
]];
Us[t_]=E^(-Chop[ToExciton2x2.Hs.FromExciton2x2]I/\[HBar] t)DiagonalMatrix[Table[1,{j,1,Length[Hs]}]];
Print["Us[t] zavedeno"];

P[ii_,jj_]=Table[Table[If[j==jj&&i==ii,1,0],{i,1,Length[Hs]}],{j,1,Length[Hs]}];
\[Delta]elta=IdentityMatrix[Length[Hs]];
ones=Table[Table[1,{i,1,Length[Hs]}],{j,1,Length[Hs]}];
Secular2x2[AA_,S_,BA_]:=(AA.(S*\[Delta]elta).BA)*\[Delta]elta+(AA \[Delta]elta).(S*(ones-\[Delta]elta)).(BA \[Delta]elta);
Print["Secular2x2 zavedeno"];


(*Constructing list for solving \[CapitalLambda]ambdas*)
(*Carefully, |g> has index 1, not 0 !!*)
exc[i_]:=If[i==1,0,If[i-1<=MoleculeNumber,1,1]];

Table[
\[CapitalLambda][i]=Table[Table[\[Lambda][i][j][k][t],{j,1,Length[Hs]}],{k,1,Length[Hs]}]Table[Table[If[(exc[i]==exc[j]&&exc[j]==exc[k])(*||(exc[i]>0&&exc[j]>0&&exc[k]>0)*),1,0],{j,1,Length[Hs]}],{k,1,Length[Hs]}];
,{i,1,Length[Hs]}];
Projector[i_]=Table[If[k==i,1,0],{k,1,Length[Hs]}];


epsilon=10^-2;
time1=AbsoluteTime[];
Equations\[CapitalLambda][s_]:=Module[{i,ret,tmp,\[CapitalLambda]t,rs},
ret={};
\[CapitalLambda]t=D[\[CapitalLambda][s],t];
(*Rewrite ToExcitonBasis explicitely for higher speed*)
rs=Sum[Cor[s,n,t+epsilon]Us[t].Ue[n].Us[-t],{n,1,Length[Hs]}];
For[i=1,i<=Length[Hs]^2,i++,
If[exc[Mod[i-1,Length[Hs]] +1]==exc[(i-1-Mod[i-1,Length[Hs]])/Length[Hs] +1]&&exc[Mod[i-1,Length[Hs]] +1]==exc[s],
ret=Join[ret,
{Projector[Mod[i-1,Length[Hs]] +1].\[CapitalLambda]t.Projector[(i-1-Mod[i-1,Length[Hs]])/Length[Hs] +1]==Projector[Mod[i-1,Length[Hs]] +1].rs.Projector[(i-1-Mod[i-1,Length[Hs]])/Length[Hs] +1]}]
];
];
Return[ret];
];
Module[{i},
Equations={};
For[i=1,i<=Length[Hs],i++,
Equations=Join[Equations,Equations\[CapitalLambda][i]];
]];
time2=AbsoluteTime[];
Print["List of equations created, ",time2-time1," s"];
Print["* Equations: ",ByteCount[Equations]];
(*Initial=Join[{\[Lambda][1][1][1][0]==0},Flatten[Table[Table[Table[\[Lambda][i][j][k][0]==0,{i,1+1,MoleculeNumber+1}],{j,1+1,MoleculeNumber+1}],{k,1+1,MoleculeNumber+1}]],Flatten[Table[Table[Table[\[Lambda][i][j][k][0]==0,{i,MoleculeNumber+1+1,Length[Hs]}],{j,MoleculeNumber+1+1,Length[Hs]}],{k,MoleculeNumber+1+1,Length[Hs]}]]];
Lambdas=Join[{\[Lambda][1][1][1][t]},Flatten[Table[Table[Table[\[Lambda][i][j][k][t],{i,1+1,MoleculeNumber+1}],{j,1+1,MoleculeNumber+1}],{k,1+1,MoleculeNumber+1}]],Flatten[Table[Table[Table[\[Lambda][i][j][k][t],{i,MoleculeNumber+1+1,Length[Hs]}],{j,MoleculeNumber+1+1,Length[Hs]}],{k,MoleculeNumber+1+1,Length[Hs]}]]];*)
(*zde pokracuj*)
Lambdas=DeleteCases[Union[Flatten[Table[\[CapitalLambda][i],{i,1,Length[Hs]}]]],0];
Initial=Table[(Lambdas[[i]]/.{t->0})==0,{i,1,Length[Lambdas]}];
Print["Other lists created"];
time1=AbsoluteTime[];
\[Lambda]ambdas=NDSolve[Join[Equations,Initial],Lambdas,{t,0,TIME}(*,WorkingPrecision->7*)][[1]];
exc[i_]:=If[i==1,0,If[i-1<=MoleculeNumber,1,2]];
time2=AbsoluteTime[];
Print["Equation for \[CapitalLambda] solved, ",time2-time1, " s"];
Print["* \[Lambda]ambdas:   ",ByteCount[\[Lambda]ambdas]];

Heff=Chop[ToExciton2x2.Hs.FromExciton2x2]-Sum[I/\[HBar] Ue[i].\[CapitalLambda][i],{i,1,Length[Hs]}];
\[Rho][t_]=Table[Table[\[Rho]\[Rho][i][j][t],{i,1,Length[Hs]}],{j,1,Length[Hs]}];
time1=AbsoluteTime[];
Print["1"];

EqnList=Flatten[Table[Table[Chop[Expand[(-\[Rho]'[t]+ Evaluate[-I/\[HBar] (Heff.\[Rho][t]-\[Rho][t].Transpose[Conjugate[Heff]])+1/\[HBar]^2 Sum[Ue[j].\[Rho][t].Transpose[Conjugate[\[CapitalLambda][j]]]+\[CapitalLambda][j].\[Rho][t].Ue[j],{j,1,Length[Hs]}]]/.\[Lambda]ambdas)[[ii]][[jj]]]==0,10^(-12)],{ii,1,Length[Hs]}],{jj,1,Length[Hs]}]];
Print["1.5"];
(*Module[{\[Mu],i,j,\[Rho]G,E0},
E0=0.4;
\[Rho]G=Table[Table[If[i==1&&j==1,1,0],{i,1,Length[Hs]}],{j,1,Length[Hs]}];
\[Mu]=Table[Table[If[(exc[i]==1&&exc[j]==0)||(exc[j]==1&&exc[i]==0),1,0],{i,1,Length[Hs]}],{j,1,Length[Hs]}];
tmp=\[Rho]G+E0 I Commutator[-\[Mu],\[Rho]G]-E0^2/2 Commutator[-\[Mu],Commutator[-\[Mu],\[Rho]G]];
Print[MatrixForm[tmp]];
];*)
(*tmp=Table[Table[If[(i==1&&j==3)||(i==3&&j==1),1,0],{i,1,Length[Hs]}],{j,1,Length[Hs]}];*)
InitList=Flatten[Table[Table[\[Rho][0][[ii]][[jj]]==(\[Rho]0(*ToExciton2x2.tmp.FromExciton2x2*))[[ii]][[jj]],{ii,1,Length[Hs]}],{jj,1,Length[Hs]}]];
(*tmp=.;*)
Print["3"];

reseni=NDSolve[{Join[EqnList,InitList]},Flatten[\[Rho][t]],{t,0,TIME},MaxSteps->10^6][[1]];
Print["4"];
time2=AbsoluteTime[];
Print["2x2 reseni calculated, time ",time2-time1, " s"];
Print["* reseni:    ",ByteCount[reseni]];


time1=AbsoluteTime[];
EqnList=Flatten[Table[Table[Chop[Expand[(-\[Rho]'[t]+ Evaluate[-I/\[HBar] (Secular2x2[Heff,\[Rho][t],\[Delta]elta]-Secular2x2[\[Delta]elta,\[Rho][t],Transpose[Conjugate[Heff]]])+1/\[HBar]^2Sum[Secular2x2[Ue[j],\[Rho][t],Transpose[Conjugate[\[CapitalLambda][j]]]]+Secular2x2[\[CapitalLambda][j],\[Rho][t],Ue[j]],{j,1,Length[Hs]}]/.\[Lambda]ambdas])[[ii]][[jj]]]==0,10^-12],{ii,1,Length[Hs]}],{jj,1,Length[Hs]}]];
(*InitList=Flatten[Table[Table[\[Rho][0][[ii]][[jj]]==(ToExciton2x2.Table[Table[If[i==1&&j==1,1,0],{i,1,Length[Hs]}],{j,1,Length[Hs]}].FromExciton2x2)[[ii]][[jj]],{ii,1,Length[Hs]}],{jj,1,Length[Hs]}]];*)
resenisec=NDSolve[{Join[EqnList,InitList]},Flatten[\[Rho][t]],{t,0,TIME},MaxSteps->10^6][[1]];
time2=AbsoluteTime[];
Print["secular 2x2 reseni calculated, time ",time2-time1, " s"];
Print["* resenisec: ",ByteCount[resenisec]];
time1=.;time2=.
epsilon=.;
Print["Memory: ",MemoryInUse[]];


Module[{tmp,s},
tmp=To4x4FormalismVector[\[Rho][t],0,1];
(*Print[Plot[Evaluate[Flatten[Re[tmp E^(-I/\[HBar] *\[CapitalOmega] t)]/.reseni]],{t,0,TIME},PlotRange->All,AxesLabel->{"t [fs]","p"}]];*)
For[s=1,s<=Length[tmp],s++,
Print["@@@Or",s," ",Re[Table[Evaluate[(tmp[[s]]E^(-I/\[HBar] *\[CapitalOmega] t)/.reseni)/.t->i/(NN/2) *TIME],{i,1,NN/2}]],Im[Table[Evaluate[(tmp[[s]]E^(-I/\[HBar] *\[CapitalOmega] t)/.reseni)/.t->i/(NN/2) *TIME],{i,1,NN/2}]]];
Print["@@@Os",s," ",Re[Table[Evaluate[(tmp[[s]]E^(-I/\[HBar] *\[CapitalOmega] t)/.resenisec)/.t->i/(NN/2) *TIME],{i,1,NN/2}]],Im[Table[Evaluate[(tmp[[s]]E^(-I/\[HBar] *\[CapitalOmega] t)/.resenisec)/.t->i/(NN/2) *TIME],{i,1,NN/2}]]];
];
tmp=To4x4FormalismVector[\[Rho][t],1,1];
For[s=1,s<=Length[tmp],s++,
Print["@@@Er",s," ",Re[Table[Evaluate[(tmp[[s]]/.reseni)/.t->i/(NN/2) *TIME],{i,1,NN/2}]],Im[Table[Evaluate[(tmp[[s]]/.reseni)/.t->i/(NN/2) *TIME],{i,1,NN/2}]]];
Print["@@@Es",s," ",Re[Table[Evaluate[(tmp[[s]]/.resenisec)/.t->i/(NN/2) *TIME],{i,1,NN/2}]],Im[Table[Evaluate[(tmp[[s]]/.resenisec)/.t->i/(NN/2) *TIME],{i,1,NN/2}]]];
];
tmp=To4x4FormalismVector[\[Rho][t],1,2];
For[s=1,s<=Length[tmp],s++,
Print["@@@2r",s," ",Re[Table[Evaluate[(tmp[[s]]E^(-I/\[HBar] *\[CapitalOmega] t)/.reseni)/.t->i/(NN/2) *TIME],{i,1,NN/2}]],Im[Table[Evaluate[(tmp[[s]]E^(-I/\[HBar] *\[CapitalOmega] t)/.reseni)/.t->i/(NN/2) *TIME],{i,1,NN/2}]]];
Print["@@@2s",s," ",Re[Table[Evaluate[(tmp[[s]]E^(-I/\[HBar] *\[CapitalOmega] t)/.resenisec)/.t->i/(NN/2) *TIME],{i,1,NN/2}]],Im[Table[Evaluate[(tmp[[s]]E^(-I/\[HBar] *\[CapitalOmega] t)/.resenisec)/.t->i/(NN/2) *TIME],{i,1,NN/2}]]];
];
];

Quit[];


(* ::Input:: *)
(*tmp=Table[Table[If[!(i==1)||!(j>=2&&j<=MoleculeNumber+1),0,1],{i,1,Length[Hs]}],{j,1,Length[Hs]}];*)
(*mezcor=Max[Flatten[Abs[(\[Rho][t]/.reseni/.t->0)tmp]]];*)
(*(*Print[Plot[Evaluate[Flatten[Abs[tmp * \[Rho][t]]/.reseni]],{t,0,TIME},PlotRange->{0,mezcor},AxesLabel->{"t [fs]","p"}]];*)
(*Print[Plot[Evaluate[Flatten[Abs[tmp * \[Rho][t]]/.resenisec]],{t,0,TIME},PlotRange->{0,mezcor},AxesLabel->{"t [fs]","p"}]];*)
(*Print[Plot[Evaluate[Flatten[Re[tmp * \[Rho][t] E^(I/\[HBar] *\[CapitalOmega] t)]/.reseni]],{t,0,TIME},PlotRange->{-mezcor,mezcor},AxesLabel->{"t [fs]","p"}]];*)
(*Print[Plot[Evaluate[Flatten[Re[tmp * \[Rho][t]E^(I/\[HBar] *\[CapitalOmega] t)]/.resenisec]],{t,0,TIME},PlotRange->{-mezcor,mezcor},AxesLabel->{"t [fs]","p"}]];*)
(*Print[Plot[Evaluate[Flatten[Re[tmp * \[Rho][t]]/.reseni]],{t,0,TIME},PlotRange->{-mezcor,mezcor},AxesLabel->{"t [fs]","p"}]];*)
(*Print[Plot[Evaluate[Flatten[Re[tmp * \[Rho][t]]/.resenisec]],{t,0,TIME},PlotRange->{-mezcor,mezcor},AxesLabel->{"t [fs]","p"}]];*)*)
(**)
(*MatrixForm[tmp]*)
(*Print[Plot[Evaluate[Flatten[Abs[tmp * \[Rho][t]]/.reseni]],{t,0,1000},PlotRange->All,AxesLabel->{"t [fs]","p"}]];*)
(*Print[Plot[Evaluate[Flatten[Abs[tmp * \[Rho][t]]/.resenisec]],{t,0,1000},PlotRange->All,AxesLabel->{"t [fs]","p"}]];*)
(**)
(*tmp=Table[Table[If[!(i>=2&&i<=MoleculeNumber+1)||!(j>=2&&j<=MoleculeNumber+1)||i==j,0,1],{i,1,Length[Hs]}],{j,1,Length[Hs]}];*)
(*Print[MatrixForm[tmp]]*)
(*Print[Plot[Evaluate[Flatten[Re[tmp * \[Rho][t] ]/.reseni]],{t,0,1000},PlotRange->All,AxesLabel->{"t [fs]","p"}]];*)
(*Print[Plot[Evaluate[Flatten[Re[tmp * \[Rho][t]]/.resenisec]],{t,0,1000},PlotRange->All,AxesLabel->{"t [fs]","p"}]];*)
(**)
(*tmp=Table[Table[If[!(i>=2&&i<=MoleculeNumber+1)||!(j>=2&&j<=MoleculeNumber+1)||!(i==j),0,1],{i,1,Length[Hs]}],{j,1,Length[Hs]}];*)
(*MatrixForm[tmp]*)
(*Print[Plot[Evaluate[Flatten[Abs[tmp * \[Rho][t]]/.reseni]],{t,0,TIME},PlotRange->All,AxesLabel->{"t [fs]","p"}]];*)
(*Print[Plot[Evaluate[Flatten[Abs[tmp * \[Rho][t]]/.resenisec]],{t,0,TIME},PlotRange->All,AxesLabel->{"t [fs]","p"}]];*)
(**)
(*tmp=Table[Table[If[!(j==1)||!(i>=2&&i<=MoleculeNumber+1)||i==j,0,1],{i,1,Length[Hs]}],{j,1,Length[Hs]}];*)
(*Print[MatrixForm[tmp]]*)
(*Print[Plot[Evaluate[Flatten[Re[tmp * \[Rho][t]E^(-I \[CapitalOmega] t/\[HBar]) ]/.reseni]],{t,0,1000},PlotRange->All,AxesLabel->{"t [fs]","p"}]];*)
(*Print[Plot[Evaluate[Flatten[Re[tmp * \[Rho][t]E^(-I \[CapitalOmega] t/\[HBar])]/.resenisec]],{t,0,1000},PlotRange->All,AxesLabel->{"t [fs]","p"}]];*)
(**)
(*tmp=Table[Table[If[!(i>MoleculeNumber+1)||!(j>=2&&j<=MoleculeNumber+1),0,1],{i,1,Length[Hs]}],{j,1,Length[Hs]}];*)
(*Print[MatrixForm[tmp]]*)
(*Print[Plot[Evaluate[Flatten[Re[tmp * \[Rho][t]E^(-I \[CapitalOmega] t/\[HBar])]/.reseni]],{t,0,1000},PlotRange->All,AxesLabel->{"t [fs]","p"}]];*)
(*Print[Plot[Evaluate[Flatten[Re[tmp * \[Rho][t]E^(-I \[CapitalOmega] t/\[HBar])]/.resenisec]],{t,0,1000},PlotRange->All,AxesLabel->{"t [fs]","p"}]];*)
(**)
