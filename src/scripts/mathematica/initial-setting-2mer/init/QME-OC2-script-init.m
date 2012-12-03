(* ::Package:: *)

MoleculeNumber=2;

MezCC=8; (* urcuje do jakeho clenu se scita Matsubarova suma, nekonecno je legalni volba, ale trva hodne dlouho *)
SECULAR=False; (*=True; - prepina QME sekularni aproximaci*)

(* nasledujici hodnoty je potreba zadavat s desetinnou teckou, jinak se Mathematica pokousi o presny vypocet, coz trva dlouho, casem lze odstranit*)
T[s_]={300,300,300}[[s]];
\[Tau]c[s_]={50,50,50}[[s]];
\[Lambda]\[Lambda][s_]={30,30,30}[[s]];
j12=-233.0+0.01;
j13=-233.0*0+0.01;
j23=-233.0*0+0.01;
J[r_,s_] :={{0.0, j12 , j13},{j12, 0.0, j23},{j13,j23,0.0}}[[r]][[s]];
J2[r_,s_]:={{0.0, j23 , j13},{j23, 0.0, j12},{j13,j12,0.0}}[[r]][[s]];
\[CapitalOmega]=10000;				(*Frekvence, ktera se odecita*)
\[Epsilon][i_]={-500.0+\[CapitalOmega],0.0+\[CapitalOmega],500.0+\[CapitalOmega]}[[i]];


NN=1024*2*8;
dt=1/4; (* dt*NN/2 = celkovy cas *)

(* \[Rho]0={1,0,0,0,0,0,0,0,0} *)
