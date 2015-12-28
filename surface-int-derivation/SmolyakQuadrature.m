(* ::Package:: *)

(* Some common definitions *)
ClearAll[iteratePrec,roundPrec,defaultPrecMult]
defaultPrecMult:=3/2
(* Rounds all real numbers in v to the nearest multiple of 10^(-prec+Ceiling[Log[10,v]]) *)
roundPrec[v_,prec_Integer]:=v/.x_Real:>Round[x,10^(-prec+Ceiling[Log[10,10.^-100+Abs[N[x]]]])];
(* Successively computes func[prec] with increasing precision until the result stops changing *)
iteratePrec[func_,precStart_Integer:20,precMult_:defaultPrecMult]:=Module[{prec=Ceiling[precStart precMult],prevRounded,curRounded,cur,round},
prevRounded=roundPrec[func[prec],precStart];
While[True,
prec=Ceiling[prec precMult];
cur=func[prec];
curRounded=roundPrec[cur,precStart];
If[curRounded===prevRounded,Break[]];
prevRounded=curRounded
];
N[cur,precStart]
]


(* Legendre polynomials *)
ClearAll[P]
P[n_Integer]:=P[n]=Expand[LegendreP[n,x]]
legendreRoots[n_]:=Table[Root[P[n]/.x->#,i],{i,n}]


(* Weights for Legendre-Gaussian quadrature *)
ClearAll[gw,x,f,calcGw,gint]
calcGw[n_Integer,prec_Integer]:=Module[{p,dp,xs,ws},
p=P[n];
dp=D[p,x];
xs=Table[N[Root[Evaluate[p/.x->#1]&,k],prec],{k,n}];
ws=N[1/((1-#^2)(dp/.x->#)^2)&/@xs,prec];
Sum[ws[[i]] f[N[(xs[[i]]+1)/2,prec]],{i,n}]
]
gw[n_Integer,prec_Integer]:=gw[n,prec]=iteratePrec[calcGw[n,#]&,prec]
gint[func_,n_,prec_Integer]:=gw[n,prec]/.f->func


(* checkPrec: checks that a quadrature rule Q[n_,prec_] is accurate up to a certain degree *)
ClearAll[accurateToOrder,accuracyOrder,checkAccuracyOrder]
(* accurateToOrder[q_,prec_,k_,precAllowance_:2]:=Abs[roundPrec[q/.f->If[k>0,#^k&,1&],prec]-1/(k+1)]<10^(-prec+precAllowance) *)
accurateToOrder[q_,prec_,k_,precAllowance_:2]:=Module[{
func=HornerForm[LegendreP[k,2#-1]],
val=If[k==0,1,0]
},
Abs[roundPrec[roundPrec[q,prec]/.f->(Evaluate[func]&),prec]-val]<10^(-prec+precAllowance)]
accuracyOrder[Q_,n_,prec_,precAllowance_:2]:=Module[{k=0},
While[accurateToOrder[Q[n,prec],prec,k,precAllowance],k++];
k-1]
checkAccuracyOrder[Q_,maxn_,prec_,orderF_]:=Do[Assert[accuracyOrder[Q,n,prec]==orderF[n]],{n,1,maxn}]


(* Smolyak quadrature *)
ClearAll[indices]
(* indices[d,k] is the list of all d-tuples of integers with the sum of them <= d+k-1 *)
indices[d_Integer,k_Integer]/;d>0:=indices[d,k]=
Join@@Table[Append[#,m]&/@indices[d-1,k+1-m],{m,k}]
indices[0,k_Integer]:={{}}


(* The Smolyak rule is defined for a dimension d and a 1d quadrature rule Q
such that Subscript[Q, i] is exact up to polynomials of order 2i-1 *)
ClearAll[smolyak,Qcoeff,f,g,multTerms,smolyakInt,smolyakW]
(* Qcoeff: extracts coefficients from a term in an integration rule *)
Qcoeff[w_Real f[x_Real]]:={w,x}
Qcoeff[x_]:=(Message[Qcoeff::eval,x];Abort[])
Qcoeff::eval:="Argument `1` not of the form w_ f[x_]"
(* multTerms: given a number of terms of the form Subscript[w, i]f[Subscript[x, i]], forms a product term Subscript[w, 1]... Subscript[w, d]f[Subscript[x, 1],...,Subscript[x, d]] *)
multTerms[xs__]:=With[{coeff=Transpose[Qcoeff/@{xs}]},Times@@coeff[[1]] f@@coeff[[2]]]
(* smolyak: the smolyak rule *)
smolyak[Q_,d_,k_,prec_]:=smolyak[Q,d,k,prec]=Module[{\[CapitalDelta]},
(* First, compute the differential rules Subscript[\[CapitalDelta], k]=Subscript[Q, Subscript[m, k]]-Subscript[Q, Subscript[m, k-1]] *)
\[CapitalDelta][i_]:=\[CapitalDelta][i]=Q[i,prec]-If[i>1,Q[i-1,prec],0];
Total[Map[
Distribute[g@@(\[CapitalDelta]/@#)]/.g->multTerms&,
indices[d,k]]]
]
smolyakW[Q_,d_,k_,prec_:MachinePrecision]:=smolyakW[Q,d,k]=N[iteratePrec[smolyak[Q,d,k,#]&,prec+1],prec]
smolyakInt[func_Function,Q_,d_,k_]:=smolyakW[Q,d,k]/.f->func
smolyakArray[Q__]:=List@@smolyakW[Q]/.a_ f[x__]->{a,x}


(* Kronrod-Patterson rules for functions defined on (-1, 1) *)
ClearAll[kpExtend,kpAbscissas,kp]
kpExtend[xs_,prec_]:=Module[{n,p,q,A,b,cs,G},
n=Length[xs];
p=n+1;
q=n-2Floor[n/2];
A=Table[P[2i-2+p+q]/.x->xs[[j]],{j,Floor[n/2]},{i,Floor[n/2]}];
(*Print[n," ",Precision[xs]," ",LinearAlgebra`MatrixConditionNumber[N[A]]];  *)
b=Table[-P[n+p]/.x->xs[[j]],{j,Floor[n/2]}];
cs=Append[LinearSolve[A,b],1];
G=SetPrecision[cs.Table[P[2i-2+p+q],{i,Floor[n/2]+1}],prec]; (*/.x->#,*)
(*Table[Root[G,k],{k,n+p}] *)
Cases[NRoots[G==0,x,PrecisionGoal->prec],x==v_->v,\[Infinity]]
]
kpAbscissas[1,prec_]:={0}
kpAbscissas[2,prec_]:=N[Table[Root[P[3]/.x->#,k],{k,3}],prec]
kpAbscissas[n_,prec_]/;n>2&&n<7:=kpAbscissas[n,prec]=iteratePrec[kpExtend[kpAbscissas[n-1,#],#]&,prec]
kpAbscissas[n_,prec_]/;n>=7:=Abort[]


ClearAll[weights,quadrature]
(* weights: given a set of abscissas on (0, 1), calculate the corresponding integration weights *)
weights[xxs_,prec_,range_:{0,1}]:=Module[{n,xs,ws,Fs,gauss},
xs=N[xxs,prec];
n=Length[xs];
Fs[x_,i_]:=Product[x-xs[[j]],{j,i-1}]Product[x-xs[[j]],{j,i+1,n}];
gauss=gw[Ceiling[n/2],prec];
Table[With[{fs=Fs[#,i]/Fs[xs[[i]],i]},gauss/.f->(fs&)],{i,n}]
]
(* quadrature: given a set of abscissas on (0, 1), construct a quadrature formula *)
quadrature[xs_,prec_]:=Total[#[[2]] f[#[[1]]]&/@({N[xs,prec],weights[xs,prec]}\[Transpose])]


(* kp: Kronrod-Patterson quadrature on (0, 1) *)
kp[n_,prec_]:=kp[n,prec]=iteratePrec[quadrature[(kpAbscissas[n,prec]+1)/2,#]&,prec]


(* Delayed Kronrod-Patterson *)
ClearAll[kpDelay]
kpDelay/:kpDelay[ds_][n_,prec_]:=kp[ds[[n]],prec]
delayFull:={1,2,2,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,6}
delayLittle:={1,2,2,3,4,5,5,6,6,6,6}


ClearAll[ecfRule, ecfSymmetric,ecfSymmetricSimplex,a,b,c,d,changeSign,changeSigns,ecfRotate,ecfRotateSigns,padCoordinates,ecfFilterResult]
padCoordinates[ecfSymmetric][{xs_,w_}]:={xs,w}
ecfFilterResult[ecfSymmetric][x_]:=x
ecfSymmetric[4]={
{0,0,0,0},
{a,0,0,0},
{a,a,0,0},
{a,b,0,0},

{a,a,a,0},
{a,a,b,0},
{a,b,c,0},
{a,a,a,a},

{a,a,a,b},
{a,a,b,b},
{a,a,b,c},
{a,b,c,d}
};
ecfFilterResult[ecfSymmetricSimplex][x_]:=Select[x,Min[#[[2;;]]]>=0&]
padCoordinates[ecfSymmetricSimplex][{xs_,w_}]:={Append[xs,1-Total[xs]],w}
ecfSymmetricSimplex[2]={
{1,0,0},
{1/2,1/2,0},
{a,b,0},
{1/3,1/3,1/3},
{a,a,b},
{a,b,c}
};
ecfSymmetricSimplex[3]={
{1,0,0,0},
{1/2,1/2,0,0},
{a,b,0,0},
{1/3,1/3,1/3,0},

{a,a,b,0},
{a,b,c,0},
{1/4,1/4,1/4,1/4},
{a,a,a,b},

{a,a,b,b},
{a,a,b,c},
{a,b,c,d}
};
changeSign[s_,w_]:=MapThread[If[#1==1,-#2,#2]&,{s,w}]
changeSigns[w_]:=Union[Table[changeSign[IntegerDigits[i,2,Length[w]],w],{i,2^Length[w]}]]
ecfRotate[p_,{xs_,w_}]:=Map[Prepend[#,w]&,Permutations[p]/.{x:(a|b|c|d):>xs[[Position[p,x][[1,1]]]]}]
ecfRotateSigns[p_,ws_]:=Union[Join@@Map[ecfRotate[#,ws]&,changeSigns[p]]]
ecfRule[{deg_,order_,points_,type_,structure_,weights__}]:=With[{
ruleStructure=Join@@MapThread[Table[#1,{i,#2}]&,{type[deg],structure}]
},
N[ecfFilterResult[type][Join@@MapThread[ecfRotateSigns,{ruleStructure,padCoordinates[type]/@{weights}}]],32]
]


(* Rules from ECF *)
ClearAll[rule4d9k145,simplex3d8k43]
rule4d9k145={
4, (* degree *)
9, (* order of accuracy *)
145, (* number of points *)
ecfSymmetric,(* rule type *)
{1,1,1,0,1,0,0,1,1,0,0,0}, (* rule structure *)
{{0.,0.,0.,0.},-1.90110453939705164190484531991390},
{{0.572452539771813136387507979953173,0.,0.,0.},0.799177722185472093345436906868952},
{{0.997379250347217930412807106259392,0.997379250347217930412807106259392,0.,0.},0.0322754507138618007859731508509594},
{{0.794567852198484877426583579329962,0.794567852198484877426583579329962,0.794567852198484877426583579329962,0.},0.0991704634972827999232052020858457},
{{0.918048087601965086301798204660102,0.918048087601965086301798204660102,0.918048087601965086301798204660102,0.918048087601965086301798204660102},0.0168278157115170721061830811272882},
{{0.465887724933933437289233730125708,0.465887724933933437289233730125708,0.465887724933933437289233730125708,0.909586231538307152694330588148460},0.113912063460676076953695291871180}
};
(* Simplex rules from ECF *)
simplex3d8k43={
3, (* degree *)
8, (* order of accuracy *)
43, (* number of points *)
ecfSymmetricSimplex,(* rule type *)
{0,0,0,0,0,0,1,3,1,2,0}, (* rule structure *)
(* points *)
{{0.25,0.25,0.25},-0.0205001886586399158405865177642941},
{{0.206829931610673204083980900024961,0.206829931610673204083980900024961,0.206829931610673204083980900024961},0.0142503058228669012484397415358704},
{{0.0821035883105467230906058078714215,0.0821035883105467230906058078714215,0.0821035883105467230906058078714215},10^-3 1.96703331313390098756280342445466},
{10^-3 {5.78195050519799725317663886414270,5.78195050519799725317663886414270,5.78195050519799725317663886414270},10^-4 1.69834109092887379837744566704016},
{{0.0505327400188942244256245285579071,0.0505327400188942244256245285579071,0.449467259981105775574375471442092},10^-3 4.57968382446728180074351446297276},
{{0.229066536116811139600408854554753,0.229066536116811139600408854554753,0.0356395827885340437169173969506114},10^-3 5.70448580868191850680255862783040},
{{0.0366077495531974236787738546327104,0.0366077495531974236787738546327104,0.190486041934633455699433285315099},10^-3 2.14051914116209259648335300092023}
};
simplex2d10k25a={
2, (* degree *)
10, (* order of accuracy *)
25, (* number of points *)
ecfSymmetricSimplex,(* rule type *)
{0,0,0,1,2,3}, (* rule structure *)
(* points *)
{{0.333333333333333333333333333333333, 0.333333333333333333333333333333333},0.0399472523706198539156235226066932},
{{0.425086210602090572969529511638044, 0.425086210602090572969529511638044},0.0355619011161886673196456436993290},
{{0.0233088675100001907144663868959796, 0.0233088675100001907144663868959796},10^-3 4.11190934523209775932331018123594},
{{0.628307400213492556420837666078834, 0.223766973576973006225686490268204},0.0227152961480850090035368146219665},
{{0.611313826181397648918755002253901, 0.358740141864431464578155300723852},0.0186799281171526384131182495009876},
{{0.821072069985629373373544413472177, 0.143295370426867145305856630617323},0.0154433284422819943912565385023143}
};
simplex2d3k4={
2, (* degree *)
3, (* order of accuracy *)
4, (* number of points *)
ecfSymmetricSimplex,(* rule type *)
{0,0,0,1,1,0}, (* rule structure *)
(* points *)
{{0.333333333333333333333333333333333, 0.333333333333333333333333333333333},-0.28125},
{{0.2, 0.2},0.260416666666666666666666666666666}
};



to01cube4d[ws_]:=Map[Prepend[N[(1+#[[2;;5]])/2,32],#[[1]]/16]&,ws]


(* Tetrahedral integration *)
ClearAll[f,rs,intTet,intTri,polyVolume,polySurface,intTetWeights,intTriWeights]
ecfTetRule=N[ecfRule[simplex3d8k43]];
ecfTriRule=N[ecfRule[simplex2d10k25a]];
polyVolume[poly_]:=-Det[#-poly[[1]]&/@poly[[2;;]]]
polySurface[tri_]:=Norm[Cross[tri[[1]]-tri[[2]],tri[[1]]-tri[[3]]]]
intTet[f_,rs_List]:=Sum[ecfTetRule[[i,1]]f@@(ecfTetRule[[i,2;;]].rs),{i,Length[ecfTetRule]}]polyVolume[rs]
intTri[f_,rs_List]:=Sum[ecfTriRule[[i,1]]f@@(ecfTriRule[[i,2;;]].rs),{i,Length[ecfTriRule]}]polySurface[rs]
intTetWeights[rs_List]:=With[{v=polyVolume[rs]}, {v #[[1]], #[[2;;]].rs}&/@ecfTetRule]
intTriWeights[rs_List]:=With[{v=polySurface[rs]}, {v #[[1]], #[[2;;]].rs}&/@ecfTriRule]


(* Weights for line integration *)
gaussianLineRule[n_]:=List@@gw[n,20]/.a_ f[x__]:>{a,x,1-x}
