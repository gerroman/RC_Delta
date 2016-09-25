(* ::Package:: *)

(* ::Section::Closed:: *)
(*Particles*)


<<Bremsstrahlung`Traces`Particles`


{me, Mp, Md, Gd}//TraditionalForm


(* ::Section::Closed:: *)
(*FormFactors*)


<<Bremsstrahlung`Traces`FormFactors`


{F1,F2,GE,GM,G1,G2,G3,GEstar,GMstar,GCstar}//TraditionalForm


(* Test conversion rules F1, F2 \[LeftRightArrow] GE, GM *)
Simplify[
	And[
		(F1[q2]/.ruleFtoG/.ruleGtoF)==F1[q2],
		(F2[q2]/.ruleFtoG/.ruleGtoF)==F2[q2]
	]
]
Simplify[
	And[
		(GM[q2]/.ruleGtoF/.ruleFtoG)==GM[q2],
		(GE[q2]/.ruleGtoF/.ruleFtoG)==GE[q2]
	]
]


(* Test conversion rules G1, G2, G3 \[LeftRightArrow] GEstar, GMstar, GCstar *)
Simplify[
	And[
		(G1[q2]/.ruleG123toGEMC/.ruleGEMCtoG123) == G1[q2],
		(G2[q2]/.ruleG123toGEMC/.ruleGEMCtoG123) == G2[q2],
		(G3[q2]/.ruleG123toGEMC/.ruleGEMCtoG123) == G3[q2]
	]
]
Simplify[
	And[
		(GMstar[q2]/.ruleGEMCtoG123/.ruleG123toGEMC) == GMstar[q2],
		(GEstar[q2]/.ruleGEMCtoG123/.ruleG123toGEMC) == GEstar[q2],
		(GCstar[q2]/.ruleGEMCtoG123/.ruleG123toGEMC) == GCstar[q2]
	]
]


(* ::Section::Closed:: *)
(*Vertices*)


<<Bremsstrahlung`Traces`Vertices`


(* ::Subsection::Closed:: *)
(*\[Gamma]p -> p*)


With[
	{
		v1 = Hold[gamma["gp->p"][\[Mu]][q]],
		v2 = Hold[gamma["p->gp"][\[Mu]][q]],
		v3 = Hold[diracConjugate[gamma["gp->p"][\[Mu]][q]]],
		v4 = Hold[diracConjugate[gamma["p->gp"][\[Mu]][q]]]
	},
	{
		Plus@@(ReleaseHold[v1 * commonFactors[v1]]),
		Plus@@(ReleaseHold[v2 * commonFactors[v2]]),
		Plus@@(ReleaseHold[v3 * commonFactors[v3]]),
		Plus@@(ReleaseHold[v4 * commonFactors[v4]]),
		Resolve[Plus@@(ReleaseHold[v1 * commonFactors[v1]]) == Plus@@(ReleaseHold[v4 * commonFactors[v4]])],
		Resolve[Plus@@(ReleaseHold[v2 * commonFactors[v2]]) == Plus@@(ReleaseHold[v3 * commonFactors[v3]])]
	}
]//TableForm


(* ::Subsection::Closed:: *)
(*\[Gamma]p -> \[CapitalDelta]*)


With[
	{
		v1 = Hold[gamma["gp->d"][\[Nu], \[Beta]][Subscript[p, \[CapitalDelta]], q]],
		v2 = Hold[gamma["d->gp"][\[Nu], \[Beta]][Subscript[p, \[CapitalDelta]], q]],
		v3 = Hold[diracConjugate[gamma["gp->d"][\[Nu], \[Beta]][Subscript[p, \[CapitalDelta]], q]]],
		v4 = Hold[diracConjugate[gamma["d->gp"][\[Nu], \[Beta]][Subscript[p, \[CapitalDelta]], q]]]
	},
	{
		Plus@@(ReleaseHold[v1 * commonFactors[v1]]),
		Plus@@(ReleaseHold[v2 * commonFactors[v2]]),
		Plus@@(ReleaseHold[v3 * commonFactors[v3]]),
		Plus@@(ReleaseHold[v4 * commonFactors[v4]]),
		Resolve[Plus@@(ReleaseHold[v1 * commonFactors[v1]]) == Plus@@(ReleaseHold[v4 * commonFactors[v4]])],
		Resolve[Plus@@(ReleaseHold[v2 * commonFactors[v2]]) == Plus@@(ReleaseHold[v3 * commonFactors[v3]])]
	}
	]//TableForm


(* ::Subsection::Closed:: *)
(*Gauge invariance and vertex transversity of \[Gamma]p->\[CapitalDelta] interaction vertex*)


(* Gauge invariance *)
With[
	{
		v1 = Hold[gamma["gp->d"][\[Nu], \[Beta]][Subscript[p, \[CapitalDelta]], q]],
		v2 = Hold[gamma["d->gp"][\[Nu], \[Beta]][Subscript[p, \[CapitalDelta]], q]],
		v3 = Hold[diracConjugate[gamma["gp->d"][\[Nu], \[Beta]][Subscript[p, \[CapitalDelta]], q]]],
		v4 = Hold[diracConjugate[gamma["d->gp"][\[Nu], \[Beta]][Subscript[p, \[CapitalDelta]], q]]]
	},
	{
		DiracReduce@Contract[ReleaseHold[v1], FV[q, \[Nu]]],
		DiracReduce@Contract[ReleaseHold[v2], FV[q, \[Nu]]],
		DiracReduce@Contract[ReleaseHold[v3], FV[q, \[Nu]]],
		DiracReduce@Contract[ReleaseHold[v4], FV[q, \[Nu]]]
	}
]//TableForm


(* Delta (1232) momentum transversity *)
With[
	{
		v1 = Hold[gamma["gp->d"][\[Nu], \[Beta]][Subscript[p, \[CapitalDelta]], q]],
		v2 = Hold[gamma["d->gp"][\[Nu], \[Beta]][Subscript[p, \[CapitalDelta]], q]],
		v3 = Hold[diracConjugate[gamma["gp->d"][\[Nu], \[Beta]][Subscript[p, \[CapitalDelta]], q]]],
		v4 = Hold[diracConjugate[gamma["d->gp"][\[Nu], \[Beta]][Subscript[p, \[CapitalDelta]], q]]]
	},
	{
		DiracReduce@Contract[ReleaseHold[v1], FV[Subscript[p, \[CapitalDelta]], \[Beta]]],
		DiracReduce@Contract[ReleaseHold[v2], FV[Subscript[p, \[CapitalDelta]], \[Beta]]],
		DiracReduce@Contract[ReleaseHold[v3], FV[Subscript[p, \[CapitalDelta]], \[Beta]]],
		DiracReduce@Contract[ReleaseHold[v4], FV[Subscript[p, \[CapitalDelta]], \[Beta]]]
	}
]//TableForm


(* ::Section::Closed:: *)
(*CurrentTensors for ep->ep and ep->e\[CapitalDelta]*)


<<Bremsstrahlung`Traces`CurrentTensors`


(* ::Subsection::Closed:: *)
(*Electron current tensors*)


ect=electronCurrentTensorElastic[\[Mu], \[Nu]]
ectd=electronCurrentTensorDelta[\[Mu], \[Nu]]


(* gauge invariance *)
{
	Contract[ect FV[qel, \[Mu]]]//FCE,
	Contract[ect FV[qel, \[Nu]]]//FCE,
	Contract[ectd FV[qd, \[Mu]]]//FCE,
	Contract[ectd FV[qd, \[Nu]]]//FCE
}


(* ::Subsection::Closed:: *)
(*Proton current tensor*)


Timing[(* ~0.5s *)
	pct=Collect[
		Plus@@(commonFactors[protonCurrentTensorElastic[\[Mu], \[Nu]]] *
			protonCurrentTensorElastic[\[Mu], \[Nu]])/.ruleFtoG,
	{GE[_],GM[_],MT[__],FV[a__]FV[b__]},
	Factor
]]


(* gauge invariance *)
{
	Contract[pct,FV[qel, \[Mu]]]//FCE,
	Contract[pct,FV[qel, \[Nu]]]//FCE
}


(* ::Subsection::Closed:: *)
(*Proton transition tensor*)


Timing[(* ~10s *)
	pdct=Collect[
		Expand[Plus@@(commonFactors[protonTransitionTensorDelta[\[Mu], \[Nu]]] *
			protonTransitionTensorDelta[\[Mu], \[Nu]])/.ruleG123toGEMC/.ruleTildePd
		],
		{GEstar[_], GMstar[_], GCstar[_], MT[__], FV[a__]FV[b__]},
		Factor
	];
]


(* gauge invariance *)
{
	Together@FCE[Contract[pdct*FV[qd,\[Mu]]]],
	Together@FCE[Contract[pdct*FV[qd,\[Nu]]]]
}


pdctSimplified = pdct/.{
	MT[a__]x_+b_+c_:>-x(-MT[a]-b/x-c/x)
}/.{
	2Mp Md + Mp^2+Md^2-SP[qd,qd]->
	Cancel[(2Mp Md + Mp^2+Md^2-SP[qd,qd])/SP[tildePd,tildePd]]*SP[HoldForm[tildePd],HoldForm[tildePd]]
}/.{
	-2Mp Md + Mp^2+Md^2-SP[qd,qd]->Simplify[-2Mp Md + Mp^2+Md^2-SP[qd,qd]]
}


pdctSimplifiedAndFactored = With[
	{mult = (Mp+Md)^2/(4Mp^2) ((Md-Mp)^2-SP[qd,qd])},
	mult (
		(#/mult)& /@ pdctSimplified
	)/.{GMstar[a_]^2 x_ + 3 GEstar[a_]^2 x_ :> x (GMstar[a]^2 + 3 GEstar[a]^2)}
]


(* ::Subsection::Closed:: *)
(*Rosenbluth cross section*)


\[Epsilon]==(SP[Kel,Pel]^2+(4 Mp^2-SP[qel,qel]) SP[qel,qel])/(SP[Kel,Pel]^2-(4 Mp^2-SP[qel,qel]) (4 me^2+SP[qel,qel]))
\[Tau]==-SP[qel,qel]/(4Mp^2)


tensorContract=Collect[FCE@Contract[ect*pct],{GE[_],GM[_]},Simplify]
tensorContractSimplified=With[
{mult = -(4Mp^2)/SP[qel,qel]*Coefficient[tensorContract,GM[SP[qel,qel]]^2]},
mult*( 
(Simplify[#/mult])&/@tensorContract
)
]//.{
-SP[qel,qel]/(4Mp^2)->\[Tau],
(SP[Kel,Pel]^2-(4 Mp^2-SP[qel,qel]) (4 me^2+SP[qel,qel]))->\[Epsilon]^-1*(SP[Kel,Pel]^2+(4 Mp^2-SP[qel,qel]) SP[qel,qel]),
-4Mp^2/(SP[qel,qel]-4Mp^2)->1/(1+\[Tau])
}



tensorContractHighElectronEnergy=((tensorContractSimplified//.{
SP[Kel,Pel]:>s-u,
SP[qel,qel]->t,
u->2me^2+2Mp^2-s-t
}/.{me->0}/.{
s->Mp^2+2Mp \[CurlyEpsilon]1,
t->-4\[CurlyEpsilon]1 \[CurlyEpsilon]3 Sin[\[Theta]/2]^2
}/.{\[CurlyEpsilon]3->\[CurlyEpsilon]1/(1+(2\[CurlyEpsilon]1)/Mp Sin[\[Theta]/2]^2)}//Simplify)/.{\[CurlyEpsilon]1^2->\[CurlyEpsilon]3 \[CurlyEpsilon]1 (1+(2\[CurlyEpsilon]1)/Mp Sin[\[Theta]/2]^2)}//Simplify)/.{-4 \[CurlyEpsilon]1 \[CurlyEpsilon]3 Sin[\[Theta]/2]^2->SP[qel,qel]}



rutherfordCrossSection=1/(4\[Pi])^2 1/(4Mp^2 \[Eta]) \[CurlyEpsilon]3/\[CurlyEpsilon]1 (Z^2 e^4)/SP[qel,qel]^2 tensorContractHighElectronEnergy/.{
e^4->(4\[Pi] \[Alpha])^2,
SP[qel,qel]^-2->(16\[CurlyEpsilon]1^2 \[CurlyEpsilon]3^2 HoldForm[Sin[\[Theta]/2]^4])^-1
}


(* ::Subsection::Closed:: *)
(*\[CapitalDelta] excitation tensor contract*)


tensorContractDelta=Collect[
	FCE@Contract[ectd*pdct],
	{GEstar[_],GMstar[_],GCstar[_]},
	FullSimplify
];
Block[{
	\[Tau] = -SP[qd, qd]/(Md+Mp)^2,
	\[Epsilon] = (SP[Kd,tildePd]^2+SP[qd,qd]SP[tildePd, tildePd])/(SP[Kd,tildePd]^2-(SP[Kd,Kd]+2SP[qd,qd])SP[tildePd,tildePd]),
	tensorContractSimplified = (SP[Kd,tildePd]^2+SP[qd,qd]SP[HoldForm[tildePd], HoldForm[tildePd]])*(Md+Mp)^2/(4Mp^2)*
			(HoldForm[\[Tau]] (GMstar[SP[qd, qd]]^2+3GEstar[SP[qd, qd]]^2+HoldForm[\[Epsilon]] (-SP[qd,qd])/Md^2 GCstar[SP[qd, qd]]^2))/(HoldForm[\[Epsilon]] (1+HoldForm[\[Tau]]))
	},
	Print[tensorContractSimplified];
	Simplify[
		tensorContractDelta == ReleaseHold[tensorContractSimplified]
	]
]


(* ::Section::Closed:: *)
(*Bremsstrahlung |iMe|^2*)


<<Bremsstrahlung`Traces`BremsstrahlungElectronSquared`


(* ::Subsection::Closed:: *)
(*electron line trace*)


Timing[(* ~0.4s *)tensorBrE[1, 1][\[Nu]1, \[Nu]2];]
Timing[(* ~0.4s *)tensorBrE[1, 2][\[Nu]1, \[Nu]2];]
Timing[(* ~0.4s *)tensorBrE[2, 1][\[Nu]1, \[Nu]2];]
Timing[(* ~0.4s *)tensorBrE[2, 2][\[Nu]1, \[Nu]2];]


SP[r1,r1]-me^2//ExpandScalarProduct//FCE//Expand//StandardForm
SP[r2,r2]-me^2//ExpandScalarProduct//FCE//Expand//StandardForm


electronTrace["iMe2"] = Collect[
	Expand[
		1/(xK-xq)^2 * Plus@@(commonFactors[tensorBrE[1, 1][\[Nu]1, \[Nu]2]]*tensorBrE[1, 1][\[Nu]1, \[Nu]2])
		 + 1/(xK-xq) * 1/(-xK-xq) * Plus@@(commonFactors[tensorBrE[2, 1][\[Nu]1, \[Nu]2]]*tensorBrE[2, 1][\[Nu]1, \[Nu]2])
		 + 1/(xK-xq) * 1/(-xK-xq) * Plus@@(commonFactors[tensorBrE[1, 2][\[Nu]1, \[Nu]2]]*tensorBrE[1, 2][\[Nu]1, \[Nu]2])
		 + 1/(-xK-xq)^2 * Plus@@(commonFactors[tensorBrE[2, 2][\[Nu]1, \[Nu]2]]*tensorBrE[2, 2][\[Nu]1, \[Nu]2])
	],
	{MT[__],FV[a__]FV[b__]},
	Factor
];


(* gauge invariance *)
FCE@Contract[
	{
		electronTrace["iMe2"] FV[qp, \[Nu]1],
		electronTrace["iMe2"] FV[qp, \[Nu]2]
	}
] /. ruleQe2toQp2 /.{
	FV[qe,\[Nu]_]:>FV[qp, \[Nu]]+FV[k, \[Nu]]
} // Factor


(* ::Subsection::Closed:: *)
(*proton line trace*)


Timing[(* ~0.5s *)tensorTp[\[Nu]1, \[Nu]2];]


protonTrace["iMe2"] = Collect[
	Plus@@(commonFactors[tensorTp[\[Nu]1, \[Nu]2]] * tensorTp[\[Nu]1, \[Nu]2]),
	{MT[__], FV[a__]FV[b__]},
	Factor
];


(* gauge invariance *)
Contract[
	{
		protonTrace["iMe2"] FV[qp, \[Nu]1],
		protonTrace["iMe2"] FV[qp, \[Nu]2]
	}
]//FCE


(* ::Subsection::Closed:: *)
(*contraction and decomposition*)


result["iMe2"] = FCE@Contract[
	protonTrace["iMe2"] * electronTrace["iMe2"]
]/.ruleQe2toQp2//.ruleFtoG//Together;


resultNumeratorSimplified["iMe2"] = With[
	{numeratorList = List@@Expand[Numerator[result["iMe2"]]]},
	Simplify[Plus@@Cases[numeratorList,#]]& /@ Table[decompositionPattern[n],{n, 2, 10}]
] //. {
	(* additional simplification *)
	-qp2 (xK^2-xq^2) (xK^2-xP^2+xq^2)+4 Mp^2 (xK^4-xq^4)->(4Mp^2-qp2)(xK^4-xq^4) + qp2 (xK^2-xq^2) xP^2,
	-qp2 (xK^2-xq^2) (xK^2+xP^2+xq^2)+4 Mp^2 (xK^4-xq^4)->(4Mp^2-qp2)(xK^4-xq^4) - qp2 (xK^2-xq^2) xP^2,
	+4 Mp^2 qp2 (xK^2-xq^2)+qp2^2 (-xK^2+xq^2)->qp2 (4Mp^2-qp2)(xK^2-xq^2),
	-qp2 xK^2+qp2 xq^2+4 Mp^2 (xK^2-xq^2)->(4Mp^2-qp2)(xK^2-xq^2),
	(4 Mp^2-qp2) qp2 (xK^2-xq^2)+4 me^2 (-nu xK xP+(4 Mp^2-qp2) (xK^2-xq^2))->
	(4 Mp^2-qp2)(xK^2-xq^2)(4me^2+qp2)-4 me^2 nu xK xP,
	nu^2-4 Mp^2 qp2+qp2^2+4 me^2 (-4 Mp^2+qp2)->FullSimplify[nu^2-4 Mp^2 qp2+qp2^2+4 me^2 (-4 Mp^2+qp2)],
	-8 xq (4 Mp^2 (4 me^2 nu xK xP+(4 Mp^2-qp2) qp2 (xK^2-xq^2)) GE[qp2]^2+qp2 (-4 me^2 nu xK xP+(4 Mp^2-qp2) (4 me^2+qp2) (xK^2-xq^2)) GM[qp2]^2)->
	8 xq (-4 Mp^2 (4 me^2 nu xK xP+(4 Mp^2-qp2) qp2 (xK^2-xq^2)) GE[qp2]^2+qp2 (4 me^2 nu xK xP-(4 Mp^2-qp2) (4 me^2+qp2) (xK^2-xq^2)) GM[qp2]^2),
	-4 (4 Mp^2 (4 me^2 xK^2 xP^2+qp2 xP^2 (xK^2-xq^2)+(4 Mp^2-qp2) (xK^4-xq^4)) GE[qp2]^2+qp2 (-4 me^2 xK^2 xP^2-qp2 xP^2 (xK^2-xq^2)+(4 Mp^2-qp2) (xK^4-xq^4)) GM[qp2]^2)->
	4 (-4 Mp^2 (4 me^2 xK^2 xP^2+qp2 xP^2 (xK^2-xq^2)+(4 Mp^2-qp2) (xK^4-xq^4)) GE[qp2]^2+qp2 (4 me^2 xK^2 xP^2+qp2 xP^2 (xK^2-xq^2)-(4 Mp^2-qp2) (xK^4-xq^4)) GM[qp2]^2)
}


(* soft photon approximation *)
resultNumeratorSimplified["iMe2"][[1]]/ Denominator[result["iMe2"]]


resultSimplified["iMe2"] = Plus@@(resultNumeratorSimplified["iMe2"]) / Denominator[result["iMe2"]]
(resultSimplified["iMe2"] - result["iMe2"])//Simplify


(* ::Subsection::Closed:: *)
(*final result*)


iMe2 = (Z^2 e^6)/qp2^2 resultSimplified["iMe2"]


With[
	{file = OpenWrite["./src/iMe2.cpp"]},
	WriteString[file, "#include \"constants.h\"\n"];
	WriteString[file, "#include \"proton_formfactors.h\"\n"];
	WriteString[file, "#include <cmath>\n"];
	WriteString[file, "\n"];
	WriteString[file, "double iMe2(double nu, double, double qp2, double kK, double kP, double kq) {\n"];
	WriteString[file, "return ("];
	Write[file, CForm[iMe2]/.{
			me^2->me2, Mp^2->Mp2,
			e^6->e6, Z^2->Z2
		} /. {
			xK->kK,
			xP->kP,
			xq->kq,
			Power->pow
		}
	];
	WriteString[file, ");\n}"];
	Close[file];
];


(* ::Section::Closed:: *)
(*Bremsstrahlung |iMp|^2*)


<<Bremsstrahlung`Traces`BremsstrahlungProtonSquared`


(* ::Subsection::Closed:: *)
(*electron line trace*)


Timing[(* ~0.1s *)tensorTe[\[Nu]1, \[Nu]2];]


electronTrace["iMp2"] = Plus@@(commonFactors[tensorTe[\[Nu]1, \[Nu]2]] * tensorTe[\[Nu]1, \[Nu]2])


(* gauge invariance *)
FCE@Contract[
	{
		electronTrace["iMp2"] FV[qe, \[Nu]1],
		electronTrace["iMp2"] FV[qe, \[Nu]2]
	}
] /. ruleQe2toQp2 /.{
	FV[qp,\[Nu]_]:>FV[qe, \[Nu]]-FV[k, \[Nu]]
} // Factor


(* ::Subsection::Closed:: *)
(*proton line trace*)


Timing[(* ~60s *)tensorBrP[1,1][\[Nu]1,\[Nu]2];]
Timing[(* ~50s *)tensorBrP[2,1][\[Nu]1,\[Nu]2];]
Timing[(* ~40s *)tensorBrP[1,2][\[Nu]1,\[Nu]2];]
Timing[(* ~17s *)tensorBrP[2,2][\[Nu]1,\[Nu]2];]


SP[t1, t1] - Mp^2//ExpandScalarProduct//FCE//Expand
SP[t2, t2] - Mp^2//ExpandScalarProduct//FCE//Expand


Timing[(* ~1s *)
	protonTrace["iMp2"] = Collect[
	Expand[
		1/(xP+xq)^2 Plus@@(commonFactors[tensorBrP[1, 1][\[Nu]1, \[Nu]2]]*tensorBrP[1, 1][\[Nu]1, \[Nu]2])+
		1/((xP+xq)(xq-xP)) Plus@@(commonFactors[tensorBrP[2, 1][\[Nu]1, \[Nu]2]]*tensorBrP[2, 1][\[Nu]1, \[Nu]2])+
		1/((xP+xq)(xq-xP)) Plus@@(commonFactors[tensorBrP[1, 2][\[Nu]1, \[Nu]2]]*tensorBrP[1, 2][\[Nu]1, \[Nu]2])+
		1/(xq-xP)^2 Plus@@(commonFactors[tensorBrP[2, 2][\[Nu]1, \[Nu]2]]*tensorBrP[2, 2][\[Nu]1, \[Nu]2])
	],
	{MT[__],FV[a__]FV[b__]},
	Factor
];]


(* gauge invariance *)
FCE@Contract[
	{
		protonTrace["iMp2"] FV[qe, \[Nu]1],
		protonTrace["iMp2"] FV[qe, \[Nu]2]
	}
]/.ruleQp2toQe2/.{FV[qp,a_]:>FV[qe,a]-FV[k,a]}//Factor


(* ::Subsection::Closed:: *)
(*contraction and decomposition*)


Timing[(* ~1s *)
	result["iMp2"] = FCE@Contract[protonTrace["iMp2"] * electronTrace["iMp2"]]/.ruleQp2toQe2//Together;
]


Timing[(* ~2s *)
	resultNumeratorSimplified["iMp2"]= With[
		{numeratorList = List@@Expand[Numerator[result["iMp2"]]]},
		Simplify[Plus@@Cases[numeratorList,#]]& /@ Table[decompositionPattern[n],{n, 2, 10}]
	];
]


(*soft photon approximation term *)
(resultNumeratorSimplified["iMp2"][[1]]/ Denominator[result["iMp2"]])/.ruleFtoG//FullSimplify


resultSimplified["iMp2"] = Plus@@(resultNumeratorSimplified["iMp2"]) / Denominator[result["iMp2"]];
(resultSimplified["iMp2"] - result["iMp2"])//Simplify


(* ::Subsection::Closed:: *)
(*final result*)


iMp2 = (Z^4 e^6)/qe2^2 resultSimplified["iMp2"]


With[{file = OpenWrite["./src/iMp2.cpp"]},
	WriteString[file, "#include \"constants.h\"\n"];
	WriteString[file, "#include \"proton_formfactors.h\"\n"];
	WriteString[file, "#include <cmath>\n"];
	WriteString[file, "\n"];
	WriteString[file, "double iMp2(double nu, double qe2, double, double kK, double kP, double kq) {\n"];
	WriteString[file, "return ("];
	Write[file, CForm[iMp2]/.{
		me^2->me2, Mp^2->Mp2, Mp^4->Mp4,
		e^6->e6, Z^2->Z2, Z^4->Z4, Z^3->Z3
	} /. {
		xK->kK,
		xP->kP,
		xq->kq,
		Power->pow
	}];
	WriteString[file, ");\n}"];
	Close[file];
]


(* ::Section::Closed:: *)
(*Bremsstrahlung |iMd|^2*)


<<Bremsstrahlung`Traces`BremsstrahlungDeltaSquared`;


(* ::Subsection::Closed:: *)
(*electron line trace*)


electronTrace["iMd2"] = Plus@@(commonFactors[tensorTe[\[Nu]1, \[Nu]2]] * tensorTe[\[Nu]1, \[Nu]2])


(* gauge invariance *)
FCE@Contract[
	{
		electronTrace["iMd2"] FV[qe, \[Nu]1],
		electronTrace["iMd2"] FV[qe, \[Nu]2]
	}
] /. ruleQe2toQp2 /.{
	FV[qp,\[Nu]_]:>FV[qe, \[Nu]]-FV[k, \[Nu]]
} // Factor


(* ::Subsection::Closed:: *)
(*proton line trace*)


Timing[tensorBrD[1,1][\[Nu]1,\[Nu]2];]
Timing[tensorBrD[2,1][\[Nu]1,\[Nu]2];]
Timing[tensorBrD[1,2][\[Nu]1,\[Nu]2];]
Timing[tensorBrD[2,2][\[Nu]1,\[Nu]2];]


SP[t1, t1] - Md^2//ExpandScalarProduct//FCE//Expand//StandardForm
SP[t2, t2] - Md^2//ExpandScalarProduct//FCE//Expand//StandardForm


protonTrace["iMd2"][1, 1] = (*1/((-Md^2+Mp^2+xP+xq)^2+Gd^2 Md^2)*) Collect[
	Plus@@(commonFactors[tensorBrD[1, 1][\[Nu]1, \[Nu]2]]*tensorBrD[1, 1][\[Nu]1, \[Nu]2]),
	{MT[__],FV[a__]FV[b__]},
	Factor
];
protonTrace["iMd2"][1, 2] = (*(-Md^2+Mp^2+xP+xq)/((-Md^2+Mp^2+xP+xq)^2+Gd^2 Md^2) 1/(-Md^2+Mp^2-xP+xq)*) Collect[
	(
		Plus@@(commonFactors[tensorBrD[1, 2][\[Nu]1, \[Nu]2]]*tensorBrD[1, 2][\[Nu]1, \[Nu]2])
	),
	{MT[__],FV[a__]FV[b__]},
	Factor
];
protonTrace["iMd2"][2, 1] = (*(-Md^2+Mp^2+xP+xq)/((-Md^2+Mp^2+xP+xq)^2+Gd^2 Md^2) 1/(-Md^2+Mp^2-xP+xq)*) Collect[
	(
		Plus@@(commonFactors[tensorBrD[2, 1][\[Nu]1, \[Nu]2]]*tensorBrD[2, 1][\[Nu]1, \[Nu]2])
	),
	{MT[__],FV[a__]FV[b__]},
	Factor
];
protonTrace["iMd2"][2, 2] = (*1/(-Md^2+Mp^2-xP+xq)^2*) Collect[
	Plus@@(commonFactors[tensorBrD[2, 2][\[Nu]1, \[Nu]2]]*tensorBrD[2, 2][\[Nu]1, \[Nu]2]),
	{MT[__],FV[a__]FV[b__]},
	Factor
];


(* gauge invariance *)
Timing[FCE@Contract[
	{
		protonTrace["iMd2"][1, 1] FV[qe, \[Nu]1],
		protonTrace["iMd2"][1, 1] FV[qe, \[Nu]2]
	}
]/.ruleQp2toQe2/.{FV[qp,a_]:>FV[qe,a]-FV[k,a]}//Factor]
Timing[FCE@Contract[
	{
		protonTrace["iMd2"][1, 2] FV[qe, \[Nu]1],
		protonTrace["iMd2"][1, 2] FV[qe, \[Nu]2]
	}
]/.ruleQp2toQe2/.{FV[qp,a_]:>FV[qe,a]-FV[k,a]}//Factor]
Timing[FCE@Contract[
	{
		protonTrace["iMd2"][2, 2] FV[qe, \[Nu]1],
		protonTrace["iMd2"][2, 2] FV[qe, \[Nu]2]
	}
]/.ruleQp2toQe2/.{FV[qp,a_]:>FV[qe,a]-FV[k,a]}//Factor]


(* ::Subsection::Closed:: *)
(*contraction and decomposition*)


result["iMd2"][i_, j_] := result["iMd2"][i, j] = (FCE@Contract[protonTrace["iMd2"][i, j] * electronTrace["iMd2"]]/.ruleQp2toQe2//Together);
Timing[result["iMd2"][1,1];]
Timing[result["iMd2"][1,2];]
Timing[result["iMd2"][2,1];]
Timing[result["iMd2"][2,2];]


resultNumeratorSimplified["iMd2"][i_, j_] := resultNumeratorSimplified["iMd2"][i, j] = 
With[
	{numeratorList = List@@Expand[Numerator[result["iMd2"][i,j]]]},
	(Simplify[Plus@@Cases[numeratorList, #]])& /@ Table[decompositionPattern[k],{k,0,10}]
]
Timing[resultNumeratorSimplified["iMd2"][1, 1];]
Timing[resultNumeratorSimplified["iMd2"][2, 1];]
Timing[resultNumeratorSimplified["iMd2"][1, 2];]
Timing[resultNumeratorSimplified["iMd2"][2, 2];]


resultSimplified["iMd2"][i_, j_] := Plus@@(resultNumeratorSimplified["iMd2"][i, j]) / Denominator[result["iMd2"][i, j]];
Timing[(resultSimplified["iMd2"][1,1] - result["iMd2"][1,1])//Simplify]
Timing[(resultSimplified["iMd2"][1,2] - result["iMd2"][1,2])//Simplify]
Timing[(resultSimplified["iMd2"][2,1] - result["iMd2"][2,1])//Simplify]
Timing[(resultSimplified["iMd2"][2,2] - result["iMd2"][2,2])//Simplify]


(* ::Subsection::Closed:: *)
(*final result*)


iMd2[1, 1] = (Z^4 e^6)/(qe2^2) * 1/((-Md^2+Mp^2+xP+xq)^2+Gd^2 Md^2) * (resultSimplified["iMd2"][1, 1]);
iMd2[1, 2] = (Z^4 e^6)/(qe2^2) * (-Md^2+Mp^2+xP+xq)/((-Md^2+Mp^2+xP+xq)^2+Gd^2 Md^2) 1/(-Md^2+Mp^2-xP+xq) * (
		resultSimplified["iMd2"][1, 2] + resultSimplified["iMd2"][2, 1]
);
iMd2[2, 2] = (Z^4 e^6)/(qe2^2) * 1/(-Md^2+Mp^2-xP+xq)^2 * (resultSimplified["iMd2"][2, 2]);


With[{file = OpenWrite["./src/iMd2.cpp"]},
	WriteString[file, "#include \"constants.h\"\n"];
	WriteString[file, "#include \"proton_formfactors.h\"\n"];
	WriteString[file, "#include \"delta_formfactors.h\"\n"];
	WriteString[file, "#include <cmath>\n"];
	WriteString[file, "\n"];
	WriteString[file, "double iMd211(double nu, double qe2, double, double kK, double kP, double kq) {\n"];
	WriteString[file, "return ("];
	Write[file, CForm[iMd2[1,1]]/.{
		me^2->me2, Mp^2->Mp2, Md^2->Md2, Gd^2->Gd2,
		e^6->e6, Z^2->Z2, Z^4->Z4
	} /. {
		xK->kK,
		xP->kP,
		xq->kq,
		Power->pow
	}];
	WriteString[file, ");\n}"];
	WriteString[file, "\n"];
	WriteString[file, "double iMd212(double nu, double qe2, double, double kK, double kP, double kq) {\n"];
	WriteString[file, "return ("];
	Write[file, CForm[iMd2[1, 2]]/.{
		me^2->me2, Mp^2->Mp2, Md^2->Md2, Gd^2->Gd2,
		e^6->e6, Z^2->Z2, Z^4->Z4
	} /. {
		xK->kK,
		xP->kP,
		xq->kq,
		Power->pow
	}];
	WriteString[file, ");\n}"];
	WriteString[file, "\n"];
	WriteString[file, "double iMd222(double nu, double qe2, double, double kK, double kP, double kq) {\n"];
	WriteString[file, "return ("];
	Write[file, CForm[iMd2[2, 2]]/.{
		me^2->me2, Mp^2->Mp2, Md^2->Md2, Gd^2->Gd2,
		e^6->e6, Z^2->Z2, Z^4->Z4
	} /. {
		xK->kK,
		xP->kP,
		xq->kq,
		Power->pow
	}];
	WriteString[file, ");\n}"];
	Close[file];
];


(* ::Section::Closed:: *)
(*Bremsstrahlung 2(iMp iMe^+)*)


<<Bremsstrahlung`Traces`BremsstrahlungElectronProtonInterference`


(* ::Subsection::Closed:: *)
(*electron line trace*)


Timing[(* 0.1s *)tensorL[1,1][\[Mu],\[Nu]1,\[Nu]2];]
Timing[(* 0.1s *)tensorL[1,2][\[Mu],\[Nu]1,\[Nu]2];]
Timing[(* 0.1s *)tensorL[2,1][\[Mu],\[Nu]1,\[Nu]2];]
Timing[(* 0.1s *)tensorL[2,2][\[Mu],\[Nu]1,\[Nu]2];]


SP[r1,r1]-me^2//ExpandScalarProduct//FCE//Expand
SP[r2,r2]-me^2//ExpandScalarProduct//FCE//Expand


electronTrace["iMpiMe"][1, 1] = 1/(xK-xq) * Plus@@(commonFactors[tensorL[1, 1][\[Mu], \[Nu]1, \[Nu]2]]*tensorL[1, 1][\[Mu], \[Nu]1, \[Nu]2]);
electronTrace["iMpiMe"][1, 2] = 1/(-xK-xq)* Plus@@(commonFactors[tensorL[1, 2][\[Mu], \[Nu]1, \[Nu]2]]*tensorL[1, 2][\[Mu], \[Nu]1, \[Nu]2]);
electronTrace["iMpiMe"][2, 1] = 1/(xK-xq) * Plus@@(commonFactors[tensorL[2, 1][\[Mu], \[Nu]1, \[Nu]2]]*tensorL[2, 1][\[Mu], \[Nu]1, \[Nu]2]);
electronTrace["iMpiMe"][2, 2] = 1/(-xK-xq)* Plus@@(commonFactors[tensorL[2, 2][\[Mu], \[Nu]1, \[Nu]2]]*tensorL[2, 2][\[Mu], \[Nu]1, \[Nu]2]);


(* gauge invariance *)
Timing[(*~0.3s*)
FCE@Contract[
	{
		(electronTrace["iMpiMe"][1, 1]+electronTrace["iMpiMe"][1, 2]) FV[qe, \[Nu]1],
		(electronTrace["iMpiMe"][1, 1]+electronTrace["iMpiMe"][1, 2]) FV[k, \[Mu]],
		(electronTrace["iMpiMe"][1, 1]+electronTrace["iMpiMe"][1, 2]) FV[qp, \[Nu]2]
	}
] /. ruleQe2toQp2 /.{
	FV[qe,\[Nu]_]:>FV[qp, \[Nu]]+FV[k, \[Nu]]
} // Factor
]


(* ::Subsection::Closed:: *)
(*proton line trace*)


Timing[(* ~4s *)tensorP[1,1][\[Mu],\[Nu]1,\[Nu]2];]
Timing[(* ~3s *)tensorP[1,2][\[Mu],\[Nu]1,\[Nu]2];]
Timing[(* ~4s *)tensorP[2,1][\[Mu],\[Nu]1,\[Nu]2];]
Timing[(* ~3s *)tensorP[2,2][\[Mu],\[Nu]1,\[Nu]2];]


SP[t1, t1] - Mp^2//ExpandScalarProduct//FCE//Expand
SP[t2, t2] - Mp^2//ExpandScalarProduct//FCE//Expand


protonTrace["iMpiMe"][1, 1] = 1/(xP+xq) Plus@@(commonFactors[tensorP[1, 1][\[Mu], \[Nu]1, \[Nu]2]]*tensorP[1, 1][\[Mu], \[Nu]1, \[Nu]2]);
protonTrace["iMpiMe"][2, 1] = 1/(xq-xP) Plus@@(commonFactors[tensorP[2, 1][\[Mu], \[Nu]1, \[Nu]2]]*tensorP[2, 1][\[Mu], \[Nu]1, \[Nu]2]);
protonTrace["iMpiMe"][1, 2] = 1/(xP+xq) Plus@@(commonFactors[tensorP[1, 2][\[Mu], \[Nu]1, \[Nu]2]]*tensorP[1, 2][\[Mu], \[Nu]1, \[Nu]2]);
protonTrace["iMpiMe"][2, 2] = 1/(xq-xP) Plus@@(commonFactors[tensorP[2, 2][\[Mu], \[Nu]1, \[Nu]2]]*tensorP[2, 2][\[Mu], \[Nu]1, \[Nu]2]);


(* gauge invariance *)
Timing[(* ~6s *)FCE@Contract[
	{
		(protonTrace["iMpiMe"][1,1]+protonTrace["iMpiMe"][2,1]) FV[k, \[Mu]],
		(protonTrace["iMpiMe"][1,1]+protonTrace["iMpiMe"][2,1]) FV[qe, \[Nu]1],
		(protonTrace["iMpiMe"][1,1]+protonTrace["iMpiMe"][2,1]) FV[qp, \[Nu]2]
	}
] /. ruleQe2toQp2 /.{
	FV[qe,\[Nu]_]:>FV[qp, \[Nu]]+FV[k, \[Nu]]
} // Factor]


(* ::Subsection::Closed:: *)
(*contraction and decomposition*)


Timing[(*~17s*)
	result["iMpiMe"] = (
		FCE@Contract[protonTrace["iMpiMe"][1,1] * electronTrace["iMpiMe"][1,1]]+
		FCE@Contract[protonTrace["iMpiMe"][2,1] * electronTrace["iMpiMe"][2,1]]+
		FCE@Contract[protonTrace["iMpiMe"][1,2] * electronTrace["iMpiMe"][1,2]]+
		FCE@Contract[protonTrace["iMpiMe"][2,2] * electronTrace["iMpiMe"][2,2]]
		//.ruleQp2toQe2//Together
	)//.(Reverse/@ruleQp2toQe2);
]


Timing[(*~4s*)
	resultNumeratorSimplified["iMpiMe"] = With[
	{numeratorList = List@@Expand[Numerator[result["iMpiMe"]]]},
	Simplify[Plus@@Cases[numeratorList,#]]& /@ Table[decompositionPattern[n],{n, 2, 10}]
];]


(*soft photon approximation term *)
((resultNumeratorSimplified["iMpiMe"][[1]]/ Denominator[result["iMpiMe"]])/.ruleFtoG/.{
	SP[qe,qe]->SP[q,q], SP[qp,qp]->SP[q,q]
})//FullSimplify


resultSimplified["iMpiMe"] = Plus@@(resultNumeratorSimplified["iMpiMe"]) / Denominator[result["iMpiMe"]]
(resultSimplified["iMpiMe"] - result["iMpiMe"])//Simplify


(* ::Subsection::Closed:: *)
(*final result*)


iMpiMe = 2(Z^3 e^6)/(qe2 qp2) resultSimplified["iMpiMe"];


With[{file = OpenWrite["./src/iMpiMe.cpp"]},
	WriteString[file, "#include \"constants.h\"\n"];
	WriteString[file, "#include \"proton_formfactors.h\"\n"];
	WriteString[file, "#include <cmath>\n"];
	WriteString[file, "\n"];
	WriteString[file, "double iMpiMe(double nu, double qe2, double qp2, double kK, double kP, double kq) {\n"];
	WriteString[file, "return ("];
	Write[file, CForm[iMpiMe]/.{
		me^2->me2, Mp^2->Mp2,
		e^6->e6, Z^2->Z2, Z^3->Z3
	} /. {
		xK->kK,
		xP->kP,
		xq->kq,
		Power->pow
	}];
	WriteString[file, ");\n}"];
	Close[file];
]


(* ::Section::Closed:: *)
(*Bremsstrahlung 2(iMd iMe^+)*)


<<Bremsstrahlung`Traces`BremsstrahlungElectronDeltaInterference`


(* ::Subsection::Closed:: *)
(*electron line trace*)


Timing[(* 0.1s *)tensorL[1,1][\[Mu],\[Nu]1,\[Nu]2];]
Timing[(* 0.1s *)tensorL[1,2][\[Mu],\[Nu]1,\[Nu]2];]
Timing[(* 0.1s *)tensorL[2,1][\[Mu],\[Nu]1,\[Nu]2];]
Timing[(* 0.1s *)tensorL[2,2][\[Mu],\[Nu]1,\[Nu]2];]


SP[r1,r1]-me^2//ExpandScalarProduct//FCE//Expand
SP[r2,r2]-me^2//ExpandScalarProduct//FCE//Expand


electronTrace["iMdiMe"][1, 1] = 1/(xK-xq) * Plus@@(commonFactors[tensorL[1, 1][\[Mu], \[Nu]1, \[Nu]2]]*tensorL[1, 1][\[Mu], \[Nu]1, \[Nu]2]);
electronTrace["iMdiMe"][1, 2] = 1/(-xK-xq)* Plus@@(commonFactors[tensorL[1, 2][\[Mu], \[Nu]1, \[Nu]2]]*tensorL[1, 2][\[Mu], \[Nu]1, \[Nu]2]);
electronTrace["iMdiMe"][2, 1] = 1/(xK-xq) * Plus@@(commonFactors[tensorL[2, 1][\[Mu], \[Nu]1, \[Nu]2]]*tensorL[2, 1][\[Mu], \[Nu]1, \[Nu]2]);
electronTrace["iMdiMe"][2, 2] = 1/(-xK-xq)* Plus@@(commonFactors[tensorL[2, 2][\[Mu], \[Nu]1, \[Nu]2]]*tensorL[2, 2][\[Mu], \[Nu]1, \[Nu]2]);


(* gauge invariance *)
Timing[(*~0.3s*)FCE@Contract[
	{
		(electronTrace["iMdiMe"][1, 1]+electronTrace["iMdiMe"][1, 2]) FV[qe, \[Nu]1],
		(electronTrace["iMdiMe"][1, 1]+electronTrace["iMdiMe"][1, 2]) FV[k, \[Mu]],
		(electronTrace["iMdiMe"][1, 1]+electronTrace["iMdiMe"][1, 2]) FV[qp, \[Nu]2]
	}
] /. ruleQe2toQp2 /.{
	FV[qe,\[Nu]_]:>FV[qp, \[Nu]]+FV[k, \[Nu]]
} // Factor]


(* ::Subsection::Closed:: *)
(*proton line trace*)


Timing[(* ~85s  *)tensorD[1,1][\[Mu],\[Nu]1,\[Nu]2];]
Timing[(* ~100s *)tensorD[2,1][\[Mu],\[Nu]1,\[Nu]2];]
Timing[(* ~80s  *)tensorD[1,2][\[Mu],\[Nu]1,\[Nu]2];]
Timing[(* ~96s  *)tensorD[2,2][\[Mu],\[Nu]1,\[Nu]2];]


SP[t1, t1] - Md^2//ExpandScalarProduct//FCE//Expand//StandardForm
SP[t2, t2] - Md^2//ExpandScalarProduct//FCE//Expand//StandardForm


protonTrace["iMdiMe"][1, 1] = (*(xP+xq-Md^2+Mp^2)/((xP+xq-Md^2+Mp^2)^2+Gd^2 Md^2) *)
	Plus@@(commonFactors[tensorD[1, 1][\[Mu], \[Nu]1, \[Nu]2]]*tensorD[1, 1][\[Mu], \[Nu]1, \[Nu]2]);
protonTrace["iMdiMe"][1, 2] = (*(xP+xq-Md^2+Mp^2)/((xP+xq-Md^2+Mp^2)^2+Gd^2 Md^2)*)
	Plus@@(commonFactors[tensorD[1, 2][\[Mu], \[Nu]1, \[Nu]2]]*tensorD[1, 2][\[Mu], \[Nu]1, \[Nu]2]);
protonTrace["iMdiMe"][2, 1] = (*1/(-xP+xq-Md^2+Mp^2)*)
	Plus@@(commonFactors[tensorD[2, 1][\[Mu], \[Nu]1, \[Nu]2]]*tensorD[2, 1][\[Mu], \[Nu]1, \[Nu]2]);
protonTrace["iMdiMe"][2, 2] =(* 1/(-xP+xq-Md^2+Mp^2)*)
	Plus@@(commonFactors[tensorD[2, 2][\[Mu], \[Nu]1, \[Nu]2]]*tensorD[2, 2][\[Mu], \[Nu]1, \[Nu]2]);


(* gauge invariance *)
Table[Print@Timing[(* ~4*20s *)FCE@Contract[
	{
		(protonTrace["iMdiMe"][i,j]) FV[k, \[Mu]],
		(protonTrace["iMdiMe"][i,j]) FV[qe, \[Nu]1],
		(protonTrace["iMdiMe"][i,j]) FV[qp, \[Nu]2]
	}
] /. ruleQe2toQp2 /.{
	FV[qe,\[Nu]_]:>FV[qp, \[Nu]]+FV[k, \[Nu]]
} // Factor],{i,1,2},{j,1,2}];


(* ::Subsection::Closed:: *)
(*contraction and decomposition*)


Timing[(*~70s*)
	result["iMdiMe"][1] = (
		FCE@Contract[protonTrace["iMdiMe"][1,1] * electronTrace["iMdiMe"][1,1]]+
		FCE@Contract[protonTrace["iMdiMe"][1,2] * electronTrace["iMdiMe"][1,2]]//.ruleQp2toQe2//Together
	)//.(Reverse/@ruleQp2toQe2);
]
Timing[(*~70s*)
	result["iMdiMe"][2] = (
		FCE@Contract[protonTrace["iMdiMe"][2,1] * electronTrace["iMdiMe"][2,1]]+
		FCE@Contract[protonTrace["iMdiMe"][2,2] * electronTrace["iMdiMe"][2,2]]//.ruleQp2toQe2//Together
	)//.(Reverse/@ruleQp2toQe2);
]


resultNumeratorSimplified["iMdiMe"][i_] := resultNumeratorSimplified["iMdiMe"][i] = With[
	{numeratorList = List@@Expand[Numerator[result["iMdiMe"][i]]]},
	Simplify[Plus@@Cases[numeratorList,#]]& /@ Table[decompositionPattern[n],{n, 2, 10}]
]
Timing[(*~170s*)resultNumeratorSimplified["iMdiMe"][1];]
Timing[(*~160s*)resultNumeratorSimplified["iMdiMe"][2];]


(*soft photon approximation term *)
((resultNumeratorSimplified["iMdiMe"][1][[1]]/ Denominator[result["iMdiMe"][1]])/.ruleFtoG/.Md->Mp/.{
	SP[qe, qe]->SP[q,q],SP[qp, qp]->SP[q,q]
}//Simplify)/.ruleG123toGEMC/.Md->Mp//Simplify
((resultNumeratorSimplified["iMdiMe"][2][[1]]/ Denominator[result["iMdiMe"][2]])/.ruleFtoG/.Md->Mp/.{
	SP[qe, qe]->SP[q,q],SP[qp, qp]->SP[q,q]
}//Simplify)/.ruleG123toGEMC/.Md->Mp//Simplify


resultSimplified["iMdiMe"][i_] := resultSimplified["iMdiMe"][i] = 
	Plus@@(resultNumeratorSimplified["iMdiMe"][i]) / Denominator[result["iMdiMe"][i]]
(resultSimplified["iMdiMe"][1] - result["iMdiMe"][1])//Simplify
(resultSimplified["iMdiMe"][2] - result["iMdiMe"][2])//Simplify


(* ::Subsection::Closed:: *)
(*final result*)


iMdiMe[1] = 2(Z^3 e^6)/(qe2 qp2) ((xP+xq-Md^2+Mp^2)/((xP+xq-Md^2+Mp^2)^2+Gd^2 Md^2) * resultSimplified["iMdiMe"][1]);
iMdiMe[2] = 2(Z^3 e^6)/(qe2 qp2) (1/(-xP+xq-Md^2+Mp^2) * resultSimplified["iMdiMe"][2]);


With[{file = OpenWrite["./src/iMdiMe.cpp"]},
	WriteString[file, "#include \"constants.h\"\n"];
	WriteString[file, "#include \"proton_formfactors.h\"\n"];
	WriteString[file, "#include \"delta_formfactors.h\"\n"];
	WriteString[file, "#include <cmath>\n"];
	WriteString[file, "\n"];
	WriteString[file, "double iMdiMe1(double nu, double qe2, double qp2, double kK, double kP, double kq) {\n"];
	WriteString[file, "return ("];
	Write[file, CForm[iMdiMe[1]]/.{
		me^2->me2, Mp^2->Mp2, Md^2->Md2, Gd^2->Gd2,
		e^6->e6, Z^2->Z2, Z^3->Z3
	} /. {
		xK->kK,
		xP->kP,
		xq->kq,
		Power->pow
	}];
	WriteString[file, ");\n}"];
	WriteString[file, "\n"];
	WriteString[file, "double iMdiMe2(double nu, double qe2, double qp2, double kK, double kP, double kq) {\n"];
	WriteString[file, "return ("];
	Write[file, CForm[iMdiMe[2]]/.{
		me^2->me2, Mp^2->Mp2, Md^2->Md2, Gd^2->Gd2,
		e^6->e6, Z^2->Z2, Z^3->Z3
	} /. {
		xK->kK,
		xP->kP,
		xq->kq,
		Power->pow
	}];
	WriteString[file, ");\n}"];
	Close[file];
]


(* ::Subsection::Closed:: *)
(*soft/hard decomposition*)


(* ::Subsubsection::Closed:: *)
(*electron line traces (soft/hard decomposition)*)


Timing[(* 0.1s *)tensorLSoft[1,1][\[Mu],\[Nu]1,\[Nu]2];]
Timing[(* 0.1s *)tensorLSoft[1,2][\[Mu],\[Nu]1,\[Nu]2];]
Timing[(* 0.1s *)tensorLSoft[2,1][\[Mu],\[Nu]1,\[Nu]2];]
Timing[(* 0.1s *)tensorLSoft[2,2][\[Mu],\[Nu]1,\[Nu]2];]
Timing[(* 0.1s *)tensorLHard[1,1][\[Mu],\[Nu]1,\[Nu]2];]
Timing[(* 0.1s *)tensorLHard[1,2][\[Mu],\[Nu]1,\[Nu]2];]
Timing[(* 0.1s *)tensorLHard[2,1][\[Mu],\[Nu]1,\[Nu]2];]
Timing[(* 0.1s *)tensorLHard[2,2][\[Mu],\[Nu]1,\[Nu]2];]


Table[Simplify[
	Plus@@(commonFactors[tensorL[i,j][\[Mu],\[Nu]1,\[Nu]2]]*tensorL[i,j][\[Mu],\[Nu]1,\[Nu]2]) == 
	Plus@@(
		commonFactors[tensorLSoft[i,j][\[Mu],\[Nu]1,\[Nu]2]]tensorLSoft[i,j][\[Mu],\[Nu]1,\[Nu]2]+
		commonFactors[tensorLHard[i,j][\[Mu],\[Nu]1,\[Nu]2]]tensorLHard[i,j][\[Mu],\[Nu]1,\[Nu]2]
	)
],{i,1,2},{j,1,2}]


electronTraceSoft["iMdiMe"][1, 1] = 1/(xK-xq) * Plus@@(commonFactors[tensorLSoft[1, 1][\[Mu], \[Nu]1, \[Nu]2]]*tensorLSoft[1, 1][\[Mu], \[Nu]1, \[Nu]2]);
electronTraceSoft["iMdiMe"][1, 2] = 1/(-xK-xq)* Plus@@(commonFactors[tensorLSoft[1, 2][\[Mu], \[Nu]1, \[Nu]2]]*tensorLSoft[1, 2][\[Mu], \[Nu]1, \[Nu]2]);
electronTraceSoft["iMdiMe"][2, 1] = 1/(xK-xq) * Plus@@(commonFactors[tensorLSoft[2, 1][\[Mu], \[Nu]1, \[Nu]2]]*tensorLSoft[2, 1][\[Mu], \[Nu]1, \[Nu]2]);
electronTraceSoft["iMdiMe"][2, 2] = 1/(-xK-xq)* Plus@@(commonFactors[tensorLSoft[2, 2][\[Mu], \[Nu]1, \[Nu]2]]*tensorLSoft[2, 2][\[Mu], \[Nu]1, \[Nu]2]);

electronTraceHard["iMdiMe"][1, 1] = 1/(xK-xq) * Plus@@(commonFactors[tensorLHard[1, 1][\[Mu], \[Nu]1, \[Nu]2]]*tensorLHard[1, 1][\[Mu], \[Nu]1, \[Nu]2]);
electronTraceHard["iMdiMe"][1, 2] = 1/(-xK-xq)* Plus@@(commonFactors[tensorLHard[1, 2][\[Mu], \[Nu]1, \[Nu]2]]*tensorLHard[1, 2][\[Mu], \[Nu]1, \[Nu]2]);
electronTraceHard["iMdiMe"][2, 1] = 1/(xK-xq) * Plus@@(commonFactors[tensorLHard[2, 1][\[Mu], \[Nu]1, \[Nu]2]]*tensorLHard[2, 1][\[Mu], \[Nu]1, \[Nu]2]);
electronTraceHard["iMdiMe"][2, 2] = 1/(-xK-xq)* Plus@@(commonFactors[tensorLHard[2, 2][\[Mu], \[Nu]1, \[Nu]2]]*tensorLHard[2, 2][\[Mu], \[Nu]1, \[Nu]2]);


(* gauge invariance *)
Timing[(*~0.3s*)FCE@Contract[
	{
		(electronTraceSoft["iMdiMe"][1, 1]+electronTraceSoft["iMdiMe"][1, 2]) FV[qe, \[Nu]1],
		(electronTraceSoft["iMdiMe"][1, 1]+electronTraceSoft["iMdiMe"][1, 2]) FV[k, \[Mu]],
		(electronTraceSoft["iMdiMe"][1, 1]+electronTraceSoft["iMdiMe"][1, 2]) FV[qp, \[Nu]2]
	}
] /. ruleQe2toQp2 /.{
	FV[qe,\[Nu]_]:>FV[qp, \[Nu]]+FV[k, \[Nu]]
} // Factor]
Timing[(*~0.3s*)FCE@Contract[
	{
		(electronTraceHard["iMdiMe"][1, 1]+electronTraceHard["iMdiMe"][1, 2]) FV[qe, \[Nu]1],
		(electronTraceHard["iMdiMe"][1, 1]+electronTraceHard["iMdiMe"][1, 2]) FV[k, \[Mu]],
		(electronTraceHard["iMdiMe"][1, 1]+electronTraceHard["iMdiMe"][1, 2]) FV[qp, \[Nu]2]
	}
] /. ruleQe2toQp2 /.{
	FV[qe,\[Nu]_]:>FV[qp, \[Nu]]+FV[k, \[Nu]]
} // Factor]


(* ::Subsubsection::Closed:: *)
(*contraction and decomposition (soft/hard decomposition)*)


Timing[(*~35s*)
	resultSoft["iMdiMe"][1] = (
		FCE@Contract[protonTrace["iMdiMe"][1,1] * electronTraceSoft["iMdiMe"][1,1]]+
		FCE@Contract[protonTrace["iMdiMe"][1,2] * electronTraceSoft["iMdiMe"][1,2]]//.ruleQp2toQe2//Together
	)//.(Reverse/@ruleQp2toQe2);
]
Timing[(*~40s*)
	resultSoft["iMdiMe"][2] = (
		FCE@Contract[protonTrace["iMdiMe"][2,1] * electronTraceSoft["iMdiMe"][2,1]]+
		FCE@Contract[protonTrace["iMdiMe"][2,2] * electronTraceSoft["iMdiMe"][2,2]]//.ruleQp2toQe2//Together
	)//.(Reverse/@ruleQp2toQe2);
]
Timing[(*~75s*)
	resultHard["iMdiMe"][1] = (
		FCE@Contract[protonTrace["iMdiMe"][1,1] * electronTraceHard["iMdiMe"][1,1]]+
		FCE@Contract[protonTrace["iMdiMe"][1,2] * electronTraceHard["iMdiMe"][1,2]]//.ruleQp2toQe2//Together
	)//.(Reverse/@ruleQp2toQe2);
]
Timing[(*~80s*)
	resultHard["iMdiMe"][2] = (
		FCE@Contract[protonTrace["iMdiMe"][2,1] * electronTraceHard["iMdiMe"][2,1]]+
		FCE@Contract[protonTrace["iMdiMe"][2,2] * electronTraceHard["iMdiMe"][2,2]]//.ruleQp2toQe2//Together
	)//.(Reverse/@ruleQp2toQe2);
]


resultNumeratorSimplifiedSoft["iMdiMe"][i_] := resultNumeratorSimplifiedSoft["iMdiMe"][i] = With[
	{numeratorList = List@@Expand[Numerator[resultSoft["iMdiMe"][i]]]},
	Simplify[Plus@@Cases[numeratorList,#]]& /@ Table[decompositionPattern[n],{n, 2, 10}]
]
resultNumeratorSimplifiedHard["iMdiMe"][i_] := resultNumeratorSimplifiedHard["iMdiMe"][i] = With[
	{numeratorList = List@@Expand[Numerator[resultHard["iMdiMe"][i]]]},
	Simplify[Plus@@Cases[numeratorList,#]]& /@ Table[decompositionPattern[n],{n, 2, 10}]
]
Timing[(*~150s*)resultNumeratorSimplifiedSoft["iMdiMe"][1];]
Timing[(*~130s*)resultNumeratorSimplifiedSoft["iMdiMe"][2];]
Timing[(*~140s*)resultNumeratorSimplifiedHard["iMdiMe"][1];]
Timing[(*~120s*)resultNumeratorSimplifiedHard["iMdiMe"][2];]


(*soft photon approximation term *)
((resultNumeratorSimplifiedSoft["iMdiMe"][1][[1]]/ Denominator[resultSoft["iMdiMe"][1]])/.ruleFtoG/.{
	Md->M, Mp->M
}/.{
	SP[qe, qe]->SP[q,q],SP[qp, qp]->SP[q,q]
}//Simplify)/.ruleG123toGEMC/.{
	Md->M, Mp->M
}//FullSimplify
((resultNumeratorSimplifiedSoft["iMdiMe"][2][[1]]/ Denominator[resultSoft["iMdiMe"][2]])/.ruleFtoG/.{
	Md->M, Mp->M
}/.{
	SP[qe, qe]->SP[q,q],SP[qp, qp]->SP[q,q]
}//Simplify)/.ruleG123toGEMC/.{
	Md->M, Mp->M
}//FullSimplify


Block[
	{tmp = (-MT[\[Lambda],\[Lambda]1])*
		LC[\[Lambda],\[Tau],\[Rho],\[Mu]]*FV[HoldForm[t1],\[Tau]]FV[k,\[Rho]](FV[HoldForm[p3],\[Mu]]/SP[HoldForm[p3],k]-FV[HoldForm[p1],\[Mu]]/SP[HoldForm[p1],k])*
		LC[\[Lambda]1,\[Tau]1,\[Sigma],\[Nu]]FV[HoldForm[t1],\[Tau]1]FV[qe,\[Sigma]]FV[K,\[Nu]],
		tmp2,tmp3
	},
	Print[tmp];
	tmp2 = Together[FCE@Contract[ReleaseHold@tmp]]/.ruleQp2toQe2//Together;
	tmp3 = (
		(Plus@@Cases[Numerator[tmp2],decompositionPattern[2]]//Expand)/.{Mp->M,SP[qe,qe]->SP[q,q]}
	)+
	(-16 M^2 me^2+nu^2) xq^2+SP[q,q] (nu xK (xP-xq)+4 me^2 xq (xP+xq)+4 M^2 (-xK^2+xq^2)+(xK^2-xq (xP+xq)) SP[q,q])//Expand
]


resultSimplifiedSoft["iMdiMe"][i_] := resultSimplifiedSoft["iMdiMe"][i] = 
	Plus@@(resultNumeratorSimplifiedSoft["iMdiMe"][i]) / Denominator[resultSoft["iMdiMe"][i]];
resultSimplifiedHard["iMdiMe"][i_] := resultSimplifiedHard["iMdiMe"][i] = 
	Plus@@(resultNumeratorSimplifiedHard["iMdiMe"][i]) / Denominator[resultHard["iMdiMe"][i]];
(resultSimplifiedSoft["iMdiMe"][1] - resultSoft["iMdiMe"][1])//Simplify
(resultSimplifiedSoft["iMdiMe"][2] - resultSoft["iMdiMe"][2])//Simplify
(resultSimplifiedHard["iMdiMe"][1] - resultHard["iMdiMe"][1])//Simplify
(resultSimplifiedHard["iMdiMe"][2] - resultHard["iMdiMe"][2])//Simplify


(* ::Subsubsection::Closed:: *)
(*final result (soft/hard decomposition)*)


iMdiMeSoft[1] = 2(Z^3 e^6)/(qe2 qp2) ((xP+xq-Md^2+Mp^2)/((xP+xq-Md^2+Mp^2)^2+Gd^2 Md^2) * resultSimplifiedSoft["iMdiMe"][1]);
iMdiMeSoft[2] = 2(Z^3 e^6)/(qe2 qp2) (1/(-xP+xq-Md^2+Mp^2) * resultSimplifiedSoft["iMdiMe"][2]);
iMdiMeHard[1] = 2(Z^3 e^6)/(qe2 qp2) ((xP+xq-Md^2+Mp^2)/((xP+xq-Md^2+Mp^2)^2+Gd^2 Md^2) * resultSimplifiedHard["iMdiMe"][1]);
iMdiMeHard[2] = 2(Z^3 e^6)/(qe2 qp2) (1/(-xP+xq-Md^2+Mp^2) * resultSimplifiedHard["iMdiMe"][2]);


With[{file = OpenWrite["./src/iMdiMeSoft.cpp"]},
	WriteString[file, "#include \"constants.h\"\n"];
	WriteString[file, "#include \"proton_formfactors.h\"\n"];
	WriteString[file, "#include \"delta_formfactors.h\"\n"];
	WriteString[file, "#include <cmath>\n"];
	WriteString[file, "\n"];
	WriteString[file, "double iMdiMeSoft1(double nu, double qe2, double qp2, double kK, double kP, double kq) {\n"];
	WriteString[file, "return ("];
	Write[file, CForm[iMdiMeSoft[1]]/.{
		me^2->me2, Mp^2->Mp2, Md^2->Md2, Gd^2->Gd2,
		e^6->e6, Z^2->Z2, Z^3->Z3
	} /. {
		xK->kK,
		xP->kP,
		xq->kq,
		Power->pow
	}];
	WriteString[file, ");\n}"];
	WriteString[file, "\n"];
	WriteString[file, "double iMdiMeSoft2(double nu, double qe2, double qp2, double kK, double kP, double kq) {\n"];
	WriteString[file, "return ("];
	Write[file, CForm[iMdiMeSoft[2]]/.{
		me^2->me2, Mp^2->Mp2, Md^2->Md2, Gd^2->Gd2,
		e^6->e6, Z^2->Z2, Z^3->Z3
	} /. {
		xK->kK,
		xP->kP,
		xq->kq,
		Power->pow
	}];
	WriteString[file, ");\n}"];
	Close[file];
]
With[{file = OpenWrite["./src/iMdiMeHard.cpp"]},
	WriteString[file, "#include \"constants.h\"\n"];
	WriteString[file, "#include \"proton_formfactors.h\"\n"];
	WriteString[file, "#include \"delta_formfactors.h\"\n"];
	WriteString[file, "#include <cmath>\n"];
	WriteString[file, "\n"];
	WriteString[file, "double iMdiMeHard1(double nu, double qe2, double qp2, double kK, double kP, double kq) {\n"];
	WriteString[file, "return ("];
	Write[file, CForm[iMdiMeHard[1]]/.{
		me^2->me2, Mp^2->Mp2, Md^2->Md2, Gd^2->Gd2,
		e^6->e6, Z^2->Z2, Z^3->Z3
	} /. {
		xK->kK,
		xP->kP,
		xq->kq,
		Power->pow
	}];
	WriteString[file, ");\n}"];
	WriteString[file, "\n"];
	WriteString[file, "double iMdiMeHard2(double nu, double qe2, double qp2, double kK, double kP, double kq) {\n"];
	WriteString[file, "return ("];
	Write[file, CForm[iMdiMeHard[2]]/.{
		me^2->me2, Mp^2->Mp2, Md^2->Md2, Gd^2->Gd2,
		e^6->e6, Z^2->Z2, Z^3->Z3
	} /. {
		xK->kK,
		xP->kP,
		xq->kq,
		Power->pow
	}];
	WriteString[file, ");\n}"];
	Close[file];
]


(* ::Section::Closed:: *)
(*Bremsstrahlung 2(iMd iMp^+)*)


<<Bremsstrahlung`Traces`BremsstrahlungProtonDeltaInterference`


(* ::Subsection::Closed:: *)
(*electron line trace*)


electronTrace["iMdiMp"] = Plus@@(commonFactors[tensorTe[\[Nu]1, \[Nu]2]] * tensorTe[\[Nu]1, \[Nu]2])


(* gauge invariance *)
FCE@Contract[
	{
		electronTrace["iMdiMp"] FV[qe, \[Nu]1],
		electronTrace["iMdiMp"] FV[qe, \[Nu]2]
	}
] /. ruleQe2toQp2 /.{
	FV[qp,\[Nu]_]:>FV[qe, \[Nu]]-FV[k, \[Nu]]
} // Factor


(* ::Subsection::Closed:: *)
(*proton line trace*)


Timing[(*~1000s*)tensorPD[1,1][\[Nu]1,\[Nu]2];]
Timing[(*~2000s*)tensorPD[2,1][\[Nu]1,\[Nu]2];]
Timing[(*~600s*)tensorPD[1,2][\[Nu]1,\[Nu]2];]
Timing[(*~400s*)tensorPD[2,2][\[Nu]1,\[Nu]2];]


SP[t1, t1] - Mp^2//ExpandScalarProduct//FCE//Expand//StandardForm
SP[t2, t2] - Mp^2//ExpandScalarProduct//FCE//Expand//StandardForm
SP[t1, t1] - Md^2//ExpandScalarProduct//FCE//Expand//StandardForm
SP[t2, t2] - Md^2//ExpandScalarProduct//FCE//Expand//StandardForm


Timing[(*~35s*)protonTrace["iMdiMp"][1] = (*(xP+xq+Mp^2-Md^2)/((xP+xq+Mp^2-Md^2)^2+Gd^2 Md^2)*) Collect[
	Expand[
		1/(xP+xq)*Plus@@(commonFactors[tensorPD[1, 1][\[Nu]1, \[Nu]2]]*tensorPD[1, 1][\[Nu]1, \[Nu]2])+
		1/(-xP+xq)*Plus@@(commonFactors[tensorPD[1, 2][\[Nu]1, \[Nu]2]]*tensorPD[1, 2][\[Nu]1, \[Nu]2])
	],
	{MT[__],FV[a__]FV[b__]},
	Factor
];]
Timing[(*~35s*)protonTrace["iMdiMp"][2] = (*1/(-Md^2+Mp^2-xP+xq)*)Collect[
	Expand[
		1/(xP+xq)*Plus@@(commonFactors[tensorPD[2, 1][\[Nu]1, \[Nu]2]]*tensorPD[2, 1][\[Nu]1, \[Nu]2])+
		1/(-xP+xq)*Plus@@(commonFactors[tensorPD[2, 2][\[Nu]1, \[Nu]2]]*tensorPD[2, 2][\[Nu]1, \[Nu]2])
	],
	{MT[__],FV[a__]FV[b__]},
	Factor
];]


(* gauge invariance *)
Timing[(*~120s*)FCE@Contract[
	{
		protonTrace["iMdiMp"][1] FV[qe, \[Nu]1],
		protonTrace["iMdiMp"][1] FV[qe, \[Nu]2]
	}
]/.ruleQp2toQe2/.{FV[qp,a_]:>FV[qe,a]-FV[k,a]}//Factor]
Timing[(*~120s*)FCE@Contract[
	{
		protonTrace["iMdiMp"][2] FV[qe, \[Nu]1],
		protonTrace["iMdiMp"][2] FV[qe, \[Nu]2]
	}
]/.ruleQp2toQe2/.{FV[qp,a_]:>FV[qe,a]-FV[k,a]}//Factor]


(* ::Subsection::Closed:: *)
(*contraction and decomposition*)


result["iMdiMp"][i_] := result["iMdiMp"][i] = FCE@Contract[protonTrace["iMdiMp"][i] * electronTrace["iMdiMp"]]/.ruleQp2toQe2//Together;
Timing[result["iMdiMp"][1];]
Timing[result["iMdiMp"][2];]


resultNumeratorSimplified["iMdiMp"][i_] := resultNumeratorSimplified["iMdiMp"][i] = With[
	{numeratorList = List@@Expand[Numerator[result["iMdiMp"][i]]]},
	Simplify[Plus@@Cases[numeratorList,#]]& /@ Table[decompositionPattern[n],{n, 2, 10}]
];
Timing[(*~400s*)resultNumeratorSimplified["iMdiMp"][1];]
Timing[(*~400s*)resultNumeratorSimplified["iMdiMp"][2];]


(*soft photon approximation term *)
((resultNumeratorSimplified["iMdiMp"][1][[1]]/ Denominator[result["iMdiMp"][1]])/.ruleFtoG/.Md->Mp//FullSimplify)/.ruleG123toGEMC/.Md->Mp//Simplify
((resultNumeratorSimplified["iMdiMp"][2][[1]]/ Denominator[result["iMdiMp"][2]])/.ruleFtoG/.Md->Mp//FullSimplify)/.ruleG123toGEMC/.Md->Mp//Simplify


resultSimplified["iMdiMp"][i_] := resultSimplified["iMdiMp"][i] = Plus@@(resultNumeratorSimplified["iMdiMp"][i]) / Denominator[result["iMdiMp"][i]]
(resultSimplified["iMdiMp"][1] - result["iMdiMp"][1])//Simplify
(resultSimplified["iMdiMp"][2] - result["iMdiMp"][2])//Simplify


(* ::Subsection::Closed:: *)
(*final result*)


iMdiMp[1] = 2(Z^4 e^6)/qe2^2 (xP+xq+Mp^2-Md^2)/((xP+xq+Mp^2-Md^2)^2+Gd^2 Md^2)(resultSimplified["iMdiMp"][1]);
iMdiMp[2] = 2(Z^4 e^6)/qe2^2 1/(-Md^2+Mp^2-xP+xq)(resultSimplified["iMdiMp"][2]);


With[{file = OpenWrite["./src/iMdiMp.cpp"]},
	WriteString[file, "#include \"constants.h\"\n"];
	WriteString[file, "#include \"proton_formfactors.h\"\n"];
	WriteString[file, "#include \"delta_formfactors.h\"\n"];
	WriteString[file, "#include <cmath>\n"];
	WriteString[file, "\n"];
	WriteString[file, "double iMdiMp1(double nu, double qe2, double, double kK, double kP, double kq) {\n"];
	WriteString[file, "return ("];
	Write[file, CForm[iMdiMp[1]]/.{
		me^2->me2, Mp^2->Mp2, Mp^4->Mp4,
		e^6->e6, Z^2->Z2, Z^4->Z4, Z^3->Z3
	} /. {
		xK->kK,
		xP->kP,
		xq->kq,
		Power->pow
	}];
	WriteString[file, ");\n}"];
	WriteString[file, "\n"];
	WriteString[file, "double iMdiMp2(double nu, double qe2, double, double kK, double kP, double kq) {\n"];
	WriteString[file, "return ("];
	Write[file, CForm[iMdiMp[2]]/.{
		me^2->me2, Mp^2->Mp2, Mp^4->Mp4,
		e^6->e6, Z^2->Z2, Z^4->Z4, Z^3->Z3
	} /. {
		xK->kK,
		xP->kP,
		xq->kq,
		Power->pow
	}];
	WriteString[file, ");\n}"];
	Close[file];
]
