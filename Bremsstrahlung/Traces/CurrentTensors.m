BeginPackage["Bremsstrahlung`Traces`CurrentTensors`",
	{
		"Bremsstrahlung`Traces`Particles`",
		"Bremsstrahlung`Traces`FormFactors`",
		"HighEnergyPhysics`FeynCalc`",
		"Bremsstrahlung`Traces`Vertices`",
		"Bremsstrahlung`Traces`MomentaElasticScattering`",
		"Bremsstrahlung`Traces`MomentaDeltaExcitation`"
	}
];


electronCurrentTensorElastic::usage = "
	electronCurrentTensorElastic[mu, nu]
";
protonCurrentTensorElastic::usage = "
	protonCurrentTensorElastic[mu, nu]
";
electronCurrentTensorDelta::usage = "
	electronCurrentTensorDelta[mu, nu]
";
protonTransitionTensorDelta::usage = "
	protonTransitionTensorDelta[mu, nu]
";


Begin["Private`"];


GA = HighEnergyPhysics`FeynCalc`GA`GA;
FV = HighEnergyPhysics`FeynCalc`FV`FV;
GS = HighEnergyPhysics`FeynCalc`GS`GS;
SP = HighEnergyPhysics`FeynCalc`SP`SP;
MT = HighEnergyPhysics`FeynCalc`MT`MT;
FCE = HighEnergyPhysics`fctools`FeynCalcExternal`FCE;
SPE = HighEnergyPhysics`FeynCalc`ScalarProductExpand`ScalarProductExpand;


electronCurrentTensorElastic[mu_, nu_] := electronCurrentTensorElastic[mu, nu] = 
	Collect[
		Expand@FCE@SPE@Block[
			{p1el, p3el},
			(1/2) * Tr[GA[mu].(GS[p1el] + me).GA[nu].(GS[p3el] + me)]
		],
		{MT[___], FV[a___]FV[b___]},
		Factor
	];


electronCurrentTensorDelta[mu_, nu_] := electronCurrentTensorDelta[mu, nu] = 
	Collect[
		Expand@FCE@SPE@Block[
			{p1d, p3d},
			(1/2) * Tr[GA[mu].(GS[p1d] + me).GA[nu].(GS[p3d] + me)]
		],
		{MT[___], FV[a___]FV[b___]},
		Factor
	];


commonFactors[protonCurrentTensorElastic[mu_, nu_]] := 
	FCE@SPE@Flatten@Outer[
		(#1 * #2)&,
		commonFactors[gamma["gp->p"][mu][qel]],
		commonFactors[diracConjugate[gamma["gp->p"][nu][qel]]]
	];
protonCurrentTensorElastic[mu_, nu_] := protonCurrentTensorElastic[mu, nu] = 
	Collect[
		Expand@FCE@SPE@Flatten@(Block[{p2el, p4el, qel},
			Outer[
				((1/2)*Tr[
					#1.
					(GS[p2el] + Mp).
					#2.
					(GS[p4el] + Mp)
				])&,
				gamma["gp->p"][mu][qel],
				diracConjugate[gamma["gp->p"][nu][qel]]
			]
		]),
		{MT[___], FV[a___]FV[b___]},
		Factor
	];


commonFactors[protonTransitionTensorDelta[mu_, nu_]] := Module[{alpha, beta},
	FCE@SPE@Flatten@Outer[
		(#1 * #2)&,
		commonFactors[gamma["gp->d"][mu, beta][p4d, qd]],
		commonFactors[diracConjugate[gamma["gp->d"][nu, alpha][p4d, qd]]]
	]
];
protonTransitionTensorDelta[mu_, nu_] := protonTransitionTensorDelta[mu, nu] = 
	Collect[
		Expand@FCE@SPE@Flatten@(Block[{p2d, p4d, qd},
			Module[{alpha, beta},
				Outer[
					((1/2)*Tr[
						(-MT[alpha, beta]+(1/3)*GA[alpha, beta]).
						#1.
						(GS[p2d] + Mp).
						#2.
						(GS[p4d] + Md)
					])&,
					gamma["gp->d"][mu, beta][p4d, qd],
					diracConjugate[gamma["gp->d"][nu, alpha][p4d, qd]]
				]
			]
		]),
		{MT[___], FV[a___]FV[b___]},
		Factor
	];

(*
protonDeltaCurrentTensor[mu_, nu_][p2_, p4_] := Collect[
	Plus@@Flatten[
		Outer[
			(#1 * #2)&,
			commonFactors[gamma["gp->d"][mu,  beta][p4, p4 - p2]],
			commonFactors[diracConjugate[gamma["gp->d"][nu, alpha][p4, p4 - p2]]]
		] * Block[{p2, p4, q}
		Outer[(
				(1/2)*FCE@Tr[
					(-MT[alpha, beta]+(1/3)*GA[alpha, beta]).
					#1.
					(GS[p2] + Mp).
					#2.
					(GS[p4] + Md)
				]
			)&,
			gamma["gp->d"][mu, beta][p4, p4 - p2],
			diracConjugate[
				gamma["gp->d"][nu, alpha][p4, p4 - p2]
			]
		]
	],
	{MT[___], FV[a___]FV[b___]},
	Factor
];
*)


End[];


EndPackage[];

