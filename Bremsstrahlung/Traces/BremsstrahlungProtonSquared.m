BeginPackage["Bremsstrahlung`Traces`BremsstrahlungProtonSquared`",
	{
		"Bremsstrahlung`Traces`Particles`",
		"Bremsstrahlung`Traces`FormFactors`",
		"HighEnergyPhysics`FeynCalc`",
		"Bremsstrahlung`Traces`Vertices`",
		"Bremsstrahlung`Traces`Momenta`"
	}
];


tensorTe::usage =
	"tensorTp[nu1, nu2]";

tensorp::usage =
	"tensorp[mu, nu][t, k, q]";

tensorBrP::usage = 
	"tensorBrP[i, j][nu1, nu2]";


Begin["Private`"];


GA = HighEnergyPhysics`FeynCalc`GA`GA;
FV = HighEnergyPhysics`FeynCalc`FV`FV;
GS = HighEnergyPhysics`FeynCalc`GS`GS;
SP = HighEnergyPhysics`FeynCalc`SP`SP;
MT = HighEnergyPhysics`FeynCalc`MT`MT;
FCE = HighEnergyPhysics`fctools`FeynCalcExternal`FeynCalcExternal;
SPE = HighEnergyPhysics`FeynCalc`ScalarProductExpand`ScalarProductExpand;


tensorTe[nu1_, nu2_] := tensorTe[nu1, nu2] = 
	Collect[Expand@FCE@SPE@Flatten@(Block[
			{p1, p2, p3, p4, k, qp, qe, r1, r2, t1, t2},
			Outer[
				((1/2) * Tr[
					#1.
					(GS[p1] + me).
					#2.
					(GS[p3] + me)
				])&,
				{GA[nu1]},
				{GA[nu2]}
			]
		]),
		{MT[__], FV[a__]FV[b__]},
		Factor
	];
commonFactors[tensorTe[nu1_, nu2_]] :=
	FCE@SPE@Flatten@Outer[
		(#1*#2)&,
		{1},
		{1}
	] // Factor;


tensorp[mu_, nu_][t_, k_, q_] :=
	Flatten@Outer[
		(#1.(GS[t] + Mp).#2)&,
		gamma["p->gp"][mu][k],
		gamma["gp->p"][nu][q]
	];
commonFactors[tensorp[mu_, nu_][t_, k_, q_]] :=
	Flatten@Outer[
		(#1 * #2)&,
		commonFactors[gamma["p->gp"][mu][k]],
		commonFactors[gamma["gp->p"][nu][q]]
	];
diracConjugate[tensorp[mu_, nu_][t_, k_, q_]] :=
	Flatten@Outer[
		(#1.(GS[t] + Mp).#2)&,
		diracConjugate[gamma["gp->p"][nu][q]],
		diracConjugate[gamma["p->gp"][mu][k]]
	];
commonFactors[diracConjugate[tensorp[mu_, nu_][t_, k_, q_]]] :=
	Flatten@Outer[
		(#1 * #2)&,
		commonFactors[diracConjugate[gamma["gp->p"][nu][q]]],
		commonFactors[diracConjugate[gamma["p->gp"][mu][k]]]
	];


tensorBrP[i_, j_][nu1_, nu2_] := tensorBrP[i, j][nu1, nu2] =
	Collect[
		Expand@FCE@SPE@Flatten@(Block[
			{p1, p2, p3, p4, k, qp, qe, r1, r2, t1, t2},
			Module[{mu},
				Outer[
					(
						(-1/2) * Tr[
							#1.
							(GS[p2] + Mp).
							#2.
							(GS[p4] + Mp)
						]
					)&,
					{
						Hold[tensorp[mu, nu1][t1, k, qe]],
						Hold[tensorp[nu1, mu][t2, -qe, -k]]
					}[[i]] // ReleaseHold,
					{
						Hold[diracConjugate[tensorp[mu, nu2][t1, k, qe]]],
						Hold[diracConjugate[tensorp[nu2, mu][t2, -qe, -k]]]
					}[[j]] // ReleaseHold
				]
			]
		]),
		{MT[__], FV[a__]FV[b__]},
		Factor
	];
commonFactors[tensorBrP[i_, j_][nu1_, nu2_]] := Module[{mu},
	FCE@SPE@Flatten@Outer[
		(#1 * #2)&,
		{
			Hold[commonFactors[tensorp[mu, nu1][t1, k, qe]]],
			Hold[commonFactors[tensorp[nu1, mu][t2, -qe, -k]]]
		}[[i]] // ReleaseHold,
		{
			Hold[commonFactors[diracConjugate[tensorp[mu, nu2][t1, k, qe]]]],
			Hold[commonFactors[diracConjugate[tensorp[nu2, mu][t2, -qe, -k]]]]
		}[[j]] // ReleaseHold
	] // Factor
];


End[];


EndPackage[];
