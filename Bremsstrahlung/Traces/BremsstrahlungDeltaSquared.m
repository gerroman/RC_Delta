BeginPackage["Bremsstrahlung`Traces`BremsstrahlungDeltaSquared`",
	{
		"Bremsstrahlung`Traces`Particles`",
		"Bremsstrahlung`Traces`FormFactors`",
		"HighEnergyPhysics`FeynCalc`",
		"Bremsstrahlung`Traces`Vertices`",
		"Bremsstrahlung`Traces`Momenta`",
		"Bremsstrahlung`Traces`BremsstrahlungProtonSquared`"
	}
];


tensord::usage =
	"tensord[mu, nu][t, k, q]";

tensorBrD::usage = 
	"tensorBrD[i, j][nu1, nu2]";


Begin["Private`"];


GA = HighEnergyPhysics`FeynCalc`GA`GA;
FV = HighEnergyPhysics`FeynCalc`FV`FV;
GS = HighEnergyPhysics`FeynCalc`GS`GS;
SP = HighEnergyPhysics`FeynCalc`SP`SP;
MT = HighEnergyPhysics`FeynCalc`MT`MT;
FCE = HighEnergyPhysics`fctools`FeynCalcExternal`FeynCalcExternal;
SPE = HighEnergyPhysics`FeynCalc`ScalarProductExpand`ScalarProductExpand;


tensord[mu_, nu_][t_, k_, q_] := Module[{alpha, beta},
	Flatten@Outer[
		(
			#1.
			(GS[t] + Md).
			(-MT[alpha, beta] + (1/3) * GA[alpha, beta]).
			#2
		)&,
		gamma["d->gp"][mu, alpha][t, k],
		gamma["gp->d"][nu, beta][t, q]
	]
];
commonFactors[tensord[mu_, nu_][t_, k_, q_]] := Module[{alpha, beta},
	Flatten@Outer[
		(#1 * #2)&,
		commonFactors[gamma["d->gp"][mu, alpha][t, k]],
		commonFactors[gamma["gp->d"][nu, beta][t, q]]
	]
];
diracConjugate[tensord[mu_, nu_][t_, k_, q_]] := Module[{alpha, beta},
	Flatten@Outer[
		(
			#1.
			(-MT[beta, alpha] + (1/3) * GA[beta, alpha]).
			(GS[t] + Md).
			#2
		)&,
		diracConjugate[gamma["gp->d"][nu, beta][t, q]],
		diracConjugate[gamma["d->gp"][mu, alpha][t, k]]
	]
];
commonFactors[diracConjugate[tensord[mu_, nu_][t_, k_, q_]]] := Module[{alpha, beta},
	Flatten@Outer[
		(#1 * #2)&,
		commonFactors[diracConjugate[gamma["gp->d"][nu, beta][t, q]]],
		commonFactors[diracConjugate[gamma["d->gp"][mu, alpha][t, k]]]
	]
];


tensorBrD[i_, j_][nu1_, nu2_] := tensorBrD[i, j][nu1, nu2] =
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
						Hold[tensord[mu, nu1][t1, k, qe]],
						Hold[tensord[nu1, mu][t2, -qe, -k]]
					}[[i]] // ReleaseHold,
					{
						Hold[diracConjugate[tensord[mu, nu2][t1, k, qe]]],
						Hold[diracConjugate[tensord[nu2, mu][t2, -qe, -k]]]
					}[[j]] // ReleaseHold
				]
			]
		]),
		{MT[__], FV[a__]FV[b__]},
		Factor
	];
commonFactors[tensorBrD[i_, j_][nu1_, nu2_]] := Module[{mu},
	FCE@SPE@Flatten@Outer[
		(#1 * #2)&,
		{
			Hold[commonFactors[tensord[mu, nu1][t1, k, qe]]],
			Hold[commonFactors[tensord[nu1, mu][t2, -qe, -k]]]
		}[[i]] // ReleaseHold,
		{
			Hold[commonFactors[diracConjugate[tensord[mu, nu2][t1, k, qe]]]],
			Hold[commonFactors[diracConjugate[tensord[nu2, mu][t2, -qe, -k]]]]
		}[[j]] // ReleaseHold
	] // Factor
];


End[];


EndPackage[];
