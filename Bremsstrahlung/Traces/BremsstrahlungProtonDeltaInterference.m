BeginPackage["Bremsstrahlung`Traces`BremsstrahlungProtonDeltaInterference`",
	{
		"Bremsstrahlung`Traces`Particles`",
		"Bremsstrahlung`Traces`FormFactors`",
		"HighEnergyPhysics`FeynCalc`",
		"Bremsstrahlung`Traces`Vertices`",
		"Bremsstrahlung`Traces`Momenta`",
		"Bremsstrahlung`Traces`BremsstrahlungProtonSquared`",
		"Bremsstrahlung`Traces`BremsstrahlungDeltaSquared`"
	}
];


tensorPD::usage = "
	tensorPD[i, j][nu1, nu2]
";


Begin["Private`"];


GA = HighEnergyPhysics`FeynCalc`GA`GA;
FV = HighEnergyPhysics`FeynCalc`FV`FV;
GS = HighEnergyPhysics`FeynCalc`GS`GS;
SP = HighEnergyPhysics`FeynCalc`SP`SP;
MT = HighEnergyPhysics`FeynCalc`MT`MT;
FCE = HighEnergyPhysics`fctools`FeynCalcExternal`FeynCalcExternal;
SPE = HighEnergyPhysics`FeynCalc`ScalarProductExpand`ScalarProductExpand;


tensorPD[i_, j_][nu1_, nu2_] := tensorPD[i, j][nu1, nu2] =
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
						Hold[diracConjugate[tensorp[mu, nu2][t1, k, qe]]],
						Hold[diracConjugate[tensorp[nu2, mu][t2, -qe, -k]]]
					}[[j]] // ReleaseHold
				]
			]
		]),
		{MT[__], FV[a__]FV[b__]},
		Factor
	];
commonFactors[tensorPD[i_, j_][nu1_, nu2_]] := Module[{mu},
	FCE@SPE@Flatten@Outer[
		(#1 * #2)&,
		{
			Hold[commonFactors[tensord[mu, nu1][t1, k, qe]]],
			Hold[commonFactors[tensord[nu1, mu][t2, -qe, -k]]]
		}[[i]] // ReleaseHold,
		{
			Hold[commonFactors[diracConjugate[tensorp[mu, nu2][t1, k, qe]]]],
			Hold[commonFactors[diracConjugate[tensorp[nu2, mu][t2, -qe, -k]]]]
		}[[j]] // ReleaseHold
	] // Factor
];


End[];


EndPackage[];
