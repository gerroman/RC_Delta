BeginPackage["Bremsstrahlung`Traces`BremsstrahlungElectronDeltaInterference`",
	{
		"Bremsstrahlung`Traces`Particles`",
		"Bremsstrahlung`Traces`FormFactors`",
		"HighEnergyPhysics`FeynCalc`",
		"Bremsstrahlung`Traces`Vertices`",
		"Bremsstrahlung`Traces`Momenta`",
		"Bremsstrahlung`Traces`BremsstrahlungElectronSquared`",
		"Bremsstrahlung`Traces`BremsstrahlungDeltaSquared`",
		"Bremsstrahlung`Traces`BremsstrahlungElectronProtonInterference`"
	}
];


tensorD::usage = "
	tensorD[i, j][mu, nu, rho]
"


Begin["Private`"];


GA = HighEnergyPhysics`FeynCalc`GA`GA;
FV = HighEnergyPhysics`FeynCalc`FV`FV;
GS = HighEnergyPhysics`FeynCalc`GS`GS;
SP = HighEnergyPhysics`FeynCalc`SP`SP;
MT = HighEnergyPhysics`FeynCalc`MT`MT;
FCE = HighEnergyPhysics`fctools`FeynCalcExternal`FeynCalcExternal;
SPE = HighEnergyPhysics`FeynCalc`ScalarProductExpand`ScalarProductExpand;


tensorD[i_, j_][mu_, nu1_, nu2_] := tensorD[i, j][mu, nu1, nu2] =
	Collect[
		Expand@FCE@SPE@Flatten@(Block[
			{p1, p2, p3, p4, k, qp, qe, r1, r2, t1, t2},
			Outer[
				(
					(1/2) * Tr[
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
				diracConjugate[gamma["gp->p"][nu2][qp]]
			]
		]),
		{MT[__]FV[__], FV[a__]FV[b__]FV[c__]},
		Factor
	];
commonFactors[tensorD[i_, j_][mu_, nu1_, nu2_]] :=
	FCE@SPE@Flatten@Outer[
		(#1 * #2)&,
		{
			Hold[commonFactors[tensord[mu, nu1][t1, k, qe]]],
			Hold[commonFactors[tensord[nu1, mu][t2, -qe, -k]]]
		}[[i]] // ReleaseHold,
		commonFactors[diracConjugate[gamma["gp->p"][nu2][qp]]]
	] // Factor;


End[];


EndPackage[];
