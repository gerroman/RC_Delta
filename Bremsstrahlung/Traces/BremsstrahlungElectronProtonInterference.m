BeginPackage["Bremsstrahlung`Traces`BremsstrahlungElectronProtonInterference`",
	{
		"Bremsstrahlung`Traces`Particles`",
		"Bremsstrahlung`Traces`FormFactors`",
		"HighEnergyPhysics`FeynCalc`",
		"Bremsstrahlung`Traces`Vertices`",
		"Bremsstrahlung`Traces`Momenta`",
		"Bremsstrahlung`Traces`BremsstrahlungElectronSquared`",
		"Bremsstrahlung`Traces`BremsstrahlungProtonSquared`"
	}
];


tensorL::usage = "
	tensorL[i, j][mu, nu, rho]
";
tensorLSoft::usage = "
	tensorLSoft[i, j][mu, nu, rho]
";
tensorLHard::usage = "
	tensorLHard[i, j][mu, nu, rho]
";
tensorP::usage = "
	tensorP[i, j][mu, nu, rho]
";


Begin["Private`"];


GA = HighEnergyPhysics`FeynCalc`GA`GA;
FV = HighEnergyPhysics`FeynCalc`FV`FV;
GS = HighEnergyPhysics`FeynCalc`GS`GS;
SP = HighEnergyPhysics`FeynCalc`SP`SP;
MT = HighEnergyPhysics`FeynCalc`MT`MT;
FCE = HighEnergyPhysics`fctools`FeynCalcExternal`FeynCalcExternal;
SPE = HighEnergyPhysics`FeynCalc`ScalarProductExpand`ScalarProductExpand;


tensorL[i_, j_][mu_, nu1_, nu2_] := tensorL[i, j][mu, nu1, nu2] = 
	Collect[
		Expand@FCE@SPE@Flatten@(Block[
			{p1, p2, p3, p4, k, qp, qe, r1, r2, t1, t2},
			Outer[
				((1/2) * Tr[
					#1.
					(GS[p1] + me).
					#2.
					(GS[p3] + me)
				])&,
				{
					GA[nu1]
				},
				{
					Hold[diracConjugate[tensorl[mu, nu2][r1, k, -qp]]],
					Hold[diracConjugate[tensorl[nu2, mu][r2, qp, -k]]]
				}[[j]] // ReleaseHold
			]
		]),
		{MT[__]FV[__], FV[a__]FV[b__]FV[c__]},
		Factor
	];
commonFactors[tensorL[i_, j_][mu_, nu1_, nu2_]] :=
	FCE@SPE@Flatten@Outer[
		(#1*#2)&,
		{
			1
		},
		{
			Hold[commonFactors[diracConjugate[tensorl[mu, nu2][r1, k, -qp]]]],
			Hold[commonFactors[diracConjugate[tensorl[nu2, mu][r2, qp, -k]]]]
		}[[j]] // ReleaseHold
	] // Factor;


tensorLSoft[i_, j_][mu_, nu1_, nu2_] := tensorLSoft[i, j][mu, nu1, nu2] = 
	Collect[
		Expand@FCE@SPE@Flatten@(Block[
			{p1, p2, p3, p4, k, qp, qe, r1, r2, t1, t2},
			Outer[
				((1/2) * Tr[
					#1.
					(GS[p1] + me).
					#2.
					(GS[p3] + me)
				])&,
				{
					GA[nu1]
				},
				{
					{2*FV[p3, mu]*GA[nu2]},
					{2*FV[p1, mu]*GA[nu2]}
				}[[j]] // ReleaseHold
			]
		]),
		{MT[__]FV[__], FV[a__]FV[b__]FV[c__]},
		Factor
	];
commonFactors[tensorLSoft[i_, j_][mu_, nu1_, nu2_]] :=
	FCE@SPE@Flatten@Outer[
		(#1*#2)&,
		{
			1
		},
		{
			{1},
			{1}
		}[[j]]
	] // Factor;


tensorLHard[i_, j_][mu_, nu1_, nu2_] := tensorLHard[i, j][mu, nu1, nu2] = 
	Collect[
		Expand@FCE@SPE@Flatten@(Block[
			{p1, p2, p3, p4, k, qp, qe, r1, r2, t1, t2},
			Outer[
				((1/2) * Tr[
					#1.
					(GS[p1] + me).
					#2.
					(GS[p3] + me)
				])&,
				{
					GA[nu1]
				},
				{
					{GA[nu2].GS[k].GA[mu]},
					{GA[mu].GS[k].GA[nu2]}
				}[[j]] // ReleaseHold
			]
		]),
		{MT[__]FV[__], FV[a__]FV[b__]FV[c__]},
		Factor
	];
commonFactors[tensorLHard[i_, j_][mu_, nu1_, nu2_]] :=
	FCE@SPE@Flatten@Outer[
		(#1*#2)&,
		{
			1
		},
		{
			{1},
			{-1}
		}[[j]]
	] // Factor;


tensorP[i_, j_][mu_, nu1_, nu2_] := tensorP[i, j][mu, nu1, nu2] =
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
					Hold[tensorp[mu, nu1][t1, k, qe]],
					Hold[tensorp[nu1, mu][t2, -qe, -k]]
				}[[i]] // ReleaseHold,
				diracConjugate[gamma["gp->p"][nu2][qp]]
			]
		]),
		{MT[__]FV[__], FV[a__]FV[b__]FV[c__]},
		Factor
	];
commonFactors[tensorP[i_, j_][mu_, nu1_, nu2_]] :=
	FCE@SPE@Flatten@Outer[
		(#1 * #2)&,
		{
			Hold[commonFactors[tensorp[mu, nu1][t1, k, qe]]],
			Hold[commonFactors[tensorp[nu1, mu][t2, -qe, -k]]]
		}[[i]] // ReleaseHold,
		commonFactors[diracConjugate[gamma["gp->p"][nu2][qp]]]
	] // Factor;


End[];


EndPackage[];
