BeginPackage["Bremsstrahlung`Traces`BremsstrahlungElectronSquared`",
	{
		"Bremsstrahlung`Traces`Particles`",
		"Bremsstrahlung`Traces`FormFactors`",
		"HighEnergyPhysics`FeynCalc`",
		"Bremsstrahlung`Traces`Vertices`",
		"Bremsstrahlung`Traces`Momenta`"
	}
];


tensorTp::usage =
	"tensorTp[nu1, nu2]";
tensorl::usage =
	"tensorl[mu, nu][r, k, q]";
tensorBrE::usage = 
	"tensorBrE[i, j][nu1, nu2]";


Begin["Private`"];


GA = HighEnergyPhysics`FeynCalc`GA`GA;
FV = HighEnergyPhysics`FeynCalc`FV`FV;
GS = HighEnergyPhysics`FeynCalc`GS`GS;
SP = HighEnergyPhysics`FeynCalc`SP`SP;
MT = HighEnergyPhysics`FeynCalc`MT`MT;
FCE = HighEnergyPhysics`fctools`FeynCalcExternal`FeynCalcExternal;
SPE = HighEnergyPhysics`FeynCalc`ScalarProductExpand`ScalarProductExpand;


tensorTp[nu1_, nu2_] := tensorTp[nu1, nu2] = 
	Collect[Expand@FCE@SPE@Flatten@(Block[
			{p1, p2, p3, p4, k, qp, qe, r1, r2, t1, t2},
			Outer[
				((1/2) * Tr[
					#1.
					(GS[p2] + Mp).
					#2.
					(GS[p4] + Mp)
				])&,
				gamma["gp->p"][nu1][qp],
				diracConjugate[gamma["gp->p"][nu2][qp]]
			]
		]),
		{MT[__], FV[a__]FV[b__]},
		Factor
	];
commonFactors[tensorTp[nu1_, nu2_]] :=
	FCE@SPE@Flatten@Outer[
		(#1*#2)&,
		commonFactors[gamma["gp->p"][nu1][qp]],
		commonFactors[diracConjugate[gamma["gp->p"][nu2][qp]]]
	] // Factor;



tensorl[mu_, nu_][r_, k_, q_] := {
	GA[mu].(GS[r] + me).GA[nu]
};
commonFactors[tensorl[mu_, nu_][r_, k_, q_]] := {
	1
};
diracConjugate[tensorl[mu_, nu_][r_, k_, q_]] := {
	GA[nu].(GS[r] + me).GA[mu]
};
commonFactors[diracConjugate[tensorl[mu_, nu_][r_, k_, q_]]] := {
	1
};


tensorBrE[i_, j_][nu1_, nu2_] := tensorBrE[i, j][nu1, nu2] = 
	Collect[
		Expand@FCE@SPE@Flatten@(Block[
			{p1, p2, p3, p4, k, qp, qe, r1, r2, t1, t2},
			Module[{mu},
				Outer[
					((-1/2) * Tr[
						#1.
						(GS[p1] + me).
						#2.
						(GS[p3] + me)
					])&,
					{
						Hold[tensorl[mu, nu1][r1, k, -qp]],
						Hold[tensorl[nu1, mu][r2, qp, -k]]
					}[[i]] // ReleaseHold,
					{
						Hold[diracConjugate[tensorl[mu, nu2][r1, k, -qp]]],
						Hold[diracConjugate[tensorl[nu2, mu][r2, qp, -k]]]
					}[[j]] // ReleaseHold
				]
			]
		]),
		{MT[__], FV[a__]FV[b__]},
		Factor
	];
commonFactors[tensorBrE[i_, j_][nu1_, nu2_]] := Module[{mu},
	FCE@SPE@Flatten@Outer[
		(#1*#2)&,
		{
			Hold[commonFactors[tensorl[mu, nu1][r1, k, -qp]]],
			Hold[commonFactors[tensorl[nu1, mu][r2, qp, -k]]]
		}[[i]] // ReleaseHold,
		{
			Hold[commonFactors[diracConjugate[tensorl[mu, nu2][r1, k, -qp]]]],
			Hold[commonFactors[diracConjugate[tensorl[nu2, mu][r2, qp, -k]]]]
		}[[j]] // ReleaseHold
	] // Factor
];


End[];


EndPackage[];
