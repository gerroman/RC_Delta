BeginPackage["Bremsstrahlung`Traces`Vertices`",
	{
		"HighEnergyPhysics`FeynCalc`",
		"Bremsstrahlung`Traces`Particles`",
		"Bremsstrahlung`Traces`FormFactors`"
	}
];


gamma::usage =
	"gamma[process][indexes][momenta]
	return list of gamma-matrix structures";
commonFactors::usage ==
	"commonFactors[gamma[process][indexes][momenta]]
	return list of form factors and numeric multiplyers";
diracConjugate::usage = 
	"diracConjugate[gamma[process][indexes][momenta]]
	return GA[0].gamma[proccess]^{+}.GA[0]";


Begin["Private`"];


GA = HighEnergyPhysics`FeynCalc`GA`GA;
FV = HighEnergyPhysics`FeynCalc`FV`FV;
GS = HighEnergyPhysics`FeynCalc`GS`GS;
SP = HighEnergyPhysics`FeynCalc`SP`SP;
MT = HighEnergyPhysics`FeynCalc`MT`MT;


SetAttributes[commonFactors, HoldAll];
SetAttributes[diracConjugate, HoldAll];


gamma["gp->p"][mu_][q_] := {
	GA[mu],
	GA[mu].GS[q] - GS[q].GA[mu]
};
commonFactors[gamma["gp->p"][mu_][q_]] := {
	F1[SP[q, q]],
	- F2[SP[q, q]] / (4 Mp)
};
diracConjugate[gamma["gp->p"][mu_][q_]] := {
	(* inverse order *)
	GA[mu],
	GS[q].GA[mu] - GA[mu].GS[q] 
};
commonFactors[diracConjugate[gamma["gp->p"][mu_][q_]]] := 
	(* same common factors *)
	commonFactors[gamma["gp->p"][mu][q]];


gamma["p->gp"][mu_][q_] := 
	diracConjugate[gamma["gp->p"][mu][q]];
commonFactors[gamma["p->gp"][mu_][q_]] := 
	commonFactors[diracConjugate[gamma["gp->p"][mu][q]]]
diracConjugate[gamma["p->gp"][mu_][q_]] := 
	gamma["gp->p"][mu][q];
commonFactors[diracConjugate[gamma["p->gp"][mu_][q_]]] :=
	commonFactors[gamma["gp->p"][mu][q]];


gamma["gp->d"][nu_, beta_][pd_, q_] := {
	MT[nu, beta] * GA[5].GS[q].GS[pd]
		- FV[pd, nu] * GA[5].GS[q].GA[beta]
		- SP[pd, q] * GA[5].GA[beta].GA[nu]
		+ FV[q,beta] * GA[5].GS[pd].GA[nu],
	FV[pd, nu] * FV[q, beta] * GA[5]
		- MT[nu, beta] * SP[pd, q] * GA[5],
	SP[q,q] * FV[pd, nu] * GA[5].GA[beta]
		- SP[q, q] * MT[nu, beta] * GA[5].GS[pd]
		+ FV[q, nu] * FV[q, beta] * GA[5].GS[pd]
		- FV[q, nu] * SP[pd, q] * GA[5].GA[beta]
};
commonFactors[gamma["gp->d"][nu_, beta_][pd_, q_]] := -Sqrt[2/3]*(1/(2*Md^2)){
	G1[SP[q,q]],
	G2[SP[q,q]],
	-G3[SP[q,q]]/(Md)
}
diracConjugate[gamma["gp->d"][nu_, beta_][pd_, q_]] := {
	(* inverse order and multiply by (-1) due to GA[5] anticommutation with GA[0] *)
	-MT[nu, beta] * GS[pd].GS[q].GA[5]
		+ FV[pd, nu] * GA[beta].GS[q].GA[5]
		+ SP[pd, q] * GA[nu].GA[beta].GA[5]
		- FV[q,beta] * GA[nu].GS[pd].GA[5],
	-FV[pd, nu] * FV[q, beta] * GA[5]
		+ MT[nu, beta] * SP[pd, q] * GA[5],
	-SP[q,q] * FV[pd, nu] * GA[beta].GA[5]
		+ SP[q, q] * MT[nu, beta] * GS[pd].GA[5]
		- FV[q, nu] * FV[q, beta] * GS[pd].GA[5]
		+ FV[q, nu] * SP[pd, q] * GA[beta].GA[5]
};
commonFactors[diracConjugate[gamma["gp->d"][nu_, beta_][pd_, q_]]] :=
	commonFactors[gamma["gp->d"][nu, beta][pd, q]];


gamma["d->gp"][nu_, beta_][pd_, q_] := 
	diracConjugate[gamma["gp->d"][nu, beta][pd, q]];
commonFactors[gamma["d->gp"][nu_,beta_][pd_, q_]] := 
	commonFactors[diracConjugate[gamma["gp->d"][nu, beta][pd, q]]];
diracConjugate[gamma["d->gp"][nu_, beta_][pd_, q_]] := 
	gamma["gp->d"][nu, beta][pd, q];
commonFactors[diracConjugate[gamma["d->gp"][nu_, beta_][pd_, q_]]] :=
	commonFactors[gamma["gp->d"][nu, beta][pd, q]];


End[];


EndPackage[];
