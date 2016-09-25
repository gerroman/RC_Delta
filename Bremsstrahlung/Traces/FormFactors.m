BeginPackage["Bremsstrahlung`Traces`FormFactors`",
	{
		"Bremsstrahlung`Traces`Particles`"
	}
];


F1::usage = "F1[q2], Dirac form factor of the proton";
F2::usage = "F2[q2], Pauli form factor of the proton";
GE::usage = "GE[q2], electric form factor of the proton";
GM::usage = "GM[q2], magnetic form factor of the proton";


G1::usage = "G1[q2], form factor Delta(1232)";
G2::usage = "G2[q2], form factor Delta(1232)";
G3::usage = "G3[q2], form factor Delta(1232)";
GEstar::usage = "GEstar[q2], electric form factor of Delta(1232)";
GMstar::usage = "GMstar[q2], magnetic form factor of Delta(1232)";
GCstar::usage = "GCstar[q2], Coulumb  form factor of Delta(1232)";


ruleFtoG::usage = "convert F1, F2 to GM, GE";
ruleGtoF::usage = "convert GM, GE to F1, F2";


ruleG123toGEMC::usage = "convert {G1, G2, G3} to {GE, GM, GC}";
ruleGEMCtoG123::usage = "convert {GE, GM, GC} to {G1, G2, G3}";


Begin["Private`"];


Format[F1, TraditionalForm] = Subscript["F", 1];
Format[F2, TraditionalForm] = Subscript["F", 2];
Format[GE, TraditionalForm] = Subscript["G", "E"];
Format[GM, TraditionalForm] = Subscript["G", "M"];
Format[G1, TraditionalForm] = Subscript["G", 1];
Format[G2, TraditionalForm] = Subscript["G", 2];
Format[G3, TraditionalForm] = Subscript["G", 3];
Format[GEstar, TraditionalForm] = Subsuperscript["G", "E", "*"];
Format[GMstar, TraditionalForm] = Subsuperscript["G", "M", "*"];
Format[GCstar, TraditionalForm] = Subsuperscript["G", "C", "*"];


ruleGtoF = {
	GE[q2_] :> F1[q2] + F2[q2] * q2 / (4 Mp^2),
	GM[q2_] :> F1[q2] + F2[q2]
};

ruleFtoG = {
	F1[q2_] :> (4*Mp^2*GE[q2] - q2*GM[q2])/(4*Mp^2 - q2),
	F2[q2_] :> 4*Mp^2*(GM[q2] - GE[q2])/(4*Mp^2 - q2)
};


ruleGEMCtoG123 = {
	GMstar[q2_] :> Mp / (3 * (Md + Mp)) * (
		((Md + Mp)^2 - q2) / Md^2 * G1[q2]
			-  (Md^2 - Mp^2 + q2) / (2*Md^2) * (G1[q2] - G2[q2])
			-  (-q2)/Md^2 * G3[q2]
	),
	GEstar[q2_] :> Mp / (3 * (Md + Mp)) * (
		-  (Md^2 - Mp^2 + q2) / (2*Md^2) * (G1[q2] - G2[q2])
			-  (-q2)/Md^2 * G3[q2]
	),
	GCstar[q2_] :> 2 * Mp / (3 * (Md + Mp)) * (
		- (G1[q2] - G2[q2])
			+  (Md^2 - Mp^2 + q2) / (2*Md^2) * G3[q2]
	)
};


ruleG123toGEMC ={
	G1[q2_] :> (3*Md^2*(Md + Mp) * (GMstar[q2] - GEstar[q2])) / (Mp*((Md + Mp)^2 - q2)),
	G2[q2_] :> (3*Md^2*(Md + Mp) * (
					((Md - Mp)*(Md + 3*Mp) + 3*q2) * GEstar[q2]
					+ ((Md - Mp)^2 - q2) * GMstar[q2]
					- 2*q2 * GCstar[q2]
				)) / (Mp * ((Md-Mp)^2-q2) * ((Md+Mp)^2-q2)),
	G3[q2_] :> (3*Md^2*(Md + Mp)*(
					(Md^2 - Mp^2 + q2)*GCstar[q2]
					- 4*Md^2*GEstar[q2]
				)) / (Mp * ((Md-Mp)^2-q2) * ((Md+Mp)^2-q2))
};


End[];


EndPackage[];
