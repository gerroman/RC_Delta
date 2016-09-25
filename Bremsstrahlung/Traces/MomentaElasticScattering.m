BeginPackage["Bremsstrahlung`Traces`MomentaElasticScattering`",
	{
		"HighEnergyPhysics`FeynCalc`",
		"Bremsstrahlung`Traces`Particles`"
	}
];


p1el::usage = "p1el, initial electron momentum";
p2el::usage = "p2el, initial proton momentum";
p3el::usage = "p3el, final electron momentum";
p4el::usage = "p4el, final proton momentum";
Kel::usage  = "Kel = p1el + p3el";
Pel::usage  = "Pel = p2el + p4el";
qel::usage  = "qel = p1el - p3el = p4el - p2el, momentum transfer";


Begin["Private`"];


SP = HighEnergyPhysics`FeynCalc`SP`SP;


Format[qel, TraditionalForm] = Subscript["\[ScriptQ]", "el"];
Format[Kel, TraditionalForm] = Subscript["\[ScriptCapitalK]", "el"];
Format[Pel, TraditionalForm] = Subscript["\[ScriptCapitalP]", "el"];


p1el = Kel/2 + qel/2;
p3el = Kel/2 - qel/2;
p2el = Pel/2 - qel/2;
p4el = Pel/2 + qel/2;


SP[Kel, Kel] = 4 * me^2 - SP[qel, qel];
SP[qel, Kel] = 0;
SP[Pel, Pel] = 4 * Mp^2 - SP[qel, qel];
SP[qel, Pel] = 0;


End[];


EndPackage[];
