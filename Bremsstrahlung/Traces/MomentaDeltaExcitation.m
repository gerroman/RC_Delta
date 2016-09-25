BeginPackage["Bremsstrahlung`Traces`MomentaDeltaExcitation`",
	{
		"HighEnergyPhysics`FeynCalc`",
		"Bremsstrahlung`Traces`Particles`"
	}
];


p1d::usage = "p1d, initial electron momentum";
p2d::usage = "p2d, initial proton momentum";
p3d::usage = "p3d, final electron momentum";
p4d::usage = "p4d, final Delta(1232) momentum in e p -> e Delta";
qd::usage  = "qd = qed = qpd momentum transfer in e p -> e Delta";
Kd::usage  = "Kd = p1d + p3d";
Pd::usage  = "Pd = p2d + p4d";
tildePd::usage = "tildePd = Pd - (qd * Pd)/(qd * qd) qd"


ruleTildePd::usage = "
	./ruleTildePd conver Pd to Tilde[Pd]
"


Begin["Private`"];


SP = HighEnergyPhysics`FeynCalc`SP`SP;
FV = HighEnergyPhysics`FeynCalc`FV`FV;
SPE = HighEnergyPhysics`FeynCalc`ScalarProductExpand`ScalarProductExpand;
FCE = HighEnergyPhysics`fctools`FeynCalcExternal`FeynCalcExternal;


Format[qd, TraditionalForm] = Subscript["\[ScriptQ]","\[CapitalDelta]"];
Format[Kd, TraditionalForm] = Subscript["\[ScriptCapitalK]", "\[CapitalDelta]"];
Format[Pd, TraditionalForm] = Subscript["\[ScriptCapitalP]", "\[CapitalDelta]"];
Format[tildePd, TraditionalForm] = Subscript[OverTilde["\[ScriptCapitalP]"], "\[CapitalDelta]"];


p1d = Kd/2 + qd/2;
p3d = Kd/2 - qd/2;
p2d = Pd/2 - qd/2;
p4d = Pd/2 + qd/2;


SP[Kd, Kd] = 4 * me^2 - SP[qd, qd];
SP[qd, Kd] = 0;
SP[Pd, Pd] = 2 * Mp^2 + 2 * Md^2 - SP[qd, qd];
SP[qd, Pd] = Md^2 - Mp^2;
SP[tildePd, qd] = 0;
SP[tildePd, tildePd] := Simplify@Factor@FCE@SPE@(
	SP[Pd,Pd] 
	- 2 SP[Pd,qd]/SP[qd,qd] SP[Pd, qd]
	+ (SP[Pd,qd]/SP[qd,qd])^2 SP[qd, qd]
)

ruleTildePd = {
	FV[Pd, mu_] :> FV[tildePd, mu] + SP[qd, Pd] / SP[qd, qd] * FV[qd, mu]
}


End[];


EndPackage[];
