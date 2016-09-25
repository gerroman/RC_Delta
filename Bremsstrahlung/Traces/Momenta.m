BeginPackage["Bremsstrahlung`Traces`Momenta`",
	{
		"HighEnergyPhysics`FeynCalc`",
		"Bremsstrahlung`Traces`Particles`"
	}
];


p1::usage = "p1, initial electron momentum";
p2::usage = "p2, initial proton momentum";
p3::usage = "p3, final electron momentum";
p4::usage = "p4, final proton momentum";
k::usage  = "k, final photon momentum";
qe::usage = "qe = p1 - p3, electron momentum transfer";
qp::usage = "qp = p4 - p2, proton momentum transfer";
K::usage  = "K = p1 + p3";
P::usage  = "P = p2 + p4";
r1::usage = "r1 = p3 + k";
r2::usage = "r2 = p1 - k";
t1::usage = "t1 = p4 + k";
t2::usage = "t2 = p2 - k";

xK::usage  = "xK = SP[k, K]";
xP::usage  = "xP = SP[k, P]";
xq::usage  = "xq = SP[k, qe] = SP[k, qp]";
qe2::usage = "qe2 = SP[qe, qe]";
qp2::usage = "qp2 = SP[qp, qp]";
Remove[Global`nu];
nu::usage  = "nu = SP[K, P]";

decompositionPattern::usage = "
	decompositionPattern[i]
"


ruleQe2toQp2::usage = "
	./ruleQe2toQp2
";
ruleQp2toQe2::usage = "
	./ruleQp2toQe2
";


Begin["Private`"];


SP = HighEnergyPhysics`FeynCalc`SP`SP;


Format[K,  TraditionalForm]  = "\[ScriptCapitalK]";
Format[P,  TraditionalForm]  = "\[ScriptCapitalP]";
Format[qe, TraditionalForm]  = Subscript["\[ScriptQ]", "e"];
Format[qp, TraditionalForm]  = Subscript["\[ScriptQ]", "p"];
Format[r1, TraditionalForm]  = Subscript["\[ScriptR]", 1];
Format[r2, TraditionalForm]  = Subscript["\[ScriptR]", 2];
Format[t1, TraditionalForm]  = Subscript["\[ScriptT]", 1];
Format[t2, TraditionalForm]  = Subscript["\[ScriptT]", 2];
Format[xK, TraditionalForm]  = Subscript["x", "\[ScriptCapitalK]"];
Format[xP, TraditionalForm]  = Subscript["x", "\[ScriptCapitalP]"];
Format[xq, TraditionalForm]  = Subscript["x", "\[ScriptQ]"];
Format[qe2, TraditionalForm] = Subsuperscript["\[ScriptQ]", "e", 2];
Format[qp2, TraditionalForm] = Subsuperscript["\[ScriptQ]", "p", 2];
Format[nu, TraditionalForm]  = "\[Nu]";



p1 = K/2 + qe/2;
p3 = K/2 - qe/2;
p2 = P/2 - qp/2;
p4 = P/2 + qp/2;
r1 = p3 + k;
r2 = p1 - k;
t1 = p4 + k;
t2 = p2 - k;


SP[k, k]  = 0;
SP[k, K]  = xK;
SP[k, P]  = xP;
SP[k, qe] = xq;
SP[k, qp] = xq;
SP[qe, qe] = qe2;
SP[qp, qp] = qp2;
SP[K, K]  = 4 * me^2 - qe2;
SP[P, P]  = 4 * Mp^2 - qp2;
SP[K, P]  = nu;
SP[qe, K] = 0;
SP[qp, P] = 0;
(* qe = qp + k *)
SP[qp, K] = - xK;
SP[qe, P] = + xP;
SP[qe, qp] = qe2/2 + qp2/2;


decompositionPattern[n_] := Alternatives[
	(x_. xP^a_. xq^b_. xK^c_. /; a + b + c == n && FreeQ[x, xP] && FreeQ[x, xq] && FreeQ[x, xK]),
	(x_. xP^a_. xq^b_. /; a + b == n && FreeQ[x, xP] && FreeQ[x, xq] && FreeQ[x, xK]),
	(x_. xK^a_. xq^b_. /; a + b == n && FreeQ[x, xP] && FreeQ[x, xq] && FreeQ[x, xK]),
	(x_. xK^a_. xP^b_. /; a + b == n && FreeQ[x, xP] && FreeQ[x, xq] && FreeQ[x, xK]),
	(x_. xP^a_. /; a == n && FreeQ[x, xP] && FreeQ[x, xq] && FreeQ[x, xK]),
	(x_. xK^a_. /; a == n && FreeQ[x, xP] && FreeQ[x, xq] && FreeQ[x, xK]),
	(x_. xq^a_. /; a == n && FreeQ[x, xP] && FreeQ[x, xq] && FreeQ[x, xK])
];


ruleQe2toQp2 = {
	qe2 -> qp2 + 2 xq
};

ruleQp2toQe2 = {
	qp2 -> qe2 - 2 xq
};


End[];


EndPackage[];
