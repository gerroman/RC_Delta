BeginPackage["Bremsstrahlung`Traces`Particles`"];


me::usage = "me, electron mass, GeV";
Mp::usage = "Mp, proton mass, GeV";
Md::usage = "Md, Delta(1232) mass, GeV";
Gd::usage = "Gd, Delta(1232) full width, GeV";


Begin["Private`"];


Format[me, TraditionalForm] = Subscript["m", "e"];
Format[Mp, TraditionalForm] = Subscript["M", "p"];
Format[Md, TraditionalForm] = Subscript["M", "\[CapitalDelta]"];
Format[Gd, TraditionalForm] = Subscript["\[CapitalGamma]", "\[CapitalDelta]"];


End[];


EndPackage[];
