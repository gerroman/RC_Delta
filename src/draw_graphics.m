(* ::Package:: *)

file = OpenRead["./output_1.txt"]
data =ReadList[file,{Real,Real,Real,Real}];
Show[
	ListLogPlot[data[[All,{1,2}]],Joined->True,PlotRange->{{0.9,2.4},{10^-7,0.01}}],
	ListLogPlot[Abs[data[[All,{1,3}]]],Joined->True],
	ListLogPlot[Abs[data[[All,{1,4}]]],Joined->True]
]


file2 = OpenRead["./output_2.txt"]
data2 =ReadList[file2,{Real,Real,Real,Real}];
Show[
	ListLogPlot[data2[[All,{1,2}]],Joined->True,PlotRange->{{0.9,2.4},{10^-7,0.01}}],
	ListLogPlot[Abs[data2[[All,{1,3}]]],Joined->True],
	ListLogPlot[Abs[data2[[All,{1,4}]]],Joined->True]
]



