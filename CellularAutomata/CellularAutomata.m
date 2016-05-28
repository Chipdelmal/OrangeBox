(* Mathematica Package *)
(* Hector Manuel Sanchez Castellanos: sanchez.hmsc@itesm.mx *)
(* Created by the Wolfram Workbench Nov 22, 2013 *)

BeginPackage["CellularAutomata`"]
ListInitialAvailableIndexes::usage
KernelCenterIndex::usage
PadCellularAutomaton::usage
ReplaceAutomatonWithKernel::usage
GetAutomatonValues::usage
ReplacePaddedAutomatonWithKernel::usage
GetPaddedAutomatonValues::usage
UpdateAutomatonStep::usage
UpdateAutomatonEpoch::usage
UpdateAutomataStep::usage
UpdateAutomataEpoch::usage
RunAutomataThroughEpochs::usage
RunAutomataThroughEpochsList::usage

Begin["`Private`"]
ListInitialAvailableIndexes[cellularAutomaton_] := Tuples[Range[cellularAutomaton // Length], 2]
KernelCenterIndex[kernelMatrix_] := Module[{i},
  	i = Floor[(kernelMatrix // Length)/2] + 1;
  	{i, i}
]
PadCellularAutomaton[cellularAutomaton_, padValue_, kernelMatrix_] := Module[{n},
  	n = Ceiling[(kernelMatrix // Length)/2] - 1;
  	ArrayPad[cellularAutomaton, n, padValue]
]
ReplaceAutomatonWithKernel[kernelCenter_, kernelValuesMatrix_, cellularAutomaton_] := Module[{kvmL, kvmLR, rawIndexes, values, replacementRawIndexesAndValues, movedReplacements, rangeL},
  	kvmL = kernelValuesMatrix // Length;
  	rangeL = Floor[kvmL/2];
  	kvmLR = Range[-rangeL, rangeL];
  	rawIndexes = Tuples[kvmLR, 2];
  	values = kernelValuesMatrix // Flatten;
  	replacementRawIndexesAndValues = Transpose[{rawIndexes, values}];
  	movedReplacements = ({#[[1]][[1]] + kernelCenter[[1]], #[[1]][[2]] + kernelCenter[[2]]} -> #[[2]]) & /@ replacementRawIndexesAndValues;
  	ReplacePart[cellularAutomaton, movedReplacements]
]
GetAutomatonValues[kernelCenter_, kernelMatrix_, cellularAutomaton_] := Module[{kvmL, rangeL, kvmLR, rawIndexes, parts},
  	kvmL = kernelMatrix // Length;
  	rangeL = Floor[kvmL/2];
  	kvmLR = Range[-rangeL, rangeL];
  	rawIndexes = Tuples[kvmLR, 2];
  	parts = {#[[1]] + kernelCenter[[1]], #[[2]] + kernelCenter[[2]]} & /@ rawIndexes;
  	Partition[cellularAutomaton[[#[[1]], #[[2]]]] & /@ parts, kvmL]
]
ReplacePaddedAutomatonWithKernel[kernelCenter_, kernelValuesMatrix_, cellularAutomaton_, padValue_] := Module[{kvmL, kvmLR, rawIndexes, values, replacementRawIndexesAndValues, movedReplacements, rangeL, padN, padCA},
  	padN = Ceiling[(kernelValuesMatrix // Length)/2] - 1;
  	(*Requires a Squared Kernel*)
  
  	kvmL = kernelValuesMatrix // Length;
  	rangeL = Floor[kvmL/2];
  	kvmLR = Range[-rangeL, rangeL];
  	rawIndexes = Tuples[kvmLR, 2];
  	values = kernelValuesMatrix // Flatten;
  
  	replacementRawIndexesAndValues = {rawIndexes, values} // Transpose;
  	movedReplacements = ({#[[1]][[1]] + kernelCenter[[1]] + padN, #[[1]][[2]] + kernelCenter[[2]] + padN} -> #[[2]]) & /@ replacementRawIndexesAndValues;
  	padCA = PadCellularAutomaton[cellularAutomaton, padValue, kernelValuesMatrix];
  
  	ReplacePart[padCA, movedReplacements][[padN + 1 ;; (padCA // Length) - padN, padN + 1 ;; (padCA // Length) - padN]]
]
GetPaddedAutomatonValues[kernelCenter_, kernelMatrix_, cellularAutomaton_, padValue_] := Module[{kvmL, rangeL, kvmLR, rawIndexes, parts, padN, paddedCellularAutomaton},
  	padN = Ceiling[(kernelMatrix // Length)/2] - 1;
  
  	kvmL = kernelMatrix // Length;
  	rangeL = Floor[kvmL/2];
  	kvmLR = Range[-rangeL, rangeL];
  	rawIndexes = Tuples[kvmLR, 2];
  	parts = {#[[1]] + kernelCenter[[1]] + padN, #[[2]] + kernelCenter[[2]] + padN} & /@ rawIndexes;
  	paddedCellularAutomaton = PadCellularAutomaton[cellularAutomaton, padValue, kernelMatrix];
  
  	Partition[paddedCellularAutomaton[[#[[1]], #[[2]]]] & /@ parts, kvmL]
]
UpdateAutomatonStep[kernelMatrix_, cellularAutomaton_, padValue_, updateFunction_, remainingIndexes_] := Module[{kernelValues, kernelCenter, updatedKernelValues},
  	kernelCenter = RandomChoice[remainingIndexes];
  	kernelValues = GetPaddedAutomatonValues[kernelCenter, kernelMatrix, cellularAutomaton, padValue];
  	updatedKernelValues = Apply[updateFunction, {kernelValues}];(*Check this step*)
  	{
  	 	kernelMatrix,
   		ReplacePaddedAutomatonWithKernel[kernelCenter, updatedKernelValues,cellularAutomaton, padValue],
   		padValue,
   		updateFunction,
   		remainingIndexes // DeleteCases[#, kernelCenter] &
   }
]
UpdateAutomatonEpoch[kernelMatrix_, cellularAutomaton_, padValue_, updateFunction_] := Module[{indexesPool},
  indexesPool = ListInitialAvailableIndexes[cellularAutomaton];
  Nest[Apply[UpdateAutomatonStep, #] &, {kernelMatrix, cellularAutomaton, padValue, updateFunction, indexesPool}, indexesPool // Length]
]
UpdateAutomataStep[kernelMatrix_, cellularAutomata_, padValue_, updateFunction_, remainingIndexes_] := Module[{kernelValues, kernelCenter, updatedKernelValues},
  	kernelCenter = RandomChoice[remainingIndexes];
  	kernelValues = GetPaddedAutomatonValues[kernelCenter, kernelMatrix, #, padValue] & /@ cellularAutomata;
  	updatedKernelValues = Apply[updateFunction, {kernelValues}];
  	{
   		kernelMatrix,
   		ReplacePaddedAutomatonWithKernel[kernelCenter, #[[1]], #[[2]], padValue] & /@ ({updatedKernelValues, cellularAutomata} // Transpose),
   		padValue,
   		updateFunction,
   		remainingIndexes // DeleteCases[#, kernelCenter] &
  	}
]
UpdateAutomataEpoch[kernelMatrix_, cellularAutomata_, padValue_, updateFunction_] := Module[{indexesPool},
  	indexesPool = ListInitialAvailableIndexes[cellularAutomata[[1]]];
  	Nest[Apply[UpdateAutomataStep, #] &, {kernelMatrix, cellularAutomata, padValue, updateFunction, indexesPool}, indexesPool // Length][[1 ;; 4]]
]
RunAutomataThroughEpochs[kernelMatrix_, cellularAutomata_, padValue_, updateFunction_, epochsNumber_] := Module[{},
  	Nest[Apply[UpdateAutomataEpoch, #] &, {kernelMatrix, cellularAutomata, padValue, updateFunction}, epochsNumber]
]
RunAutomataThroughEpochsList[kernelMatrix_, cellularAutomata_, padValue_, updateFunction_, epochsNumber_] := Module[{},
  	NestList[Apply[UpdateAutomataEpoch, #] &, {kernelMatrix, cellularAutomata, padValue, updateFunction}, epochsNumber]
]

End[]

EndPackage[]

