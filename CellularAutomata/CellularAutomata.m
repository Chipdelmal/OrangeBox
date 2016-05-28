(* Mathematica Package *)

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
SeparateAutomaton::usage = "Separates the different values of a CA matrix creating different separate automata"
MergeAutomata::usage = "Merges the values of CA matrix into a full automata one"
EvaluateCellularAutomaton::usage = "Evaluates one epoch of a multilayered CA in which each layer has its own function and after those functions are evluated the system is evaluated by another function as a whole"
FilterVariablesWithRanges::usage = "Filters a list of samples with respect with a restriction pattern that allows ranges of values"
NestListEvaluateCellularAutomaton::usage = "Evaluates a CA iteratively returning every epoch's value"
NestEvaluateCellularAutomaton::usage = "Evaluates a CA iteratively returning only the last epoch's value"
ExtractKernelValuesAndPositions::usage = "Gets the values of the CA matrix filtered by a given kernel and returns them along with their positions"
ReplaceKernelValuesWithPositions::usage = "Replaces the values of a CA matrix after being evaluated so that the new CA contains the processed outputs of the kernel" 
GenerateProcessedCellsMatrix::usage = "Generates a Boolean matrix that will hold the cells state (processed or unprocessed) in a detemrinate epoch"
UpdateProcessedMatrix::usage = "Updates the value of a given cell of the processed matrix so that the cells that have been processed are updated"
RandomUnprocessedPosition::usage = "Select a random position with a TRUE value within a matrix of boolean values"
FixSize::usage = "Adds padding to a CA so that the kernel may be applied even if it goes beyond its scope"
EvaluateRandomCell::usage = "Evaluates a random cell (according to the values of a matrix that holds if a cell has been updated or not) and applies a given function to it"
EvaluateOneEpoch::usage = "Evaluates the whole CA through one epoch"
EvaluateCAThroughEpochs::usage = "Evaluates the CA trhough a finite number of epochs"

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
                                                                                            
                                                                                            
(*----------------CA Second Package--------------*)
SeparateAutomaton[automata_] := Table[#[[All, i]] & /@ automata, {i, 1, automata[[1]][[1]] // Length}]
EvaluateCellularAutomaton[automata_, neighborhoodType_, functions_, cellFunction_] := Module[{vonNewmanNeighborhood, mooreNeighborhood, automataSeparated, automataNeighbours, automataNeighboursReleased, appliedSeparateAutomata, mergedMatrix},
  	(*Define the type of neighbourhood to use*)
  	vonNewmanNeighborhood = {{{a_, b_, c_}, {d_, e_, f_}, {g_, h_, i_}} -> Hold[{b, d, e, f, h}]};
  	mooreNeighborhood = {{{a_, b_, c_}, {d_, e_, f_}, {g_, h_, i_}} -> Hold[{a, b, c, d, e, f, g, h, i}]};
  	(*Separate different automata*)
	automataSeparated = SeparateAutomaton[automata];
  	(*Obtain the neighbourhoods of the automata*)
  	automataNeighbours = CellularAutomaton[vonNewmanNeighborhood, #, 1] & /@ automataSeparated;
  	automataNeighboursReleased = (Map[ReleaseHold, automataNeighbours[[#]]] // Last) & /@ Range[automataNeighbours // Length];
  	(*Evaluate Selected Functions for each Separate Automata*)
  	appliedSeparateAutomata = Map[functions[[#]], automataNeighboursReleased[[#]], {2}] & /@ Range[automataNeighboursReleased // Length];
  	(*Merge the Separate Cellular Automata*)
  	mergedMatrix = MergeAutomata[appliedSeparateAutomata];
  	{Map[Apply[cellFunction, #] &, mergedMatrix, {2}],neighborhoodType,functions,cellFunction}
]
MergeAutomata[automata_]:=(automata[[All, #]] // Transpose) & /@Range[automata[[1]] // Length];
NestListEvaluateCellularAutomaton[automata_,neighborhoodType_,functions_,cellFunction_,epochs_]:=NestList[(Apply[EvaluateCellularAutomaton[#1, #2, #3, #4] &, #]) &, {automata,neighborhoodType,functions,cellFunction}, epochs]
NestEvaluateCellularAutomaton[automata_,neighborhoodType_,functions_,cellFunction_,epochs_]:=Nest[(Apply[EvaluateCellularAutomaton[#1, #2, #3, #4] &, #]) &, {automata,neighborhoodType,functions,cellFunction}, epochs]
ExtractKernelValuesAndPositions[automaton_, kernelCenter_, kernel_] :=Module[{tempR, tempC, kLength, sLength, rRange, cRange, rValues, cValues, listOfPossibleElements, submatrix, selectedWithKernel, kernelValues, kernelPositions},
  	tempR = kernelCenter[[1]];
  	tempC = kernelCenter[[2]];
  	kLength = kernel // Length; sLength = (kLength - 1)/2;
  	(*Ranges of values to extract from matrix*) 
  	rRange = {tempR - sLength, tempR + sLength};
  	cRange = {tempC - sLength, tempC + sLength};
  	(*Specific positions to extract from matrix*)  
  	rValues = Apply[Range, rRange];
 	cValues = Apply[Range, cRange];
  	listOfPossibleElements = Partition[Tuples[{rValues, cValues}], kLength];
  	(*Values of positions of matrix*)
  	submatrix = automaton[[rRange[[1]] ;; rRange[[2]], cRange[[1]] ;; cRange[[2]]]];
  	(*One to one reference of elements and positions*)
  	selectedWithKernel = Position[kernel, 1];
  	kernelValues = Map[Apply[submatrix[[##]] &, #] &, selectedWithKernel];
  	kernelPositions = Map[Apply[listOfPossibleElements[[##]] &, #] &, selectedWithKernel];
  	{kernelValues, kernelPositions}
]
ReplaceKernelValuesWithPositions[automaton_, kernelValuesAndPositions_] := Module[{replaceRules},
  	replaceRules = #[[2]] -> #[[1]] & /@ (kernelValuesAndPositions // Transpose);
  	replaceRules = DeleteCases[replaceRules, {0, _} -> _];
  	replaceRules = DeleteCases[replaceRules, {_, 0} -> _];
  	ReplacePart[automaton, replaceRules]
]
GenerateProcessedCellsMatrix[automaton_] := ConstantArray[False, {automaton // Length, automaton[[1]] // Length}]
UpdateProcessedMatrix[processedMatrix_, index_] := ReplacePart[processedMatrix, index -> True]
RandomUnprocessedPosition[processedMatrix_] := (Position[processedMatrix, False] // RandomSample)[[1]]
FixSize[automaton_, kernel_, paddingValue_] := ArrayPad[automaton,((kernel // Length) - 1)/2, paddingValue]
(*Under Review*)
EvaluateRandomCell[automaton_, processedMatrix_, kernel_, kernelFunction_, matrixFixValue_] := Module[{valuesAndPositions, evaluatedKernel, updateProcessedMatrix, updatedAutomata, fixedSizeMatrix, sKernel, randomCell},
  	fixedSizeMatrix = FixSize[automaton, kernel, matrixFixValue];
  	sKernel = ((kernel // Length) - 1)/2;
  	randomCell = RandomUnprocessedPosition[processedMatrix] + sKernel;
  	valuesAndPositions = ExtractKernelValuesAndPositions[fixedSizeMatrix, randomCell, kernel];
  	(*Apply Function*)
	evaluatedKernel = {Apply[kernelFunction, {valuesAndPositions[[1]]}],valuesAndPositions[[2]] - sKernel};
  	(*Apply Function*)
	updateProcessedMatrix = UpdateProcessedMatrix[processedMatrix, randomCell - sKernel];
	updatedAutomata = ReplaceKernelValuesWithPositions[automaton, evaluatedKernel];
	(*Print[{randomCell-sKernel,updatedAutomata//MatrixForm}];*)
  	{updatedAutomata, updateProcessedMatrix, kernel, kernelFunction, matrixFixValue}
]
EvaluateOneEpoch[automaton_, kernel_, kernelFunction_, matrixFixValue_] := Module[{processedMatrix,output},
  	processedMatrix = GenerateProcessedCellsMatrix[automaton];
  	output = Nest[(Apply[EvaluateRandomCell[#1, #2, #3, #4, #5] &, #]) &, {automaton, processedMatrix, kernel, kernelFunction, matrixFixValue}, (automaton // Length)*(automaton[[1]] // Length)];
	Delete[output,2]
]
EvaluateCAThroughEpochs[initialAutomaton_, kernel_, kernelFunction_, padValue_, epochsNumber_, historicBool_] :=
 	If[historicBool,
 		NestList[(Apply[EvaluateOneEpoch[#1, #2, #3, #4] &, #]) &, {initialAutomaton, kernel, kernelFunction, padValue}, epochsNumber][[All, 1]],
  		Nest[(Apply[EvaluateOneEpoch[#1, #2, #3, #4] &, #]) &, {initialAutomaton, kernel, kernelFunction, padValue}, epochsNumber][[1]]
]                                                                                            

End[]

EndPackage[]

