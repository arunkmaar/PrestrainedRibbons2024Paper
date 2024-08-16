(* ::Package:: *)

(*
   A parser for AUTO's solution files
   by Basile Audoly (baudoly@gmail.com), September 2012 
   based on reverse engineering of Solution::parse() and Solution::read()
      found in auto/07p/plaut04/src/Solution.c++
*)

(*
 Structure of AUTO's solution file:
    Snapshot1
    Snapshot2
    ...

Where each snapshot reads
    Header
    Data
    [UnusedData (optional)]
    Parameters

The *Header* is made up of at least 12 integer digits:
 1 - ibr    + branch ID (kBranchID)
 2 - ntot   - ? used to assign a label related to stability
 3 - itp    - ?
 4 - lab    + solution label (kSolLabel)
 5 - nfpr   + number of extra data given in UnusedData section
 6 - isw    - unused
 7 - ntpl   + number of vertices (kNVert)
 8 - nar    + data per vertex (kNDataPerVertex)
 9 - nrowpr + line count for Data+[UnusedData]+Parameters (kNLinesAfterHeader)
10 - ntst   - tested for !=0, which is always true as far as we tested
11 - ncol   - unused
12 - npar1  + number of parameters (kNPar)
[...]

The *Data* is made up of nVertices records of DataForVertex

The *DataForVertex* is made up of nDataPerVertex floating numbers
ending up with a newline, with a maximum of 7 records per line

The presence of *UnusedData* can be determined by comparing nrowpr
to nrowprsmall, which is the line count of (Data+Parameters) based
on the values of nar, ntpl and npar1

The optional *UnusedData* has the following structure:
  kNVert records of (DataPerVertex-1) floating numbers terminated
     by a newline, each record having a maximum of 7 numbers per line
  nfpr integers terminated by a newline, with a maximum of 20 ints
     per line
  nfpr float

The *Parameters* section lists the npar1 parameters with a maximum of
    7 per line, special values being par[1]=mass and par[10]=period
*)


(* 7 floating numbers per line is a universal constant that is hardcoded into AUTO *)
nLinesForData[nItems_Integer,nItemsPerLine_:7]:=Ceiling[nItems/nItemsPerLine];

ReadHeaderSection[l_List] := {
 kBranchID -> l[[1]], 
 kSolLabel -> l[[4]],
 kNFPR -> l[[5]],
 kNVert -> l[[7]],
 kNDataPerVertex -> l[[8]],
 kNLinesAfterHeader -> l[[9]],
 kNPar -> l[[12]]
};

ParseAUTOSnapshot[fort8l_List, startLine_Integer] := 
Module[{
    head, startLineForNextSnapshot, lcData(* line count of data section *),
    data,lcParams,lcUnused,params}, Assert[startLine <= Length[fort8l]];
  head = ReadHeaderSection[fort8l[[startLine]]];
  startLineForNextSnapshot = startLine + 1(*header*) + kNLinesAfterHeader /.head;
  lcData  = {kNVert,kNDataPerVertex}/.head // #[[1]] * nLinesForData[#[[2]]]&;
  data = Partition[Flatten[fort8l[[startLine+1;;startLine +lcData]]],kNDataPerVertex/.head];
  Assert[Dimensions[data]===({kNVert,kNDataPerVertex}/.head)];
  lcParams = nLinesForData[kNPar/.head];
  lcUnused = If[ (kNLinesAfterHeader/.head) > lcData+lcParams,
     (* UnusedData section detected, skip lines *)
     (kNVert/.head)*nLinesForData[kNDataPerVertex-1/.head] + nLinesForData[kNFPR/.head] + nLinesForData[kNFPR/.head,20],
     0
  ];
  params=Flatten[fort8l[[startLine+lcData+lcUnused+1;;startLine+lcData+lcUnused+lcParams]]];
  {
    startLineForNextSnapshot,
    Join[head, {kData -> data,kContinuationParams -> params}]
   }
  ];

(*
Parse the entire fort .8 file
*)
ParseAUTOSolutionFileAsTable[l_List] := Module[{cl, ret, t, \[Lambda]},
   cl = 1;(* current line *)
   ret = {};
   While[
    cl <= Length[l],
    t = ParseAUTOSnapshot[l, cl];
    cl = t[[1]];
    ret = Append[ret, t[[2]]];
    ];
   ret
   ];

(*
   loading
*)
reloadSolutionFile[solutionFilePath_String] := Module[{autoSols},
   autoSols = 
    ReadList[solutionFilePath, Number, RecordLists -> True] // ParseAUTOSolutionFileAsTable;
   Print["** AUTO solution file has been parsed, with ", autoSols // Length, " snapshots found"];
   autoSols
];
