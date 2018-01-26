Mesh.RemeshAlgorithm=1;
Mesh.CharacteristicLengthFactor=0.05;
Merge "/home/bkj/sintering/cubeWithHole.msh";
CreateTopology;

Compound Surface(23)={8};
Compound Surface(24)={7};
Compound Surface(25)={6};
Compound Surface(26)={5};
Compound Surface(27)={4};
Compound Surface(28)={3};
Compound Surface(29)={2};

Surface Loop(30) = {23,24,25,26,27,28,29};
Volume(31) = {30};

Physical Surface(32) = {23,24,25,26,27,28,29};
Physical Volume(33) = {31};

DefineConstant[
  hide = {Geometry.HideCompounds, Choices{0,1},
    AutoCheck 0, GmshOption "Geometry.HideCompounds",
    Name "Parameters/Hide compound sub-entities"}
];
//+
Recombine Surface {29};
//+
Recombine Surface {29, 26, 25, 23, 24, 27, 28};
//+
Coherence;
//+
Coherence;
//+
Coherence;
//+
Coherence;
//+
Coherence;
//+
Coherence;
//+
Coherence;
