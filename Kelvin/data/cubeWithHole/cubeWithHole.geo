Merge "/home/bkj/sintering/cubeWithHole.msh";
CreateTopology;

Compound Surface(23)={8};
Compound Surface(24)={7};
Compound Surface(25)={6};
Compound Surface(26)={5};
Compound Surface(27)={4};
Compound Surface(28)={3};
Compound Surface(29)={2};

Mesh.RemeshAlgorithm=1;
Mesh.CharacteristicLengthFactor=0.05;
Mesh.RemeshParametrization = 7; // conformal finite element
Mesh.Algorithm = 6; // Frontal
