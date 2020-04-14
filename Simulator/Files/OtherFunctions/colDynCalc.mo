within Simulator.Files.OtherFunctions;

function colDynCalc
  extends Modelica.Icons.Function;
  input Integer noOfStages, dynamic;
  output Real d[noOfStages];
algorithm
if dynamic == 1 then
  d := fill(1, noOfStages);
else 
  d := fill(0, noOfStages);
end if;
end colDynCalc;
