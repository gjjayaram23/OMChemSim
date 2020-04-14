within Simulator.UnitOperations;

model Cooler "Model of a cooler to heat a material stream"
  extends Simulator.Files.Icons.Cooler; 
  
    parameter Simulator.Files.ChemsepDatabase.GeneralProperties C[Nc] "Component instances array" annotation(
    Dialog(tab = "Cooler Specifications", group = "Component Parameters"));
    parameter Integer Nc "number of components" annotation(
    Dialog(tab = "Cooler Specifications", group = "Component Parameters"));
  //====================================================================================
  Real Fin(unit = "mol/s", min = 0, start = Fg) "Inlet stream molar flow rate";
  Real Pin(unit = "Pa", min = 0, start =Pg) "Inlet stream pressure";
  Real Tin(unit = "K", min = 0, start = Tg) "Inlet stream temperature";
  Real Hin(unit = "kJ/kmol") "Inlet stream molar enthalpy";
  Real Sin(unit = "kJ/[kmol.K]") "Inlet stream molar entropy";
  Real xvapin(unit = "-", min = 0, max = 1, start = xvapg) "Inlet stream vapor phase mole fraction";
  
  Real Q(unit = "W") "Heat removed";
  Real Tdel(unit = "K") "Temperature drop";
   
  Real Fout(unit = "mol/s", min = 0, start = Fg) "Outlet stream molar flow rate";
  Real Pout(unit = "Pa", min = 0, start = Pg) "Outlet stream pressure";
  Real Tout(unit = "K", min = 0, start = Tg) "Outlet stream temperature";
  Real xvapout(unit = "-", min = 0, max = 1, start = xvapg) "Outlet stream vapor phase mole fraction";
  Real Hout(unit = "kJ/kmol") "Outlet stream molar enthalpy";
  Real Sout(unit = "kJ/[kmol.K]") "Outlet stream molar entropy"; 
  Real x_c[Nc](each unit = "-", each min = 0, each max = 1, start=xg) "Component mole fraction";
  Real MT(unit = "mol", start = 1000) "Total Molar Holdup";
  Real ML(unit = "mol") "Liquid Holdup", MV(unit = "mol") "Vapor Holdup";
  Real VG(unit = "m^3") "Vapor Volume", VL(unit = "m^3") "Liquid Volume";
  Real PG(unit = "Pa") "Gas Pressure", h(unit = "m") "Level inside heater", PT(unit = "Pa") "Total Pressure";
  Real rholiq_c[Nc](each unit = "mol/m^3") "Component Density", rholiq(unit = "mol/m^3") "Mixture Density";
  
  //========================================================================================
  parameter Real Pdel(unit = "Pa") "Pressure drop" annotation(
    Dialog(tab = "Cooler Specifications", group = "Calculation Parameters"));
  parameter Real Eff(unit = "-") "Efficiency" annotation(
    Dialog(tab = "Cooler Specifications", group = "Calculation Parameters"));
  parameter Boolean Dynamics "Type of Operation" annotation(
    Dialog(tab = "Physical Specifications"));
  parameter Real Cd(unit = "-") = 0.5 "Discharge Coefficient of Valve" annotation(
    Dialog(tab = "Physical Specifications", group = "Vessel Geometry"));
  parameter Real VT(unit = "m^3") = 5 "Volume of Vessel" annotation(
    Dialog(tab = "Physical Specifications", group = "Vessel Geometry"));
  parameter Real A(unit = "m^2") = 2 "Area of Vessel" annotation(
    Dialog(tab = "Physical Specifications", group = "Vessel Geometry"));
  parameter Real xtl(fixed = false), xtv(fixed = false);
  parameter Real Pset(unit = "Pa") = 202650 "Initial Pressure"annotation(
    Dialog(tab = "Physical Specifications", group = "Initial Conditions"));
  parameter Real hset(unit = "Pa") = 0.5 "Initial Level" annotation(
    Dialog(tab = "Physical Specifications", group = "Initial Conditions"));
//========================================================================================
  Simulator.Files.Interfaces.matConn In(Nc = Nc) annotation(
    Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Simulator.Files.Interfaces.matConn Out(Nc = Nc) annotation(
    Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Simulator.Files.Interfaces.enConn En annotation(
    Placement(visible = true, transformation(origin = {0, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, -100}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  //========================================================================================
  extends GuessModels.InitialGuess;
  
  initial equation 
  PG = Pset;
  h = hset;
  if Dynamics == true then
  der(MT) = 0;
  der(Hout) = 0;
  end if;
  equation
//connector equations
  In.P = Pin;
  In.T = Tin;
  In.F = Fin;
  In.H = Hin;
  In.S = Sin;
  In.x_pc[1, :] = x_c[:];
  In.xvap = xvapin;
  Out.P = Pout;
  Out.T = Tout;
  Out.F = Fout;
  Out.H = Hout;
  Out.S = Sout;
  Out.x_pc[1, :] = x_c[:];
  Out.xvap = xvapout;
  En.Q = Q;
//=============================================================================================
//Mole balance
  Fin - Fout = if Dynamics == true then der(MT) else 0;
  MT = ML + MV ;
  ML = rholiq * VL ;
  PG * VG = MV * 8.314 * Tout;
  VT = VL + VG;
  VL = A * h;
  PT = PG + rholiq * 9.81 * h;
  Fout * (1 - xvapout) = xtl * Cd * sqrt(PT - 1e5) ;
  Fout * xvapout = xtv * Cd * sqrt(PG - 1e5);
  for i in 1:Nc loop
  rholiq_c[i] = Simulator.Files.ThermodynamicFunctions.Dens(C[i].LiqDen, C[i].Tc, Tout, Pout);
 end for;    
  rholiq = sum(x_c[:] .* rholiq_c) ;
//=============================================================================================
//Energy balance
  Hin + Eff * Q / Fin - Hout = if Dynamics == true then der(MT*Hout) else 0;
//Pressure calculation
  Pin - Pdel = Pout;
//=============================================================================================
//Temperature calculation
  Tin + Tdel = Tout;
  annotation(
    Documentation(info = "<html><head></head><body><span style=\"font-size: 12px;\">The <b>Cooler</b> is used to simulate the cooling process of a material stream.</span><div><br></div><div><div><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13px; orphans: 2; widows: 2;\">The cooler model have following connection ports:</span></div><div><div style=\"orphans: 2; widows: 2;\"><ol><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Two Material Streams:</span></font></li><ul><li><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13px;\">feed stream</span></li><li><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13px;\">outlet stream</span></li></ul><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">One Energy Stream:</span></font></li><ul><li><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13px;\">heat removed</span></li></ul></ol></div></div><div style=\"font-size: 12px;\"><br></div><div style=\"font-size: 12px;\">Following calculation parameters must be provided to the cooler:</div><div style=\"font-size: 12px;\"><ol><li>Pressure Drop (<b>Pdel</b>)</li><li>Efficiency (<b>Eff</b>)</li></ol><div><div>The above variables have been declared of type&nbsp;<i>parameter Real.&nbsp;</i></div><div>During simulation, their values can specified directly under&nbsp;<b>Cooler Specifications</b>&nbsp;by double clicking on the cooler model instance.</div></div><div><br></div><div><br></div><div>In addition to the above parameters, any one additional variable from the below list must be provided for the model to simulate successfully:</div><div><ol><li>Outlet Temperature (<b>Tout</b>)</li><li>Temperature Drop (<b>Tdel</b>)</li><li>Heat Removed (<b>Q</b>)</li><li>Outlet Stream Vapor Phase Mole Fraction (<b>xvapout</b>)</li></ol><div><div>These variables are declared of type&nbsp;<i>Real.</i></div><div>During simulation, value of one of these variables need to be defined in the equation section.</div></div><div><br></div><div><br></div></div><div>For detailed explaination on how to use this model to simulate a Cooler, go to <a href=\"modelica://Simulator.Examples.Cooler\">Cooler Example</a></div></div></div></body></html>"));
    
    end Cooler;
