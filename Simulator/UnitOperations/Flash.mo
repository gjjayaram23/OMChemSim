within Simulator.UnitOperations;

model Flash "Model of a flash column to separate vapor and liquid phases from a mixed phase material stream"
  //==============================================================================
  //Header Files and Parameters
  extends Simulator.Files.Icons.Flash;
  import Simulator.Files.*;
  parameter ChemsepDatabase.GeneralProperties C[Nc] "Component instances array" annotation(
    Dialog(tab = "Flash Specifications", group = "Component Parameters"));
  parameter Integer Nc "Number of components" annotation(
    Dialog(tab = "Flash Specifications", group = "Component Parameters"));
  parameter Boolean BTdef = false "True if flash is operated at temperature other than feed temp else false" annotation(
    Dialog(tab = "Flash Specifications", group = "Calculation Parameters"));
  parameter Boolean BPdef = false "True if flash is operated at pressure other than feed pressure else false" annotation(
    Dialog(tab = "Flash Specifications", group = "Calculation Parameters"));
  parameter Boolean Dynamics "True if flash is operated in Unsteady mode" annotation(
    Dialog(tab = "Physical Specifications"));
  parameter Real Tdef(unit = "K") = 298.15 "Separation temperature if BTdef is true" annotation(
    Dialog(tab = "Flash Specifications", group = "Calculation Parameters"));
  parameter Real Pdef(unit = "Pa") = 101325 "Separation pressure if BPdef is true" annotation(
    Dialog(tab = "Flash Specifications", group = "Calculation Parameters"));
  parameter Real A(unit = "m^2") = 4 "Wetted area of the flash tank" annotation(
    Dialog(tab = "Physical Specifications", group = "Calculation Parameters"));
  parameter Real VT(unit = "m^3") = 8 "Total Volume of the flash tank" annotation(
    Dialog(tab = "Physical Specifications", group = "Calculation Parameters"));
  parameter Real Cd = 0.4 "Discharge coefficient of valve"annotation(
    Dialog(tab = "Physical Specifications", group = "Calculation Parameters"));
  parameter Real xtl(min = 0, fixed = false) "Stem Position of the Valve", xtv(min = 0, fixed = false) "Stem position of the valve";
  parameter Real Pset(unit = "Pa") = 101325 "Initial Column Pressure" annotation(
    Dialog(tab = "Physical Specifications", group = "Initial Condition"));
  parameter Real hset(unit = "m") = 1 "Initial Column level" annotation(
    Dialog(tab = "Physical Specifications", group = "Initial Condition"));
  
  //==============================================================================
//Model Variables
  Real T(unit = "K", min = 0, start = Tg) "Flash column temperature";
  Real Pin(unit = "Pa", min = 0, start = Pg) "Flash column pressure";
  Real F_p[3](each unit = "mol/s", each min = 0, each start = Fg)"Feed stream mole flow";
  Real x_pc[3, Nc](each unit = "-", each min = 0, each max = 1, start = {xg,yg,xguess}) "Component mole fraction";
  Real Hin(unit = "kJ/kmol") "Molar enthalpy in phase";
  Real Hliq[Nc](each unit = "kJ/kmol") "Comopent molar enthalpy";
  Real Hvap[Nc](each unit = "kJ/kmol") "Comopent molar enthalpy";
  Real xliq(unit = "-", min = 0, max = 1, start = xliqg)"Liquid phase mole fraction";
  Real xvap(unit = "-", min = 0, max = 1, start = xvapg) "Vapor phase mole fraction";
  Real M[Nc](each unit = "mol") "Component Molar Holdup", ML(unit = "mol", start = 2000) "Liquid Holdup", MV(unit = "mol") "Vapor Holdup", MT(unit = "mol");
  Real VL(unit = "m^3") "Volume of Liquid holdup", VG(unit = "m^3") "Volume of Gas Holdup";
  Real Q(unit = "W") "Heat supplied or Removed", rholiq_c[Nc](each unit = "mol/m^3") "Liquid Density", h(unit = "m") "Height of liquid holdup", rholiq(unit = "mol/m^3");
  Real P(unit = "Pa", start = 101325) "Pressure of Gas Holdup", Pt "Total Column Pressure";
  Real Pbubl(unit = "Pa", min = 0, start = Pmin) "Bubble point pressure";
  Real Pdew(unit = "Pa", min = 0, start = Pmax) "Dew point pressure";
  
//===============================================================================
//Instantiation of Connectors
  Simulator.Files.Interfaces.matConn In(Nc = Nc) annotation(
    Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Simulator.Files.Interfaces.matConn Out1(Nc = Nc) annotation(
    Placement(visible = true, transformation(origin = {102, 72}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Simulator.Files.Interfaces.matConn Out2(Nc = Nc) annotation(
    Placement(visible = true, transformation(origin = {100, -72}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  
  extends GuessModels.InitialGuess;
initial equation
P = Pset;
h = hset;
if Dynamics == true then
der(MT) = 0;
for i in 1:Nc-1 loop
der(M[i]) = 0;
end for;
end if;  
equation
//================================================================================
//Connector equation
  if BTdef then
    Tdef = T;
  else
    In.T = T;
  end if;
  if BPdef then
    Pdef = P;
  else
    In.P = Pin;
  end if;
  In.F = F_p[1];
  In.x_pc[1, :] = x_pc[1, :];
  In.H = Hin ;
  Out2.T = T;
  Out2.P = Pin;
  Out2.F = F_p[2];
  Out2.x_pc[1, :] = x_pc[2, :];
  Out1.T = T;
  Out1.P = Pin;
  Out1.F = F_p[3];
  Out1.x_pc[1, :] = x_pc[3, :];
  
//=================================================================================
//Thermodynamic Equations
for i in 1:Nc loop
rholiq_c[i] = Simulator.Files.ThermodynamicFunctions.Dens(C[i].LiqDen, C[i].Tc, T, P);
end for;
//=================================================================================
//Mole Balance
  F_p[1] - F_p[2] - F_p[3] = if Dynamics == true then der(MT) else 0;
for i in 1:Nc-1 loop
  x_pc[1, i] .* F_p[1] - x_pc[2, i] .* F_p[2] - x_pc[3, i] .* F_p[3] = if Dynamics == true then der(M[i]) else 0;
  M[i] = ML .* x_pc[2, i] + MV .* x_pc[3, i];
end for;
  sum(M) = MT;
  MT = ML + MV;
  rholiq =  sum(x_pc[2, :] .* rholiq_c) ;
  VL = ML /sum(x_pc[2, :] .* rholiq_c) ;
  Pt = P + rholiq * 9.81 * h;
  F_p[2] = xtl * Cd * sqrt(Pt - 1e5);
  F_p[3] = xtv * Cd * sqrt(P - 1e5) ;
  P = (MV * 8.314 * T) / VG;
  VL = A * h;
  VT = VL + VG;
//===================================================================================
//Bubble point calculation
  Pbubl = sum(gmabubl_c[:] .* x_pc[1, :] .* exp(C[:].VP[2] + C[:].VP[3] / T + C[:].VP[4] * log(T) + C[:].VP[5] .* T .^ C[:].VP[6]) ./ philiqbubl_c[:]);
//Dew point calculation
  Pdew = 1 / sum(x_pc[1, :] ./ (gmadew_c[:] .* exp(C[:].VP[2] + C[:].VP[3] / T + C[:].VP[4] * log(T) + C[:].VP[5] .* T .^ C[:].VP[6])) .* phivapdew_c[:]);
  if P >= Pbubl then
//below bubble point region
    F_p[3] = 1e-15;
    F_p[2] = F_p[1];
    x_pc[2, :] = x_pc[1, :];
//VLE region
  elseif P >= Pdew then
    for i in 1:Nc loop
    x_pc[3, i] = K_c[i] * x_pc[2, i];
    end for;
    sum(x_pc[3, :]) = 1;
    sum(x_pc[2, :]) = 1;
  else
//above dew point region
    F_p[3] = F_p[1];
    F_p[2] = 1e-15;
    x_pc[3, :] = x_pc[1, :];
  end if;
//===================================================================================
//Energy Balance 
  F_p[1] * Hin + Q - F_p[2] * sum(Hliq[:] .* x_pc[2, :]) - F_p[3] * sum(Hvap[:] .* x_pc[3, :]) = 0;
  
//===================================================================================
//Enthalpy calculation from Thermodynamic Functions
   for i in 1:Nc loop 
    Hliq[i] = ThermodynamicFunctions.HLiqId(C[i].SH, C[i].VapCp, C[i].HOV, C[i].Tc, T);
    Hvap[i] = ThermodynamicFunctions.HVapId(C[i].SH, C[i].VapCp, C[i].HOV, C[i].Tc, T);
  end for;

//=======================================================================================
//phase molar fractions
  xliq = F_p[2] / F_p[1];
  xvap = F_p[3] / F_p[1];
annotation(
    Icon(coordinateSystem(extent = {{-100, -200}, {100, 200}})),
    Diagram(coordinateSystem(extent = {{-100, -200}, {100, 200}})),
    __OpenModelica_commandLineOptions = "",
  Documentation(info = "<html><head></head><body>The <b>Flash</b> column is used to calculate the vapor and liquid phase distribution for a mixed phase material stream.<div><br></div><div><div><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13px; orphans: 2; widows: 2;\">The flash column model have three Material Stream connection ports as:</span></div><div><div style=\"orphans: 2; widows: 2;\"><ol><li><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13px;\">a feed stream</span></li><li><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13px;\">two outlet streams as vapor and liquid</span></li></ol></div></div><div><br></div><div>Following calculation parameters may be provided for the flash column:</div><div><ol><li>Boolean for Separation Temperature (<b>BTdef</b>)</li><li>Boolean for Separation Pressure (<b>BPdef</b>)</li><li>Separation Temperature (<b>Tdef</b>)</li><li>Separation Pressure (<b>Pdef</b>)</li></ol><div>The first two variables (<b>BTdef</b> and <b>BPdef</b>) are of type boolean to indicate if the flash is operated at temperature and pressure other than the temperature and pressure at which feed stream is entering. If so, the variables are to be entered as <b>True</b> else can be left to be blank to take the default value as False.</div><div><br></div><div>The other two variables (<b>Tdef</b> and <b>Pdef</b>) are of type <i>parameter Real</i> and are to be specified only if the corresponding boolean variables, <b>BTdef</b> and <b>BPdef</b> respectively<b>,</b>&nbsp;are <b>True</b>. In some cases, any one of the boolean variables can be True only and therefore its value at corresponding variable needs to specified.</div><div><br></div><div>In case either both or any one of the boolean variables among BTdef and BPdef are defined as True and its corresponding value at Tdef and Pdef are not specified, the flash will be operated at temperature of 298.15 K and pressure of 101325 Pa.</div><div><br></div><div><span style=\"font-size: 12px;\">During simulation, their values can be specified directly under&nbsp;</span><b style=\"font-size: 12px;\">Flash Specifications&nbsp;</b><span style=\"font-size: 12px;\">by double clicking on the flash model instance.</span></div><div><span style=\"font-size: 12px;\"><br></span></div><div><br></div></div><div><span style=\"font-size: 12px;\">For detailed explaination on how to use this model to simulate a Flash, go to&nbsp;</span><a href=\"modelica://Simulator.Examples.Flash\" style=\"font-size: 12px;\">Flash Example</a><span style=\"font-size: 12px;\">.</span></div></div></body></html>"));
  end Flash;
