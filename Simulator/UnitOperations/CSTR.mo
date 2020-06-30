within Simulator.UnitOperations;

model CSTR
  //===============================================================================
  //Header Files and Parameters
  import data = Simulator.Files.ChemsepDatabase;
  parameter data.GeneralProperties C[Nc] "Component instances array" annotation(
    Dialog(tab = "Reactor Specifications", group = "Component Parameters"));
  import Simulator.Files.*;
  import Simulator.Files.ThermodynamicFunctions.*;
  extends Simulator.Files.Icons.CSTR;
  extends Simulator.GuessModels.InitialGuess;
  parameter String Mode "Required mode of operation: Isothermal, Define_Out_Temperature, Adiabatic" annotation(
    Dialog(tab = "Reactor Specifications", group = "Calculation Parameters"));
  parameter Integer Nc "Number of Components" annotation(
    Dialog(tab = "Reactor Specifications", group = "Component Parameters"));
  parameter Real Pdel(unit = "Pa") "Pressure Drop" annotation(
    Dialog(tab = "Reactor Specifications", group = "Calculation Parameters"));
  parameter Real Tdef(unit = "K") "Define outlet Temperature , applicable if Define_Out_Temperature mode is chosen" annotation(
    Dialog(tab = "Reactor Specifications", group = "Calculation Parameters"));
  parameter String Phase "Reaction Phase: Mixture, Liquid, Vapour" annotation(
    Dialog(tab = "Reactions", group = "Kinetic Reaction Parameters"));
  parameter Real V(unit = "m^3") "Reactor Volume" annotation(
    Dialog(tab = "Reactor Specifications", group = "Calculation Parameters"));
  parameter Boolean Dynamics "To operate the reactor in unsteady mode" annotation(
    Dialog(tab = "Reactor Specifications"));
  parameter Real A(unit = "m^2") = 1 "Area of the Reactor" annotation(
    Dialog(tab = "Physical Specifications", group = "Vessel Geometry"));
  parameter Real Cdl(unit = "-") = 0.5 "Discharge coefficient of valve" annotation(
    Dialog(tab = "Physical Specifications", group = "Vessel Geometry"));
  parameter Real Cdv(unit = "-") = 0.5;
  parameter Real xtl(fixed = false), xtv(fixed = false) "Stem Position of the valve";
  parameter Real Pset = 101325 "Initial Pressure" annotation(
    Dialog(tab = "Physical Specifications", group = "Vessel Geometry")); 
//===============================================================================
//Model Variables
  Real Pin(unit = "Pa", start = Pg) "Inlet Pressure";
  Real Tin(unit = "K", start = Tg) "Inlet Temperature";
  Real Pout(unit = "Pa", start = Pg) "Outlet Pressure";
  Real Tout(unit = "K", start = Tg) "Outlet Temperature";
  Real Hin(unit = "kJ/kmol", start = Htotg) "Mixture Mol Enthalpy in the inlet";
  Real Hout(unit = "kJ/kmol", start = Htotg) "Mixture Mol Enthalpy in the outlet";
  Real Fin(unit = "mol/s") "Inlet total mole flow";
  Real Fout_p[3](each unit = "mol/s", start={Fg,Fliqg,Fvapg}) "Outlet total mole flow";
  Real Fout_pc[3, Nc](each unit = "mol/s", each start = Fg) "Component mole flow at outlet";
  Real xin_pc[1, Nc](each unit = "-") "Component mole fraction at inlet";
  Real xout_pc[3, Nc](each unit = "-", start = {xguess,xg,yg}) "Component mol-frac at outlet";
  Real Q(unit = "W") "Energy supplied/removed";
  Real Hr(unit = "kJ/kmol") "Heat of Reaction";
  Real C_c[Nc](each unit = "mol/m^3") "Concentration of Component";
  Real rholiq_c[Nc](each unit = "mol/m^3");
  Real rholiq(unit = "mol/m^3");
  Real Fvout_p[3](each unit = "m^3/s") "Mixture Volumetric Flow rate";
  Real FO_c[Nc](each unit = "mol/s") "Inlet Flow";
  Real F_c[Nc](each unit = "mol/s") "Outlet Flow";
  Real X(unit = "-") "Conversion of Base component";
  Real k_r[Nr] "Rate constant";
  Real r_rc[Nc, Nr](each unit = "mol/m^3.S") "Rate of Reaction";
  Real r_c[Nc](each unit = "mol/m^3.S") "Overall Rate of Reaction";
  Real xvap(unit = "-", start = xvapg) "Vapor Phase Mole Fraction";
  Real Vliq(unit = "m^3") "Volume of Liquid in Reactor";
  Real Vvap(unit = "m^3") "Volume of Liquid in Reactor";
  Real tau(unit = "s") "Residence Time";
  Real Vrxn(start = V) "Actual Reaction Volume";
  Real M[Nc](each unit = "mol", each start = 2000) "Component Molar Holdup";
  Real MT(unit = "mol") "Total Holdup";
  Real ML(unit = "mol") "Liquid Holdup", MV(unit = "mol") "Vapour Holdup";
  Real VL(unit = "m^3") "Volume of Liquid holdup";
  Real VG(unit = "m^3") "Volume of Vapour holdup";
  Real h(unit = "m", start = 0.25) "Height of liquid holdup";
  Real Pt(unit = "Pa") "Total Reactor Pressure", P(start = Pset);
  //===============================================================================
  //Instatiation of Connectors
  Simulator.Files.Interfaces.matConn In(Nc = Nc) annotation(
    Placement(visible = true, transformation(origin = {-94, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-98, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Simulator.Files.Interfaces.matConn Out(Nc = Nc) annotation(
    Placement(visible = true, transformation(origin = {94, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Simulator.Files.Interfaces.enConn En annotation(
    Placement(visible = true, transformation(origin = {0, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -142}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  extends Simulator.Files.Models.ReactionManager.KineticReaction(Nr = 1, BC_r = {1}, Coef_cr = {{-1}, {-1}, {1}}, DO_cr = {{1}, {0}, {0}}, Af_r = {0.5}, Ef_r = {0});
  
  initial equation
  P = Pset;
  h = 0.3;
  if Dynamics == true then
  for i in 1:Nc loop
  der(M[i]) = 0 ;
  end for;
  end if;
  
  equation
//===============================================================================  
//Connector-Equations
  In.P = Pin;
  In.T = Tin;
  In.F = Fin;
  In.H = Hin;
  In.x_pc[1, :] = xin_pc[1, :];
  Out.P = Pout;
  Out.T = Tout;
  Out.F = Fout_p[1];
  Out.H = Hout;
  Out.x_pc[:, :] = xout_pc[:, :];
  Out.xvap = xvap;
  En.Q = Q;
//===============================================================================
//Mole Balance
  FO_c[:] = xin_pc[1, :] * Fin;
 for i in 1:Nc loop
    FO_c[i] - F_c[i] + r_c[i] * Vrxn  = if Dynamics == true then der(M[i]) else 0;
 end for;
  X = (FO_c[BC_r[1]] - F_c[BC_r[1]]) / FO_c[BC_r[1]] ;
  sum(F_c[:]) = Fout_p[1];
  for i in 1:Nc loop
    xout_pc[1, i] = F_c[i] / Fout_p[1];
    Fout_pc[:, i] = Fout_p[:] .* xout_pc[:, i];
  end for;
  Fout_p[1] = Fout_p[2] + Fout_p[3];
  xvap = Fout_p[3] / Fout_p[1];
  k_r[:] = Simulator.Files.Models.ReactionManager.Arhenious(Nr, Af_r[:], Ef_r[:], Tout);
//===============================================================================
//Holdup Calculations and Valve equation  
  M = MT .* xout_pc[1, :] ;
  ML = rholiq * VL ;
  MV = (P * VG) / (8.314 * Tout) ;
  VL = A * h;
  VG = V - VL ;
  MT = ML + MV;
  Pt = P + rholiq * 9.81 * h;
  if xvap >= 1 then
  h = 0;
  else
  Fout_p[2] = xtl * Cdl * sqrt(Pt - 1e5);
  end if;
  if xvap <= 0 then
  P = 101325;
  else
  Fout_p[3] = xtv * Cdv * sqrt(P - 1e5);
  end if;
//===============================================================================
//Calculation of Mixer Density and Volumetric Flowrate
  for i in 1:Nc loop
    rholiq_c[i] = Simulator.Files.ThermodynamicFunctions.Dens(C[i].LiqDen, C[i].Tc, Tout, P);
  end for;
  rholiq = 1 / sum(xout_pc[2, :] ./ rholiq_c[:]);
  Fvout_p[1] = Fvout_p[2] + Fvout_p[3];
  Fvout_p[2] = Fout_p[2] / rholiq;
  Fvout_p[3] = Fout_p[3] * 8.314 * Tout / Pout;
//===============================================================================
//Volume of Reaction Mixture
  tau = V / Fvout_p[1] ;
  tau = Vliq / Fvout_p[2] ;
  tau = Vvap / Fvout_p[3] ;
//===============================================================================
//Concentration of the Components
  if Phase == "Mixture" then
    C_c[:] = Fout_pc[1, :] ./ Fvout_p[1];
    Vrxn = V;
  elseif Phase == "Liquid" then
    C_c[:] = Fout_pc[2, :] ./ Fvout_p[2];
    Vrxn = Vliq;
  elseif Phase == "Vapour" then
    C_c[:] = P .* xout_pc[3, :] / (8.314 * Tout);
    Vrxn = Vvap;
  end if;
//===============================================================================
//Rate Equation
  for i in 1:Nr loop
    r_rc[:, i] = (Coef_cr[:, i] / abs(Coef_cr[BC_r[i], i])) * k_r[i] * product(C_c[:] .^ DO_cr[:, i]);
  end for;
  for i in 1:Nc loop
    r_c[i] = sum(r_rc[i, :]);
  end for;
  Pout = Pin - Pdel;
//===============================================================================
//Energy Balance
//===============================================================================
//Isothermal Operation
  if Mode == "Isothermal" then
    Hr = Hr_r[1] * FO_c[BC_r[Nr]] .* X;
    Fin * Hin - Fout_p[1] * Hout + Q - Hr = 0;
    Tin = Tout;
//===============================================================================
//Outlet Temperature Defined
  elseif Mode == "Define_Out_Temperature" then
    Tout = Tdef;
    Hr = Hr_r[1] * FO_c[BC_r[Nr]] .* X;
    Fin * Hin - Fout_p[1] * Hout + Q - Hr = 0;
//===============================================================================
//Adiabatic Operation
  elseif Mode == "Adiabatic" then
    Q = 0;
    Hr = Hr_r[1] * FO_c[BC_r[Nr]] .* X;
    Fin * Hin + Q = Fout_p[1] * Hout + Hr;
  end if;
//===============================================================================
  annotation(
    Documentation(info = "<html><head></head><body><div style=\"font-size: 12px;\">The <b>Continuous Stirred Tank&nbsp;Reactor (CSTR)</b>&nbsp;is used to calculate the mole fraction of components at outlet stream when the reaction kinetics and reactor volume is defined.</div><div style=\"font-size: 12px;\"><br></div><div style=\"font-size: 12px;\"><div><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13px; orphans: 2; widows: 2;\">The CSTR model have following connection ports:</span></div><div><div style=\"orphans: 2; widows: 2;\"><ol><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Two Material Streams:</span></font></li><ul><li><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13px;\">feed stream</span></li><li><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13px;\">outlet stream</span></li></ul><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">One Energy Stream:</span></font></li><ul><li><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13px;\">heat added</span></li></ul></ol></div></div></div><div style=\"font-size: 12px;\"><br></div><span style=\"font-size: 12px;\">To simulate a CSTR, following calculation parameters must be provided:</span><div><ol style=\"font-size: 12px;\"><li>Calculation Mode (<b>Mode</b>)</li><li>Reaction Phase (<b>Phase</b>)</li><li>Outlet Temperature&nbsp;(<b>Tdef</b>)&nbsp;(If calculation mode is Define_Out_Temperature)</li><li>Number of Reactions (<b>Nr</b>)</li><li>Base Component (<b>Base_C</b>)</li><li>Stoichiometric Coefficient of Components in Reaction (<b>Coef_cr</b>)</li><li>Reaction Order (<b>DO_cr</b>)</li><li>Pre-exponential Factor (<b>Af_r</b>)</li><li>Activation Energy (<b>Ef_r</b>)</li><li>Pressure Drop (<b>Pdel</b>)</li><li>Reactor Volume (<b>V</b>)</li><li>Reactor Headspace (<b>VHspace</b>)</li></ol><div><div style=\"font-size: 12px; orphans: 2; widows: 2;\"><span style=\"orphans: auto; widows: auto;\">Among the above variables, only the first variable is of type&nbsp;<i>parameter String</i>. Calculation Mode (<b>Mode</b>) can have either of the sting values among the following:</span></div><div style=\"orphans: 2; widows: 2;\"><ol style=\"font-size: 12px;\"><li><b>Isothermal</b>: If the reactor is operated isothermally</li><li><b>Define_Out_Temperature</b>: If the reactor is operated at specified outlet temperature</li><li><b>Adiabatic</b>: If the reactor is operated adiabatically</li></ol><div style=\"font-size: 12px;\"><br></div></div></div><div style=\"font-size: 12px;\">The other variables are of type&nbsp;<i>parameter Real.</i></div></div><div style=\"font-size: 12px;\"><div>During simulation, their values can specified directly under&nbsp;<b>Reactor Specifications and Reactions&nbsp;</b>by double clicking on the CSTR model instance.</div><div><br></div></div><div style=\"font-size: 12px;\"><br></div><div style=\"font-size: 12px;\">For detailed explaination on how to use this model to simulate a Continuous Stirred Tank Reactor, go to&nbsp;<a href=\"modelica://Simulator.Examples.CSTR\">CSTR Example</a>.</div><div><br></div></body></html>"));
    end CSTR;
