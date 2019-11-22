within Simulator.UnitOperations;

model ShortcutColumn
//==============================================================================
//Header Files and Parameters
  extends Simulator.Files.Icons.DistillationColumn;
  import data = Simulator.Files.ChemsepDatabase;
  parameter data.GeneralProperties C[Nc];
  parameter Integer Nc "Number of Components";
  parameter Integer HKey "Heavy Key Component", LKey "Light Key Component";
  parameter String Ctype = "Total";
  
//==============================================================================
//Model Variables
  Real F_p[3](each min = 0, each start = 100) "Inlet Molar Flow";
  Real x_pc[3, Nc](each start = 1 / (Nc + 1), each min = 0, each max = 1) "Inlet Mole Fraction";
  Real H_p[3]"Inlet Molar Enthalpy ";
  Real S_p[3] "Inlet Molar Entropy";
  Real Pin(min = 0, start = 101325) "Feed Pressure";
  Real Tin(min = 0, start = 273.15)"Feed Temperature";
  Real xin_pc[3, Nc](each min = 0, each max = 1, each start = 1 / (Nc + 1)) "Adjustment for Thermodynamics";
  
  Real Ntmin(min = 0, start = 10) "Minimum Number of trays";
  Real RRmin(start = 1) "Minimum Reflux Ratio";
  Real alpha_c[Nc] "Relative Volatility";
  Real theta(start = 1) "Fraction";
  Real T "Thermodynamic Adjustment", P "Thermodynamic Adjustment";
  Real Tcond(start = max(C[:].Tb), min = 0)"Condenser Temperature";
  Real Pcond(min = 0, start = 101325) "Condenser Pressure";
  Real Preb(min = 0, start = 101325)"Re-boiler Pressure";
  Real Treb(start = min(C[:].Tb), min = 0) "Re-boiler Temperature";
  Real xvap_p[3](each min = 0, each max = 1, each start = 0.5) "Vapor Phase Mole Fraction";
  Real Hliqcond "Total Enthalpy of Liquid in Condenser";
  Real Hvapcond "Total Enthalpy of Vapor in Condenser";
  Real Hvapcond_c[Nc] "Component Enthalpy of Vapor in Condenser";
  Real Hliqcond_c[Nc] "Component Enthalpy of Vapor in Condenser";
  Real xliqcond_c[Nc](each min = 0, each max = 1, each start = 1 / (Nc + 1))"Liquid Mole Fraction in Condenser";
  Real xvapcond_c[Nc](each min = 0, each max = 1, each start = 1 / (Nc + 1))"Vapor Mole Fraction in Condenser";
  
  Real Pdew(min = 0, start = sum(C[:].Pc) / Nc)"Dew Point Pressure";
  Real Pbubl(min = 0, start = sum(C[:].Pc) / Nc)"Bubble Point Pressure";
  Real RR "Reflux Ratio";
  Real Nt "Number of Trays";
  Real x "Intermediate 1";
  Real y "Intermediate 1";
  Real Intray "Feed Tray";
  Real Fliqrec(min = 0, start = 100) "Liquid Flow in Rectification Section";
  Real Fvaprec(min = 0, start = 100)"Vapor Flow in Rectification Section";
  Real Fliqstrip(min = 0, start = 100) "Liquid Flow in Stripping Section";
  Real Fvapstrip(min = 0, start = 100)"Vapor Flow in Stripping Section";
  Real Qr "Re-boiler Duty";
  Real Qc "Condenser Duty";

//==============================================================================
//Instantiation of Connections
  Simulator.Files.Interfaces.matConn In(Nc = Nc) annotation(
    Placement(visible = true, transformation(origin = {-250, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-250, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Simulator.Files.Interfaces.matConn Out1(Nc = Nc) annotation(
    Placement(visible = true, transformation(origin = {250, 336}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {250, 300}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Simulator.Files.Interfaces.matConn Out2(Nc = Nc) annotation(
    Placement(visible = true, transformation(origin = {250, -266}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {250, -300}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Simulator.Files.Interfaces.enConn En1 annotation(
    Placement(visible = true, transformation(origin = {248, 594}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {250, 600}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Simulator.Files.Interfaces.enConn En2 annotation(
    Placement(visible = true, transformation(origin = {254, -592}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {250, -600}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

equation
//==============================================================================
// Connector equations
  In.P = Pin;
  In.T = Tin;
  In.F = F_p[1];
  In.x_pc[1, :] = x_pc[1, :];
  In.H = H_p[1];
  In.S = S_p[1];
  In.xvap = xvap_p[1];
  Out2.P = Preb;
  Out2.T = Treb;
  Out2.F = F_p[2];
  Out2.x_pc[1, :] = x_pc[2, :];
  Out2.H = H_p[2];
  Out2.S = S_p[2];
  Out2.xvap = xvap_p[2];
  Out1.P = Pcond;
  Out1.T = Tcond;
  Out1.F = F_p[3];
  Out1.x_pc[1, :] = x_pc[3, :];
  Out1.H = H_p[3];
  Out1.S = S_p[3];
  Out1.xvap = xvap_p[3];
  En2.Q = Qr;
  En1.Q = Qc;

//==============================================================================
//Adjustment for thermodynamic packages
  xin_pc[1, :] = x_pc[1, :];
  xin_pc[2, :] = xin_pc[1, :] ./ (1 .+ xvap_p[1] .* (K_c[:] .- 1));
  for i in 1:Nc loop
    xin_pc[3, i] = K_c[i] * xin_pc[2, i];
  end for;
  T = Tin;
  P = Pin;
//==============================================================================
//Bubble point calculation
  Pbubl = sum(gmabubl_c[:] .* xin_pc[1, :] .* exp(C[:].VP[2] + C[:].VP[3] / Tin + C[:].VP[4] * log(Tin) + C[:].VP[5] .* Tin .^ C[:].VP[6]) ./ philiqbubl_c[:]);
  
//==============================================================================
//Dew point calculation
  Pdew = 1 / sum(xin_pc[1, :] ./ (gmadew_c[:] .* exp(C[:].VP[2] + C[:].VP[3] / Tin + C[:].VP[4] * log(Tin) + C[:].VP[5] .* Tin .^ C[:].VP[6])) .* phivapdew_c[:]);
  for i in 1:Nc loop
    if x_pc[1, i] == 0 then
      x_pc[3, i] = 0;
    else
      F_p[1] .* x_pc[1, i] = F_p[2] .* x_pc[2, i] + F_p[3] .* x_pc[3, i];
    end if;
  end for;
  sum(x_pc[3, :]) = 1;
  sum(x_pc[2, :]) = 1;
  
//==============================================================================
//Distillate and Bottom composition
  for i in 1:Nc loop
    if i <> HKey then
      if Ctype == "Total" then
        x_pc[3, i] / x_pc[3, HKey] = alpha_c[i] ^ Ntmin * (x_pc[2, i] / x_pc[2, HKey]);
      elseif Ctype == "Partial" then
        x_pc[3, i] / x_pc[3, HKey] = alpha_c[i] ^ (Ntmin + 1) * (x_pc[2, i] / x_pc[2, HKey]);
      end if;
    end if;
  end for;
//==============================================================================
//Relative Volatility
  alpha_c[:] = K_c[:] / K_c[HKey];
  
//==============================================================================
//Calculation of temperature at Distillate and Bottoms
  if Tcond <= 0 and Treb <= 0 then
  if Ctype == "Partial" then
    1 / Pcond = sum(x_pc[3, :] ./ (gma_c[:] .* exp(C[:].VP[2] + C[:].VP[3] / 1 + C[:].VP[4] * 1 + C[:].VP[5] .* Tcond .^ C[:].VP[6])));
    
    Preb = sum(gma_c[:] .* x_pc[2, :] .* exp(C[:].VP[2] + C[:].VP[3] / 1 + C[:].VP[4] * 1 + C[:].VP[5] .* Treb .^ C[:].VP[6]));
    
  elseif Ctype == "Total" then
    Pcond = sum(gma_c[:] .* x_pc[3, :] .* exp(C[:].VP[2] + C[:].VP[3] / 1 + C[:].VP[4] * 1 + C[:].VP[5] .* Tcond .^ C[:].VP[6]));
    
    Preb = sum(gma_c[:] .* x_pc[2, :] .* exp(C[:].VP[2] + C[:].VP[3] / 1 + C[:].VP[4] * 1 + C[:].VP[5] .* Treb .^ C[:].VP[6]));
  end if;
 //============================================================================== 
elseif Tcond<=0 then
  if Ctype == "Partial" then
    1 / Pcond = sum(x_pc[3, :] ./ (gma_c[:] .* exp(C[:].VP[2] + C[:].VP[3] / 1 + C[:].VP[4] * 1 + C[:].VP[5] .* Tcond .^ C[:].VP[6])));
    
    Preb = sum(gma_c[:] .* x_pc[2, :] .* exp(C[:].VP[2] + C[:].VP[3] / Treb + C[:].VP[4] * log(Treb) + C[:].VP[5] .* Treb .^ C[:].VP[6]));
    
  elseif Ctype == "Total" then
    Pcond = sum(gma_c[:] .* x_pc[3, :] .* exp(C[:].VP[2] + C[:].VP[3] / 1 + C[:].VP[4] * 1 + C[:].VP[5] .* Tcond .^ C[:].VP[6]));
    
    Preb = sum(gma_c[:] .* x_pc[2, :] .* exp(C[:].VP[2] + C[:].VP[3] / Treb + C[:].VP[4] * log(Treb) + C[:].VP[5] .* Treb .^ C[:].VP[6]));
  end if;
//==============================================================================
elseif Treb<=0 then
  if Ctype == "Partial" then
    1 / Pcond = sum(x_pc[3, :] ./ (gma_c[:] .* exp(C[:].VP[2] + C[:].VP[3] / Tcond + C[:].VP[4] * log(Tcond) + C[:].VP[5] .* Tcond .^ C[:].VP[6])));
    
    Preb = sum(gma_c[:] .* x_pc[2, :] .* exp(C[:].VP[2] + C[:].VP[3] / 1 + C[:].VP[4] * 1 + C[:].VP[5] .* Treb .^ C[:].VP[6]));
    
  elseif Ctype == "Total" then
    Pcond = sum(gma_c[:] .* x_pc[3, :] .* exp(C[:].VP[2] + C[:].VP[3] / Tcond + C[:].VP[4] * log(Tcond) + C[:].VP[5] .* Tcond .^ C[:].VP[6]));
    
    Preb = sum(gma_c[:] .* x_pc[2, :] .* exp(C[:].VP[2] + C[:].VP[3] / 1 + C[:].VP[4] * 1 + C[:].VP[5] .* Treb .^ C[:].VP[6]));
  end if;
//==============================================================================
else
  if Ctype == "Partial" then
    1 / Pcond = sum(x_pc[3, :] ./ (gma_c[:] .* exp(C[:].VP[2] + C[:].VP[3] / Tcond + C[:].VP[4] * log(Tcond) + C[:].VP[5] .* Tcond .^ C[:].VP[6])));
    
    Preb = sum(gma_c[:] .* x_pc[2, :] .* exp(C[:].VP[2] + C[:].VP[3] / Treb + C[:].VP[4] * log(Treb) + C[:].VP[5] .* Treb .^ C[:].VP[6]));
    
  elseif Ctype == "Total" then
    Pcond = sum(gma_c[:] .* x_pc[3, :] .* exp(C[:].VP[2] + C[:].VP[3] / Tcond + C[:].VP[4] * log(Tcond) + C[:].VP[5] .* Tcond .^ C[:].VP[6]));
    
    Preb = sum(gma_c[:] .* x_pc[2, :] .* exp(C[:].VP[2] + C[:].VP[3] / Treb + C[:].VP[4] * log(Treb) + C[:].VP[5] .* Treb .^ C[:].VP[6]));
  end if;
end if;

//==============================================================================
//Minimum Reflux, Underwoods method
if theta > alpha_c[LKey] or theta < alpha_c[HKey] then
  0= sum((alpha_c[:] .* x_pc[1, :])./ (alpha_c[:] .- theta));
else
  xvap_p[1] = sum((alpha_c[:] .* x_pc[1, :])./ (alpha_c[:] .- theta));
end if;
  RRmin + 1 = sum(alpha_c[:] .* x_pc[3, :] ./ (alpha_c[:] .- theta));
  
//==============================================================================
//Actual number of trays,Gillilands method
  x = (RR - RRmin) / (RR + 1);
  y = (Nt - Ntmin) / (Nt + 1);
  if x >= 0 then
    y = 0.75 * (1 - x ^ 0.5668);
  else
    y = -1;
  end if;
//==============================================================================
// Feed location, Fenske equation
  Intray = Nt / Ntmin * log(x_pc[1, LKey] / x_pc[1, HKey] * (x_pc[2, HKey] / x_pc[2, LKey])) / log(K_c[LKey] / K_c[HKey]);
  
//==============================================================================
//Rectifying and Stripping flows Calculation
  Fliqrec = RR * F_p[3];
  Fliqstrip = (1 - xvap_p[1]) * F_p[1] + Fliqrec;
  Fvapstrip = Fliqstrip - F_p[2];
  Fvaprec = xvap_p[1] * F_p[1] + Fvapstrip;
  for i in 1:Nc loop
    Hvapcond_c[i] = Simulator.Files.ThermodynamicFunctions.HVapId(C[i].SH, C[i].VapCp, C[i].HOV, C[i].Tc, Tcond);
    Hliqcond_c[i] = Simulator.Files.ThermodynamicFunctions.HLiqId(C[i].SH, C[i].VapCp, C[i].HOV, C[i].Tc, Tcond);
  end for;
  if Ctype == "Total" then
    Hliqcond = H_p[3];
  elseif Ctype == "Partial" then
    Hliqcond = sum(xliqcond_c[:] .* Hliqcond_c[:]);
  end if;
  Hvapcond = sum(xvapcond_c[:] .* Hvapcond_c[:]);
  Fvaprec .* xvapcond_c[:] = Fliqrec .* xliqcond_c[:] + F_p[3] .* x_pc[3, :];
  if Ctype == "Partial" then
    x_pc[3, :] = K[:] .* xliqcond_c[:];
  elseif Ctype == "Total" then
    x_pc[3, :] = xliqcond_c[:];
  end if;
//==============================================================================
//Energy Balance
  F_p[1] * H_p[1] + Qr - Qc = F_p[2] * H_p[2] + F_p[3] * H_p[3];
  Fvaprec * Hvapcond = Qc + F_p[3] * H_p[3] + Fliqrec * Hliqcond;
annotation(
    Icon(coordinateSystem(extent = {{-250, -600}, {250, 600}})),
    Diagram(coordinateSystem(extent = {{-250, -600}, {250, 600}})),
    __OpenModelica_commandLineOptions = "");end ShortcutColumn;