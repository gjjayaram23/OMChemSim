within Simulator.UnitOperations;

package PackedColumn
extends Modelica.Icons.Package;

  model Cond "Model of a condenser used in distillation column"
    import Simulator.Files.*;
    parameter ChemsepDatabase.GeneralProperties C[Nc];
    parameter Integer Nc = 2 "Number of components";
    parameter Boolean Bin = false;
    Real P(unit = "K", min = 0, start = Pg) "Pressure";
    Real T(unit = "Pa", min = 0, start = Tg) "Temperature";
    Real Fin(unit = "mol/s", min = 0, start = Fg) "Feed molar flow rate";
    Real xin_c[Nc](each unit = "-", each min = 0, each max = 1, start = xg) "Feed components mole fraction";
    Real xvapin_c[Nc](each unit = "-", each min = 0, each max = 1, start = xvapg) "Inlet components vapor molar fraction";
    Real Hin(unit = "kJ/kmol", start = Htotg) "Feed inlet molar enthalpy";
    Real Fout(unit = "mol/s", min = 0, start = Fg) "Side draw molar flow";
    Real Fvapin(unit = "mol/s", min = 0, start = Fg) "Inlet vapor molar flow";
    Real Fliqout(unit = "mol/s", min = 0, start = Fg) "Outlet liquid molar flow";
    Real xout_c[Nc](each unit = "-", each min = 0, each max = 1, start = xg) "Side draw components mole fraction";
    Real xliqout_c[Nc](each unit = "-", each min = 0, each max = 1, start = xliqg) "Outlet components liquid mole fraction";
    Real Hvapin(unit = "kJ/kmol", start = Hvapg) "Inlet vapor molar enthalpy";
    Real Hliqout(unit = "kJ/kmol", start = Hliqg) "Outlet liquid molar enthalpy";
    Real Q(unit = "W") "Heat load";
    Real Hout(unit = "kJ/kmol", start = Htotg) "Side draw molar enthalpy";
    Real Hliqout_c[Nc](each unit = "kJ/kmol") "Outlet liquid components molar enthalpy";
    Real x_pc[3, Nc](each unit = "-", each min = 0, each max = 1, start = {xguess, xguess, xguess}) "Component mole fraction";
    Real Pdew(unit = "Pa", min = 0, start = Pmax) "Dew point pressure";
    Real Pbubl(unit = "Pa", min = 0, start = Pmin) "Bubble point pressure";
    //String sideDrawType(start = "Null");
    //L or V
    parameter String Ctype "Condenser type: Partial or Total";
    replaceable Simulator.Files.Interfaces.matConn In(Nc = Nc) if Bin annotation(
      Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.matConn In_Dmy(Nc = Nc, P = 0, T = 0, x_pc = zeros(3, Nc), F = 0, H = 0, S = 0, xvap = 0) if not Bin annotation(
      Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.matConn Out(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.trayConn Out_Liq(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {-50, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-50, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.trayConn In_Vap(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {50, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {50, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.enConn En annotation(
      Placement(visible = true, transformation(origin = {100, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    extends Simulator.GuessModels.InitialGuess;
  equation
//connector equation
    if Bin then
      In.x_pc[1, :] = xin_c[:];
      In.H = Hin;
      In.F = Fin;
    else
      In_Dmy.x_pc[1, :] = xin_c[:];
      In_Dmy.H = Hin;
      In_Dmy.F = Fin;
    end if;
    Out.P = P;
    Out.T = T;
    Out.x_pc[1, :] = xout_c[:];
    Out.F = Fout;
    Out.H = Hout;
    Out_Liq.F = Fliqout;
    Out_Liq.H = Hliqout;
    Out_Liq.x_c[:] = xliqout_c[:];
    In_Vap.F = Fvapin;
    In_Vap.H = Hvapin;
    In_Vap.x_c[:] = xvapin_c[:];
    En.Q = Q;
//Adjustment for thermodynamic packages
    x_pc[1, :] = (Fout .* xout_c[:] + Fliqout .* xliqout_c[:]) ./ (Fout + Fliqout);
    x_pc[2, :] = xliqout_c[:];
    x_pc[3, :] = K_c[:] .* x_pc[2, :];
//Bubble point calculation
    Pbubl = sum(gmabubl_c[:] .* x_pc[1, :] .* Pvap_c[:] ./ philiqbubl_c[:]);
//Dew point calculation
    Pdew = 1 / sum(x_pc[1, :] ./ (gmadew_c[:] .* Pvap_c[:]) .* phivapdew_c[:]);
//molar balance
//Fin + Fvapin = Fout + Fliqout;
    Fin .* xin_c[:] + Fvapin .* xvapin_c[:] = Fout .* xout_c[:] + Fliqout .* xliqout_c[:];
//equillibrium
    if Ctype == "Partial" then
      xout_c[:] = K_c[:] .* xliqout_c[:];
    elseif Ctype == "Total" then
      xout_c[:] = xliqout_c[:];
    end if;
//summation equation
//  sum(xliqout_c[:]) = 1;
    sum(xout_c[:]) = 1;
// Enthalpy balance
    Fin * Hin + Fvapin * Hvapin = Fout * Hout + Fliqout * Hliqout + Q;
//Temperature calculation
    if Ctype == "Total" then
      P = sum(xout_c[:] .* Pvap_c[:]);
    elseif Ctype == "Partial" then
      1 / P = sum(xout_c[:] ./ Pvap_c[:]);
    end if;
// outlet liquid molar enthalpy calculation
    for i in 1:Nc loop
      Hliqout_c[i] = Simulator.Files.ThermodynamicFunctions.HLiqId(C[i].SH, C[i].VapCp, C[i].HOV, C[i].Tc, T);
    end for;
    Hliqout = sum(xliqout_c[:] .* Hliqout_c[:]) + Hres_p[2];
    annotation(
      Diagram(coordinateSystem(extent = {{-100, -40}, {100, 40}})),
      Icon(coordinateSystem(extent = {{-100, -40}, {100, 40}})),
      __OpenModelica_commandLineOptions = "");
  end Cond;

  model PackingSection "Model of a tray used in distillation column"
    import Simulator.Files.*;
    parameter ChemsepDatabase.GeneralProperties C[Nc];
    parameter Integer Nc = 2 "Number of components";
    parameter Boolean Bin = true;
    parameter Real DL = 5.6e-09 "Liquid Diffusivity";
    parameter Real DV = 4.5e-06 "Vapour Diffusivity";
    parameter Real sigma = 20e-03 "Liquid Surface Tension";
    parameter Real z = 4 "Height of Packing";
    parameter Real dsec(unit = "m") = 1 "Diameter of the packing section";
    parameter Real dp(unit = "m") = 13e-03 "Diameter of Raschig Ring packing" ;
    parameter Real epslon = 0.64 "Voidage";
    parameter String Corrln "Packing Correlation: Shulman or Onda";
    Real P(unit = "Pa", min = 0, start = Pg) "Pressure";
    Real T(unit = "K", min = 0, start = Tg) "Temperature";
    Real Fin(unit = "mol/s", min = 0, start = Fg) "Feed molar flow";
    Real xin_c[Nc](each unit = "-", each min = 0, each max = 1, start = xg) "Feed components mole fraction";
    Real Hin(unit = "kJ/kmol", start = Htotg) "Feed molar enthalpy";
    Real Fout(unit = "mol/s", min = 0, start = Fg) "Sidedraw molar flow";
    Real Fvap_s[2](each unit = "mol/s", each min = 0, start = {Fg, Fg}) "Vapor molar flow";
    Real Fliq_s[2](each unit = "mol/s", each min = 0, start = {Fg, Fg}) "Liquid molar flow";
    Real xout_c[Nc](each unit = "-", each min = 0, each max = 1, start = xg) "Components mole fraction at sidedraw";
    Real xvap_sc[2, Nc](each unit = "-", each min = 0, each max = 1, start = {yg, yg}) "Components vapor mole fraction";
    Real xliq_sc[2, Nc](each unit = "-", each min = 0, each max = 1, start = {xg, xg}) "Components liquid mole fraction";
    Real Hvap_s[2](unit = "kJ/kmol", start = Hvapg) "Vapor molar enthalpy";
    Real Hliq_s[2](unit = "kJ/kmol", start = Hliqg) "Liquid molar enthalpy";
    Real Q(unit = "W") "Heat load";
    Real Hout(unit = "kJ/kmol", start = Htotg) "Side draw molar enthalpy";
    Real Hvapout_c[Nc](unit = "kJ/kmol", start = Hvapg) "Outlet components vapor molar enthalpy";
    Real Hliqout_c[Nc](unit = "kJ/kmol", start = Hliqg) "Outlet components liquid molar enthalpy";
    Real xI[Nc](start = xg) "Interface Liquid Composition", yI[Nc](start = yg) "Interface Vapor Composition";
    Real x_pc[3, Nc](each min =0, each max = 0,start={xguess,xguess,xguess});
    Real Nliq[Nc](each unit = "mol/m^2.s") "Component Molar Flux of liquid", Nvap[Nc](each unit = "mol/m^2.s") "Component Molar Flux of liquid", NV(unit = "mol/m^2.s") "Total Molar Flux of vapor", NL(unit = "mol/m^2.s") "Total Molar Flux of liquid";
    Real Pdmy1, Tdmy1, xdmy1_pc[3, Nc], Fdmy1, Hdmy1, Sdmy1, xvapdmy1;
    Real rholiq_c[Nc] "Molar Liquid Density", rhovap_c[Nc] "Vapor Density", vv(start = 10) "Vapor Velocity", vl(start = 1) "Liquid Velocity";
    Real muL[Nc] "Component Viscosity of Liquid", ReL "Reynolds Number of liquid", ScL "Schmid number of Liquid", muliq "Mixture Viscosity", rholiq "Mixture density" ;
    Real muV[Nc] "Component Viscosity of Vapor", ReV "Reynolds Number of Vapor", ScV "Schmid number of Vapor", muvap "Mixture Viscosity", rhovap "Mixture density" ;
    Real MW_p[3] "Molecular Weight of Phases";
    Real KL "Liquid Mass Transfer Coefficient", KV "Vapor Mass Transfer Coefficient";
    Real ap "Specific Surface area of packing section";
    Real ae "Effective area of packing section";
    //this is adjustment done since OpenModelica 1.11 is not handling array modification properly
    String OutType(start = "Null");
    //L or V
    replaceable Simulator.Files.Interfaces.matConn In(Nc = Nc) if Bin annotation(
      Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable Simulator.Files.Interfaces.matConn In_Dmy(Nc = Nc, P = 0, T = 0, F = 0, x_pc = zeros(3, Nc), H = 0, S = 0, xvap = 0) if not Bin annotation(
      Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.matConn Out(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.trayConn In_Liq(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {-50, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-50, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.trayConn Out_Liq(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {-50, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-50, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.trayConn Out_Vap(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {50, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {50, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.trayConn In_Vap(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {50, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {50, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.enConn En annotation(
      Placement(visible = true, transformation(origin = {100, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    extends Simulator.GuessModels.InitialGuess;
  equation
  //connector equation
    if Bin then
      In.P = Pdmy1;
  //this is adjustment done since OpenModelica 1.11 is not handling array modification properly
      In.T = Tdmy1;
      In.x_pc = xdmy1_pc;
      In.F = Fdmy1;
      In.H = Hdmy1;
      In.S = Sdmy1;
      In.xvap = xvapdmy1;
    else
      In_Dmy.P = Pdmy1;
      In_Dmy.T = Tdmy1;
      In_Dmy.x_pc = xdmy1_pc;
      In_Dmy.F = Fdmy1;
      In_Dmy.H = Hdmy1;
      In_Dmy.S = Sdmy1;
      In_Dmy.xvap = xvapdmy1;
    end if;
  //this is adjustment done since OpenModelica 1.11 is not handling array modification properly
    xdmy1_pc[1, :] = xin_c[:];
    Hdmy1 = Hin;
    Fdmy1 = Fin;
    Out.P = P;
    Out.T = T;
    Out.F = Fout;
    Out.H = Hout;
    In_Liq.F = Fliq_s[1];
    In_Liq.H = Hliq_s[1];
    In_Liq.x_c[:] = xliq_sc[1, :];
    Out_Liq.F = Fliq_s[2];
    Out_Liq.H = Hliq_s[2];
    Out_Liq.x_c[:] = xliq_sc[2, :];
    In_Vap.F = Fvap_s[1];
    In_Vap.H = Hvap_s[1];
    In_Vap.x_c[:] = xvap_sc[1, :];
    Out_Vap.F = Fvap_s[2];
    Out_Vap.H = Hvap_s[2];
    Out_Vap.x_c[:] = xvap_sc[2, :];
    En.Q = Q;
    
  // Adjustment to find the molecular weight
    x_pc[1, :] = (Fout .* xout_c[:] + Fvap_s[2] .* xvap_sc[2, :] + Fliq_s[2] .* xliq_sc[2, :]) / (Fliq_s[2] + Fvap_s[2] + Fout);
     x_pc[2, :] = xliq_sc[2,:];
     x_pc[3, :] = xvap_sc[2,:];
  
  //Density and Viscosity Calculation
    for i in 1:Nc loop
      muL[i] = Simulator.Files.TransportProperties.LiqVis(C[i].LiqVis, T);
      muV[i] = Simulator.Files.TransportProperties.VapVisc(C[i].VapVis, T);
      rholiq_c[i] = Simulator.Files.ThermodynamicFunctions.Dens(C[i].LiqDen, C[i].Tc, T, P);
      rhovap_c[i] = P * C[i].MW / (8.314 * T) ;
    end for;
    
  //Mixture Density and Viscosity
    muliq = sum(muL .* xliq_sc[2, :]);
    muvap = sum(muV .* xvap_sc[2, :]) ;
    rholiq = sum(rholiq_c .* xliq_sc[2, :]) * MW_p[2] * 1e-03;
    rhovap = sum(rhovap_c .* xvap_sc[2, :]) * 1e-03 ;
    
  //Calculation of Reynolds and Schmid Numbers
    ReL =  dp * rholiq * vl / (muliq * (1 - epslon)) ;
    ScL =  muliq / (rholiq * DL) ;
    ReV =  dp * rhovap * vv / (muvap * (1 - epslon)) ;
    ScV =  muvap / (rhovap * DV) ;
  //Specific and effective area of packing
    ap = 6 * (1 - epslon) / dp ;
    ae / ap = (1 - exp((-1.45 * (rholiq * vl / (ap * muliq))^0.1 * (vl^2 * ap / 9.81)^(-0.05) * (rholiq * vl^2 / (ap * sigma))^0.2)));
    
  //Correlation to find Mass Transfer coefficients
    if Corrln == "Onda" then
     KL = (0.0051 / (ap * dp)^(-0.4)) * (muliq * 9.81 / rholiq)^0.33 * (rholiq * vl / (ap * muliq))^0.66 * ScL^(-0.5);
     KV = 2 * (DV / (ap * dp^2)) * (rhovap * vv / (ap * muvap))^0.7 * ScV^0.33;
   elseif Corrln == "Shulman" then
     KL = (DL / dp) * 25.1 * (ReL)^0.45 * (ScL)^0.5 ;
     KV = 1.195 * vv * (ReV)^(-0.36) * (ScV)^(-0.66);
   end if;
  //Finding the velocities
    Fliq_s[1] = sum(rholiq_c .* xliq_sc[1, :]) * vl * (3.14 / 4) * dsec ^ 2;
    Fvap_s[1] = (rhovap * 1e3 / MW_p[3]) * vv * (3.14 / 4) * dsec ^ 2;
    
  //Molar Component Balance
  for i in 1:Nc loop
    Nliq[i] = Nvap[i];
  //Liquid Flux
    Fin * xin_c[i] + Fliq_s[1] * xliq_sc[1, i] - Fliq_s[2] * xliq_sc[2, i] = - Nliq[i] * (3.14/4) * dsec^2 * z * ap;
  //Vapour Flux
    Fvap_s[1] * xvap_sc[1, i] - Fvap_s[2] * xvap_sc[2, i] =  Nvap[i] * (3.14/4) * dsec^2 * z * ap;
  end for;
  
  //Flux Equations
    Nliq =   KL * (xI[:] - xliq_sc[2, :]) * (rholiq * 1e3 / MW_p[2])  + xliq_sc[2, :] * NL;
    Nvap =   KV * (xvap_sc[2, :] - yI[:]) * (rhovap * 1e3 / MW_p[3]) + xvap_sc[2, :] * NV;
    
  //Equillibrium
    yI = K_c[:] .* xI;
  //summation equation
    sum(xliq_sc[2, :]) = 1;
    sum(xvap_sc[2, :]) = 1;
    sum(yI) = 1;
    sum(xI) = 1;
    
  // Enthalpy balance
    Fin * Hin + Fvap_s[1] * Hvap_s[1] + Fliq_s[1] * Hliq_s[1] = Fout * Hout + Fvap_s[2] * Hvap_s[2] + Fliq_s[2] * Hliq_s[2] + Q;
  //enthalpy calculation
    for i in 1:Nc loop
      Hliqout_c[i] = ThermodynamicFunctions.HLiqId(C[i].SH, C[i].VapCp, C[i].HOV, C[i].Tc, T);
      Hvapout_c[i] = ThermodynamicFunctions.HVapId(C[i].SH, C[i].VapCp, C[i].HOV, C[i].Tc, T);
    end for;
    Hliq_s[2] = sum(xliq_sc[2, :] .* Hliqout_c[:]) + Hres_p[2];
    Hvap_s[2] = sum(xvap_sc[2, :] .* Hvapout_c[:]) + Hres_p[3];
  //sidedraw calculation
    if OutType == "L" then
      xout_c[:] = xliq_sc[2, :];
    elseif OutType == "V" then
      xout_c[:] = xvap_sc[2, :];
    else
      xout_c[:] = zeros(Nc);
    end if;
    algorithm
    for i in 1:Nc loop
      MW_p[:] := MW_p[:] + C[i].MW * x_pc[:, i];
    end for;
    annotation(
      Diagram(coordinateSystem(extent = {{-100, -40}, {100, 40}})),
      Icon(coordinateSystem(extent = {{-100, -40}, {100, 40}})),
      __OpenModelica_commandLineOptions = "");
  end PackingSection;

  model Reb "Model of a reboiler used in distillation column"
    import Simulator.Files.*;
    parameter Integer Nc = 2 "Number of components";
    parameter ChemsepDatabase.GeneralProperties C[Nc];
    parameter Boolean Bin = false;
    Real P(unit = "Pa", min = 0, start = Pg) "Pressure";
    Real T(unit = "K", min = 0, start = Tg) "Temperature";
    Real Fin(unit = "mol/s", min = 0, start = Fg) "Feed molar flow";
    Real Hin(unit = "kJ/kmol", start = Htotg) "Feed molar enthalpy";
    Real Hliqin(unit = "kJ/kmol", start = Hliqg) "Inlet liquid molar enthalpy";
    Real xin_c[Nc](each unit = "-", each min = 0, each max = 1, start = xguess) "Feed components mole fraction";
    Real Fliqin(unit = "mol/s", min = 0, start = Fg) "Inlet liquid molar flow";
    Real xliqin_c[Nc](each unit = "-", each min = 0, each max = 1, start = xg) "Inlet liquid component mole fraction";
    Real Fout(unit = "mol/s", min = 0, start = Fg) "Side draw molar flow";
    Real Fvapout(unit = "mol/s", min = 0, start = Fvapg) "Outlet vapor molar flow";
    Real xout_c[Nc](each unit = "-", each min = 0, each max = 1, start = xg) "Side draw mole fraction";
    Real xvapout_c[Nc](each unit = "-", each min = 0, each max = 1, start = xvapg) "Outlet vapor component mole fraction";
    Real Hvapout(unit = "kJ/kmol", start = Hvapg) "Outlet vapor molar enthalpy";
    Real Hvapout_c[Nc](each unit = "kJ/kmol") "Outlet vapor component molar enthalpy";
    Real Q(unit = "W") "Heat load";
    Real Hout(unit = "kJ/kmol") "Side draw molar enthalpy";
    Real x_pc[3, Nc](each unit = "-", each min = 0, each max = 1, each start = 1 / (Nc + 1)) "Component mole fraction";
    Real Pdew(unit = "Pa", min = 0, start = sum(C[:].Pc) / Nc) "Dew point pressure";
    Real Pbubl(unit = "Pa", min = 0, start = sum(C[:].Pc) / Nc) "Bubble point pressure";
    replaceable Simulator.Files.Interfaces.matConn In(Nc = Nc) if Bin annotation(
      Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    replaceable Simulator.Files.Interfaces.matConn In_Dmy(Nc = Nc, P = 0, T = 0, x_pc = zeros(3, Nc), F = 0, H = 0, S = 0, xvap = 0) if not Bin annotation(
      Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.matConn Out(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.trayConn In_Liq(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {-50, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-50, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.trayConn Out_Vap(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {50, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {50, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.enConn En annotation(
      Placement(visible = true, transformation(origin = {100, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    extends Simulator.GuessModels.InitialGuess;
  equation
//connector equation
    if Bin then
      In.x_pc[1, :] = xin_c[:];
      In.H = Hin;
      In.F = Fin;
    else
      In_Dmy.x_pc[1, :] = xin_c[:];
      In_Dmy.H = Hin;
      In_Dmy.F = Fin;
    end if;
    Out.P = P;
    Out.T = T;
    Out.x_pc[1, :] = xout_c;
    Out.F = Fout;
    Out.H = Hout;
    In_Liq.F = Fliqin;
    In_Liq.H = Hliqin;
    In_Liq.x_c[:] = xliqin_c[:];
    Out_Vap.F = Fvapout;
    Out_Vap.H = Hvapout;
    Out_Vap.x_c[:] = xvapout_c[:];
    En.Q = Q;
//Adjustment for thermodynamic packages
    x_pc[1, :] = (Fout .* xout_c[:] + Fvapout .* xvapout_c[:]) ./ (Fout + Fvapout);
    x_pc[2, :] = xout_c[:];
//This equation is temporarily valid since this is only "partial" reboiler. Rewrite equation when "total" reboiler functionality is added
    x_pc[3, :] = xvapout_c[:];
//Bubble point calculation
    Pbubl = sum(gmabubl_c[:] .* x_pc[1, :] .* Pvap_c[:] ./ philiqbubl_c[:]);
//Dew point calculation
    Pdew = 1 / sum(x_pc[1, :] ./ (gmadew_c[:] .* Pvap_c[:]) .* phivapdew_c[:]);
//molar balance
//  Fin + Fliqin = Fout + Fvapout;
    for i in 1:Nc loop
      Fin .* xin_c[i] + Fliqin .* xliqin_c[i] = Fout .* xout_c[i] + Fvapout .* xvapout_c[i];
    end for;
//equillibrium
    xvapout_c[:] = K_c[:] .* xout_c[:];
//summation equation
//    sum(xvapout_c[:]) = 1;
    sum(xout_c[:]) = 1;
    for i in 1:Nc loop
      Hvapout_c[i] = Simulator.Files.ThermodynamicFunctions.HVapId(C[i].SH, C[i].VapCp, C[i].HOV, C[i].Tc, T);
    end for;
    Hvapout = sum(xvapout_c[:] .* Hvapout_c[:]) + Hres_p[3];
// bubble point calculations
    P = sum(xout_c[:] .* exp(C[:].VP[2] + C[:].VP[3] / T + C[:].VP[4] * log(T) + C[:].VP[5] .* T .^ C[:].VP[6]));
//    Fout = 10;
    Fin * Hin + Fliqin * Hliqin = Fout * Hout + Fvapout * Hvapout + Q;
    annotation(
      Diagram(coordinateSystem(extent = {{-100, -40}, {100, 40}})),
      Icon(coordinateSystem(extent = {{-100, -40}, {100, 40}})),
      __OpenModelica_commandLineOptions = "");
  end Reb;

  model DistCol "Model of a distillation column representing fractionating towers where mixture is separated in equilibrium stages"
    extends Simulator.Files.Icons.DistillationColumn;
    parameter Simulator.Files.ChemsepDatabase.GeneralProperties C[Nc] "Component instances array" annotation(
      Dialog(tab = "Column Specifications", group = "Component Parameters"));
    parameter Integer Nc "Number of components" annotation(
      Dialog(tab = "Column Specifications", group = "Component Parameters"));
    import data = Simulator.Files.ChemsepDatabase;
    parameter Boolean Bin_t[Nt] = Simulator.Files.OtherFunctions.colBoolCalc(Nt, Ni, InT_s) "Stream stage associations" annotation(
      Dialog(tab = "Column Specifications", group = "Component Parameters"));
    parameter Integer Nt = 4 "Number of stages" annotation(
      Dialog(tab = "Column Specifications", group = "Calculation Parameters"));
    parameter Integer Nout = 0 "Number of side draws" annotation(
      Dialog(tab = "Column Specifications", group = "Calculation Parameters"));
    parameter Integer NQ = 0 "Number of heat load" annotation(
      Dialog(tab = "Column Specifications", group = "Calculation Parameters"));
    parameter Integer Ni = 1 "Number of feed streams" annotation(
      Dialog(tab = "Column Specifications", group = "Calculation Parameters"));
    parameter Integer InT_s[Ni] "Feed stage location" annotation(
      Dialog(tab = "Column Specifications", group = "Calculation Parameters"));
    parameter String Ctype = "Total" "Condenser type: Total or Partial" annotation(
      Dialog(tab = "Column Specifications", group = "Calculation Parameters"));
    parameter String Corrln = "Onda" "Packing Correlation Shulman or Onda" annotation(
      Dialog(tab = "Column Specifications", group = "Calculation Parameters"));
    parameter Real DL = 5.6e-09 "Liquid Diffusivity" annotation(
      Dialog(tab = "Column Specifications", group = "Calculation Parameters"));
    parameter Real DV = 4.5e-06 "Vapour Diffusivity"annotation(
      Dialog(tab = "Column Specifications", group = "Calculation Parameters"));
    parameter Real sigma = 20e-03 "Liquid Surface Tension"annotation(
      Dialog(tab = "Column Specifications", group = "Calculation Parameters"));
    parameter Real z = 4 "Height of Packing"annotation(
      Dialog(tab = "Column Specifications", group = "Calculation Parameters"));
    parameter Real dsec(unit = "m") = 1 "Diameter of the packing section"annotation(
      Dialog(tab = "Column Specifications", group = "Calculation Parameters"));
    parameter Real dp(unit = "m") = 13e-03 "Diameter of Raschig Ring packing" annotation(
      Dialog(tab = "Column Specifications", group = "Calculation Parameters"));
    parameter Real epslon = 0.64 "Voidage"annotation(
      Dialog(tab = "Column Specifications", group = "Calculation Parameters"));
    Real RR(min = 0);
    Simulator.Files.Interfaces.matConn In_s[Ni](each Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {-248, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-250, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.matConn Dist(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {250, 316}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {250, 298}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.matConn Bot(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {250, -296}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {252, -300}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.enConn Cduty annotation(
      Placement(visible = true, transformation(origin = {246, 590}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {250, 600}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.enConn Rduty annotation(
      Placement(visible = true, transformation(origin = {252, -588}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {250, -598}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.matConn Out_s[Nout](each Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {-36, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-70, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.enConn En[NQ](each Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {-34, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-70, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    for i in 1:Ni loop
      if InT_s[i] == 1 then
        connect(In_s[i], condenser.In);
      elseif InT_s[i] == Nt then
        connect(In_s[i], reboiler.In);
      elseif InT_s[i] > 1 and InT_s[i] < Nt then
  //this is adjustment done since OpenModelica 1.11 is not handling array modification properly
        In_s[i].P = tray[InT_s[i] - 1].Pdmy1;
        In_s[i].T = tray[InT_s[i] - 1].Tdmy1;
        In_s[i].F = tray[InT_s[i] - 1].Fdmy1;
        In_s[i].x_pc = tray[InT_s[i] - 1].xdmy1_pc;
        In_s[i].H = tray[InT_s[i] - 1].Hdmy1;
        In_s[i].S = tray[InT_s[i] - 1].Sdmy1;
        In_s[i].xvap = tray[InT_s[i] - 1].xvapdmy1;
      end if;
    end for;
    connect(condenser.Out, Dist);
    connect(reboiler.Out, Bot);
    connect(condenser.En, Cduty);
    connect(reboiler.En, Rduty);
    for i in 1:Nt-3 loop
      connect(tray[i].Out_Liq, tray[i+1].In_Liq);
      connect(tray[i].In_Vap, tray[i+1].Out_Vap);
    end for;
    connect(tray[1].Out_Vap, condenser.In_Vap);
    connect(condenser.Out_Liq, tray[1].In_Liq);
    connect(tray[Nt-2].Out_Liq, reboiler.In_Liq);
    connect(reboiler.Out_Vap, tray[Nt-2].In_Vap);
  //tray pressures
    for i in 1:Nt-2 loop
      tray[i].P = condenser.P + i * (reboiler.P - condenser.P) / (Nt - 1);
    end for;
    for i in 2:Nt - 1 loop
      tray[i - 1].OutType = "Null";
      tray[i - 1].Out.x_pc = zeros(3, Nc);
      tray[i - 1].Out.F = 0;
      tray[i - 1].Out.H = 0;
      tray[i - 1].Out.S = 0;
      tray[i - 1].Out.xvap = 0;
      tray[i - 1].Q = 0;
    end for;
    RR = condenser.Fliqout / condenser.Out.F;
    annotation(
      Icon(coordinateSystem(extent = {{-250, -600}, {250, 600}})),
      Diagram(coordinateSystem(extent = {{-250, -600}, {250, 600}})),
      __OpenModelica_commandLineOptions = "",
      Documentation(info = "<html><head></head><body><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13.3333px; orphans: 2; widows: 2;\">The&nbsp;<b>Distillation Column</b>&nbsp;is used to separate the component mixture into component parts or fraction based on difference in volatalities.</span><div style=\"font-size: 12px;\"><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13.3333px; orphans: 2; widows: 2;\"><br></span></div><div style=\"font-size: 12px;\"><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13.3333px; orphans: 2; widows: 2;\">The distillation column model have following connection ports:</span></div><div><ol style=\"font-size: 12px;\"><li><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13.3333px; orphans: 2; widows: 2;\">Material Streams</span></li><ul><li><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13.3333px; orphans: 2; widows: 2;\">any number of feed stage</span></li><li><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13.3333px; orphans: 2; widows: 2;\">two products (distillate and bottom)</span></li></ul><li><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13.3333px; orphans: 2; widows: 2;\">Two Energy Streams</span></li><ul><li><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13.3333px; orphans: 2; widows: 2;\">condenser (total or partial)</span></li><li><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13.3333px; orphans: 2; widows: 2;\">reboiler</span></li></ul></ol><div style=\"font-size: 12px; orphans: 2; widows: 2;\"><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\"><br></span></font></div><div style=\"font-size: 12px; orphans: 2; widows: 2;\"><div><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">The results are:</span></font></div><div><ol><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Molar flow rate of Distillate and Bottoms</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Composition of Components in Distillate and Bottoms</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Condenser and Reboiler Duty</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Stagewise Liquid and Vapor Flows</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Temperature Profile</span></font></li></ol><div><br></div></div><div><br></div></div><div style=\"font-size: 12px; orphans: 2; widows: 2;\"><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">To simulate a distillation column, following calculation parameters must be provided:</span></font></div><div style=\"font-size: 12px; orphans: 2; widows: 2;\"><ol><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Number of Stages (<b>Nt</b>)</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Number of Feed Streams (<b>Ni</b>)</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Feed Tray Location (<b>InT_s</b>)</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Condenser Type (<b>Ctype</b>)</span></font></li></ol><div><span style=\"orphans: auto; widows: auto;\">All the variables are of type <i>parameter Real</i> except the last one (<b>Ctype</b>) which is of type&nbsp;<i>parameter String</i>. It can have either of the sting values among following:</span></div><div><ol><li><span style=\"orphans: auto; widows: auto;\"><b>Total</b>: To indicate that the condenser is Total Condenser</span></li><li><span style=\"orphans: auto; widows: auto;\"><b>Partial</b>: To indicate that the condenser is Partial Condenser</span></li></ol><span style=\"orphans: auto; widows: auto;\">During simulation, their values can specified directly under&nbsp;</span><b style=\"orphans: auto; widows: auto;\">Column Specifications</b><span style=\"orphans: auto; widows: auto;\">&nbsp;by double clicking on the column instance.</span></div><div><span style=\"orphans: auto; widows: auto;\"><br></span></div><div><span style=\"orphans: auto; widows: auto;\"><br></span></div></div><div style=\"font-size: 12px; orphans: 2; widows: 2;\"><span style=\"font-family: Arial, Helvetica, sans-serif; font-size: 13px;\">Additionally, following input for following variables must be provided:</span></div><div style=\"orphans: 2; widows: 2;\"><ol style=\"font-size: 12px;\"><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Condenser Pressure (<b>Pcond)</b></span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Reboiler Pressure (<b>Preb</b>)</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Top Specification</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Bottom Specification</span></font></li></ol><div><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Any one of the following variables can be considered as Top Specification:</span></font></div><div><ol><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Reflux Ratio (<b>RR</b>)</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Molar Flow rate (<b>F_p[1]</b>)</span></font></li></ol><div><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Any one of the following can be considered as Bottoms Specification:</span></font></div></div><div><ol><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Molar Flow rate (<b>F_p[1]</b>)</span></font></li><li><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">Mole Fraction of Component (<b>x_pc[1,:]</b>)</span></font></li></ol></div></div><div style=\"font-size: 12px; orphans: 2; widows: 2;\"><div><span style=\"orphans: auto; widows: auto;\">These variables are declared of type&nbsp;</span><i style=\"orphans: auto; widows: auto;\">Real</i><span style=\"orphans: auto; widows: auto;\">&nbsp;and therefore all these variables need to be declared in the equation section while simulating the model.</span></div><div><span style=\"orphans: auto; widows: auto;\"><br></span></div><div><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\"><br></span></font></div><div><span style=\"orphans: auto; widows: auto;\">For detailed explaination on how to use this model to simulate a Distillation Column,</span><font face=\"Arial, Helvetica, sans-serif\"><span style=\"font-size: 13px;\">&nbsp;go to&nbsp;<a href=\"modelica://Simulator.Examples.Distillation\">Distillation Column Example</a></span></font></div></div></div></body></html>"));
  end DistCol;
  annotation(
    Documentation(info = "<html><head></head><body>This is a package of models that have been developed to perfrom simulation of a Distillation Column. It contains following models:<div><br></div><div><ol><li><a href=\"modelica://Simulator.UnitOperations.DistillationColumn.Cond\">Cond</a>: Model of a Condenser used in Distillation Column</li><li><a href=\"modelica://Simulator.UnitOperations.DistillationColumn.DistTray\">DistTray</a>: Model of a Tray used in Distillation Column</li><li><a href=\"modelica://Simulator.UnitOperations.DistillationColumn.Reb\">Reb</a>: Model of a Reboiler used in Distillation Column</li><li><a href=\"modelica://Simulator.UnitOperations.DistillationColumn.DistCol\">DistCol</a>:&nbsp;Model of a distillation column representing fractionating towers where mixture is separated in equilibrium stages</li></ol></div></body></html>"));

end PackedColumn;
