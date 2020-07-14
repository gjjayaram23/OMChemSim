within Simulator.UnitOperations;

package BatchDistillation
extends Modelica.Icons.Package;

  model Cond "Model of a condenser used in distillation column"
  //==============================================================================
    import Simulator.Files.*;
    parameter ChemsepDatabase.GeneralProperties C[Nc];
    parameter Integer Nc = 2 "Number of components";
    parameter String Ctype "Condenser type: Partial or Total";
    parameter Real Cd(unit = "-") = 0.5 "Discharge Coefficient of Valve", A(unit = "m^2") = 3 "Reflux Drum Area";
    parameter Real RR(unit = "-") = 5 "Reflux Ratio";
  //==============================================================================
    Real P(unit = "K", min = 0, start = Pg) "Pressure";
    Real T(unit = "Pa", min = 0, start = Tg) "Temperature";
    Real xvapin_c[Nc](each unit = "-", each min = 0, each max = 1, start = xvapg) "Inlet components vapor molar fraction";
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
    Real M[Nc](each start = 500), ML(start = 2000), rholiq_c[Nc], h(start = 0.1);
  //===============================================================================
  
    Simulator.Files.Interfaces.matConn Out(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.trayConn Out_Liq(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {-50, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-50, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.trayConn In_Vap(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {50, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {50, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.enConn En annotation(
      Placement(visible = true, transformation(origin = {100, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, 14}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    extends Simulator.GuessModels.InitialGuess;

  equation
  //=====================================================================================
//Connector Equation
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

    for i in 1:Nc loop
    rholiq_c[i] = Simulator.Files.ThermodynamicFunctions.Dens(C[i].LiqDen, C[i].Tc, T, P);
    end for;
  //========================================================================================
  //Mole Balance
    Fvapin = Fout + RR*Fout + der(ML);
    for i in 1:Nc-1 loop
    Fvapin .* xvapin_c[i] - Fout .* xout_c[i] - Fliqout .* xliqout_c[i] = der(M[i]);
    end for;
    M = ML .* xliqout_c ;
    ML = sum(rholiq_c .* xliqout_c) * A * h;
    Fout = Cd * sqrt(2*9.81*h);
  //=========================================================================================
//Equillibrium
    if Ctype == "Partial" then
      xout_c[:] = K_c[:] .* xliqout_c[:];
    elseif Ctype == "Total" then
      xout_c[:] = xliqout_c[:];
    end if;
  //==========================================================================================
//Summation Equation
    sum(xout_c[:]) = 1;
  //========================================================================================
// Energy Balance
     Fvapin * Hvapin = Fout * Hout + Fliqout * Hliqout + Q ;
//Temperature calculation
    if Ctype == "Total" then
      P = sum(xout_c[:] .* Pvap_c[:]);
    elseif Ctype == "Partial" then
      1 / P = sum(xout_c[:] ./ Pvap_c[:]);
    end if;
  //========================================================================================
// Outlet liquid molar enthalpy calculation
    for i in 1:Nc loop
      Hliqout_c[i] = Simulator.Files.ThermodynamicFunctions.HLiqId(C[i].SH, C[i].VapCp, C[i].HOV, C[i].Tc, T);
    end for;
    Hliqout = sum(xliqout_c[:] .* Hliqout_c[:]) + Hres_p[2];
    annotation(
      Diagram(coordinateSystem(extent = {{-100, -40}, {100, 40}})),
      Icon(coordinateSystem(extent = {{-100, -40}, {100, 40}})),
      __OpenModelica_commandLineOptions = "");
  end Cond;

  model DistTray "Model of a tray used in distillation column"
  //=========================================================================================
    import Simulator.Files.*;
    parameter ChemsepDatabase.GeneralProperties C[Nc];
    parameter Integer Nc = 2 "Number of components";
    parameter Real A_Active(unit = "m^2") = 0.75 "Active Area of the tray";
    parameter Real hweir(unit = "m") = 0.1 "Weir Height";
    
  //=========================================================================================
    Real P(unit = "Pa", min = 0, start = Pg) "Pressure";
    Real T(unit = "K", min = 0, start = Tg) "Temperature";
    Real Fvap_s[2](each unit = "mol/s", each min = 0, start = {Fg, Fg}) "Vapor molar flow";
    Real Fliq_s[2](each unit = "mol/s", each min = 0, start = {Fg, Fg}) "Liquid molar flow";
    Real xvap_sc[2, Nc](each unit = "-", each min = 0, each max = 1, start = {yg, yg}) "Components vapor mole fraction";
    Real xliq_sc[2, Nc](each unit = "-", each min = 0, each max = 1, start = {xg, xg}) "Components liquid mole fraction";
    Real Hvap_s[2](unit = "kJ/kmol", start = Hvapg) "Vapor molar enthalpy";
    Real Hliq_s[2](unit = "kJ/kmol", start = Hliqg) "Liquid molar enthalpy";
    parameter Real Q(unit = "W") = 0  "Heat load";
    Real Hvapout_c[Nc](unit = "kJ/kmol", start = Hvapg) "Outlet components vapor molar enthalpy";
    Real Hliqout_c[Nc](unit = "kJ/kmol", start = Hliqg) "Outlet components liquid molar enthalpy";
    Real lweir(unit = "m") "Weir Length";
    Real M[Nc](each unit = "mol", each start = 1000) "Component Molar Holdup";
    Real ML(unit = "mol", start = 1000) "Molar Holdup of Liquid";
    Real rho_c[Nc](each unit = "mol/m^3") "Component Density";
    Real rho(unit = "mol/m^3") "Mixture Density";
    //=========================================================================================
    Simulator.Files.Interfaces.trayConn In_Liq(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {-50, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-50, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.trayConn Out_Liq(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {-50, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-50, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.trayConn Out_Vap(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {50, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {50, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.trayConn In_Vap(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {50, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {50, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    extends Simulator.GuessModels.InitialGuess;
    
  
  equation
  //=========================================================================================
  //Connector Equation
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
  
    for i in 1:Nc loop
  rho_c[i] = Simulator.Files.ThermodynamicFunctions.Dens(C[i].LiqDen, C[i].Tc, T, P);
    end for;
  //=========================================================================================
  //Molar balance
    Fvap_s[1] + Fliq_s[1] = Fvap_s[2] + Fliq_s[2]  ;
    for i in 1:Nc-1 loop
    Fvap_s[1] .* xvap_sc[1, i] + Fliq_s[1] .* xliq_sc[1, i] - Fvap_s[2] .* xvap_sc[2, i] - Fliq_s[2] .* xliq_sc[2, i] = der(M[i]);
    end for;
    rho =  sum(xliq_sc[2, :] .* rho_c) ;
    M = ML .* xliq_sc[2, :];
    lweir = (A_Active * 4 / 3.14)^0.5;
    ML = rho * A_Active * (hweir + 1.41*(Fliq_s[2] / (rho * lweir * 9.81))^0.66);
  //=========================================================================================
  //Equillibrium
    xvap_sc[2, :] = K_c[:] .* xliq_sc[2, :];
  //Summation Equation
  sum(xliq_sc[2, :]) = 1;
  sum(xvap_sc[2, :]) = 1;
  //=========================================================================================
  // Energy Balance
    Fvap_s[1] * Hvap_s[1] + Fliq_s[1] * Hliq_s[1] - Fvap_s[2] * Hvap_s[2] - Fliq_s[2] * Hliq_s[2] + Q  = 0;
  //Enthalpy calculation
    for i in 1:Nc loop
      Hliqout_c[i] = ThermodynamicFunctions.HLiqId(C[i].SH, C[i].VapCp, C[i].HOV, C[i].Tc, T);
      Hvapout_c[i] = ThermodynamicFunctions.HVapId(C[i].SH, C[i].VapCp, C[i].HOV, C[i].Tc, T);
    end for;
    Hliq_s[2] = sum(xliq_sc[2, :] .* Hliqout_c[:]) + Hres_p[2];
    Hvap_s[2] = sum(xvap_sc[2, :] .* Hvapout_c[:]) + Hres_p[3];
    annotation(
      Diagram(coordinateSystem(extent = {{-100, -40}, {100, 40}})),
      Icon(coordinateSystem(extent = {{-100, -40}, {100, 40}})),
      __OpenModelica_commandLineOptions = "");
  end DistTray;

  model Reb "Model of a reboiler used in distillation column"
  //=====================================================================================
    import Simulator.Files.*;
    parameter Integer Nc = 2 "Number of components";
    parameter ChemsepDatabase.GeneralProperties C[Nc];
    parameter Real VT(unit = "m^3") = 5 "Volume of the Holdup Tank";
  //=====================================================================================
    Real P(unit = "Pa", min = 0, start = Pg) "Pressure";
    Real T(unit = "K", min = 0, start = Tg) "Temperature";
    Real Hliqin(unit = "kJ/kmol", start = Hliqg) "Inlet liquid molar enthalpy";
    Real Fliqin(unit = "mol/s", min = 0, start = Fg) "Inlet liquid molar flow";
    Real xliqin_c[Nc](each unit = "-", each min = 0, each max = 1, start = xg) "Inlet liquid component mole fraction";
    Real Fvapout(unit = "mol/s", min = 0, start = Fvapg) "Outlet vapor molar flow";
    Real xb[Nc](each unit = "-", each min = 0, each max = 1, each start = 1/Nc) "Mole fraction at any instant";
    Real xvapout_c[Nc](each unit = "-", each min = 0, each max = 1, start = xvapg) "Outlet vapor component mole fraction";
    Real Hvapout(unit = "kJ/kmol", start = Hvapg) "Outlet vapor molar enthalpy";
    Real Hvapout_c[Nc](each unit = "kJ/kmol") "Outlet vapor component molar enthalpy";
    Real Q(unit = "W") "Heat load";
    Real M[Nc](each unit = "mol") "Component Molar Holdup";
    Real rholiq_c[Nc](each unit = "mol/m^3") "Component Molar Density";
    Real ML(unit = "mol", start = 100) "Molar Holdup of Liquid";
    
  //=====================================================================================
    Simulator.Files.Interfaces.trayConn In_Liq(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {-50, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-50, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.trayConn Out_Vap(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {50, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {50, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.enConn En annotation(
      Placement(visible = true, transformation(origin = {-4, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-2, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    extends Simulator.GuessModels.InitialGuess;
   
   initial equation
   der(xb) = zeros(Nc);
  equation
  
  //=====================================================================================
//Connector Equation
    In_Liq.F = Fliqin;
    In_Liq.H = Hliqin;
    In_Liq.x_c[:] = xliqin_c[:];
    Out_Vap.F = Fvapout;
    Out_Vap.H = Hvapout;
    Out_Vap.x_c[:] = xvapout_c[:];
    En.Q = Q;
  //=====================================================================================
//Mole Balance
    for i in 1:Nc loop
    rholiq_c[i] = Simulator.Files.ThermodynamicFunctions.Dens(C[i].LiqDen, C[i].Tc, T, P);
    end for;
    Fliqin - Fvapout = der(ML);
    for i in 1:Nc-1 loop
    Fliqin .* xliqin_c[i] - Fvapout .* xvapout_c[i] = der(M[i]) ;
    end for;
    M = ML .* xb;
    ML = sum(rholiq_c .* xb) * VT;
  //  Fvapout = 120;
//Equillibrium
    xvapout_c =  K_c .* xb;
//Summation Equation
    sum(xb) - sum(xvapout_c) = 0;
    for i in 1:Nc loop
      Hvapout_c[i] = Simulator.Files.ThermodynamicFunctions.HVapId(C[i].SH, C[i].VapCp, C[i].HOV, C[i].Tc, T);
    end for;
    Hvapout = sum(xvapout_c[:] .* Hvapout_c[:]) + Hres_p[3];
  //=====================================================================================
// Energy Balance
    Fliqin * Hliqin - Fvapout * Hvapout + Q = 0;
    P = sum(xb.* exp(C[:].VP[2] + C[:].VP[3] / T + C[:].VP[4] * log(T) + C[:].VP[5] .* T .^ C[:].VP[6]));
    annotation(
      Diagram(coordinateSystem(extent = {{-100, -40}, {100, 40}})),
      Icon(coordinateSystem(extent = {{-100, -40}, {100, 40}})),
      __OpenModelica_commandLineOptions = "");
  end Reb;

  model DistCol
  //=============================================================================
    extends Simulator.Files.Icons.DistillationColumn;
    parameter Simulator.Files.ChemsepDatabase.GeneralProperties C[Nc];
    parameter Integer Nc "Number of components";
    import data = Simulator.Files.ChemsepDatabase;
    parameter Integer Nt = 4 "Number of stages";
    parameter Integer Nout = 0 "Number of side draws";
    parameter Integer NQ = 0 "Number of heat load";
    parameter String Ctype = "Total" "Condenser type: Total or Partial";
  //=============================================================================
    Simulator.Files.Interfaces.matConn Dist(Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {250, 316}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {250, 298}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.enConn Cduty annotation(
      Placement(visible = true, transformation(origin = {246, 590}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {250, 600}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.enConn Rduty annotation(
      Placement(visible = true, transformation(origin = {252, -588}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {250, -598}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Files.Interfaces.enConn En[NQ](each Nc = Nc) annotation(
     Placement(visible = true, transformation(origin = {-34, -54}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-70, -60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  
  equation
  //=============================================================================
  //Connecting Equations
    connect(condenser.Out, Dist);
    connect(condenser.En, Cduty);
    connect(reboiler.En, Rduty);
    for i in 2:Nt loop
      connect(tray[i-1].Out_Liq, tray[i].In_Liq);
      connect(tray[i-1].In_Vap, tray[i].Out_Vap);
    end for;
    connect(tray[1].Out_Vap, condenser.In_Vap);
    connect(condenser.Out_Liq, tray[1].In_Liq);
    connect(tray[Nt].Out_Liq, reboiler.In_Liq);
    connect(reboiler.Out_Vap, tray[Nt].In_Vap);
  //=============================================================================
//Tray pressures
    for i in 1:Nt loop
      tray[i].P = condenser.P + i * (reboiler.P - condenser.P) / (Nt - 1);
    end for;
    annotation(
      Icon(coordinateSystem(extent = {{-250, -600}, {250, 600}})),
      Diagram(coordinateSystem(extent = {{-250, -600}, {250, 600}})),
      __OpenModelica_commandLineOptions = "");
  end DistCol;
end BatchDistillation;
