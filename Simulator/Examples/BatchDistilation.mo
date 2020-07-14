within Simulator.Examples;

package BatchDistilation

extends Modelica.Icons.ExamplesPackage;

  model Condenser
    extends Simulator.UnitOperations.BatchDistillation.Cond;
    extends Simulator.Files.ThermodynamicPackages.RaoultsLaw;
  end Condenser;

  model Tray
    extends Simulator.UnitOperations.BatchDistillation.DistTray;
    extends Simulator.Files.ThermodynamicPackages.RaoultsLaw;
  end Tray;

  model Reboiler
    extends Simulator.UnitOperations.BatchDistillation.Reb;
    extends Simulator.Files.ThermodynamicPackages.RaoultsLaw;
  end Reboiler;

  model DistColumn
    extends Simulator.UnitOperations.BatchDistillation.DistCol;
    Condenser condenser(Nc = Nc, C = C, Ctype = Ctype);
    Reboiler reboiler(Nc = Nc, C = C);
    Tray tray[Nt](each Nc = Nc, each C = C);
  end DistColumn;

  model ms
    extends Simulator.Streams.MaterialStream;
    extends Simulator.Files.ThermodynamicPackages.RaoultsLaw;
  end ms;

  model Test
    extends Modelica.Icons.Example;
    parameter Integer Nc = 2;
    import data = Simulator.Files.ChemsepDatabase;
    parameter data.Toluene tol;
    parameter data.Benzene ben;
    parameter Simulator.Files.ChemsepDatabase.GeneralProperties C[Nc] = {tol,ben};
    Simulator.Examples.BatchDistilation.DistColumn distCol( C = C, Ctype = "Total", Nc = Nc, Nt = 3, condenser.RR = 2, reboiler.xb(start = {0.5, 0.5})) annotation(
      Placement(visible = true, transformation(origin = {-22, 8}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Distillation.ms distillate( C = C, Nc = Nc) annotation(
      Placement(visible = true, transformation(origin = {48, 38}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Streams.EnergyStream cond_duty annotation(
      Placement(visible = true, transformation(origin = {36, 68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Simulator.Streams.EnergyStream reb_duty annotation(
      Placement(visible = true, transformation(origin = {48, -52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(distCol.Cduty, cond_duty.In) annotation(
      Line(points = {{3, 68}, {26, 68}}));
    connect(distCol.Dist, distillate.In) annotation(
      Line(points = {{3, 38}, {38, 38}}));
    connect(distCol.Rduty, reb_duty.In) annotation(
      Line(points = {{3, -52}, {38, -52}}));
  
    distCol.condenser.P = 1e5;
    distCol.reboiler.P = 1.5e5;
    distCol.reboiler.Fvapout = 100;
  end Test;
end BatchDistilation;
