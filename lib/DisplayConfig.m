function DisplayConfig

fprintf(1, 'displaying configuration\n');
figure(1);
A = gca;
pos = get(A, 'Position');
vue = get(A, 'View');

clf;
hold on;

grid on;
DisplayEndoscope;

DisplayRobot;

DisplayInstrument;

%DisplayPatient;
DisplayPatient_v2;

DisplayOrgan;

DisplayLocalizer;

drawnow;

GetTargetPosition;
GetInstrumentPosition;

A = gca;
set(A, 'View', vue);
set(A, 'Position', pos);

axis equal;
drawnow;

figure(2);
clf;
hold on;
DisplayImage;
drawnow;

fprintf(1, 'configuration displayed\n');
