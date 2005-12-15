% Print coordinates to file "mcoords" that can be read by 
% "vasputil importcoords". E.g. usage scenario is:
%
% First dump coordinates from POSCAR into stdout, redirect to file
% coords:
% % vasputil dumpcoords POSCAR > coords
% Enter octave/matlab, load coordinates:
% octave:1> load coords
% Manipulate coordinates
% print coords with this function:
% octave:42> printcoords(coords)
% Import coordinates into POSCAR file:
% % vasputil importcoords POSCAR < mcoords
function printcoords (coords)

[fd,msg] = fopen ("mcoords", "w");

if fd == -1,
  printf (msg);
end;

for i=1:size (coords, 1),
  fprintf(fd, "    %18.14f    %18.14f    %18.14f\n", coords(i,:));
end;

fclose (fd);
