function [X0,Y0,Vx,Vy,radii] =InputData(Bfile)
% Read the data from the input text file
formatSpec = '%s';
N = 4;
C_text = textscan(Bfile,formatSpec,6,'Delimiter',',');
C = textscan(Bfile,'%f %f %f %f %f %f','Delimiter',',');
fclose(Bfile);
%celldisp(C)
%
X0 = C{2};
Y0 = C{3};
Vx = C{4};
Vy = C{5};
radii  = C{6};