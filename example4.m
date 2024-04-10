 
clc; clear; close;

obj = Mesh();

obj.nen = 4;
obj.ndm = 2;

s = sind(45);

% Nodes
obj.x(1,:)  = [ s,  s];  
obj.x(2,:)  = [-s,  s];  
obj.x(3,:)  = [-s, -s];
obj.x(4,:)  = [ s, -s];  
obj.x(5,:)  = [ 1,  0];   
obj.x(6,:)  = [ 0,  1];                   
obj.x(7,:)  = [-1,  0];                  
obj.x(8,:)  = [ 0, -1];                  
obj.x(9,:)  = [-33/25 * s,    s/2];                   
obj.x(10,:) = [ 33/25 * s,   -s/2];                    
obj.x(11,:) = [      -s/2,    s/2]; 
obj.x(12,:) = [       s/2,   -s/2]; 


% Elements
obj.ix(1,:) = [ 8    4   10   12];
obj.ix(2,:) = [12   10    5    1]; 
obj.ix(3,:) = [ 3    8   12   11]; 
obj.ix(4,:) = [11   12    1    6]; 
obj.ix(5,:) = [11    6    2    9]; 
obj.ix(6,:) = [ 7    3   11    9];

obj.numnp = length(obj.x);
obj.numel = length(obj.ix);

nodes = [1   2   3   4   5   6   7   8   9    10];
obj.setCoordinateSystem(nodes, 'POLAR', [ 0,  0]);
 
obj.split( 1,  5, 25);
obj.split( 2,  6, 25);
obj.split( 6,  1, 25);
obj.split( 2,  9, 25);
obj.split( 9,  7, 25);


obj.smoothInternal(50);

obj.plotMesh();
