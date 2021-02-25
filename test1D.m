%test simple 1-d random walk
clc
close
clear 
dbstop if error

N=8;
A=diag([0,repmat(0.5,1,N-1)],1)...
+diag([1,repmat(-1,1,N-1),1],0)...
+diag([repmat(0.5,1,N-1),0],-1);

B=[0,repmat(0,1,N-1),1]';

p=A\B;