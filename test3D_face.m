%test simple 1-d random walk
clc
close
clear 
dbstop if error

%%
%遍历每个面的相邻面，计算概率值，并归一化
%如何获取每个面片的相邻面

[V,F]=readOBJ('潘羽萍_261711_lowertooth2.obj');
%drawMesh(V,F, 'facecolor','y','edgecolor','w');
%hold on
N=size(F,1);
P=zeros(N,N);

fring = compute_face_ring(F);
[normal,normalf] = compute_normal(V,F);
[A,C] = adjacency_dihedral_angle_matrix(V,F);
P = zeros(N,N);
%%
%每个面片的概率方程构成最终方程组的一行，4个非零元素
for face_i=1:N
    Nring=3;%length(fring{face_i});%==3
    dinternal=zeros(Nring,1);
    for i=1:Nring
        dihedral = A(face_i,fring{face_i}(i));
        if A > pi
            eta = 1;
        else
            eta =0.2;
        end
        dinternal(i)=eta*(1-cos(dihedral));
    end
    dinternal=dinternal./sum(dinternal);
    ptemp=zeros(Nring,1);
    for i=1:Nring
        commonedge = intersect(F(face_i,:),F(fring{face_i}(i),:));
        edge = V(commonedge(1),:)-V(commonedge(2),:);
        length = norm(edge);
        ptemp(i)=length*exp(-dinternal(i));
    end
    ptemp=ptemp./sum(ptemp);
    P(face_i,fring{face_i})=ptemp;
end

A=P+diag(-1*ones(N,1),0);

B=zeros(N,1);
%人工指定 前景背景的面的id
%%
background =[7100 8837 1244 892 3158 12854];
foreground =[ 4300 3762 4836 6038 3758 937 906 8273 11296 5492];
B(background)=zeros(size(background,2),1);
A(background,:)=zeros(size(background,2),N);
A(background,background)=diag (ones(size(background,2),1),0);

B(foreground)=ones(size(foreground,2),1);
A(foreground,:)=zeros(size(foreground,2),N);
A(foreground,foreground)=diag(ones(size(foreground,2),1),0);

S = sparse(A) ;
B1 = sparse(B) ;
pseg=S\B1;

%%
%show
k0=find(pseg>0.5);
k1=find(pseg<0.5);

drawMesh(V,F(k0,:), 'facecolor','g','edgecolor','y');
drawMesh(V,F(k1,:), 'facecolor','r','edgecolor','m');

     view(3);
     axis equal
     axis off
     camlight
     lighting gouraud
     set(gca, 'Position',[0 0 1 1]);
     cameratoolbar
     cameratoolbar('SetCoordSys','none');    