%test simple 3-d random walk
clc
close all
clear 
dbstop if error

%%
%遍历每个顶点的一邻域点，计算概率值，并归一化
[V,F]=readOBJ('潘羽萍_261711_lowertooth2.obj');
%drawMesh(V,F, 'facecolor','y','edgecolor','w');
%hold on
N=size(V,1);
P=zeros(N,N);

vring = compute_vertex_ring(F);
Nor = per_vertex_normals(V,F);
[A,C] = adjacency_dihedral_angle_matrix(V,F);

sumofring = 0;
sumofP=0;
Eta_total =cell(N,1);
b=cell(N,1);
n=cell(N,1);
for vertex_i=1:N
    Nring=length(vring{vertex_i});
    
    b1=repmat(V(vertex_i,:),Nring,1)-V(vring{vertex_i},:);
    b{vertex_i}=b1;
    
    n1=repmat(Nor(vertex_i,:),Nring,1)-Nor(vring{vertex_i},:);
    n{vertex_i}=n1;
   
    %%
    %判断边所在二面角是钝角还是锐角
    %method1 利用顶点法线求面法线
    if 1
    normface1 = normalizerow(Nor(circshift(vring{vertex_i},-1,2),:)+Nor(vring{vertex_i},:)+repmat(Nor(vertex_i,:),Nring,1));
    normface2 = normalizerow(Nor(circshift(vring{vertex_i},1,2),:)+Nor(vring{vertex_i},:)+repmat(Nor(vertex_i,:),Nring,1));
    D = dot(normface1,normface2,2);
    Eta=5*ones(Nring,1).*(D>0)+ 1*ones(Nring,1).*(D<0);
    Eta_total{vertex_i}=Eta;
    end
    
    %%
    %利用adjacency_dihedral_angle_matrix，需要到F中寻找当前面的索引， unfinished
   if 0
    Nf=size(F,1);
    Eta =ones(Nring,1);
    for j=1:Nring
    result=zeros(2,1);
    for i=1:Nf
        p=intersect(F(i,:),[vertex_i,vring{vertex_i}(j)]);
        if length(p)==2 && result(1)==0
            result(1)=i;
        end
        if length(p)==2 && result(1)>0 && i ~= result(1)
            result(2)=i;
        end
        if result(1)>0 && result(2)>0%必须要水密网格才行
            break;
        end
    end
    if A(result(1),result(2)) > pi
        Eta(j)=5;
    end
    end
   end
    Eta_total{vertex_i}=Eta;
    
    %%
    %利用混合积计算体积
    if 0
    Tet=[circshift(vring{vertex_i},-1,2)',vring{vertex_i}',repmat(vertex_i,Nring,1),circshift(vring{vertex_i},1,2)'];
    %Vol=calcVol(V,Tet(1,:));
    Vol=calcVol(V,Tet);
    Eta=1*ones(Nring,1).*(Vol>0)+ 5*ones(Nring,1).*(Vol<0);
    end
    
    originP=Eta.*dot(b1,n1,2);
    sumofP = sumofP+sum(originP);
    sumofring = sumofring+Nring;
end

bard = sumofP/sumofring;

for vertex_i=1:N    
    avgP=Eta_total{vertex_i}.*dot(b{vertex_i},n{vertex_i},2)./bard;
    w=exp(-avgP.*avgP);    
    kk=vring{vertex_i};
    P(vertex_i,kk)=w./sum(w);
end
%%
%每个顶点的概率方程构成最终方程组的一行

A=P+diag(-1*ones(N,1),0);

B=zeros(N,1);
%人工指定 前景背景的点的id
background =[674 2869 780 5285 2091 6607 6698 6534 6568 6420 6452 3200 6554 6328 6401 6092 6355 134];
foreground =[2136 2140 2696 1874 ];
B(background)=zeros(size(background,2),1);
A(background,:)=zeros(size(background,2),N);
A(background,background)=diag (ones(size(background,2),1),0);

B(foreground)=ones(size(foreground,2),1);
A(foreground,:)=zeros(size(foreground,2),N);
A(foreground,foreground)=diag(ones(size(foreground,2),1),0);

S = sparse(A) ;
pseg=S\B;

%%
%show
k0=double(pseg>0.5);
k1=double(pseg<0.5);

drawMesh(V,F, 'facecolor','g','edgecolor','flat','FaceVertexCData',k0);

     view(3);
     axis equal
     axis off
     camlight
     lighting gouraud
     set(gca, 'Position',[0 0 1 1]);
     cameratoolbar
     cameratoolbar('SetCoordSys','none');    
     
