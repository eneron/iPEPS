function [ tensor_C,tensor_T ] =...
    CTMRG_vertical( tensor_Cell,tensor_C,tensor_T )
% CTMRG following 2010.Coroboz.PRB.82.245119
% 20141125.KIAS

BOND_DIM=size(tensor_T{1},1);
CHI=size(tensor_T{1},3);

tensor_a=ncon({tensor_Cell{1},conj(tensor_Cell{1})},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_b=ncon({tensor_Cell{2},conj(tensor_Cell{2})},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_c=ncon({tensor_Cell{3},conj(tensor_Cell{3})},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_d=ncon({tensor_Cell{4},conj(tensor_Cell{4})},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});

tensor_left=ncon({tensor_C{1},tensor_C{4},...
    tensor_T{1},tensor_T{6},tensor_T{7},tensor_T{8},...
    tensor_a,tensor_c},...
    {[2,1],[15,12],[3,4,1,-1],[13,14,-4,15],[10,11,12,7],[5,6,7,2],...
    [3,5,8,-2,4,6,9,-3],[8,10,13,-5,9,11,14,-6]});
tensor_right=ncon({tensor_C{2},tensor_C{3},...
    tensor_T{2},tensor_T{3},tensor_T{4},tensor_T{5},...
    tensor_b,tensor_d},...
    {[1,2],[12,15],[3,4,-1,1],[5,6,2,7],[10,11,7,12],[13,14,15,-4],...
    [3,-2,8,5,4,-3,9,6],[8,-5,13,10,9,-6,14,11]});

% Isometry u1, u2 construction
tensor_Q1=ncon({tensor_left,tensor_right},...
    {[-1,-2,-3,1,2,3],[-4,-5,-6,1,2,3]});
matrix_Q1=reshape(tensor_Q1,[CHI*BOND_DIM^2,CHI*BOND_DIM^2]);
[u1,s1,~]=svd(matrix_Q1);
u1=u1(:,1:CHI);
u1_tensor=reshape(u1,[CHI,BOND_DIM,BOND_DIM,CHI]);
u1_tensor_dagger=reshape(u1',[CHI,CHI,BOND_DIM,BOND_DIM]);
projection_1=reshape(u1*u1',[CHI,BOND_DIM,BOND_DIM,CHI,BOND_DIM,BOND_DIM]);

tensor_Q2=ncon({tensor_left,tensor_right},...
    {[1,2,3,-1,-2,-3],[1,2,3,-4,-5,-6]});
matrix_Q2=reshape(tensor_Q2,[CHI*BOND_DIM^2,CHI*BOND_DIM^2]);
[u2,s2,~]=svd(matrix_Q2);
u2=u2(:,1:CHI);
u2_tensor=reshape(u2,[CHI,BOND_DIM,BOND_DIM,CHI]);
u2_tensor_dagger=reshape(u2',[CHI,CHI,BOND_DIM,BOND_DIM]);
projection_2=reshape(u2*u2',[CHI,BOND_DIM,BOND_DIM,CHI,BOND_DIM,BOND_DIM]);

% Isometry u3, u4 construction
tensor_left2=ncon({tensor_left,projection_1,projection_2,...
    tensor_T{2},tensor_T{5},tensor_b,tensor_d},...
    {[1,2,3,4,5,6],[9,10,11,1,2,3],[14,15,16,4,5,6],...
    [12,13,-1,9],[17,18,-4,14],...
    [12,10,19,-2,13,11,20,-3],[19,15,17,-5,20,16,18,-6]});
tensor_right2=ncon({tensor_right,projection_1,projection_2,...
    tensor_T{1},tensor_T{6},tensor_a,tensor_c},...
    {[1,2,3,4,5,6],[1,2,3,9,10,11],[4,5,6,14,15,16],...
    [12,13,-1,9],[17,18,14,-4],...
    [12,-2,19,10,13,-3,20,11],[19,-5,17,15,20,-6,18,16]});

tensor_Q3=ncon({tensor_left2,tensor_right2},...
    {[-1,-2,-3,1,2,3],[-4,-5,-6,1,2,3]});
matrix_Q3=reshape(tensor_Q3,[CHI*BOND_DIM^2,CHI*BOND_DIM^2]);
[u3,s3,~]=svd(matrix_Q3);
u3=u3(:,1:CHI);
u3_tensor=reshape(u3,[CHI,BOND_DIM,BOND_DIM,CHI]);
u3_tensor_dagger=reshape(u3',[CHI,CHI,BOND_DIM,BOND_DIM]);

tensor_Q4=ncon({tensor_left2,tensor_right2},...
    {[1,2,3,-1,-2,-3],[1,2,3,-4,-5,-6]});
matrix_Q4=reshape(tensor_Q4,[CHI*BOND_DIM^2,CHI*BOND_DIM^2]);
[u4,s4,~]=svd(matrix_Q4);
u4=u4(:,1:CHI);
u4_tensor=reshape(u4,[CHI,BOND_DIM,BOND_DIM,CHI]);
u4_tensor_dagger=reshape(u4',[CHI,CHI,BOND_DIM,BOND_DIM]);

% Renormalization
tensor_C{1}=ncon({tensor_C{1},tensor_T{8},u3_tensor_dagger},...
    {[4,1],[2,3,-1,4],[-2,1,2,3]});
tensor_C{1}=normalize_CTM(tensor_C{1});

tensor_C{2}=ncon({tensor_C{2},tensor_T{3},u3_tensor},...
    {[1,4],[2,3,4,-2],[1,2,3,-1]});
tensor_C{2}=normalize_CTM(tensor_C{2});

tensor_C{3}=ncon({tensor_C{3},tensor_T{4},u4_tensor},...
    {[4,1],[2,3,-1,4],[1,2,3,-2]});
tensor_C{3}=normalize_CTM(tensor_C{3});

tensor_C{4}=ncon({tensor_C{4},tensor_T{7},u4_tensor_dagger},...
    {[1,4],[2,3,4,-2],[-1,1,2,3]});
tensor_C{4}=normalize_CTM(tensor_C{4});

tensor_T{1}=ncon({tensor_T{1},tensor_a,u3_tensor,u1_tensor_dagger},...
    {[4,5,1,6],[4,2,-1,7,5,3,-2,8],[1,2,3,-3],[-4,6,7,8]});
tensor_T{1}=normalize_CTM(tensor_T{1});

tensor_T{2}=ncon({tensor_T{2},tensor_b,u1_tensor,u3_tensor_dagger},...
    {[4,5,1,6],[4,2,-1,7,5,3,-2,8],[1,2,3,-3],[-4,6,7,8]});
tensor_T{2}=normalize_CTM(tensor_T{2});

tensor_T{6}=ncon({tensor_T{6},tensor_c,u4_tensor,u2_tensor_dagger},...
    {[4,5,1,6],[-1,7,4,2,-2,8,5,3],[6,7,8,-4],[-3,1,2,3]});
tensor_T{6}=normalize_CTM(tensor_T{6});

tensor_T{5}=ncon({tensor_T{5},tensor_d,u2_tensor,u4_tensor_dagger},...
    {[4,5,1,6],[-1,7,4,2,-2,8,5,3],[6,7,8,-4],[-3,1,2,3]});
tensor_T{5}=normalize_CTM(tensor_T{5});

tensor_temp=tensor_T{3};
tensor_T{3}=tensor_T{4};
tensor_T{4}=tensor_temp;

tensor_temp=tensor_T{8};
tensor_T{8}=tensor_T{7};
tensor_T{7}=tensor_temp;
end

