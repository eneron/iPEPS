function [ tensor_Env ] = tensor_Env_left...
    (tensor_A,tensor_B,tensor_C,tensor_T)

CHI=size(tensor_C{1},1);
BOND_DIM=size(tensor_T{1},1);

tensor_a=ncon({tensor_A,conj(tensor_A)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_b=ncon({tensor_B,conj(tensor_B)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});

%% RG with the 1st inserted column
% Absorption
c1_tilde=ncon({tensor_C{1},tensor_T{8}},{[1,-2],[-3,-4,-1,1]});
c2_tilde=ncon({tensor_C{2},tensor_T{3}},{[-1,1],[-2,-3,1,-4]});
t1_tilde=ncon({tensor_T{1},tensor_a},...
    {[1,2,-3,-6],[1,-4,-1,-7,2,-5,-2,-8]});
t2_tilde=ncon({tensor_T{2},tensor_b},...
    {[1,2,-3,-6],[1,-4,-1,-7,2,-5,-2,-8]});

% Renormalization applying isometries

% Calculate Z-tensor
c1_tilde_matrix=reshape(c1_tilde,[CHI,CHI*BOND_DIM^2]);
c2_tilde_matrix=reshape(c2_tilde,[CHI*BOND_DIM^2,CHI]);
RG_matrix_Z=...
    c2_tilde_matrix*c2_tilde_matrix'+...
    c1_tilde_matrix.'*conj(c1_tilde_matrix);
[z_matrix,~]=eig_descend(RG_matrix_Z,CHI); % "~" for redundunt outputs
z_tensor=reshape(z_matrix,[CHI,BOND_DIM,BOND_DIM,CHI]);
z_tensor_dagger=permute(conj(z_tensor),[4,1,2,3]);
% Check isometry of z_tensor:
%   ncon({z_tensor,z_tensor_dagger},{[1,2,3,-1],[-2,1,2,3]})-eye(CHI)
%   = z_matrix'*z_matrix-eye(CHI)
%   = ~O(1.0e-15)
% Check the validity of z_tensor_dagger:
%   isequal(permute(conj(z_tensor),[4,1,2,3]),z_tensor_dagger)
%   ans = 1

% Obtain renormalized tensors C1', C2'
c2_prime=normalize_CTM(z_matrix'*c2_tilde_matrix);
c1_prime=normalize_CTM(c1_tilde_matrix*z_matrix);

% Calculate W-tensor
q1_tilde=ncon({c1_tilde,t1_tilde},{[-1,1,2,3],[-2,-3,1,2,3,-4,-5,-6]});
q2_tilde=ncon({c2_tilde,t2_tilde},{[1,2,3,-4],[-5,-6,-1,-2,-3,1,2,3]});

q1_matrix=reshape(q1_tilde,[CHI*BOND_DIM^2,CHI*BOND_DIM^2]);
q2_matrix=reshape(q2_tilde,[CHI*BOND_DIM^2,CHI*BOND_DIM^2]);

RG_matrix_Q=q2_matrix*q2_matrix'+q1_matrix.'*conj(q1_matrix);
[w_matrix,~]=eig_descend(RG_matrix_Q,CHI);
w_tensor=reshape(w_matrix,[CHI,BOND_DIM,BOND_DIM,CHI]);
w_tensor_dagger=reshape(w_matrix',[CHI,CHI,BOND_DIM,BOND_DIM]);

% Obtain renormalized tensors T1',T2'
t1_prime=normalize_CTM(ncon({z_tensor_dagger,t1_tilde,w_tensor},...
    {[-3,1,2,3],[-1,-2,1,2,3,4,5,6],[4,5,6,-4]}));
t2_prime=normalize_CTM(ncon({w_tensor_dagger,t2_tilde,z_tensor},...
    {[-3,1,2,3],[-1,-2,1,2,3,4,5,6],[4,5,6,-4]}));

% Calculating the environment tensor
tensor_Env=ncon({c1_prime,c2_prime,tensor_C{3},tensor_C{4},...
    t1_prime,t2_prime,tensor_T{4},...
    tensor_T{5},tensor_T{6},tensor_T{7}},...
    {[1,2],[4,5],[6,7],[9,10],...
    [-1,-2,2,3],[-7,-8,3,4],[-11,-12,5,6],...
    [-9,-10,7,8],[-5,-6,8,9],[-3,-4,10,1]});
end
