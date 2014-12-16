function [ tensor_Env ] = tensor_Env_up...
    (tensor_A,tensor_B,tensor_C,tensor_T)

CHI=size(tensor_C{1},1);
BOND_DIM=size(tensor_T{7},1);

tensor_a=ncon({tensor_A,conj(tensor_A)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_b=ncon({tensor_B,conj(tensor_B)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});

%% RG with the 1st inserted column
% Absorption
c1_tilde=ncon({tensor_C{1},tensor_T{1}},{[-1,1],[-2,-3,1,-4]});
c4_tilde=ncon({tensor_C{4},tensor_T{6}},{[1,-2],[-3,-4,-1,1]});
t8_tilde=ncon({tensor_T{8},tensor_a},...
    {[1,2,-3,-6],[-7,1,-4,-1,-8,2,-5,-2]});
t7_tilde=ncon({tensor_T{7},tensor_b},...
    {[1,2,-3,-6],[-7,1,-4,-1,-8,2,-5,-2]});

% Renormalization applying isometries

% Calculate Z-tensor
c1_tilde_matrix=reshape(c1_tilde,[CHI*BOND_DIM^2,CHI]);
c4_tilde_matrix=reshape(c4_tilde,[CHI,CHI*BOND_DIM^2]);
RG_matrix_Z=...
    c1_tilde_matrix*c1_tilde_matrix'+...
    c4_tilde_matrix.'*conj(c4_tilde_matrix);
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

% Obtain renormalized tensors C1', C4'
c1_prime=normalize_CTM(z_matrix'*c1_tilde_matrix);
c4_prime=normalize_CTM(c4_tilde_matrix*z_matrix);

% Calculate W-tensor
q1_tilde=ncon({c1_tilde,t8_tilde},{[1,2,3,-4],[-5,-6,-1,-2,-3,1,2,3]});
q4_tilde=ncon({c4_tilde,t7_tilde},{[-1,1,2,3],[-2,-3,1,2,3,-4,-5,-6]});

q1_matrix=reshape(q1_tilde,[CHI*BOND_DIM^2,CHI*BOND_DIM^2]);
q4_matrix=reshape(q4_tilde,[CHI*BOND_DIM^2,CHI*BOND_DIM^2]);

RG_matrix_Q=q1_matrix*q1_matrix'+q4_matrix.'*conj(q4_matrix);
[w_matrix,~]=eig_descend(RG_matrix_Q,CHI);
w_tensor=reshape(w_matrix,[CHI,BOND_DIM,BOND_DIM,CHI]);
w_tensor_dagger=reshape(w_matrix',[CHI,CHI,BOND_DIM,BOND_DIM]);
%     Check isometry of w_tensor:
%     ncon({w_tensor,w_tensor_dagger},{[1,2,3,-1],[-2,1,2,3]})-eye(CHI)
%     = O(1.0e-13)

% Obtain renormalized tensors T7',T8'
t7_prime=normalize_CTM(ncon({w_tensor,t7_tilde,z_tensor_dagger},...
    {[4,5,6,-4],[-1,-2,1,2,3,4,5,6],[-3,1,2,3]}));
t8_prime=normalize_CTM(ncon({z_tensor,t8_tilde,w_tensor_dagger},...
    {[4,5,6,-4],[-1,-2,1,2,3,4,5,6],[-3,1,2,3]}));

% Calculating the environment tensor
tensor_Env=ncon({c1_prime,tensor_C{2},tensor_C{3},c4_prime,...
    tensor_T{2},tensor_T{3},tensor_T{4},...
    tensor_T{5},t7_prime,t8_prime},...
    {[1,2],[3,4],[6,7],[8,9],...
    [-1,-2,2,3],[-5,-6,4,5],[-11,-12,5,6],...
    [-9,-10,7,8],[-7,-8,9,10],[-3,-4,10,1]});
end
