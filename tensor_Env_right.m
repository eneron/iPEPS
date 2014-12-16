function [ tensor_Env ] = tensor_Env_right...
    (tensor_A,tensor_B,tensor_C,tensor_T)

CHI=size(tensor_C{3},1);
BOND_DIM=size(tensor_T{5},1);

tensor_a=ncon({tensor_A,conj(tensor_A)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
tensor_b=ncon({tensor_B,conj(tensor_B)},...
    {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});

%% RG with the 1st inserted column
% Absorption
c3_tilde=ncon({tensor_C{3},tensor_T{4}},{[1,-2],[-3,-4,-1,1]});
c4_tilde=ncon({tensor_C{4},tensor_T{7}},{[-1,1],[-2,-3,1,-4]});
t5_tilde=ncon({tensor_T{5},tensor_a},...
    {[1,2,-3,-6],[-1,-7,1,-4,-2,-8,2,-5]});
t6_tilde=ncon({tensor_T{6},tensor_b},...
    {[1,2,-3,-6],[-1,-7,1,-4,-2,-8,2,-5]});

% Renormalization applying isometries

% Calculate Z-tensor
c3_tilde_matrix=reshape(c3_tilde,[CHI,CHI*BOND_DIM^2]);
c4_tilde_matrix=reshape(c4_tilde,[CHI*BOND_DIM^2,CHI]);
RG_matrix_Z=...
    c4_tilde_matrix*c4_tilde_matrix'+...
    c3_tilde_matrix.'*conj(c3_tilde_matrix);
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

% Obtain renormalized tensors C3', C4'
c4_prime=normalize_CTM(z_matrix'*c4_tilde_matrix);
c3_prime=normalize_CTM(c3_tilde_matrix*z_matrix);

% Calculate W-tensor
q3_tilde=ncon({c3_tilde,t5_tilde},{[-1,1,2,3],[-2,-3,1,2,3,-4,-5,-6]});
q4_tilde=ncon({c4_tilde,t6_tilde},{[1,2,3,-4],[-5,-6,-1,-2,-3,1,2,3]});

q3_matrix=reshape(q3_tilde,[CHI*BOND_DIM^2,CHI*BOND_DIM^2]);
q4_matrix=reshape(q4_tilde,[CHI*BOND_DIM^2,CHI*BOND_DIM^2]);

RG_matrix_Q=q4_matrix*q4_matrix'+q3_matrix.'*conj(q3_matrix);
[w_matrix,~]=eig_descend(RG_matrix_Q,CHI);
w_tensor=reshape(w_matrix,[CHI,BOND_DIM,BOND_DIM,CHI]);
w_tensor_dagger=reshape(w_matrix',[CHI,CHI,BOND_DIM,BOND_DIM]);

% Obtain renormalized tensors T5',T6'
t5_prime=normalize_CTM(ncon({w_tensor,t5_tilde,z_tensor_dagger},...
    {[4,5,6,-4],[-1,-2,1,2,3,4,5,6],[-3,1,2,3]}));
t6_prime=normalize_CTM(ncon({z_tensor,t6_tilde,w_tensor_dagger},...
    {[4,5,6,-4],[-1,-2,1,2,3,4,5,6],[-3,1,2,3]}));

% Calculating the environment tensor
tensor_Env=ncon({tensor_C{1},tensor_C{2},c3_prime,c4_prime,...
    tensor_T{1},tensor_T{2},tensor_T{3},...
    t5_prime,t6_prime,tensor_T{8}},...
    {[1,2],[4,5],[6,7],[9,10],...
    [-1,-2,2,3],[-7,-8,3,4],[-11,-12,5,6],...
    [-9,-10,7,8],[-5,-6,8,9],[-3,-4,10,1]});
end
