%% 1) Matrix A * Matrix B => Matrix A, Matrix B ? (Ans.: No)
a=randn(5);
b=randn(5);

[u,s,v]=svd(a*b);
a_re=u*sqrt(s);
b_re=sqrt(s)*v';

a-a_re
b-b_re
%% 2) A, B: SVD + reshape
PHYS_DIM=2;
BOND_DIM=3;

% initial preparation of (random) tensors A, B, a, b
[tensor_A,tensor_B,tensor_a,tensor_b]=...
    unit_cell_generation(PHYS_DIM,BOND_DIM);
matrix_A=reshape(tensor_A,[PHYS_DIM*BOND_DIM^3,BOND_DIM]);
% matrix_B=reshape(permute(tensor_B,[,[

tensor_AB=ncon({tensor_A,tensor_B},...
    {[-1,-2,-3,-4,1],[-5,-6,1,-8,-7]});

% Reshape into tensors A, B, lambda
matrix_AB=reshape(tensor_AB,[PHYS_DIM*BOND_DIM^3,PHYS_DIM*BOND_DIM^3]);
[temp_u,temp_s,temp_v]=svd(matrix_AB);
temp_v=temp_v';

rank(temp_s)
% temp_s(abs(temp_s)<10^(-14))=0;
temp_s=temp_s(1:rank(temp_s),1:rank(temp_s));

temp_u=temp_u(:,1:BOND_DIM);
temp_u_tensor=...
    reshape(temp_u,[PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM]);
tensor_A_re=ncon({temp_u_tensor,sqrt(temp_s)},...
    {[-1,-2,-3,-4,1],[1,-5]});

% check_sum=0;
tensor_A-tensor_A_re
sum(abs(tensor_A(:)-tensor_A_re(:)))
%% 3) A, B, lambda: SVD + reshape
PHYS_DIM=2;
BOND_DIM=3;

% initial preparation of (random) tensors A, B, a, b
[tensor_A,tensor_B,tensor_a,tensor_b]=...
    unit_cell_generation(PHYS_DIM,BOND_DIM);
matrix_A=reshape(tensor_A,[PHYS_DIM*BOND_DIM^3,BOND_DIM]);
% matrix_B=reshape(permute(tensor_B,[,[

lambda=sort(abs(randn(BOND_DIM,1)),'descend');
lambda=diag(lambda/norm(lambda));

tensor_AB=ncon({tensor_A,tensor_B,lambda},...
    {[-1,-2,-3,-4,1],[-5,-6,2,-8,-7],[1,2]});

% Reshape into tensors A, B, lambda
matrix_AB=reshape(tensor_AB,[PHYS_DIM*BOND_DIM^3,PHYS_DIM*BOND_DIM^3]);
[temp_u,temp_s,temp_v]=svd(matrix_AB);
temp_v=temp_v';
temp_u=temp_u(:,1:BOND_DIM);

tensor_A_re=reshape(temp_u,[PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM]);

rank(temp_s)
temp_s(abs(temp_s)<10^(-14))=0;
temp_s=temp_s(1:rank(temp_s),1:rank(temp_s));
lambda-temp_s
% check_sum=0;
sum(abs(tensor_A(:)-tensor_A_re(:)))
%% Conclusion
% As in the simplest case 1), once a product is built,
% (some) information about the original matrices is lost.
% So, one can't reconstruct the original matrices.
% For instance, consider that "SVD + reshape" strategy for A, B=I.