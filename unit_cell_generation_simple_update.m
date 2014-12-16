function [Gamma_A,Gamma_B,lambda] =...
    unit_cell_generation_simple_update( PHYS_DIM,BOND_DIM,varargin )
%   Generates and returns initial random tensors A and B, with the reduced tensors a
%   and b.
% tensor_A for the upper tensor, conj(tensor_A) for the lower tensor;
% (u,l,d,r)=(1,2,3,4), (5,6,7,8)

switch nargin
    case 2
        for i=1:4
            lambda{i}=sort(abs(randn(BOND_DIM,1)),'descend');
            lambda{i}=diag(lambda{i}/norm(lambda{i}));
        end
        Gamma_A_unnormalized=...
            randn(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM)+...
            1j*randn(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);
        Gamma_B_unnormalized=...
            randn(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM)+...
            1j*randn(PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM);
        % normalization by a largest modulus element to avoid divergence during RG
        Gamma_A=Gamma_A_unnormalized/max(abs(Gamma_A_unnormalized(:)));
        Gamma_B=Gamma_B_unnormalized/max(abs(Gamma_B_unnormalized(:)));
end
end

