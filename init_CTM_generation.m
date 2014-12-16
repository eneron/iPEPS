function [ tensor_C,tensor_T ] = init_CTM_generation( BOND_DIM,CHI )
% generates returns two cell arrays of the initial random corner transfer
% matrices C's and T's, respectively (for each time step).
% Built C{i} real matrices and T{i} with special form to ensure Hermiticity
% of the resulting reduced density matrix.

tensor_C=cell(4);
tensor_T=cell(8);

for i=1:4
    tensor_C{i}=randn(CHI,CHI);
    %         tensor_C{i}=randn(CHI,CHI)+1j*randn(CHI,CHI);
    tensor_C{i}=normalize_CTM(tensor_C{i});
end
for i=1:8
    tensor_T{i}=randn(BOND_DIM,BOND_DIM,CHI,CHI)+...
        1j*randn(BOND_DIM,BOND_DIM,CHI,CHI);
    % For the reduced matrix to be Hermitian,
    % T{i}(j,k,:,:)=T{i}(j,k,:,:)^* for j!=k
    % T{i}(j,j,:,:)=T{i}(j,j,:,:)^*
    % What about positive semidefiniteness?
    for j=1:BOND_DIM
        for k=1:j-1
            tensor_T{i}(j,k,:,:)=conj(tensor_T{i}(k,j,:,:));
        end
        tensor_T{i}(j,j,:,:)=real(tensor_T{i}(j,j,:,:));
    end
    tensor_T{i}=normalize_CTM(tensor_T{i});
end
end

