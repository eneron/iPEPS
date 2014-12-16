function [ check_sum ] = check_CTM( input_tensor )
% Check if CTM's satisfy the conditions; return '0' for proper cases.
check_sum=0;
if ismatrix(input_tensor)
    % Check that C{i}^*=C{i}
    temp_C=input_tensor-conj(input_tensor);
    temp_C(abs(temp_C)>10^(-12))=1;
    check_sum=sum(temp_C(:));
elseif ndims(input_tensor)==4
    BOND_DIM=size(input_tensor,1);
    CHI=size(input_tensor,3);
    % Check that T{i}(j,k,:,:)=T{i}(k,j,:,:)^* for j!=k
    % Check that T{i}(j,j,:,:)=T{i}(j,j,:,:)^*
    temp_T=zeros(BOND_DIM,BOND_DIM,CHI,CHI);
    for j=1:BOND_DIM
        for k=1:j-1
            temp_T(j,k,:,:)=...
                input_tensor(j,k,:,:)-conj(input_tensor(k,j,:,:));
        end
        temp_T(j,j,:,:)=...
            input_tensor(j,j,:,:)-conj(input_tensor(j,j,:,:));
    end
    temp_T(abs(temp_T)>10^(-12))=1;
    check_sum=sum(temp_T(:));
else
    error('Incorrect input tensor.')
end
end

