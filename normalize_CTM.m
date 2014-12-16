function [ output_tensor ] = normalize_CTM( input_tensor )
% Certain types of normalization seem necessary,
% since if not, CTM's either diverge or diminish.
% It should be similar to the canonical condition for iMPS tensors;
% Related to the normalization of the states.
% Dominant eigenvalue set equal to one.

% input_tensor(abs(input_tensor)<10^(-12))=0;

if max(abs(input_tensor(:)))<10^(-12)
    output_tensor=input_tensor;
else
    output_tensor=input_tensor/max(abs(input_tensor(:)));

% output_tensor=input_tensor/max(abs(input_tensor(:)));
    
%     if ismatrix(input_tensor)
%         % matrix normalization by SVD
%         output_tensor=input_tensor/svds(input_tensor,1,'L');
%         % tensor normalization by a largest modulus element
%         output_tensor=input_tensor/max(abs(input_tensor(:)));
%     elseif ndims(input_tensor)==4
%         % tensor normalization by a largest modulus singular value
%         % input_matrix <=> bonds to center bonds to C matrices
%         input_matrix=reshape(input_tensor,...
%             [size(input_tensor,1)^2,size(input_tensor,3)^2]);
%         output_tensor=input_tensor/svds(input_matrix,1,'L');
%         % tensor normalization by a largest modulus element
%         output_tensor=input_tensor/max(abs(input_tensor(:)));
%     else
%         error('Incorrect input tensor.')
%     end
end
