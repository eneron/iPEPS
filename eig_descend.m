function [ eigenvectors,eigenvalues ] = eig_descend( matrix,varargin )
% Return sorted eigen system (descending order in eigenvalues by default)
[temp_vec,temp_eig]=eig(matrix);
if ~isequal(diag(temp_eig),sort(diag(temp_eig),'descend'))
    % 'issorted' gives 1 if sorted in ascending order.
    % 'sort' sorts in ascending order by default.
    [eig_value_sorted,indices]=sort(diag(temp_eig),'descend');
    temp_vec=temp_vec(:,indices);
    temp_eig=diag(eig_value_sorted);
end
if nargin==2
    eigenvectors=temp_vec(:,1:varargin{1});
    eigenvalues=temp_eig(1:varargin{1},1:varargin{1});
else
    eigenvectors=temp_vec;
    eigenvalues=temp_eig;
end

end
