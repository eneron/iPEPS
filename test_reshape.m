%% 'reshape' vs 'kron': vectors
% Correct index-ordering: reversed order to the mathematical convention
v1=randn(2,1)+1.0j*randn(2,1);
v2=randn(3,1)+1.0j*randn(3,1);
v3=randn(5,1)+1.0j*randn(5,1);
v=kron(kron(v1,v2),v3);
v_tensor=reshape(v,[5,3,2]); % index-order is reversed from math. conv.
check_sum_1=0;
check_sum_2=0;
for i=1:2 % index for v1
    for j=1:3 % index for v2
        for k=1:5 % index for v3
            % mathematical convention
            if v_tensor(k,j,i)~=v1(i)*v2(j)*v3(k)
                check_sum_1=check_sum_1+1;
            end
            if v_tensor(k,j,i)~=v(k+5*(j-1)+5*3*(i-1))
                % RHS is the same way of counting to the math. conv.
                check_sum_2=check_sum_2+1;
            end
            
        end
    end
end
check_sum_1
check_sum_2

% Wrong index-ordering sticking to the mathematical convention
v1=randn(2,1)+1.0j*randn(2,1);
v2=randn(3,1)+1.0j*randn(3,1);
v3=randn(5,1)+1.0j*randn(5,1);
v=kron(kron(v1,v2),v3);
v_tensor=reshape(v,[2,3,5]);
check_sum_1=0;
check_sum_2=0;
for i=1:2 % index for v1
    for j=1:3 % index for v2
        for k=1:5 % index for v3
            if v_tensor(i,j,k)~=v1(i)*v2(j)*v3(k)
                check_sum_1=check_sum_1+1;
            end
            if v_tensor(i,j,k)~=v(3*5*(i-1)+5*(j-1)+k)
                % RHS is the same way of counting to the math. conv.
                check_sum_2=check_sum_2+1;
            end
            
        end
    end
end
check_sum_1
check_sum_2
%% Indexing order in 'reshape' for the Kroneck Product of matrices
% Matlab ordering rule: increase the first index, then the second index, etc...
% ; Fortran convention.
% In Python Numpy, np.reshape works in the opposite way by default('C' convention)
% ; it can be changed the otherway by option 'F'.
% a=[1,2;3,4];
% b=[1,2,3;4,5,6;7,8,9];
a=randn(2,2);
b=randn(3,3);
c=kron(a,b);
c_tensor=reshape(c,[3,2,3,2]);
% % Note that the index-ordering in 'reshape' is different from the
% % conventional mathematical ordering as in the following:
% c_tensor=reshape(c,[2,3,2,3]);

% Correct indexing order in the 'reshape' function
check=0;
for i=1:2
    for j=1:2
        for k=1:3
            for l=1:3
                % A_{i,j}, B_{k,l}
                if c_tensor(k,i,l,j)~=a(i,j)*b(k,l)
                    check=check+1;
                end
            end
        end
    end
    
end
check

% Wrong indexing order in the 'reshape' function
a=randn(2,2);
b=randn(3,3);
c=kron(a,b);
c_tensor_2=reshape(c,[2,3,2,3]);
check=0;
for i=1:2
    for j=1:2
        for k=1:3
            for l=1:3
                if c_tensor_2(i,k,j,l)~=a(i,j)*b(k,l)
                    check=check+1;
                end
            end
        end
    end
end
check
%% Index-Ordering in Matlab and 'reshape'
% Matlab ordering rule: increase the first index, then the second index, etc...
% ; Fortran convention.
% In Python Numpy, np.reshape works in the opposite way by default('C' convention)
% ; it can be changed the otherway by option 'F'.

% tensor to matrix 1
a=randn(3,4,3,4);
m_a=reshape(a,[3*4,3*4]);
test_value_1=0;
for i=1:3
    for j=1:4
        for k=1:3
            for l=1:4
                if a(i,j,k,l)~=m_a(3*(j-1)+i,3*(l-1)+k)
                    test_value_1=test_value_1+1;
                end
            end
        end
    end
end
test_value_1
% matrix to tensor 1
m_b=randn(3*4,3*4);
b=reshape(m_b,[3,4,3,4]);
test_value_2=0;
for i=1:3
    for j=1:4
        for k=1:3
            for l=1:4
                if b(i,j,k,l)~=m_b(3*(j-1)+i,3*(l-1)+k)
                    test_value_2=test_value_2+1;
                end
            end
        end
    end
end
test_value_2
% matrix to tensor 2
m_b=randn(3*4,5);
b=reshape(m_b,[3,4,5]);
test_value_3=0;
for i=1:3
    for j=1:4
        for k=1:5
            if b(i,j,k)~=m_b(3*(j-1)+i,k)
                test_value_3=test_value_2+1;
            end
        end
    end
end
test_value_3
% tensor to matrix 2
a=randn(3,4,3,4,5);
m_a=reshape(a,[3*4,3*4*5]);
test_value_4=0;
for i=1:3
    for j=1:4
        for k=1:3
            for l=1:4
                for m=1:5
                    if a(i,j,k,l,m)~=m_a(3*(j-1)+i,3*4*(m-1)+3*(l-1)+k)
                        test_value_4=test_value_4+1;
                    end
                end
            end
        end
    end
end
test_value_4
%% Conclusion
% When you 'reshape' a vector or a matrix as a kronecker product of several
% subspace-elements, i.e. as a multi-index array, one should interpret the
% index-order in the reverse order as in the mathematical convention.