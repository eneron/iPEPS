%% Conclusion
% 20141210: It seems the order in 'reshape' during the construction of the isometry
% and the RG processes does not include mistakes...
%% test
BOND_DIM_1=2;
BOND_DIM_2=3;
CHI_1=5;
CHI_2=7;
c1_tilde=randn([CHI_1,BOND_DIM_1,BOND_DIM_2,CHI_2]);
c4_tilde=randn([CHI_2,CHI_1,BOND_DIM_1,BOND_DIM_2]);

c1_matrix=reshape(c1_tilde,[CHI_1*BOND_DIM_1*BOND_DIM_2,CHI_2]);
c4_matrix=reshape(c4_tilde,[CHI_2,CHI_1*BOND_DIM_1*BOND_DIM_2]);

c1_prime=reshape(c1_matrix,[CHI_1,BOND_DIM_1,BOND_DIM_2,CHI_2]);
c4_prime=reshape(c4_matrix,[CHI_2,CHI_1,BOND_DIM_1,BOND_DIM_2]);

sum(sum(sum(sum(abs(c1_tilde-c1_prime)))))
sum(sum(sum(sum(abs(c4_tilde-c4_prime)))))

check_sum=0;
for i=1:CHI_1
    for j=1:BOND_DIM_1
        for k=1:BOND_DIM_2
            for l=1:CHI_2
                if c1_tilde(i,j,k,l)~=...
                        c1_matrix(i+CHI_1*(j-1)+CHI_1*BOND_DIM_1*(k-1),l)
                    check_sum=check_sum+1;
                end
            end
        end
    end
end
check_sum

check_sum=0;
for i=1:CHI_2
    for j=1:CHI_1
        for k=1:BOND_DIM_1
            for l=1:BOND_DIM_2
                if c4_tilde(i,j,k,l)~=...
                        c4_matrix(i,j+CHI_1*(k-1)+CHI_1*BOND_DIM_1*(l-1))
                    check_sum=check_sum+1;
                end
            end
        end
    end
end
check_sum
