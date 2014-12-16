temp_m=randn(5);
temp_m=temp_m*temp_m';
[v0,d0]=eig(temp_m)
[v1,d1]=eigdescend(temp_m)
[v2,d2]=eigdescend(temp_m,3)
