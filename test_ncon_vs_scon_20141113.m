check_sum=0;
for i=1:100
    index_list=randi([2,10],1,5);
    tensor_A=randn(index_list(1:4));
    tensor_B=randn(index_list(2:5));
    difference_tensor=...
        ncon({tensor_A,tensor_B},{[-1,1,2,3],[1,2,3,-2]})-...
        scon({tensor_A,tensor_B},{[-1,1,2,3],[1,2,3,-2]});
    check_sum=check_sum+sum(abs(difference_tensor(:)));
end
check_sum

check_sum=0;
for i=1:100
    index_list=randi([2,10],1,5);
    tensor_A=randn(index_list(1:4));
    tensor_B=randn([randi(10),index_list(2:3),randi(10)]);
    difference_tensor=...
        ncon({tensor_A,tensor_B},{[-1,1,2,-3],[-2,1,2,-4]})-...
        scon({tensor_A,tensor_B},{[-1,1,2,-3],[-2,1,2,-4]});
    check_sum=check_sum+sum(abs(difference_tensor(:)));
end
check_sum