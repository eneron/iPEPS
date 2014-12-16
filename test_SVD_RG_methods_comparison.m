% Two RG methods through SVD + Cutting are equal.
% Choose the efficient one that cuts singular values and unitary matrices,
% and then builds up the new tensors.
for i=1:10
    m=magic(5);
    [u,s,v]=svd(m);
    v=v';
    
    s_cut=s(1:3,1:3);
    u_cut=u(:,1:3);
    v_cut=v(1:3,:);
    m_cut=u_cut*s_cut*v_cut;
    
    u_prime=u*sqrt(s);
    v_prime=sqrt(s)*v;
    u_cut2=u_prime(:,1:3);
    v_cut2=v_prime(1:3,:);
    m_cut2=u_cut2*v_cut2;
    
    temp=m_cut-m_cut;
    norm(abs(temp(:)))
end