function [ c2_prime,c3_prime,t3_prime,t4_prime ] =...
    CTMRG_right( tensor_C,tensor_T,tensor_a,tensor_b )
    % perform a right-ward CTMRG and return the updated tensors.
    % Keep separate indices separte rather than merging in order to adjust the order
    % of tensor contraction for computational complexity.
    
    CHI=size(tensor_C{2},1);
    BOND_DIM=size(tensor_T{3},1);
    
    %% RG with the 1st inserted column
    % Absorption
    c2_tilde=ncon({tensor_C{2},tensor_T{2}},{[1,-2],[-3,-4,-1,1]});
    c3_tilde=ncon({tensor_C{3},tensor_T{5}},{[-1,1],[-2,-3,1,-4]});
    t3_tilde=ncon({tensor_T{3},tensor_b},...
        {[1,2,-3,-6],[-4,-1,-7,1,-5,-2,-8,2]});
    t4_tilde=ncon({tensor_T{4},tensor_a},...
        {[1,2,-3,-6],[-4,-1,-7,1,-5,-2,-8,2]});
    
    % Renormalization applying isometries
    
    % Calculate Z-tensor
    c2_tilde_matrix=reshape(c2_tilde,[CHI,CHI*BOND_DIM^2]);
    c3_tilde_matrix=reshape(c3_tilde,[CHI*BOND_DIM^2,CHI]);
    RG_matrix_Z=...
        c3_tilde_matrix*c3_tilde_matrix'+...
        c2_tilde_matrix.'*conj(c2_tilde_matrix);
    RG_matrix_Z=cutoff(RG_matrix_Z,10^(-12));
    [z_matrix,~]=eig_descend(RG_matrix_Z,CHI); % "~" for redundunt outputs
    z_tensor=reshape(z_matrix,[CHI,BOND_DIM,BOND_DIM,CHI]);
    z_tensor_dagger=permute(conj(z_tensor),[4,1,2,3]);
    % Check isometry of z_tensor:
    %   ncon({z_tensor,z_tensor_dagger},{[1,2,3,-1],[-2,1,2,3]})-eye(CHI)
    %   = z_matrix'*z_matrix-eye(CHI)
    %   = ~O(1.0e-15)
    % Check the validity of z_tensor_dagger:
    %   isequal(permute(conj(z_tensor),[4,1,2,3]),z_tensor_dagger)
    %   ans = 1
    
    
    % Obtain renormalized tensors C2', C3'
    c3_intermediate=normalize_CTM(z_matrix'*c3_tilde_matrix);
    c2_intermediate=normalize_CTM(c2_tilde_matrix*z_matrix);
    
    % Calculate W-tensor
    q2_tilde=ncon({c2_tilde,t3_tilde},{[-1,1,2,3],[-2,-3,1,2,3,-4,-5,-6]});
    q3_tilde=ncon({c3_tilde,t4_tilde},{[1,2,3,-4],[-5,-6,-1,-2,-3,1,2,3]});

    q2_matrix=reshape(q2_tilde,[CHI*BOND_DIM^2,CHI*BOND_DIM^2]);
    q3_matrix=reshape(q3_tilde,[CHI*BOND_DIM^2,CHI*BOND_DIM^2]);
    
    RG_matrix_Q=q3_matrix*q3_matrix'+q2_matrix.'*conj(q2_matrix);
    RG_matrix_Q=cutoff(RG_matrix_Q,10^(-12));
    [w_matrix,~]=eig_descend(RG_matrix_Q,CHI);
    w_tensor=reshape(w_matrix,[CHI,BOND_DIM,BOND_DIM,CHI]);
    w_tensor_dagger=reshape(w_matrix',[CHI,CHI,BOND_DIM,BOND_DIM]);
    
    % Obtain renormalized tensors T3',T4'
    t3_intermediate=normalize_CTM(ncon({z_tensor_dagger,t3_tilde,w_tensor},...
        {[-3,1,2,3],[-1,-2,1,2,3,4,5,6],[4,5,6,-4]}));
    t4_intermediate=normalize_CTM(ncon({w_tensor_dagger,t4_tilde,z_tensor},...
        {[-3,1,2,3],[-1,-2,1,2,3,4,5,6],[4,5,6,-4]}));
    
    %% RG with the 2nd inserted column
    % Absorption
    c2_tilde=ncon({c2_intermediate,tensor_T{1}},{[1,-2],[-3,-4,-1,1]});
    c3_tilde=ncon({c3_intermediate,tensor_T{6}},{[-1,1],[-2,-3,1,-4]});
    t3_tilde=ncon({t3_intermediate,tensor_a},...
        {[1,2,-3,-6],[-4,-1,-7,1,-5,-2,-8,2]});
    t4_tilde=ncon({t4_intermediate,tensor_b},...
        {[1,2,-3,-6],[-4,-1,-7,1,-5,-2,-8,2]});
    
    % Renormalization applying isometries
    
    % Calculate Z-tensor
    c2_tilde_matrix=reshape(c2_tilde,[CHI,CHI*BOND_DIM^2]);
    c3_tilde_matrix=reshape(c3_tilde,[CHI*BOND_DIM^2,CHI]);
    RG_matrix_Z=...
        c3_tilde_matrix*c3_tilde_matrix'+...
        c2_tilde_matrix.'*conj(c2_tilde_matrix);
    RG_matrix_Z=cutoff(RG_matrix_Z,10^(-12));
    [z_matrix,~]=eig_descend(RG_matrix_Z,CHI); % "~" for redundunt outputs
    z_tensor=reshape(z_matrix,[CHI,BOND_DIM,BOND_DIM,CHI]);
    z_tensor_dagger=permute(conj(z_tensor),[4,1,2,3]);
    % Check isometry of z_tensor:
    %   ncon({z_tensor,z_tensor_dagger},{[1,2,3,-1],[-2,1,2,3]})-eye(CHI)
    %   = z_matrix'*z_matrix-eye(CHI)
    %   = ~O(1.0e-15)
    % Check the validity of z_tensor_dagger:
    %   isequal(permute(conj(z_tensor),[4,1,2,3]),z_tensor_dagger)
    %   ans = 1
    
    
    % Obtain renormalized tensors C2', C3'
    c3_prime=normalize_CTM(z_matrix'*c3_tilde_matrix);
    c2_prime=normalize_CTM(c2_tilde_matrix*z_matrix);
    
    % Calculate W-tensor
    q2_tilde=ncon({c2_tilde,t3_tilde},{[-1,1,2,3],[-2,-3,1,2,3,-4,-5,-6]});
    q3_tilde=ncon({c3_tilde,t4_tilde},{[1,2,3,-4],[-5,-6,-1,-2,-3,1,2,3]});

    q2_matrix=reshape(q2_tilde,[CHI*BOND_DIM^2,CHI*BOND_DIM^2]);
    q3_matrix=reshape(q3_tilde,[CHI*BOND_DIM^2,CHI*BOND_DIM^2]);
    
    RG_matrix_Q=q3_matrix*q3_matrix'+q2_matrix.'*conj(q2_matrix);
    RG_matrix_Q=cutoff(RG_matrix_Q,10^(-12));
    [w_matrix,~]=eig_descend(RG_matrix_Q,CHI);
    w_tensor=reshape(w_matrix,[CHI,BOND_DIM,BOND_DIM,CHI]);
    w_tensor_dagger=reshape(w_matrix',[CHI,CHI,BOND_DIM,BOND_DIM]);
    
    % Obtain renormalized tensors T3',T4'
    t3_prime=normalize_CTM(ncon({z_tensor_dagger,t3_tilde,w_tensor},...
        {[-3,1,2,3],[-1,-2,1,2,3,4,5,6],[4,5,6,-4]}));
    t4_prime=normalize_CTM(ncon({w_tensor_dagger,t4_tilde,z_tensor},...
        {[-3,1,2,3],[-1,-2,1,2,3,4,5,6],[4,5,6,-4]}));
end

