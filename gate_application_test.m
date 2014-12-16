function [tensor_A,tensor_B] = gate_application...
    (tensor_A,tensor_B,gate_tensor,direction)
    % Applies gates on tensors and update them by SVD and truncation.
    
    % Warning! Check that tensors reshaped from SVD as in the case of
    % simple_update; V' <=> V.' (20141031)
    
    PHYS_DIM=size(tensor_A,1);
    BOND_DIM=size(tensor_A,2);
%     CUTOFF=10^(-12);
    
    if strcmpi(direction,'rightward')
        % A rightward B
        tensor_temp=ncon({tensor_A,tensor_B,gate_tensor},...
            {[1,-2,-3,-4,3],[2,-8,3,-6,-7],[-1,-5,1,2]});
        tensor_temp_matrix=...
            reshape(tensor_temp,[PHYS_DIM*BOND_DIM^3,PHYS_DIM*BOND_DIM^3]);
        [temp_u,temp_s,temp_v]=svd(tensor_temp_matrix);
        temp_v=temp_v';
%         temp_s(abs(temp_s)<CUTOFF)=0;
        temp_u=temp_u*sqrt(temp_s);
        temp_v=sqrt(temp_s)*temp_v;

        % Cutting the bonds connecting the tensor A, B
        temp_u=temp_u(:,1:BOND_DIM);
        temp_v=temp_v(1:BOND_DIM,:);

        tensor_A=...
            reshape(temp_u,[PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM]);
        tensor_A=tensor_A/max(abs(tensor_A(:)));
        tensor_B=...
            reshape(temp_v,[BOND_DIM,PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM]);
        tensor_B=permute(tensor_B,[2,5,1,3,4]);
        tensor_B=tensor_B/max(abs(tensor_B(:)));

%         tensor_a=...
%             ncon({tensor_A,conj(tensor_A)},...
%             {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
%         tensor_b=...
%             ncon({tensor_B,conj(tensor_B)},...
%             {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
    elseif strcmpi(direction,'leftward')
        % A leftward B
        tensor_temp=ncon({tensor_B,tensor_A,gate_tensor},...
            {[1,-2,-3,-4,3],[2,-8,3,-6,-7],[-1,-5,1,2]});
        tensor_temp_matrix=...
            reshape(tensor_temp,[PHYS_DIM*BOND_DIM^3,PHYS_DIM*BOND_DIM^3]);
        [temp_u,temp_s,temp_v]=svd(tensor_temp_matrix);
        temp_v=temp_v';
%         temp_s(abs(temp_s)<CUTOFF)=0;
        temp_u=temp_u*sqrt(temp_s);
        temp_v=sqrt(temp_s)*temp_v;

        % Cutting the bonds connecting the tensor A, B
        temp_u=temp_u(:,1:BOND_DIM);
        temp_v=temp_v(1:BOND_DIM,:);

        tensor_B=...
            reshape(temp_u,[PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM]);
        tensor_B=tensor_B/max(abs(tensor_B(:)));
        tensor_A=...
            reshape(temp_v,[BOND_DIM,PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM]);
        tensor_A=permute(tensor_A,[2,5,1,3,4]);
        tensor_A=tensor_A/max(abs(tensor_A(:)));

%         tensor_a=...
%             ncon({tensor_A,conj(tensor_A)},...
%             {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
%         tensor_b=...
%             ncon({tensor_B,conj(tensor_B)},...
%             {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
    elseif strcmpi(direction,'upward')
        % A upward B
        tensor_temp=ncon({tensor_A,tensor_B,gate_tensor},...
            {[1,3,-2,-3,-4],[2,-7,-8,3,-6],[-1,-5,1,2]});
        tensor_temp_matrix=...
            reshape(tensor_temp,[PHYS_DIM*BOND_DIM^3,PHYS_DIM*BOND_DIM^3]);
        [temp_u,temp_s,temp_v]=svd(tensor_temp_matrix);
        temp_v=temp_v';
%         temp_s(abs(temp_s)<CUTOFF)=0;
        temp_u=temp_u*sqrt(temp_s);
        temp_v=sqrt(temp_s)*temp_v;

        % Cutting the bonds connecting the tensor A, B
        temp_u=temp_u(:,1:BOND_DIM);
        temp_v=temp_v(1:BOND_DIM,:);

        tensor_A=...
            reshape(temp_u,[PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM]);
        tensor_A=permute(tensor_A,[1,5,2,3,4]);
        tensor_A=tensor_A/max(abs(tensor_A(:)));
        tensor_B=...
            reshape(temp_v,[BOND_DIM,PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM]);
        tensor_B=permute(tensor_B,[2,4,5,1,3]);
        tensor_B=tensor_B/max(abs(tensor_B(:)));

%         tensor_a=...
%             ncon({tensor_A,conj(tensor_A)},...
%             {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
%         tensor_b=...
%             ncon({tensor_B,conj(tensor_B)},...
%             {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
    elseif strcmpi(direction,'downward')
        % A downward B
        tensor_temp=ncon({tensor_B,tensor_A,gate_tensor},...
            {[1,3,-2,-3,-4],[2,-7,-8,3,-6],[-1,-5,1,2]});
        tensor_temp_matrix=...
            reshape(tensor_temp,[PHYS_DIM*BOND_DIM^3,PHYS_DIM*BOND_DIM^3]);
        [temp_u,temp_s,temp_v]=svd(tensor_temp_matrix);
        temp_v=temp_v';
%         temp_s(abs(temp_s)<CUTOFF)=0;
        temp_u=temp_u*sqrt(temp_s);
        temp_v=sqrt(temp_s)*temp_v;

        % Cutting the bonds connecting the tensor A, B
        temp_u=temp_u(:,1:BOND_DIM);
        temp_v=temp_v(1:BOND_DIM,:);

        tensor_B=...
            reshape(temp_u,[PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM,BOND_DIM]);
        tensor_B=permute(tensor_B,[1,5,2,3,4]);
        tensor_B=tensor_B/max(abs(tensor_B(:)));
        tensor_A=...
            reshape(temp_v,[BOND_DIM,PHYS_DIM,BOND_DIM,BOND_DIM,BOND_DIM]);
        tensor_A=permute(tensor_A,[2,4,5,1,3]);
        tensor_A=tensor_A/max(abs(tensor_A(:)));

%         tensor_a=...
%             ncon({tensor_A,conj(tensor_A)},...
%             {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
%         tensor_b=...
%             ncon({tensor_B,conj(tensor_B)},...
%             {[1,-1,-2,-3,-4],[1,-5,-6,-7,-8]});
    else
        error('Incorrect direction.')
    end
end