%% Gate tensor verification
% Ising interaction only: Z
gate=gate_TFIM(-1.0,0.0,100.0,'Z');
state=randn(2,2)+1.0j*randn(2,2);
state=state/norm(state(:))
for i=1:100
    state=ncon({gate,state},{[-1,-2,1,2],[1,2]});
    state=state/norm(state);
end
state
% note that only the components a|00>+b|11> survive;
% a, b change according to the (random) state.

% Ising interaction only: X
gate=gate_TFIM(1.0,0.0,100.0,'X');
input_state=randn(2,2)+1.0j*randn(2,2);
input_state=input_state/norm(input_state(:))
output_state=input_state;
for i=1:100
    output_state=ncon({gate,output_state},{[-1,-2,1,2],[1,2]});
    output_state=output_state/norm(output_state);
end
% basis change to {|+>,|->}_X
output_state=ncon({[1,1;1,-1]/sqrt(2),[1,1;1,-1]/sqrt(2),output_state},...
     {[1,-1],[2,-2],[1,2]})

% Ising interaction only: Y
gate=gate_TFIM(-1.0,0.0,100.0,'Y');
input_state=randn(2,2)+1.0j*randn(2,2);
input_state=input_state/norm(input_state(:))
output_state=input_state;
for i=1:100
    output_state=ncon({gate,output_state},{[-1,-2,1,2],[1,2]});
    output_state=output_state/norm(output_state);
end
% basis change to {|+>,|->}_Y
output_state=ncon({[1,1;-1j,1j]/sqrt(2),[1,1;-1j,1j]/sqrt(2),output_state},...
    {[1,-1],[2,-2],[1,2]}) % Caution on the index order of transformation matrix.
%% Conclusion
% gate_TFIM(i,k,j,l) works as |i><j|x|k><l|=|ik><jl|.