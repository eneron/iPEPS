function [ gate_tensor_TFIM ] =...
    gate_TFIM( J_ISING,H_FIELD,DELTA_T,varargin)
% [gate tensor generation: Transverse Field Ising model]
pauli_X=[0.0,1.0;1.0,0.0];
pauli_Y=[0.0,-1.0j;1.0j,0.0];
pauli_Z=[1.0,0.0;0.0,-1.0];

% ferromagnetic Ising Hamiltonian gate
% h_ij,kl: for the 1st particle, input 3, output 1
%          for the 2nd particle, input 4, output 2

% gate_temp=expm(-DELTA_T*...
%     (J_ISING*kron(pauli_Z,pauli_Z)+...
%     H_FIELD*kron(eye(2),pauli_Z)));

switch nargin
    case 3        
        gate_temp=expm(-DELTA_T*...
            (J_ISING*kron(pauli_X,pauli_X)+...
            H_FIELD/4*(kron(eye(2),pauli_Z)+kron(pauli_Z,eye(2)))));
        gate_tensor_TFIM=reshape(gate_temp,[2,2,2,2]);
        
    case 4
        if varargin{1}=='X'
            gate_temp=expm(-DELTA_T*...
                (J_ISING*kron(pauli_X,pauli_X)+...
                H_FIELD/4*(kron(eye(2),pauli_Z)+kron(pauli_Z,eye(2)))));
            gate_tensor_TFIM=reshape(gate_temp,[2,2,2,2]);
        elseif varargin{1}=='Y'
            gate_temp=expm(-DELTA_T*...
                (J_ISING*kron(pauli_Y,pauli_Y)+...
                H_FIELD/4*(kron(eye(2),pauli_Z)+kron(pauli_Z,eye(2)))));
            gate_tensor_TFIM=reshape(gate_temp,[2,2,2,2]);
        elseif varargin{1}=='Z'
            gate_temp=expm(-DELTA_T*...
                (J_ISING*kron(pauli_Z,pauli_Z)+...
                H_FIELD/4*(kron(eye(2),pauli_Z)+kron(pauli_Z,eye(2)))));
            gate_tensor_TFIM=reshape(gate_temp,[2,2,2,2]);
        else
            error('Incorrect Ising-interaction direction.')
        end
        
end