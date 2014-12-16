function [ output ] = cutoff( input,EPS )
% Cut-off elements smaller than the EPS when normalized.
% 20141202: to deal with 'EIG did not converge.' during the isometry
% construction.
input_norm=norm(input);
output=input/input_norm;
output(abs(output)<EPS)=0.0;
output=output*input_norm;
end

