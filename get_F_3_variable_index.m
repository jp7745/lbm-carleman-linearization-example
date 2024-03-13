%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LICENSE:
%
% Copyright 2024 L3Harris Technologies, Inc.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

function index = get_F_3_variable_index(node_x, ...
    lattice_vector_i, ...
    node_y, ...
    lattice_vector_j,...
    node_z, ...
    lattice_vector_k,...
    n,...
    Q)

%   NOTE: first-order variables are of the form:  f_i(x)
%   second-order variables are of the form:  f_i(x)*f_j(y)
%   third-order variables are of the form:  f_i(x)*f_j(y)*f_k(z)
%   WARNING!  MATLAB is Base1, so this will return `index` 
%   and `index` will be 1,2,3..., (nQ)^3
%   THIS MAPPING IS NOT UNIQUE.  JUST ONE EXAMPLE.
%   See equations 83a,b,c in the Applications of incompressible fluid
%   dynamics paper.


    index = (n*Q)^2*(Q*node_x + lattice_vector_i - Q) ...
        + n*Q*(Q*node_y + lattice_vector_j - Q) + (Q*node_z + lattice_vector_k - Q) - n*Q...
        - (n*Q)^2; 

        % have to subtract some Q, nQ, and (nQ)^2 because MATLAB base1.

end

