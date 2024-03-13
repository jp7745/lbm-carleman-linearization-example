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

function index = get_F_2_variable_index(node_x, ...
    lattice_vector_i, ...
    node_y, ...
    lattice_vector_j,...
    n,...
    Q)

%   NOTE: first-order variables are of the form:  f_i(x)
%   second-order variables are of the form:  f_i(x)*f_j(y)
%   WARNING!  MATLAB is Base1, so this will return `index` 
%   and `index` will be 1,2,3..., nQ*nQ
%   THIS MAPPING IS NOT UNIQUE.  JUST ONE EXAMPLE.
%   See equations 81a,b,c,d in the Applications of incompressible fluid
%   dynamics paper.


    index = n*Q*(Q*node_x + lattice_vector_i - Q) + (Q*node_y + lattice_vector_j - Q) - n*Q;

    % have to subtract some Q and nQ because MATLAB is base1.

end

