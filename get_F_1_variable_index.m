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


function index = get_F_1_variable_index(node_x, ...
    lattice_vector_i,...
    Q)

%   first-order variables are of the form:  f_i(x)
%   i = 1 , ..., Q
%   x = 1 , ..., n
%   index will be between 1,..., nQ
%   THIS MAPPING IS NOT UNIQUE.  JUST ONE EXAMPLE.
% 
    index = node_x*Q + lattice_vector_i - (Q);  % -(Q) because MATLAB base1.
end

