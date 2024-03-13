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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE:  
% 
%
% Considering the approximation the Lattice Boltzmann Equation:
%
% df/dt = (S)f + (F_1)f + (F_2)fxf + (F_3)fxfxf
%
% where F_hat = S + F_1...
%
% calculate the Carleman Linearized matrix A where:
% 
% A =   [ F_hat         F_2         F_3         ]
%       [ 0             F_hat^(2)   F_2^(2)     ]
%       [ 0             0           F_hat^(3)   ]
% 
% This uses a simple topology of n grid nodes in a line.  Node 1 wraps
% around to node n.  
%
% This uses the Lattice Boltzmann Method topology of D1Q3: the
% dimension is 1 and the number of lattice vectors for each node is 3.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
rng(7);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% color map for plots
clim([-20,20])
my_colormap = [ ones(100,1),zeros(100,2) ; 1,1,1 ; zeros(100,2),ones(100,1)];







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters

% number of grid nodes
n = 5;

% number of lattice edges (per grid node)
Q = 3;

% number of f_i(x) population distributions
nQ = n*Q;

% tau is a relaxation parameter chosen to match the physical properties of
% the system including viscosity and time/space discretization.
% tau should be 0.5 < tau <= 1.0
tau = 0.75;

% D1Q3 vectors (fixed for D1Q3 topology)
c = [ 0 , 1.0 , -1.0];

% D1Q3 weights (fixed for D1Q3 topology)
w = [ 2.0/3 , 1.0/6 , 1.0/6 ];







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the variable indexes.


% monomial variables of the form:  f_i(x)
F_1_var_ids = zeros(nQ, 2);

for node_x = 1 : n
    for i = 1 : Q
        index = get_F_1_variable_index(node_x, i, Q);
        F_1_var_ids(index, :) = [ node_x, i ];
    end
end

% monomial variables of the form:  f_i(x)*f_j(x)
F_2_var_ids = zeros(nQ*nQ, 4);
for node_x = 1 : n
    for i = 1 : Q
        for node_y = 1 : n
            for j = 1 : Q
                index = get_F_2_variable_index(node_x, ...
                    i, ...
                    node_y, ...
                    j, ...
                    n, ...
                    Q);
                F_2_var_ids(index,:) = [ node_x, i, node_y, j ];
            end
        end
    end
end

% monomial variables of the form:  f_j(x)*f_k(x)*f_ell(x)           
F_3_var_ids = zeros(nQ^3, 6);
for node_x = 1 : n
    for i = 1 : Q
        for node_y = 1 : n
            for j = 1 : Q
                for node_z = 1 : n
                    for k = 1 : Q
                        index = get_F_3_variable_index(node_x,...
                            i, ...
                            node_y,...
                            j,...
                            node_z,...
                            k,...
                            n,...
                            Q);
                        F_3_var_ids(index,:) = [node_x, i, node_y, j , node_z, k];
                    end
                end
            end
        end
    end
end










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S matrix
% S is effectively adding inbound streaming population and
% subtracting the current population (subtraction of the identity matrix)

S_inbound_streaming = zeros(nQ,nQ);
for node_x = 1 : n
    for i = 1 : Q
    
        inbound_node_x = node_x - c(i);

        % hacky modulo arithemetic with MATLAB since MATLAB is base1.
        if inbound_node_x == 0
            inbound_node_x = n;
        end
        if inbound_node_x == n+1
            inbound_node_x = 1;
        end

        row = node_x*Q + i - (Q); % -(Q) because MATLAB base1.
        col = inbound_node_x*Q + i - (Q); % -(Q) because MATLAB base1.

        S_inbound_streaming(row,col) = 1;
    end
end
S = S_inbound_streaming - eye(nQ, nQ);
figure(1);
imagesc(S);
title(sprintf("S matrix. size is (nQ,nQ)=(%d,%d), where n=%d, Q=%d",nQ,nQ,n,Q));
clim([-20,20])
colormap(my_colormap)
colorbar('Ticks',[-10,0,10],...
         'TickLabels',{'negative','0','positive'})
saveas(gcf, "figure--S_matrix.png");



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F_1 matrix

% see equations (67), (68) and (69) from the "Applications of
% Incompressible Computational Fluid Dynamics"

F_1 = zeros(nQ, nQ);

for node_x = 1 : n
    for i = 1 : Q
        for j = 1 : Q
            row = node_x*Q + i - (Q); % -(Q) because MATLAB base1.
            col = node_x*Q + j - (Q); % -(Q) because MATLAB base1.

            % init
            T = 0.0;
            
            % contribution from (67)
            if i == j
                T = T + (-1/tau);
            end

            % contribution from (68)
            T = T + w(i)/tau;

            % contribution from (69)
            T = T + (3*w(i)/tau)*c(i)*c(j);


            F_1(row,col) = T;
        end
    end
end



figure(2);
imagesc(F_1);
title(sprintf("F_1 matrix.. size is (nQ,nQ)=(%d,%d), where n=%d, Q=%d",nQ,nQ,n,Q));
clim([-20,20])
colormap(my_colormap)
colorbar
saveas(gcf, "figure--F_1_matrix.png");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F_hat is an nQ by nQ matrix that operates on the first-order
% variables.

F_hat = S + F_1;

figure(3);
imagesc(F_hat);
title(sprintf("F^{hat} = S + F_1 matrix. size is (nQ,nQ)=(%d,%d), where n=%d, Q=%d",nQ,nQ,n,Q));
clim([-20,20])
colormap(my_colormap)
colorbar('Ticks',[-10,0,10],...
         'TickLabels',{'negative','0','positive'})
saveas(gcf, "figure--F_hat_matrix.png");





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F_2 is an nQ by (nQ)^2 matrix that operates on the second order 
% variables.  I.e., variables of the form f_j(x_1)*f_k(x_2)


F_2 = zeros(nQ, (nQ)^2);

for node_x = 1 : n
    for i = 1 : Q
        row = get_F_1_variable_index(node_x, i, Q);
        % row corresponds to f_i(x)

        for j = 1 : Q
            for k = 1 : Q
        
                % collision terms only happen on the same node_x... 
                % so anywhere where node_x != node_y will be zero.
                % so we only lookup and calculate elements where node_x ==
                % node_x.
                % however!  lattice velocity vector indexes i,j can be
                % different.

                col = get_F_2_variable_index(node_x, j, node_x, k, n, Q);
                % col corresponds to f_j(x)*f_k(x)
    
                % using the "DENSE" version of the formula in equation (108)
                T_3_2_dense = (9.0*w(i)/tau)*(c(i)*c(j))*(c(i)*c(k));
                T_4_2_dense = (-3.0*w(i)/tau)*(c(j)*c(k));
    
                % T_3_2_dense = 1;
                % T_4_2_dense = 1;

                
                F_2(row,col) = T_3_2_dense + T_4_2_dense;

                % when c vector is zero, these elements are zeroed out.
                % just highlighting them to confirm they are there.
                % if F_2(row,col) == 0
                %     F_2(row,col) = -5;
                % end

            end
        end

    end
end


figure(4);
imagesc(F_2);
title(sprintf("F_2 matrix. size is (nQ,(nQ)^2)=(%d,%d)",nQ,(nQ)^2));
clim([-20,20])
colormap(my_colormap)
colorbar('Ticks',[-10,0,10],...
         'TickLabels',{'negative','0','positive'})
saveas(gcf, "figure--F_2_matrix.png");







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F_3 is an nQ by (nQ)^3 matrix that operates on the second order 
% variables.  I.e., variables of the form f_j(x_1)*f_k(x_2)*f_ell(x_3)


F_3 = zeros(nQ, (nQ)^3);

for node_x = 1 : n
    for i = 1 : Q
        row = get_F_1_variable_index(node_x, i, Q);
        % row corresponds to f_i(x)

        for j = 1 : Q
            for k = 1 : Q
                for ell = 1 : Q
        
                    col = get_F_3_variable_index(node_x, j, node_x, k, node_x, ell, n, Q);
                    % col corresponds to f_j(x)*f_k(x)*f_ell(x)
        
                    % using the "DENSE" version of the formula in equation (111)
                    T_3_3_dense = (-9.0*w(i)/(2*tau))*(c(i)*c(k))*(c(i)*c(ell)); % (99b)
                    T_4_3_dense = (-3.0*w(i)/(2*tau))*(c(k)*c(ell)); % (106b)
                    
                    F_3(row,col) = T_3_3_dense + T_4_3_dense;
                    
                    % visually inspecting zero terms when vectors c_i are
                    % zero:
                    % if F_3(row,col) == 0
                    %     F_3(row,col) = -10;
                    % end
                
                end

            end
        end

    end
end


figure(5);
imagesc(F_3);
title("F_3");
title(sprintf("F_3 matrix. size is (nQ,(nQ)^3)=(%dx%d)",nQ,(nQ)^3));
clim([-20,20])
colormap(my_colormap)
colorbar('Ticks',[-10,0,10],...
         'TickLabels',{'negative','0','positive'})
saveas(gcf, "figure--F_3_matrix.png");



























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the Carleman Linearized Matrix A

% A =   [ F_hat         F_2         F_3         ]
%       [ 0             F_hat^(2)   F_2^(2)     ]
%       [ 0             0           F_hat^(3)   ]

I = eye(nQ);

dim_A = nQ + (nQ)^2 + (nQ)^3;

A = zeros(dim_A, dim_A);


A(1:nQ, 1:nQ) = F_hat;
A(1:nQ, nQ + 1 : nQ + (nQ)^2) = F_2;
A(1:nQ, nQ + (nQ)^2 + 1 : end) = F_3;

% F_hat^(2)
fprintf("calculating F_hat^(2)...\n")
tic
A(nQ+1 : nQ + (nQ)^2 , nQ+1 : nQ + (nQ)^2) = kron(I,F_hat) + kron(F_hat,I);
elapsed_time = toc;
fprintf("done.  elapsed time: %f.2\n\n\n",elapsed_time);

% F_2^(2)
fprintf("calculating F_2^(2)...\n")
tic
A(nQ+1 : nQ + (nQ)^2 , nQ + (nQ)^2 + 1 : end) = kron(I,F_2) + kron(F_2,I);
elapsed_time = toc;
fprintf("done.  elapsed time: %f.2\n\n\n",elapsed_time);


% F_hat^(3)
fprintf("calculating F_hat^(3)...\n")
tic
A(nQ + (nQ)^2 + 1 : end , nQ + (nQ)^2 + 1 : end) = ...
    + kron(I,kron(I,F_hat)) ...
    + kron(I,kron(F_hat,I)) ...
    + kron(F_hat,kron(I,I));
elapsed_time = toc;
fprintf("done.  elapsed time: %f.2\n\n\n",elapsed_time);



hFig = figure(6);
set(hFig,'Position',[100,100,1000,900])
imagesc(A);
clim([-20,20])
colormap(my_colormap)
colorbar('Ticks',[-10,0,10],...
         'TickLabels',{'negative','0','positive'})
title(sprintf("Carleman Matrix A. size is (nQ + (nQ)^2 + (nQ)^3,nQ + (nQ)^2 + (nQ)^3)=(%dx%d)",dim_A,dim_A));
saveas(gcf, "figure--A_matrix.png");



A_block = zeros(dim_A, dim_A);
A_block(1:nQ, 1:nQ) = ones(nQ,nQ);
A_block(1:nQ, nQ + 1 : nQ + (nQ)^2) = 0.9*ones(nQ,(nQ)^2);
A_block(1:nQ, nQ + (nQ)^2 + 1 : end) = 0.8*ones(nQ,(nQ)^3);

A_block(nQ+1 : nQ + (nQ)^2 , nQ+1 : nQ + (nQ)^2) = 0.7*ones((nQ)^2,(nQ)^2);
A_block(nQ+1 : nQ + (nQ)^2 , nQ + (nQ)^2 + 1 : end) = 0.6*ones((nQ)^2,(nQ)^3);

A_block(nQ + (nQ)^2 + 1 : end , nQ + (nQ)^2 + 1 : end) = 0.5*ones((nQ)^3,(nQ)^3);

hFig = figure(7);
set(hFig,'Position',[100,100,1000,900])
imagesc(A_block);
title("Carleman Matrix A: highlighting possible nonzero sections")
saveas(gcf, "figure--A_matrix_sections.png");



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some basic analysis, printed to the console
% 


n
Q
fprintf("nQ = %d\n\n", nQ);
fprintf("(nQ)^2 = %d\n\n", nQ^2);
fprintf("(nQ)^3 = %d\n\n", nQ^3);

size_S = size(S)
rank_S = rank(S)

size_F_1 = size(F_1)
rank_F_1 = rank(F_1)

rank_F_hat = rank(F_hat)

size_F_2 = size(F_2)
size_F_3 = size(F_3)

fprintf("nQ + (nQ)^2 + (nQ)^3 = %d\n\n", nQ+nQ^2+nQ^3);
size_A = size(A)

A_max_value = max(max(A))
A_min_value = min(min(A))
A_abs_max_value = max(max(abs(A)))

S_eigenvalues = eig(S)
F_1_eigenvalues = eig(F_1)
F_hat_eigenvalues = eig(F_hat)
% A_eigenvalues = eig(A) % doesn't take that long... but probably don't
% want the output printed to console.









% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % forward Euler method to see what happens after some time T
% % % 


% initial_values = zeros(nQ+nQ^2+nQ^3,1);
% 
% fprintf("initializing values for forward Euler method...\n")
% tic
% for node_x = 1 : n
%     for i = 1 : Q
% 
%         f_i_x_index = get_F_1_variable_index(node_x,i,Q);
% 
%         % initialize to equilibrium weight plus some noise.
%         initial_values(f_i_x_index) = w(i) + 0.1*randn(1); 
% 
%     end
% end
% 
% for node_x = 1 : n
%     for i = 1 : Q
% 
%         f_i_x_index = get_F_1_variable_index(node_x,i,Q);
% 
%         for node_y = 1 : n
%             for j = 1 : Q
% 
%                 f_j_y_index = get_F_1_variable_index(node_y,j,Q);
% 
%                 second_order_val = initial_values(f_i_x_index)*initial_values(f_j_y_index);
% 
%                 second_order_variable_index = nQ + get_F_2_variable_index( ...
%                     node_x, ...
%                     i, ...
%                     node_y, ...
%                     j, ...
%                     n, ...
%                     Q);
% 
%                 initial_values(second_order_variable_index) = second_order_val;
% 
% 
%                 for node_z = 1 : n
%                     for k = 1 : Q
% 
%                         f_k_z_index = get_F_1_variable_index(node_z,k,Q);
% 
%                         third_order_val = initial_values(f_i_x_index)*initial_values(f_j_y_index)*initial_values(f_k_z_index);
% 
%                         third_order_variable_index = nQ + (nQ)^2 + get_F_3_variable_index( ...
%                             node_x, ...
%                             i, ...
%                             node_y, ...
%                             j, ...
%                             node_z, ...
%                             k, ...
%                             n, ...
%                             Q);
% 
%                         initial_values(third_order_variable_index) = third_order_val;
%                     end
%                 end
%             end
%         end
%     end
% end
% elapsed_time = toc;
% fprintf("done.  elapsed time: %f.2\n\n\n",elapsed_time);
% 
% 
% 
% f = initial_values;
% num_steps = 10000;
% h = 0.1;
% 
% figure(10)
% stem(f(1:nQ));% only plot the values of the original first-order variables.
% title("starting population distributions")
% 
% 
% 
% fprintf("executing forward Euler method...\n");
% tic
% for step = 1 : num_steps
%     f = f + h*A*f;
%     if mod(step,50)==0
%         fprintf("step %d/%d\n",step,num_steps)
%     end
% end
% elapsed_time = toc;
% fprintf("done.  elapsed time: %f.2\n\n\n",elapsed_time);
% 
% 
% figure(11)
% stem(f(1:nQ));% only plot the values of the original first-order variables.
% title("final population distributions")
% 
% 


