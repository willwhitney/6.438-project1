%% 6.438 FALL 2015
%% MAIN FILE FOR PART I OF PROJECT

% Will Whitney

% make sure you have unzipped all the given files in the same dir...

clc; close all; clear;

%
compressed_file_name = 'mcoli_rate_high'; % compressed file name (can modify)
%   - mcoli_rate_high
%   - mcoli_rate_moderate
%   - mcoli_rate_low

%
fprintf('Start testing compression method ... \n');
disp(compressed_file_name);
load(compressed_file_name);
load('mcoli'); % ground truth (fixed)
load('mcoli_code_dope'); % doping parameters (fixed)
addpath(genpath([pwd '/supplementary_functions_part_I'])); % add functions in path
% fetch some dimensions
m = length(s);
[k,n] = size(H);

% some initializations for variables you will use (no modification necessary)
M_from_code = ones(n,2); % binary msgs coming out of code graph
M_to_source = ones(m,4); % M_from_code after bit-to-alphabet conversion to enter source graph
M_from_source = ones(m,4); % msgs coming out of source graph
M_to_code = ones(n,2); % M_from_source after alphabet-to-bit converstion to enter code graph
vector_error = []; % a vector storing error for every ite
s_hat = zeros(m,1); % estimate of source data
s_hat_old = zeros(m,1); % previous estimate of source data

% using your favorate data structure,
% keep track of updated msgs inside the code and source graph
o_code = struct('node_to_factor', zeros(n, k), 'factor_to_node', zeros(k, n));
o_source = ones(m, m, 4)/ 2;

H = full(H);

% start BP
l = 0;
while(1)
    l = l+1;
    errs = 0;
    fprintf(['Ite num = ' num2str(l) '\n']); % print iteration number
    
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    % [1] CODE GRAPH BP ITERATION
    
    % ************************************************************
    % ************ write your code graph BP code here ************
    % ************************************************************
    % (either directly in here or as function call)
    % input variables:     
    %   M_to_code - msg passed from source to code graph
    %               [size n x 2]
    %   o_code - struct of msgs in code graph that is being updated
    %   x - the compressed data
    %       [size k x 1] (given by us)
    %   H - LDPC matrix 
    %       [size k x n] (given by us)
    %   phi_code - doped node potentials
    %       [size n x 2] (given by us)
    % output variable:
    %   M_from_code - product of incoming code graph msgs at every node
    %                 to be passed to source graph (all after msg update)
    %                 [size n x 2] (convert the LLR message to standard
    %                 message)
    %   o_code - struct of updated msgs in code graph
    display('Step 1')
    
    [M_from_code, o_code] = code_graph_BP(M_to_code, o_code, x, H, phi_code);
    
%     display(size(M_from_code))
    
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % [2] BIT-TO-ALPHABET CONVERSION FOR SOURCE GRAPH BP
    % (no modification necessary)
    display('Step 2')
    M_to_source = msgs_1to8_gray(M_from_code,1,m); % [size m x 4]

    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % [3] SOURCE GRAPH BP ITERATION
    
    % **************************************************************
    % ************ write your source graph BP code here ************
    % **************************************************************
    % (either directly in here or as function call)
    % input variables:
    %   M_to_source - msg passed from code to source graph
    %                 [size m x 4]
    %   o_source - struct of msgs in source graph that is being updated
    %   phi_source - doped node potentials 
    %                [size m x 4] (given by us)
    %   psi_source - the representative edge potential matrix 
    %                representing probabilistic transitions in markov chain
    %                from i to i+1 where the row sum is 1
    %                [size 4 x 4] (given by us)
    % output variable:
    %   M_from_source - product of incoming source graph msgs at every node
    %                   to be passed to code graph (all after msg update)
    %                   [size m x 4]
    %   o_source - struct of updated msgs in source graph
    
    o_source_new = ones(size(o_source));
    
    display('Step 3')
    % do one round of bidirectional sum-product
    for i = 1:m
%         fprintf(['node' num2str(i) '\n']);
        for j = [i-1, i+1]
           if j > m || j < 1
               continue
           end
           
           msg = zeros(4);
           for sj = 1:4
                submsg = 0;
                for si = 1:4
                    prod = 1;
                    if i > j && i+1 <= m
                        prod = prod * o_source(i+1, i);
                    elseif i < j && i-1 >= 1
                        prod = prod * o_source(i-1, i);
                    end
                    submsg = submsg + M_to_source(1, i, si) * phi_source(i, si) * psi_source(si, sj) * prod;
                end
                msg(sj) = submsg;
           end
           
           msg = msg / sum(msg);
           o_source_new(i, j, :) = msg;
        end
    end
    
    % write out the M_from_source values
    for i = 1:m
        result = 1;
        for j = [i-1, i+1]
            if j > m || j < 1
                continue
            end
            result = result .* o_source_new(j, i);
        end
        M_from_source(i, :) = result;
    end
%     display(sum(o_source_new - o_source));
    
    % save our new messages for use next time
    o_source = o_source_new;
    
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % [4] ALPHABET-TO-BIT CONVERSION FOR CODE GRAPH BP
    % (no modification necessary)
    display('Step 4')
    M_to_code = msgs_8to1_gray(M_from_source); % [size n x 2]
    
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        
    % ***************************************************************
    % **** write your code here for decoding from source msgs *******
    % ***************************************************************
    % (either directly in here or as function call)
    % input variables:
    %   M_from_source - msgs coming into nodes within source graph
    %                   [size m x 4]
    %   M_to_source - msgs coming into nodes from code graph
    %                 [size m x 4]
    %   phi_source - doped source graph node potentials
    %                [size m x 4]
    % output variable:
    %   s_hat - decoded solution using marginal mode
    %          [size m x 1]
    
    display('Step 5')
    for i = 1:m
        options = zeros(4, 1);
        for j = 1:4 
            options(j) = phi_source(i, j) * M_from_source(i, j) * M_to_source(1, i, j); 
        end
        [thing, argm] = max(options);
       
        s_hat(i) = argm-1;
    end
    

    % ************************************************************
    % ****** write your code here for computing error ************
    % ************************************************************
    % (either directly in here or as function call)
    % input variables:
    %   s_hat - decoded solution
    %           [size m x 1]
    %   s - true source data
    %       [size m x 1] (given by us)    
    % output variable:
    %   errs - abs difference error
    %          [scalar variable]
    
    display('Step 6')

    for i = 1:m
        if s_hat(i) ~= s(i)
            errs = errs + 1;
        end
    end
    
%     if l > 1
%        display( sum(abs(M_from_code - M_from_code_old)))
%     end
%     display(M_from_code(1:100, :))
%     M_from_code_old = M_from_code;
    
    vector_error = [vector_error errs];
    fprintf(['... Error = ' num2str(errs) '\n']);
    % terminate if BP gradient doesn't change
    if(l>1 && sum(abs(s_hat - s_hat_old))<0.5)
        break;  % exit BP loop
    end
    s_hat_old = s_hat; % update solution
end


% ******************************************************************************
% ****** write your code here for plotting vector_error vs. l (iteration) ******
% ******************************************************************************


% note: if vector_error converges to 0 exactly, 
%       then you have sucessfully achieved lossless compression 

