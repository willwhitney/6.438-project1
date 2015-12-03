function [s_hat] = project_part_II_decoder(...
    x,H,phi_source,phi_code,psi_source,temp)
    % [do not modify arguments (input/output) of function]
    % INPUT
    %   x           - k x 1 compressed data vector 
    %   H           - k x n code matrix
    %   phi_source  - m x 4 source node potentials from doping
    %   phi_code    - n x 2 code node potentials from doping
    %   psi_source  - 4 x 4 chain transition matrix (each column sums to 1)
    %   temp        - [optional] 1 x L array storing L nonneg integer values
    % OUTPUT
    %   s_hat       - m x 1 decoded/decompressed result

    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    % Implement your decoder (going from x to s_hat) here 

    H = full(H);
    psi_source = psi_source';
    m = size(phi_source, 1);
    n = size(phi_code, 1);
    k = size(H, 1);
    num_shotgun = temp(end);
    temp = reshape(temp(1:end-1), num_shotgun, 3);
    temp(:, 1) = temp(:, 1) / 100;
    disp(temp)
    
    l_shotgun = m / num_shotgun;
    
    
    
    messages = ones(m, m, 4, 2) / 4;
    M_from_code = ones(n,2); % binary msgs coming out of code graph
    M_to_source = ones(m,4); % M_from_code after bit-to-alphabet conversion to enter source graph
    M_from_source = ones(m,4); % msgs coming out of source graph
    M_to_code = ones(n,2); % M_from_source after alphabet-to-bit converstion to enter code graph
    vector_error = []; % a vector storing error for every ite
    s_hat = zeros(m,1); % estimate of source data
    s_hat_old = zeros(m,1); % previous estimate of source data
    o_code = struct('node_to_factor', zeros(n, k), 'factor_to_node', zeros(k, n));

    aligned_psi_source = construct_alignment_psi(temp, psi_source, l_shotgun, m);
    r = 1;
    w = 2;
    
    % start BP
    l = 0;
    for t = 1:30
        l = l+1;
%         errs = 0;
        fprintf(['Ite num = ' num2str(l) '\n']); % print iteration number
        
%         if r == 1
%             w = 1;
%             r = 2;
%         else
%             w = 2;
%             r = 1;
%         end
            
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
        
        display('code graph')
        [M_from_code, o_code] = code_graph_BP(M_to_code, o_code, x, H, phi_code);


        % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        % [2] BIT-TO-ALPHABET CONVERSION FOR SOURCE GRAPH BP
        % (no modification necessary)
        
        M_to_source = msgs_1to8_gray(M_from_code,1,m); % [size m x 4]
        M_to_source = reshape(M_to_source, m, 4);

        % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        % [3] SOURCE GRAPH BP ITERATION
        % **************************************************************
        % ************ write your source graph BP code here ************
        % **************************************************************
        % (either directly in here or as function call)
        % input variables:
        %   M_to_source - msg passed from code to source graph
        %                 [size m x 4]
        %   messages - msgs in source graph that is being updated 
        %                 [size m x m x 4 x 2]
        %   phi_source - doped node potentials 
        %                [size m x 4] (given by us)
        %   aligned_psi_source - the representative edge potential matrix 
        %                representing probabilistic transitions in markov chain
        %                from i to j [size m x m x 4 x 4] 
        % output variable:
        %   M_from_source - product of incoming source graph msgs at every node
        %                   to be passed to code graph (all after msg update)
        %                   [size m x 4]
        %   messages - updated msgs in source graph
        %                 [size m x m x 4 x 2]

        % NOTE: psi_source for chain is not symmetrical! Special case for
        % backward msgs?
        
        display('source graph')
        for i = 1:2
            if r == 1
                w = 1;
                r = 2;
            else
                w = 2;
                r = 1;
            end
            [M_from_source, messages] = source_graph_lbp(M_to_source, messages,...
                    phi_source, aligned_psi_source, r, w);
            
        end
            
        % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        % [4] ALPHABET-TO-BIT CONVERSION FOR CODE GRAPH BP
        % (no modification necessary)
        
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

        display('decoding s_hat')
        for first_node = 1:m
            options = zeros(4, 1);
            for j = 1:4 
                options(j) = phi_source(first_node, j) * M_from_source(first_node, j) * M_to_source(first_node, j); 
            end
            [~, argm] = max(options);

            s_hat(first_node) = argm-1;
        end
        
        if exist('s')
            display(['score: ' num2str(sum(s == s_hat))])
        end
        
        if s_hat == s_hat_old
            break
        else
            s_hat_old = s_hat;
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

%         display('Step 6')
% 
%         for i = 1:m
%             if s_hat(i) ~= s(i)
%                 errs = errs + 1;
%             end
%         end

    %     if l > 1
    %        display( sum(abs(M_from_code - M_from_code_old)))
    %     end
    %     display(M_from_code(1:100, :))
    %     M_from_code_old = M_from_code;

%         vector_error = [vector_error errs];
%         fprintf(['... Error = ' num2str(errs) '\n']);
%         % terminate if BP gradient doesn't change
%         if(l>1 && sum(abs(s_hat - s_hat_old))<0.5)
%             break;  % exit BP loop
%         end
%         s_hat_old = s_hat; % update solution
        
    end
    
%     s_hat = zeros(size(H,2)/2,1); % output decoded solution as column vector


end

