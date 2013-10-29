addpath heuristics/

trials = 50;
test_er_p_range = 0;
test_gamma_range = 0;
test_alpha_choice = 0;
test_scale_free = 0;
test_multi_fennel = 1;


if (test_er_p_range)
    fprintf('ER p-Range Test:\n');
    er_edge_prob = (1:20)./1000;
    
    quality = zeros(trials, length(er_edge_prob));
    edges = zeros(trials, length(er_edge_prob));
    c = 1;
    for t = er_edge_prob
        fprintf('%f\n', t)
        for tr = 1:trials
            A = sparse(rand(1000) < t);
            A = A + A';
            A = A - diag(diag(A));
            A(A>1) = 1;
            
            B = fennel(A, 1.5);
            z = nnz(B(1:500, 501:1000)) + nnz(B(501:1000, 1:500));
            edges(tr, c) = nnz(B);
            quality(tr, c) = z/nnz(B);
        end
        c = c+1;
    end
    
    ave_qual = sum(quality, 1)./trials;
    ave_edges = sum(edges, 1)./trials;
    
    plot(ave_edges, ave_qual, 'o');
    xlabel('Average Number of Edges per 1000-Nodes');
    ylabel('Quality of Partition');
end


if (test_gamma_range)
    fprintf('Gamma Range Test:\n');
    
    gamma_range = linspace(1.1,3.1,20);
    
    quality = zeros(trials, length(gamma_range));
    c = 1;
    for t = gamma_range
        fprintf(' %f\n', t)
        for tr = 1:trials
            A = sparse(rand(1000) < .001);
            A = A + A';
            A = A - diag(diag(A));
            A(A>1) = 1;
            
            B = fennel(A, t);
            z = nnz(B(1:500, 501:1000)) + nnz(B(501:1000, 1:500));
            quality(tr, c) = z/nnz(B);
        end
        c = c+1;
    end
    ave_qual = sum(quality, 1)./trials;
    
    plot(gamma_range, ave_qual, 'o');
    title('(Cost Fcn) Gamma-Range Tests')
    xlabel('Value of Gamma');
    ylabel('Quality of Partition');
end

if (test_alpha_choice)
    fprintf('Alpha Test:\n');
    quality = zeros(trials, 1);
    for t = 1:trials
        B = fennel(A, 1.5);
        z = nnz(B(1:500, 501:1000)) + nnz(B(501:1000, 1:500));
        quality(t, 1) = z/nnz(B); 
    end
    ave_qual = sum(quality, 1)./trials
end

if (test_scale_free)
    sfi = kronecker_generator(13, 2);
    adj = sparse(sfi(1, :), sfi(2,:), ones(1,length(sfi(1,:))));
    B = fennel(adj, 1.5);
    q = cutsize(B, 2)
    
end

if (test_multi_fennel)
  A = sparse(rand(1000) < .001);
  A = A + A';
  A = A - diag(diag(A));
  A(A>1) = 1;
  
  for t = 2:15
    B = fennel(A, 1.5, t);
    qualedge(t-1)=cutsize(B, t);
  end
  plot(qualedge)
end


