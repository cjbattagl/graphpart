addpath heuristics/

trials = 20;
test_er_p_range = 0;
test_gamma_range_er = 0;
test_alpha_choice_er = 0;

test_multi_fennel = 0;

test_gamma_range_kron = 0;
test_kron_param = 1;

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


if (test_gamma_range_er)
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

if (test_alpha_choice_er)
  fprintf('Alpha Test:\n');
  quality = zeros(trials, 1);
  for t = 1:trials
    B = fennel(A, 1.5);
    z = nnz(B(1:500, 501:1000)) + nnz(B(501:1000, 1:500));
    quality(t, 1) = z/nnz(B);
  end
  ave_qual = sum(quality, 1)./trials
end

% if (test_scale_free)
%   sfi = kronecker_generator(13, 2);
%   adj = sparse(sfi(1, :), sfi(2,:), ones(1,length(sfi(1,:))));
%   B = fennel(adj, 1.5);
%   q = cutsize(B, 2)
%   
% end

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

if (test_gamma_range_kron)
  fprintf('Gamma Range Test:\n');
  kron_size = 10;
  parts = 2;
  gamma_range = linspace(1.1,3.1,20);
  quality = zeros(trials, length(gamma_range));
  c = 1;
  for t = gamma_range
    fprintf(' %2.3f\n', t)
    for tr = 1:trials
      kron_list = kronecker_generator(kron_size, 1, 0.57, 0.19, 0.19);
      kron_list = kron_list + 1;
      adj = sparse(kron_list(1, :), kron_list(2,:), ones(1,length(kron_list(1,:))), 2^kron_size, 2^kron_size);
      B = fennel(adj, t, parts);
      quality(tr, c) = cutsize(B, parts);
    end
    c = c+1;
  end
  ave_qual = sum(quality, 1)./trials;
  
  plot(gamma_range, ave_qual, 'o');
  title('(Cost Fcn) Gamma-Range Tests')
  xlabel('Value of Gamma');
  ylabel('Quality of Partition');
%   set(gca,'FontSize',fsize)
%   set(findall(gcf,'type','text'),'FontSize',fsize)
%   saveas(fig,'figures/gamma_range_kron.png', 'png')

end

if (test_kron_param)
  fprintf('Kronecker Parameter Range Test:\n');
  kron_size = 8;
  parts = 2;
  p1 = linspace(0, 1, 10);
  p2 = linspace(0, 1, 10);
  quality = zeros(trials, length(gamma_range));
  c1 = 1;
  for t = p1
    fprintf(' %2.3f\n', t)
    c2 = 1;
    for u = p2
      qual_acc = 0;
      for tr = 1:trials
        kron_list = kronecker_generator(kron_size, 1, 0.57, 0.19, 0.19);
        kron_list = kron_list + 1;
        adj = sparse(kron_list(1, :), kron_list(2,:), ones(1,length(kron_list(1,:))), 2^kron_size, 2^kron_size);
        B = fennel(adj, t, parts);
        qual_acc = qual_acc + cutsize(B, parts);
      end
      quality(c2, c1) = qual_acc/trials; 
      c2 = c2+1;
    end
    c1 = c1+1;
  end
%   ave_qual = sum(quality, 1)./trials;
  
%   plot(gamma_range, ave_qual, 'o');
  plot(quality)
  title('(Cost Fcn) Gamma-Range Tests')
  xlabel('Value of Gamma');
  ylabel('Quality of Partition');
%   set(gca,'FontSize',fsize)
%   set(findall(gcf,'type','text'),'FontSize',fsize)
%   saveas(fig,'figures/gamma_range_kron.png', 'png')

end
























