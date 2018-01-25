function [Pj, Vj] = getPjVj(H,j)
% for a given Hamiltonian matrix and index j, it returns Vj and Pj
% see the paper for explanations
N = size(H,1);
n = log2(N);
k = n-1; % at kth step of the recursion we have 2x2 matrices
K = 2^k;

I = eye(2);
X = [0 1;1 0];

W = cell(k,1);
for i = 1: K
    W{i} = strrep(dec2bin(i-1, k), '1', '3');
end

%compute permuation matrix Pj
w = dec2bin(j,k);
Pj = 1;
for i = 1:k
    if (w(i) == '1')
        Pj = kron(Pj,X);
    else
        Pj = kron(Pj,I);
    end
end
Pj = kron(Pj, eye(N/K));

% % % % % % % % % % % % % % % %
Vj = zeros(N);
if (j == 0)
    for i = 1:K
        i2 = 2*i-1;
        Vj(i2:i2+1, i2:i2+1) = getAw(H, W{i});
    end
    return;
end

%j is not zero
%form the matrix Vj
for i = 1: K
    wi =  W{i};
    for ii = 1:k
        if (w(ii) == '1')
            %since thre is an X on ith,
            if( wi(ii) == '0')
                wi(ii) = '1';
            else % it is '3'
                wi(ii) = '2';
            end
        end
    end
    %update Vj
    i2 = 2*i-1;
    Vj(i2:i2+1, i2:i2+1) = getAw(H, wi);
    % % % % % % %
end