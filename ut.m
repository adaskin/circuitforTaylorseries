% function [B, V, U] = ut(H, t)
% forms Ut for a random matrix.
% This code is not optimizedfor large system sizes. 
t = 1/10;
fid = 1; %debug file
threshold = 0.000001;% error threshold
H = rand(16);
N = size(H,1);
if (N ~= size(H,2))
    error('H must be square');
end
n = log2(N);
k = n-1; % at kth step of the recursion we have 2x2 matrices
K = 2^k;
% % % % % % % % % %
I = eye(2);
In = eye(N);
X = [0 1;1 0];
warning('Hadamard gate is not normalized');
Hadamard = [1 1; 1 -1]%/sqrt(2);
fprintf(fid, 'forming BH..................\n');
UHadamard = 1;
for j = 1: k
    UHadamard = kron(UHadamard, Hadamard);
end
BH = kron(UHadamard, eye(N));
fprintf(fid, 'forming VH..................\n');
VH = zeros(N^2/2);
%     H1 = H;
for j = 0:K-1
    [Pj, Vj] = getPjVj(H,j);
    j2 = N*j;
    VH(j2+1:j2+N, j2+1:j2+N) = Vj*Pj;
    %         H1 = H1-Vj*Pj;
end
fprintf(fid, 'forming UH..................\n');
UH = BH'*VH*BH;
if(norm(H-UH(1:N,1:N)) > threshold)
    fprintf('something is wrong with BVB');
end
% % % % % % % % % % % % % % % % % %

%     circuit U(t) requires two additional qubit.
% the total system size is 2n+1
qtotal = 2*n+1;
N2 = N^2;
Ntotal = 2^qtotal;
IN2 = eye(N2);
C0UH = kron(I,blkdiag(UH, 1i*eye(N2/2))); %operates when 2nd qubit is 0
C1UH = blkdiag(eye(N2), kron(I,UH));%operates when 1st qubit is 1

Pi = blkdiag(eye(N), kron(X, eye((N2-2*N)/2)),eye(N));
CPi = blkdiag(eye(N2),Pi);
V = C1UH*CPi*C0UH;

cf = sqrt([t, 1,t^2/2,0])';
ncf = norm(cf);
cf = cf/ncf;
B = mykron(householderMatrix(cf), eye(Ntotal/4));
Ut = B'*V*B;
Utexpected = (t*H + t^2*H^2/2 + 1i*eye(N))/ncf^2;
if(norm(Ut(1:N,1:N)-Utexpected) > threshold)
    warning('Utexpected is different than Ut');
else
    fprintf('Ut is succesfully formed\n');
end

