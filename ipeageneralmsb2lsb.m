
% this estimates the phase starting from the most significant bit
% toward the least significant bit.
m = 30;                             %number of bits,accuracy
N = 16;                             %system size

[Q, ~] = qr(randn(N));              %a random unitary
l = rand(N,1);                      %random eigenvalues
H = Q*diag(l)*Q';                   %random Hamiltonian
I = eye(N);
eindx = 8;                          %which eigenvalue
Hadamard = [1 1;1 -1]/sqrt(2);      %hadarmad matrix
Uh = kron(Hadamard, I);

q = zeros(m,1);                     %estimated bit values
CU = eye(2*N);
U = expm(1i*pi*H);
in = kron([1;0],Q(:,eindx));        %input to the circuit
P0 = zeros(2,m);                    %probabilities
for j = 0: m-1
    j1 = j+1;
    % to remove the part b1b2b3.b4b5
    %after anglepi, we have b3.b4b5
    %then apply -pi/2, so that we have cos(b3.b4b5*pi-pi/2)
    %if b3=1, then cos in [pi/2, 3pi/2]
    %if b3=0, then cos in [-pi/2, pi/2]
    angle = 0;
    for jj = 1:j
        angle = angle + q(jj)*2^(j1-jj);
    end
    angle = angle*pi;
    Z = [1 0;0 exp(-1i*angle)*exp(-1i*pi/2)];
    Uz = kron(Z,I);
    
    p = 2^(j1);
    
    U1 = U^p;
    CU = blkdiag(I,U1);
    
    out = Uh*Uz*CU*Uh*in;
    
    P0(:,j1) = calculate_block_prob(abs(out).^2, 1);
    if(P0(2,j1) > P0(1,j1))
        q(j1) = 1;
    end
    
end
%compute the estimated phase
phase = 0;
for j= 1:m
    phase = phase + q(j)*2^(-j);
end
fprintf('the actual and the estimated phases are:\n%1.8f\n%1.8f\n',...
    l(eindx), phase)


