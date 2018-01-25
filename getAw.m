function A = getAw(H,w)
% finds A_w on matrix H
% e.g. getA(H,'01')
% w can only include 0,1,2,3
N = size(H,1);

% start indices of the matrix
rstart = 1;
cstart = 1;
rend = N;
cend = N;
for j = 1:length(w)
    r2 = (rstart + rend-1)/2;
    c2 = (cstart + cend-1)/2;
%  reminder switch does not have "break;" in MATLAB
    switch( w(j))
        case '0'
            rend = r2;
            cend = c2;
        case '2'
            rstart = r2+1;
            cend = c2;
        case '1'
            cstart = c2+1;
            rend = r2;
        case '3'
            rstart= r2+1;
            cstart = c2+1;
        otherwise
            error('unknown letter in w');
    end
end
A = H(rstart:rend,cstart:cend);
