function H = householderMatrix (x)
% forms a Householder matrix for the vector u, 
%the first row column of the matrix is u.
% u is a column vector;

m = length(x);

if x(1) == 0
    sx = 1;
else
    sx = sign(x(1));
end
u = x + sx*eye(m,1);
u = u/norm(u);
I = eye(m);
H = I - (2*u)*u';
H = H./(-sx);
end
