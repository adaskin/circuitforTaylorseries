function f = calculate_block_prob(psi,qubits)
% computes probabilities for the states of a given qubit group.
% psi is a probability vector for each state.
%  qubits is the set of qubits such as [1 3 5]
% % the order of the qubits 1 2 3 4;
% % @Anmer Daskin
m = size(qubits,2);
[n,kk] = size(psi);
logn = log2(n)+1;

qubits = sort(qubits,  'descend');
qubits = logn-qubits;%reverse the order of the qubits


f = zeros(2^m,kk);

for k = 1:kk %the columns of the input
    for j = 0:n-1
        ind = 1.0;
        index = bitget(j,qubits(:));
        for i = 1:m
            %            index(i)=index(i)*(2^(i-1));%in octave this is being converted to logical.
            ind = ind+index(i)*(2^(i-1));
            
        end
        
        f(ind,k) = f(ind,k)+psi(j+1,k);
    end
    
end
end
% end
