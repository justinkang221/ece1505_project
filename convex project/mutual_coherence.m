function x = mutual_coherence(A)
[~,n] = size(A);
x=0;
for i=1:n
    for j=1:i
        if i==j
            continue
        end
        x_t = abs(A(:,i)'*A(:,j))/(norm(A(:,i)*norm(A(:,j))));
        x = max([x, x_t]);
    end
end