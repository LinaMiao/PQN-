function y = A_handle_multiply(A,x,mode)

if mode == 1
    x = vec(x);
    y = A*x;
else
    m = size(A,1);
    n = numel(A)/m;
    x = reshape(x,m,n);
    y = A'*x;
end

