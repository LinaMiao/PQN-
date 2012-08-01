function op = opFFT2(dim)

% fft and inverse fft on the last two dimensios of input

m = prod(dim);

% Construct the operator
fh = @(data,mode) opFFT3_intrnl(data,dim,mode);
op = opFunction(m, m, fh);


function y = opFFT3_intrnl(u,dim,mode)

    if (mode == 0)
        
        y = {m,m,[0,0,0,0],{'opDFT'}};
        
    elseif (mode == 1)
        
        if size(size(u))~=3
            u = reshape(u,dim);
        end
        spec = u;
        spec = fft(spec,[],2)/sqrt(size(u,2));
        spec = fft(spec,[],3)/sqrt(size(u,3));
        
        y = vec(spec);
        
    else % mode = -1
        
        if size(size(u))~=3
            u = reshape(u,dim);
        end
        v = u;
        v = ifft(v,[],2)*sqrt(size(u,2));
        v = ifft(v,[],3)*sqrt(size(u,3));
        y = vec(v);
        
    end % end of if mode == 0




end

end