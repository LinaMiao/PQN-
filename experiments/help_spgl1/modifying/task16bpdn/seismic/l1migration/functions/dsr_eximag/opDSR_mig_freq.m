function op = opDSR_mig_freq(dim,f,x,y,z,v,option)

% Construct the operator
fh = @(data,mode) opDSR_freq_intrnl(data,f,x,y,z,v,option,mode);
m = prod(dim);

% by defaut extended migration, option = 1
if not(exist('option','var'))
    option = 1;
end

if option == 0 % migration
    n = length(z)*length(x);
else
    n = length(z)*length(x)*length(y);
end

op = opFunction(n, m, fh);


function output = opDSR_freq_intrnl(data,f,x,y,z,v,option,mode)

    if (mode == 0)
        
        output = {m,m,[0,0,0,0],{'opDFT'}};
        
    elseif (mode == 1)
        
        output = DSR_mig_freq(data,f,x,y,z,v,option);
        output = vec(output);
        
    else % mode = -1
        
        output = DSR_mig_freq_inv(data,f,x,y,z,v,option);
        output = vec(output);
        
    end % end of if mode == 0




end

end