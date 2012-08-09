function op = opDSR_mig_time(dim,t,x,y,z,v,option)

% function op = opTEST11(dim,t,x,y,z,v,option)

% by defaut extended migration, option = 1
if not(exist('option','var'))
    option = 1;
end

m = prod(dim);
if option == 0 % migration
    n = length(z)*length(x);
else
    n = length(z)*length(x)*length(y);
end
% Construct the operator
fh = @(data,mode) opDSR_mig_time_intrnl(data,t,x,y,z,v,option,mode);
op = opFunction(n, m, fh);


function y = opDSR_mig_time_intrnl(data,t,x,y,z,v,option,mode)

    if (mode == 0)
        
        y = {m,m,[0,0,0,0],{'opDFT'}};
        
    elseif (mode == 1)
        
        y = DSR_mig_time(data,t,x,y,z,v,option);
        y = vec(y);
        
    else % mode = -1
        
        y = DSR_mig_time_inv(data,t,x,y,z,v,option);
        y = vec(y);
       
    end % end of if mode == 0

    
   

end

end
