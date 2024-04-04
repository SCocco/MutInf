function [f,g] = cross_entropy_c1(x, pa, ga_h, qi_tot)

    
    
    ha=x(1 : qi_tot)';
    
    
    
    
    % Partition Function
    zz1=1;
    zz2=sum(exp(ha));
    
    zz=zz1+zz2;
    % energy 
    
   
    ehi=-sum(ha.*pa');
    
    
    %regularization
    reg = ga_h * norm(ha)^2 ;
    
    f=log(zz)+ehi+reg;
    

    if nargout > 1 % gradient required
   
        g = zeros(length(x),1);
        
        % gradient of h_a
        
        
        for i=1:qi_tot
            g(i) = exp(ha(i)) ;
            g(i) = g(i)/zz - pa(i) + 2 * ga_h * ha(i);
        end
        
        
    end
    
end


