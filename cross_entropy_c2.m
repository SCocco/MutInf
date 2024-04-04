function [f,g] = cross_entropy_c2(x, pa, pb, pab, ga_h, ga_J, qi_tot, qj_tot)

    
    Jab = reshape(x((qi_tot + qj_tot + 1) : (qi_tot * qj_tot + qi_tot + qj_tot)), qj_tot, qi_tot)';
    ha=x(1 : qi_tot)';
    hb=x(qi_tot + 1 : qi_tot + qj_tot)';
    
    ham=repmat(ha,[1,qj_tot]);
    hbm=repmat(hb,[1,qi_tot])';
    Jm=ham+hbm+Jab;
    
    % Partition Function
    zz1=1;
    zz2=sum(exp(ha))+sum(exp(hb));
    zz3=sum(sum(exp(Jm)));
    zz=zz1+zz2+zz3;
    % energy 
    
    ej=-sum(sum(Jab.*pab));
    ehi=-sum(ha.*pa');
    ehj=-sum(hb.*pb');
    ee=ej+ehi+ehj;
    
    %regularization
    reg = ga_h * (norm(ha)^2 + norm(hb)^2) + ga_J * norm(Jab)^2;
    
    f=log(zz)+ee+reg;
    

    if nargout > 1 % gradient required
   
        g = zeros(length(x),1);
        
        % gradient of h_a
        
        Jmb = Jab + hbm;
        for i=1:qi_tot
            g(i) = exp(ha(i)) * (1 + sum(exp(Jmb(i,:))));
            g(i) = g(i)/zz - pa(i) + 2 * ga_h * ha(i);
        end
        
        % gradient of h_b
        Jma = Jab + ham;
        for i=1:qj_tot
            g(i + qi_tot) = exp(hb(i)) * (1 + sum(exp(Jma(:,i))));
            g(i + qi_tot) = g(i + qi_tot)/zz - pb(i) + 2 * ga_h * hb(i);
        end
         
        for i=1:qi_tot
            for j=1:qj_tot
                index = (i-1)*qj_tot + j;
                g(qi_tot + qj_tot + index) = exp(Jm(i,j))/zz - pab(i,j) + 2 * ga_J * Jab(i,j);
            end
        end  
    end
    
end


