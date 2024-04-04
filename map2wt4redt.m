%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function change the couplings and fields' gauge into consensus.
%
% List of inputs:
%
% - align = matlab array. real (=not reduced) alignment (load from *.msa)
% - Jold = cell array. J{i,j} is a rectangular couplings matrix of size 
%          q_tot(i) x q_tot(j).
% - hold = cell array. h{i} is a fields vector of size q_tot(i)
% - q = maximum number of colors per site
%
% List of outputs:
% - Jnew = cell array. J{i,j} is a rectangular couplings matrix of size 
%          q_tot(i) x q_tot(j) in gauge consensus (1 line and 1 column 
%          filled with zeros)
% - hnew = cell array. h{i} is a fields vector of size q_tot(i) in gauge 
%          consensus (1 entry zero)
%
% Usage : [J_real_cons,h_real_cons]=map2cons4real(align_real,J_real,H_real,10);
%load('ER05_sJ1_sH5_Q10_maxC7_r5_B1000.msa');
%align=load('ER05_sJ1_sH5_Q10_maxC7_r5_B1000.msa');
%[J_real_cons,h_real_cons]=map2cons4real(align,J_real,H_real,10);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Jnew,hnew]=map2wt4redt(wtcut,Jold,hold,q_tot,N,contacts)
    
 %wtcut=disscut   
 
    
    %Jold=J_diss;
    %hold=h_diss;
    %N=27;
    %q_tot=q_keptg+1;
    % Change fields
    hnew=cell(1,N);
    nbc=size(contacts,1);

    Jnew=cell(1,nbc);
    
    for i=1:N
      hnew{i} = (hold{i} - hold{i}(wtcut(i)));  
     
    end      
        
    %nbc=length(contacts);
    for cc = 1:nbc
        ii=contacts(cc,1);
        jj=contacts(cc,2);
        JJ = Jold{cc}';
        hnew{ii} = hnew{ii} + sum( Jold{cc}(:,wtcut(jj)) - repmat( Jold{cc}(wtcut(ii),wtcut(jj)), q_tot(ii), 1 ), 2 )';
        hnew{jj} = hnew{jj} + sum( JJ(:,wtcut(ii)) - repmat( JJ(wtcut(jj),wtcut(ii)), q_tot(jj), 1 ), 2 )';    
        Jnew{cc} = Jold{cc} - repmat( Jold{cc}(:,wtcut(jj)), 1, q_tot(jj) ) - repmat( Jold{cc}(wtcut(ii),:), q_tot(ii), 1 ) + repmat( Jold{cc}(wtcut(ii),wtcut(jj)) , q_tot(ii), q_tot(jj) );         
    end
  
end