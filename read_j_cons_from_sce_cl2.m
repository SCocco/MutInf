%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function change the couplings and fields' gauge into consensus.
%
% List of inputs:
%
% - filename = output of the ACE code in standard format (*.j)
% - q_tot = Nx1 vector, q_tot(i) is the effective number of colors at 
%           site i.
% - conscut = Nx1 vector, conscut(i) is the consensus color at site i
%             in the reduced alphabet /!\ (in *.mat)
% - N = size of the alignment
%
% List of outputs:
% - J = cell array. J{i,j} is a rectangular couplings matrix of size 
%          q_tot(i) x q_tot(j) in gauge consensus (1 line and 1 column 
%          filled with zeros)
% - h = cell array. h{i} is a fields vector of size q_tot(i) in gauge 
%          consensus (1 entry zero)
%
% Usage : [J,h]=map2cons4real(ER05_sJ1_sH5_Q10_maxC7_r5_B100_ps0_cons.j,q_tot,conscut,50);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [J,h]=read_j_cons_from_sce_cl2(Jold,hold,q_tot,conscut,N,contacts)
	
    % a modifier
%Jold=J2
%hold=h
%contacts=contacts_sel
%conscut=disscut

	h=cell(1,N);
	for ii=1:N
        %foo=hold{i};
        h{ii}=zeros(1,q_tot(ii));
        h{ii}(1:conscut(ii)-1)=hold{ii}(1:conscut(ii)-1);
        h{ii}(conscut(ii)+1:q_tot(ii))=hold{ii}(conscut(ii):q_tot(ii)-1);
		%h{i}=[hh(i,1:q_tot(i)-1) 0];
	end

	%J=cell(N-1,N);
    ncl=size(Jold,1);
    J=cell(ncl,1);

    for kc=1:ncl
            ii=contacts(kc,1);
            jj=contacts(kc,2);
            J{kc}=zeros(q_tot(ii),q_tot(jj));
            J{kc}(1:conscut(ii)-1,1:conscut(jj)-1)=Jold{kc}(1:conscut(ii)-1,1:conscut(jj)-1); % en haut à gauche
            J{kc}(conscut(ii)+1:q_tot(ii),1:conscut(jj)-1)=Jold{kc}(conscut(ii):q_tot(ii)-1,1:conscut(jj)-1); % en bas à gauche
            J{kc}(1:conscut(ii)-1,conscut(jj)+1:q_tot(jj))=Jold{kc}(1:conscut(ii)-1,conscut(jj):q_tot(jj)-1); % en haut à droite
            J{kc}(conscut(ii)+1:q_tot(ii),conscut(jj)+1:q_tot(jj))=Jold{kc}(conscut(ii):q_tot(ii)-1,conscut(jj):q_tot(jj)-1); % en bas à droite
    end

    
  end  
    