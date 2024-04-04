%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function rewrites a general alignment with a reduced number of
% colors given a threshold on their frequencies.
%
% Reduction scheme: keep all colors which frequency is above the threshold
% and regroup the other colors (making sure q>=2). 
%
% List of input:
% - q = number of colors in the input alignment (q=21 for proteins)
% - filename = file containing a MxN matrix of numbers from 1 to q
% - pcut = alphabet reduction threshold
%
% List of output:
% - align_cut = MxN matrix of numbers from 1 to q_tot(i) in column i. 
%               Reduced alignment
% - msacut = M*sum(q_tot) matrix of spins of size q_tot(i) with 0 and 1.
%            Reduced alignment
% - q_tot = Nx1 vector, q_tot(i) is the effective number of colors at 
%           site i.
% - q_reg = Nx1 vector, q_reg(i) is the number of regrouped colors at 
%           site i.
% - q_kept = Nx1 vector, q_kept(i) is the number of kept colors at site i. 
% - index = Nxq matrix, each line is the list of colors in the following
%           order: [kept, regrouped, not observed]
% - conscut = Nx1 vector, conscut(i) is the consensus color at site i
%
% Usage: [aligncut,q_tot,q_reg,q_kept,index,conscut]=cut_msa_p_simple(align,0.01,21)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [aligncut,msacut,q_tot,q_reg,q_kept,index,conscut,disscut,pmat]=cut_msa_p_simple_d(align,pcut,q)
    
    [M,N]=size(align);
    

    % alignment msa21 on N*Q sites
    msa21=zeros(M,N*q);
    for i=1:N
        for a=1:q   
            % a.a. are labeled from 1 to q
            msa21(:,(i-1)*q+a)=(align(:,i)==a);
        end
    end

    p=sum(msa21)/M; % p = frequence of amino acid on site i

    % matrix pf frequencies QxN
    pmat=reshape(p,q,N);

    q_tot=zeros(1,N); % total # of colors
    q_reg=zeros(1,N); % # regrouped colors
    q_kept=zeros(1,N); % # kept colors
    index=zeros(N,q); % list of indices: [ kept, regrouped, not observed ]

    % special case i=1
    ind_kept=find((pmat(:,1)>pcut));
    ind_out=find((pmat(:,1))&(pmat(:,1)<=pcut));
    ind_zero=find(pmat(:,1)==0);
    q_kept(1)=size(ind_kept,1);
    q_reg(1)=size(ind_out,1);

    % if all observed colors have frequency below pcut, still make sure q_i>=2
    if(q_kept(1)==0) 
        [~,maxout]=max(pmat(ind_out,1));
        ind_kept=ind_out(maxout);
        ind_out(maxout)=[];
        q_kept(1)=size(ind_kept,1);
        q_reg(1)=size(ind_out,1);
    end 
    
%     % 2) if fi_regrouped > fi_kept_a, keep one more aminoacid
%     while(sum(repmat(sum(pmat(ind_out,1)),size(ind_kept,1),1)>pmat(ind_kept,1))>0) 
%         [~,maxout]=max(pmat(ind_out,1));
%         ind_kept=[ind_kept;ind_out(maxout)];
%         ind_out(maxout)=[];
%         q_kept(1)=size(ind_kept,1);
%         q_out(1)=size(ind_out,1);
%     end
%    ind_kept=sort(ind_kept);
%    ind_out=sort(ind_out);
    
    % list of indices
    index(1,:)=[ind_kept' ind_out' ind_zero'];

    % real indices in the N*Q spins alignment
    ind_kept_tot=ind_kept([1:q_kept(1)]);
    ind_out_tot=ind_out([1:q_reg(1)]);

    % compute total # of colors
    if(q_reg(1)>0) % if some colors are regrouped
        q_tot(1)=q_kept(1)+1;
    else 
        q_tot(1)=q_kept(1);
    end
    
    % other sites i>1
    for i=2:N
        ind_kept=find((pmat(:,i)>pcut));
        ind_out=find((pmat(:,i))&(pmat(:,i)<=pcut));
        ind_zero=find(pmat(:,i)==0);
        q_kept(i)=size(ind_kept,1);
        q_reg(i)=size(ind_out,1);

        if(q_kept(i)==0)
            [~,maxout]=max(pmat(ind_out,i));
            ind_kept=ind_out(maxout);
            ind_out(maxout)=[];
            q_kept(i)=size(ind_kept,1);
            q_reg(i)=size(ind_out,1);
        end 
%         while(sum(repmat(sum(pmat(ind_out,i)),size(ind_kept,1),1)>pmat(ind_kept,i))>0) 
%             [~,maxout]=max(pmat(ind_out,i));
%             ind_kept=[ind_kept;ind_out(maxout)];
%             ind_out(maxout)=[];
%             q_kept(i)=size(ind_kept,1);
%             q_out(i)=size(ind_out,1);
%         end
%         ind_kept=sort(ind_kept);
%         ind_out=sort(ind_out);
        
        % list of indices
        index(i,:)=[ind_kept' ind_out' ind_zero'];
        
        % real indices in the N*Q spins alignment
        ind_kept_new=ind_kept([1:q_kept(i)])+(i-1)*q;
        ind_out_new=ind_out([1:q_reg(i)])+(i-1)*q;
        
        % update and compute total # of colors 
        ind_kept_tot=cat(1,ind_kept_tot,ind_kept_new);
        if(q_reg(i)>0);
            q_tot(i)=q_kept(i)+1;
            ind_out_tot=cat(1,ind_out_tot,ind_out_new);
        else
            q_tot(i)=q_kept(i);
        end
        
    end

    % consensus colors in new alphabet
    [~,cons]=max(pmat);
    conscut=zeros(1,N);
    for i=1:N
        conscut(i)=find(index(i,1:q_kept(i))==cons(i));
    end
    %dissensus colors in new alphabet
    %[~,diss]=min(pmat);
    disscut=zeros(1,N);
     for i=1:N
      
       [pmin(i),disscut(i)]=min(pmat(index(i,1:q_kept(i)),i));

      %disscut(i)=index(i,indpmin(i));
      end
    
    
    % block indices 
    qt=zeros(N+1,1); % for the real alignment
    qtt=zeros(N+1,1); % for the new cut alignment
    qtout=zeros(N+1,1); % taking into account the regrouped color
    
    qtout(1)=0;
    for i=2:N+1
        qt(i)=qt(i-1)+q_kept(i-1);
        qtt(i)=qtt(i-1)+q_tot(i-1);
        qtout(i)=qtout(i-1)+q_reg(i-1);
    end
    
    % cutting the alignment
    msacut=zeros(M,qtt(N+1));  
    for i=1:N
        if (q_reg(i)>0) 
            msacut(:,(qtt(i)+1:qtt(i)+q_kept(i)))=msa21(:,ind_kept_tot(qt(i)+1:qt(i)+q_kept(i))');
            msacut(:,qtt(i)+q_tot(i))=sum(msa21(:,ind_out_tot(qtout(i)+1:qtout(i)+q_reg(i))),2)>0;
        else
           msacut(:,qtt(i)+1:qtt(i)+q_kept(i))=msa21(:,ind_kept_tot(qt(i)+1:qt(i)+q_kept(i))');
        end
    end
    
    % rewriting with numbers
    aligncut=zeros(M,N);
    for i=1:N
        for m=1:M
            aligncut(m,i)=find(msacut(m,qtt(i)+1:qtt(i)+q_tot(i))==1);
        end
    end
    
    %%% TO UNCOMMENT ONLY FOR ACE (writes *.p file)
    %% fix the gauge and compute correlations and frequencies for ACE
     %msacut_gauge=msacut;
     %disscut_ind=zeros(1,N);
     %regcut_ind=zeros(1,N);
     % for i=1:N
    %     conscut_ind(i)=conscut(i)+qtt(i);
    %     regcut_ind(i)=qtt(i)+q_tot(i);
    % end
    % msacut_gauge(:,regcut_ind)=[];    
    % freq=sum(msacut_gauge)/M;
    % corr=cov(msacut_gauge,1)+(freq'*freq);    
    
    
    % conscut_ind=zeros(1,N);
    % regcut_ind=zeros(1,N);
    % for i=1:N
    %     conscut_ind(i)=conscut(i)+qtt(i);
    %     regcut_ind(i)=qtt(i)+q_tot(i);
    % end
    % msacut_gauge(:,regcut_ind)=[];    
    % freq=sum(msacut_gauge)/M;
    % corr=cov(msacut_gauge,1)+(freq'*freq);    
    
     %    % write frequencies and correlations for ACE code
     %    pfile=[filename(1:find(filename=='.',1,'last')-1),'_ps',num2str(pcut),'.p'];    
     %    pfile_reg=[filename(1:find(filename=='.',1,'last')-1),'_ps',num2str(pcut),'_reg.p'];
     %    fid=fopen(pfile,'w');
     %    fid_reg=fopen(pfile_reg,'w');
     %    la=zeros(N,q);
     %    l=0;
    	% for i=1:N   
    	%     for a=1:q_tot(i)-1
     %            l=l+1;
     %            la(i,a)=l;
     %            fprintf(fid,' %e',freq(l));
     %            fprintf(fid_reg,' %e',freq(l));
     %        end
     %        fprintf(fid,'\n');
     %        fprintf(fid_reg,'\n');
     %    end
     %    for i=1:N
     %        for j=i+1:N
     %            for a=1:q_tot(i)-1
     %                for b=1:q_tot(j)-1
     %                    fprintf(fid_reg,' %e',max(corr(la(i,a),la(j,b)),1e-10));
     %                    fprintf(fid,' %e',max(corr(la(i,a),la(j,b)),1e-10));
     %                end
     %            end
     %            fprintf(fid,'\n');
     %            fprintf(fid_reg,'\n');
     %        end
     %    end
     %    fclose(fid);
     %    fclose(fid_reg);
     
    %% msacut needs to delete consensus color     
    %     % write down the MSA in compact format
    %     pfile=[filename(1:find(filename=='.',1,'last')-1),'_ps',num2str(pcut),'.cmsa'];
    %     fid = fopen(pfile,'w');
    %     for i=1:size(msacut,2)
    %         fprintf(fid,'%d\n',-1);
    %         for j=1:size(msacut,1)
    %             if (msacut(j,i)==1)
    %             fprintf(fid,'%d\n',j);
    %             end
    %         end
    %     end
    %     fclose(fid);
    

end
