
%Program to infer the field and couplings  with the 2-site clusters
%approximation on a list of pairs of sites given  in input (file contacts)
%it compares the fitness cost of single mutations given in mutations
%with the energetic cost in the inferred model  plot and calculate the spearman 
%coefficient and the correlation coefficient r2. 
% Simona Cocco 16/05/2019
%Modification here 
%change output name linea 
%non ordino secondo le frobenious ma secondo ordine naturale  vedi linea 188
%prendo un graph con 1 accoppiamento solo 

%%%%%%%%%%%%%%%%%%%%%%%List of analyzed Proteins%%%%%%%%%%%%%%%


% ncore fixed to 3 for the laptop ( to be changed if more powerful computer eg to 20)
n_core=3;
    for prot=1:1
 if (prot==1)    
    protname='RNA-bind'
 elseif (prot==2) 
    protname='WW'
 elseif(prot==3)
    protname='RL401'
 elseif(prot==4)
    protname='DNA-bind'
 elseif(prot==5)
    protname='BRCA1'
 elseif(prot==6)   
    protname='PDZ'
 elseif(prot==7)  
    protname='UBOX'
 elseif(prot==8)    
    protname='Betalactamase'
 elseif(prot==9) 
    protname='Bgluco'
 end
%to avoid paralelization comment line:290

%default='true';
%to change  default parameters put default='false' 
default='false'
%last modif K2max=''
%K2max fixed al massimo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Input DATA %%%%%%%%%%%%%%%%%%%%% 
%input directory
inputsMSA=['./',protname,'/Dati_inputs_MSAMarks/'];
%only for B-lactanmase at different concentrations

%output directory
dirnameo = ['./',protname,'/Outputs_MSAMarks/1Link_p/'];
%c='c2'
%if onelink
%    switch(c)
%        case 'c0'
%    dirnameo = ['./',protname,'/Outputs_MSAMarks/1Link_full_p_Amp0/'];
%        case 'c1'
%    dirnameo = ['./',protname,'/Outputs_MSAMarks/1Link_full_p_Amp39/'];
%        case 'c2'
%    dirnameo = ['./',protname,'/Outputs_MSAMarks/1Link_full_p_Amp156/'];
%        case 'c3'
%    dirnameo = ['./',protname,'/Outputs_MSAMarks/1Link_full_p_Amp625/'];
%        case 'c4'
%    dirnameo = ['./',protname,'/Outputs_MSAMarks/1Link_full_p_Amp2500/'];
%        case 'c5'
%    dirnameo = ['./',protname,'/Outputs_MSAMarks/1Link_full_p_R_Cefotaxine/'];
%    end
      


%oly for Betalactamase and Herve data  
%inputsMSA=['./',protname,'/Dati_input_MSAHerve_enriched/'];
%inputsMSA=['./',protname,'/Dati_input_MSAHerve/'];
%dirnameo = ['./',protname,'/Outputs_MSAHerve_enriched/'];
%dirnameo = ['./',protname,'/Output_MSAHerve/'];

% Frobenious Norms calculated with Plm  and in consensus gaugewith APC (Average Product Corrections)
%Fapc=importdata([inputsMSA,'frobenious/Fapc_PLM_lambda001_theta02_pcut0_cons.mat']);
%Fapc=importdata([inputsMSA,'frobenious/F_PLM_lambda001_theta02_pcut0_cons.mat']);

% Input Alignment
align = conversion_align([inputsMSA,'sequences.faa']);
[M,N] = size(align);

%reweighting 
theta=0.2;
w = importdata([inputsMSA,'weights/reweighting_0.200.dat']);

% Wild-Type sequence
wt=conversion_align(['./',protname,'/wt.faa']);

%read the Mutational data set: sequences and fitness fitness
X=fastaread(['./',protname,'/fitness_mutations.faa']);
[M_test] = size(X,1);

align_test=zeros(M_test,N);

fitness=zeros(M_test,1);

for i=1:M_test
    fitness(i)=sscanf(X(i).Header,'%f');
    for j=1:N
        align_test(i,j)=letter2number(X(i).Sequence(j));
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters in the loop over nb of couplings in the model K2 %%%%%%%%

%Parameters to fix the iterations in the number of couplings in the inferred model
% at each iteration the spearman between the theoretical prediction and the
% experimental ones is given in output.

% to Change only if default is different from true
default=false
if (default)
    %K2min: minimal number of couplings in the inferred model
    K2min=0;
    %K2max: maximal number of couplings in the inferred model 
    K2max=min(6000,N*(N-1)/2);
    
   
    %%%%%%%%%%%%%%%%%%%%%%% Parameter of the model%%%%%%%%%%%%%%%%%%%%%%%%
    %number of Potts states
    q=21;
    %pcut (inference done on colors with p>pcut)
    pcut=0;
    % threshold of number of states at which sites can be taken as coupled
    % sites if zero it does not play any role
    qcut=0;

    %regolarization strengths : bayesian
    lambdah=0.1/M;
    lambdaJ=10./M;
    sigmah=1;
    sigmaj=1;
    
    
else
     K2min=0;
     %K2max=min(4000,N*(N-1)/2);
     K2max=N*(N-1)/2
     %K2max=1000;
     dK2=1
     
     q=21;
     pcut=0;
     %qcut=10
     qcut=0;
     lambdah=0.1/M;
     lambdaJ=10./M;
     sigmah=1;
     sigmaj=1;
     %no reweighting
     %theta=0;
     %w=ones(M,1);

end    


%cut the alignment
[aligncut,msacut,q_tot,q_reg,q_kept,index,conscut,disscut,freqmsa]=cut_msa_p_simple_d(align,pcut,q);

%wt in the cut alignement
for i=1:N
    wtcut(i)=find(index(i,:)==wt(i));
end

%Msacut with remotion of the gauge symbols
msacut_gauge=msacut;
qkeptsum=zeros(1,N+1);
qtotsum=zeros(1,N+1);
pmin=zeros(1,N);
qkeptsum(1,2:N+1)=cumsum(q_kept);
qtotsum(1,2:N+1)=cumsum(q_tot);

%removing gauge symbols in number of symbols (one parameter less)
q_kept_gauge=q_kept-1;
qtotsum_gauge=zeros(1,N+1);
qtotsum_gauge(1,2:N+1)=cumsum(q_kept_gauge);
%finding (least probable) aa
for i=1:N
    disscut_ind(i)=disscut(i)+qtotsum(i);
    pmin(i)=freqmsa(index(i,disscut(i)),i);
end
% remove the gauge (least probable) aa from the msa
msacut_gauge(:,disscut_ind)=[];
[M_cut,Nq_cut]=size(msacut);
%unseens for an independent model
n_un=q*N-Nq_cut;
%approx of variance related to unseens  in the independent model
sigma_in_un= n_un*sigmah/((N*q));

indl=zeros(N*(N-1)/2,2);
indll=zeros(N,N);
%indexes of the pairs of couplings and FAPC
l=0;
for i=1:N-1
    for j=i+1:N
        l=l+1;
        indl(l,1)=i;
        indl(l,2)=j;
        indll(i,j)=l;
    Fv(l)=Fapc(i,j);
    %Fv(l)=Fapc(l,3);
    end
end
    




%%%%%%%%%% order couplings based on their Frobenious%%%%%%%
%[Fsort,IdxF]=sort(Fv,'descend');
%remove pairs with to much variance due to the unseens
%IdxF2=IdxF; 
if(default==false)
 %ordering by cardinal number (not taking frobenious):
 IdxF2=[1:K2max];
 
%%%In case we want to put couplings on Structural Contacts%%%%%%%%%%%%%
%contacts=importdata('list_contacts625_A.dat');
    
% %For Betalactamase Olivier double mutations on alpha helix
% IdxFo=zeros(1,55);
% l=0;
% for i=88:98;
%     for j=i+1:98
%         l=l+1;
%     IdxFo(l)=indll(i,j);
% end
% end
%IdxF2=[IdxFo,IdxF2];
%indcutpairs=find((q_kept(indl(IdxF,1))<qcut)|(q_kept(indl(IdxF,2))<qcut));
%IdxF2(indcutpairs)=[];
%K2max=size(IdxF2,2);
end

%Extracting 1 and 2 points correlations

NN=size(msacut_gauge,2);

wm=repmat(w,[1,NN]);

%reweighted probabilities
msacut_gauge_w=msacut_gauge.*wm;
p=sum(msacut_gauge_w)/sum(w);
%reweighted pairwise probabilities
p2=msacut_gauge'*msacut_gauge_w/sum(w);

%if(default==false)
%without reweighting
%p=sum(msacutg)/M;
%end

%  1 cluster Inference of Field contribution
h0=cell(1,N);
h=cell(1,N);

for i=1:N
    pa=p(qtotsum_gauge(i)+1:qtotsum_gauge(i)+q_kept_gauge(i));
    x0=zeros(1,q_kept_gauge(i));
     f = @(x)cross_entropy_c1(x,pa,lambdah,q_kept_gauge(i));
     %for matlab versions after 2017
    options = optimoptions(@fminunc,'SpecifyObjectiveGradient',true);
    %options=optimoptions(@fminunc,'GradObj','on');
    x = fminunc(f, x0, options);
    h{i}=x;
    h0{i}=x;
end

%2 cluster Inference of Field and Coupling Contribution on maximal size 
K2=K2max;
contacts=zeros(K2max,2);
contacts(1:K2max,1)=indl(IdxF2(1:K2max),1);
contacts(1:K2max,2)=indl(IdxF2(1:K2max),2);


J=cell(K2max,1);


xh=cell(K2max,1);
% only when (n_core>1)
tic
parpool('local',n_core);

parfor cc=1:K2max
%else
%for cc=1:K2max
   ii=contacts(cc,1);
   jj=contacts(cc,2);
   
    Jab=zeros(q_kept_gauge(ii),q_kept_gauge(jj));
    ha=h0{ii};
    hb=h0{jj};
   
    pa=p(qtotsum_gauge(ii)+1:qtotsum_gauge(ii)+q_kept_gauge(ii));
    pb=p(qtotsum_gauge(jj)+1:qtotsum_gauge(jj)+q_kept_gauge(jj));
    pab=p2(qtotsum_gauge(ii)+1:qtotsum_gauge(ii)+ q_kept_gauge(ii),qtotsum_gauge(jj)+1:qtotsum_gauge(jj)+q_kept_gauge(jj));
    
   
    x0=[ha,hb,reshape(Jab',[q_kept_gauge(ii)*q_kept_gauge(jj),1])'];
    
   
    f = @(x)cross_entropy_c2(x,pa,pb,pab,lambdah,lambdaJ,q_kept_gauge(ii),q_kept_gauge(jj));
    %for matlab versions after 2017
    %options = optimoptions(@fminunc,'SpecifyObjectiveGradient',true);
    %for parallelization see cluster2_crossmin_diss_w_par 
    options=optimoptions(@fminunc,'GradObj','on');
    x = fminunc(f, x0, options);
    
    
   
  J{cc}=reshape(x(q_kept_gauge(jj)+q_kept_gauge(ii)+1:q_kept_gauge(jj)+q_kept_gauge(ii)+q_kept_gauge(ii)*q_kept_gauge(jj)),[q_kept_gauge(jj),q_kept_gauge(ii)])';
  xh{cc}=x(1:q_kept_gauge(ii)+q_kept_gauge(jj));
  
 %h{ii}=h{ii}+xh{cc}(1:q_keptg(ii))-h0{ii};
  
  
end
%when (n_core>1)
delete(gcp)
toc

%Iterations in the number  K  of Couplings included  in the model %
% Potts Energy vs fitness  comparisons for each K
%J2=[];
%cc=0;

spvsk=zeros(K2max,2);
r2vsk=zeros(K2max,2);
%file di punti spearman versus k
sigma_K=zeros(K2max+1,2);
%file di punti degli unseens versus k
n_unK=zeros(K2max,2);

    
for cl=1:K2max
    J2=[];
    cc=cl;
    h=h0;
%one at the time the first iteration cl=1,cc=0 corresponds to the independent model 
%note that I have put J2=[] h=0 inside the loop so at each tilme only one edge 
if (cc>0)   
        
    ii=contacts(cc,1);
    jj=contacts(cc,2);
   
    h{ii}=h{ii}+xh{cc}(1:q_kept_gauge(ii))-h0{ii};  
    h{jj}=h{jj}+xh{cc}(q_kept_gauge(ii)+1:q_kept_gauge(ii)+q_kept_gauge(jj))-h0{jj};
    J2=J(cc);
    end
    
    % single site variance sigma1, variance due to pairs sigma2, number of
    % unseens due to pairs n_unK, number of pairs of aa unseen for each site
    
     contacts_sel=contacts(cc,:);
    [sigma1,sigma2,sigma_Kd,n_unKd,nun2]=sigmak(freqmsa,contacts_sel,wt,M_cut,q_kept,sigmah,sigmaj);
   
    
    %insert  the 0 corresponding to the fields and couplings of the diss in the
    %cut alignement
    [J_diss,h_diss]=read_j_cons_from_sce_cl2(J2,h,q_tot,disscut,N,contacts_sel);
    
    %put the parameters  in wt gauge
    
    [J_wt,h_wt]=map2wt4redt(wtcut,J_diss,h_diss,q_tot,N,contacts_sel);
    
    
    
    %Evaluate the energy of the mutants msa_test coupled and grouped model
    % this can be put in a function
    
    is_mutated = zeros(N*q, 1);
    for k=1:M_test
        %mutated site
        idd=find(align_test(k,:)~=wt);
        %index of colors present in the reduced alphabet in the mutated site
        indi=index(idd,1:q_kept(idd));
        
        % index of mutated aa in the reduced alphabet
        indi_mut=find(indi==align_test(k,idd));
        
        is_mutated((idd-1)*q + align_test(k,idd)) = 1;
        
        %if the mutated aa is not in the explicitely modeled states
        if(isempty(indi_mut))
            % se c'e' il ragruppato su questo sito si prende come riferimento
            if((index(idd,q_kept(idd)+1:q_kept(idd)+q_reg(idd)))>0)
                
                refi=q_kept(idd)+1; %%%% che e' anche q_tot(i) verificarlo
                indexreg=index(idd,q_kept(idd)+1:q_kept(idd)+q_reg(idd));
                
                frefi=  sum(freqmsa(indexreg,idd));%%% la frequenza del ragruppato
            else
                refi=disscut(idd);
                
                frefi=freqmsa(index(idd,disscut(idd)),idd);
                
            end
            indexreg=index(idd,q_kept(idd)+1:q_kept(idd)+q_reg(idd));
            indiar=find(indexreg==align_test(k,idd));
            %if the mutated aa is not in the regrouped symbols
            if(isempty(indiar))
                enh(k)=h_wt{idd}(refi)+log(lambdah/frefi) ;
                % se e' nel regrouped symbol
            else
                enh(k)=h_wt{idd}(refi)+log(freqmsa(align_test(k,idd),idd)/frefi);
            end
            
        else
            enh(k)=h_wt{idd}(indi_mut);
        end
    end
    variance = mean(mean(sigma2(is_mutated>0)));
    n_unseens=n_un;
    if (cc>0)
    n_unseens=n_unKd+n_un;
    end
    %comparison with experiments
    fitnesscomp=zeros(M_test,2);
    fitnesscomp(:,1)=fitness; %experimental fitness
    fitnesscomp(:,2)=enh'; %minus the energy of Potts model
    r2=corr(fitness,enh')^2;
    [rsp,pvalsp] = corr(fitness,enh','Type','Spearman');
    %the first point for cc=0 is the independent model
    %spearman index sp
    spvsk(cl,1)=cc;
    spvsk(cl,2)=rsp;
    %correlation index r2
    r2vsk(cl,1)=cc;
    r2vsk(cl,2)=r2;
    %variance
    sigma_K(cl,1)=cc;
    sigma_K(cl,2)=variance;
    %number of unseens
    n_unK(cl,1)=cc;
    n_unK(cl,2)=n_unseens;
   
    
    %cc=cc+1;
end

%%%%%% Output names changing  the model  and iterarions parameters %%%%%%%%%%%%%%

parname=[];
pariter=[];
% Parameters choiches are written in file name only if they are not default value
%if(default==false)
%parname=['_w',num2str(theta),'_pcut',num2str(pcut),'_qcut',num2str(qcut)];
%pariter=['_Kmin',num2str(K2min),'_Kmax',num2str(K2max),'dK',num2str(dK2)];
%parname=['FnoapC'];
parname=[]
%pariter=['_Kmax',num2str(K2max)];
pariter=['_sorted']
%end

%Recording Sperman vs Number of Couplings in the model
%outfile_rs=[dirnameo,'rs_k',parname,pariter];
%save(outfile_rs,'spvskf','-ascii');
%save(outfile_rs,'spvsk','-ascii');

%Recording Sperman vs Number of Couplings in the model
%outfile_r2=[dirnameo,'r2_k',parname,pariter];
%save(outfile_rs,'spvskf','-ascii');
%save(outfile_r2,'r2vsk','-ascii');

%%plots spearman vs k
%fig3=figure(3);
%hold on
%plot(spvsk(:,1),spvsk(:,2),'b.','MarkerSize',20);
%%txt=[name,'r2=',num2str(r2),'rs=',num2str(rsp)]
%%title(txt)
%xlabel('Number of couplings K')
%ylabel('Spearmann')
%%figname=[dirnameo,'spvsk-',protname,'_dk10_k1055_withalphah.fig']
%%saveas(fig3,figname,'fig');
%%savefig(figname)
%figname=[dirnameo,'rs_k',parname,pariter,'.pdf'];
%saveas(fig3,figname,'pdf');
%hold off



%Recording Unseens vs Number of Couplings in the model
%outfile_un=[dirnameo,'unseens_k',parname,pariter];
%save(outfile_un,'n_unK','-ascii')
 
%Recording Variance vs Number of Couplings in the model
%outfile_var=[dirnameo,'variance_k',parname,pariter];
% save(outfile_var,'sigma_K','-ascii');

 %MutI: Selection of couplings  based on spearman contribution 
 
%  %contribution to spearman of a coupling
%  diffsp=diff(spvsk(:,2));
% 
%  [diffsp_sort,Idx_siffsp_sort]=sort(diffsp,'descend');
%  %pair index
%  IdxF2_sort=IdxF2(Idx_siffsp_sort);
[diffsp_sort,Idx_siffsp_sort]=sort(spvsk(:,2),'descend');
%pair index
  IdxF2_sort=IdxF2(Idx_siffsp_sort);
 % paired sites
 Jp_sort=zeros(K2max,2);
 Jp_sort(:,1)=indl(IdxF2_sort,1);
 Jp_sort(:,2)=indl(IdxF2_sort,2);
   

% recording coupling indexes in order of spearman contribution 
name=[dirnameo,'J',pariter];
fid = fopen(name,'w');

  for i=1:size(IdxF2_sort,2)
 
      %fprintf(fid,'%5d %5d %5d %5.4f',IdxF2_sort(i), pair_sort(i,1),pair_sort(i,2),diffsp_sort(i));
      fprintf(fid,'%5d %5d %5d %9.6f %9.3f',IdxF2_sort(i), Jp_sort(i,1),Jp_sort(i,2),diffsp_sort(i),Fapc(Jp_sort(i,1),Jp_sort(i,2)));
             fprintf(fid,'\n');
              
  end
   
  
% %plot energy vs fitness 
%figure(3)
%sl = scatter(fitnesscomp(:,1),fitnesscomp(:,2),'filled','SizeData',200);
%alpha(sl,.5)
%ylabel('-DH_{Potts}')
%xlabel('Fitness')
% figname=[dirnameo,'Scatter-EPottsvsFit',parname,pariter,'.fig']
 %savefig(figname)
 %figname=[dirnameo,'Scatter-EPottsvsFit',parname,pariter,'.pdf']
%saveas(figure(3),figname,'pdf');
%hold off  
  
   
%%%%%file in outputs to plot the energy vs fitness for the independent model one round %%%%%%%%%%%   


end