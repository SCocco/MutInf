%MuTINF
% Copyright Simona Cocco 3/4/2024 
%please cite the paper 
%Functional effects of mutations in proteins can be predicted and interpreted by guided selection of
%sequence covariation information. by Simona Cocco, Lorenzo Posani, Rémi Monasson. 

%Program to infer the field and couplings  with the Mutinf Pipeline and 2-site clusters
%approximation , starting on a list of linkd given  in input 
%it compares the fitness cost of single mutations given in mutations
%with the energetic cost in the inferred model  plot and calculate the spearman coefficient. 
%olso the case of one link

%%%%%%%%%%%%%%%%%%%%%%%List of analyzed Proteins%%%%%%%%%%%%%%%
%different initial list and selection procedure
%by default plm_DCA or random or by mutationally guided
%to obtain the list according to mutational scan onelink=true run
%spACE_only_1J.m
onelink=true;
%onelink=false

%random initial list
%randomlist=true;
%randomlist=false;
%no Kmax=6000 (to speed up only in the case of blacta and Bgluco)
Kcut=false;
%Kcut=true
%from a round to the next one only dgt>0
noprune=false;
%you need to have parpool, if not comment parpool command (line 325)
n_core=3;
%n_core=20


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
    protname='Bgluco'
 elseif(prot==9) 
    protname='Betalactamase'
 end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Input DATA %%%%%%%%%%%%%%%%%%%%% 
%input directory
inputsMSA=['./',protname,'/Dati_inputs_MSAMarks/'];
%output directory

if onelink
dirnameo = ['./',protname,'/Outputs_MSAMarks/1Link_p/'];
elseif randomlist

dirnameo = ['./',protname,'/Outputs_MSAMarks/Random_p/']; 
elseif Kcut

dirnameo = ['./',protname,'/Outputs_MSAMarks/Fapc_cut_p/']; 
else

dirnameo = ['./',protname,'/Outputs_MSAMarks/Fapc_p/'];  
end

if(onelink) 
IdxF= importdata([dirnameo,'J_sorted']);
end


% Frobenious Norms calculated with Plm with APC (Average Product Corrections)
Fapc=importdata([inputsMSA,'frobenious/Fapc_PLM_lambda001_theta02_pcut0_cons.mat']);


% Input Alignment
align = conversion_align([inputsMSA,'sequences.faa']);
[M,N] = size(align);

%reweighting 
theta=0.2;
w = importdata([inputsMSA,'weights/reweighting_0.200.dat']);
%theta=0;
%w=ones(M,1);


% Wild-Type
wt=conversion_align(['./',protname,'/wt.faa']);

%read the mutational data set
X=fastaread(['./',protname,'/fitness_mutations.faa']);
%RL401 function
%X=fastaread(['./',protname,'/fitness_mutations_function.faa']);

[M_test] = size(X,1);

align_test=zeros(M_test,N);

fitness=zeros(M_test,1);

for i=1:M_test
    fitness(i)=sscanf(X(i).Header,'%f');
    for j=1:N
        align_test(i,j)=letter2number(X(i).Sequence(j));
        
    end
end
%random permutation of data set to select only 80 for train
%indperm=randperm(M_test);
%ind_train=indperm(1,1:per_train*M_test);
%ind_test=indperm(1,per_train*M_test+1:M_test);

%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters in the loop over nb of couplings in the model K2 %%%%%%%%

%Parameters to fix the iterations in the number of couplings in the inferred model
% at each iteration the spearman between the theoretical prediction and the
% experimental ones is given in output.
%K2min: minimal number of couplings in the inferred model
K2min=0;
%K2max: maximal number of couplings in the inferred model
%K2max=3000
if Kcut
K2max=min(6000,N*(N-1)/2);
else
K2max=N*(N-1)/2;
end



%%%%%%%%%%%%%%%%%%%%%%% Parameter of the model%%%%%%%%%%%%%%%%%%%%%%%%

%number of Potts states
q=21;
%pcut (inference done on colors with p>pcut)
pcut=0;
% threshold of number of states at which sites are taken for coupling parameters 
%qcut=10;
qcut=0;
%regolarization strengths 
lambdah=0.1/M;
lambdaJ=10./M;
sigmah=1;
sigmaj=1;

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
q_kept_gauge=q_kept-1
qtotsum_gauge=zeros(1,N+1);
qtotsum_gauge(1,2:N+1)=cumsum(q_kept_gauge);

for i=1:N
    disscut_ind(i)=disscut(i)+qtotsum(i);
    %disscut_ind(i)=disscut(i)+qtotsum(i);
    %   regcut_ind(i)=qtt(i)+q_tot(i);
    pmin(i)=freqmsa(index(i,disscut(i)),i);
end
%estimation  of the variance of independent model related to the unseens
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
      %  if (onelink==false)
    Fv(l)=Fapc(i,j);
       % end
    %Fv(l)=Fapc(l,3);
    end
end
    



     
  

for i=1:N
    for a=1:q
        %freq of the aa in the test data set
        freq_test(i,a)=length(find(align_test(:,i)==a))/M;
    end
end

for i=1:N
    %frequence of the wild type in the msa
    freqwtmsa(i)=freqmsa(wt(i),i);
end
%repeted matrix of Nxq for the wt frequences
freqwtmsam=repmat(freqwtmsa,[q,1]);


%%%%%%%%%%SELECT a vector of  selected and ordered LINKS %%%%%%%

%redo for 3 cycles at each time keeping only >0

%%%In case we want to put couplings on Structural Contacts%%%%%%%%%%%%%
%contacts=importdata('list_contacts625_A.dat');
if(onelink==false) 
    if(randomlist)
    IdxF=randperm(K2max);
    else
    %according to frobenius    
    [Fsort,IdxF]=sort(Fv,'descend');
    end
end
Idx_s=IdxF(1:K2max);


%K2max=size(IdxF2,2)
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

%Extracting 1 and 2 points correlations

NN=size(msacut_gauge,2);

wm=repmat(w,[1,NN]);

%reweighted probabilities
msacut_gauge_w=msacut_gauge.*wm;
p=sum(msacut_gauge_w)/sum(w);
%reweighted pairwise probabilities
p2=msacut_gauge'*msacut_gauge_w/sum(w);
%senza pesi
%p=sum(msacutg)/M;

% Fields
h0=cell(1,N);
h=cell(1,N);
%ga_h=ga_J/100;
for i=1:N
    pa=p(qtotsum_gauge(i)+1:qtotsum_gauge(i)+q_kept_gauge(i));
    x0=zeros(1,q_kept_gauge(i));
     f = @(x)cross_entropy_c1(x,pa,lambdah,q_kept_gauge(i));
     %for matlab versions after 2017
    %options = optimoptions(@fminunc,'SpecifyObjectiveGradient',true);
    options=optimoptions(@fminunc,'GradObj','on');
    x = fminunc(f, x0, options);
    h{i}=x;
    h0{i}=x;
end
for nbr=1:6

%Calculation of parameters x for the maximal size 
K2=K2max;
contacts=zeros(K2max,2);
contacts(1:K2max,1)=indl(Idx_s(1:K2max),1);
contacts(1:K2max,2)=indl(Idx_s(1:K2max),2);

%nbc=size(contacts,1);
%J=cell(N-1,N);
J=cell(K2max,1);

xh=cell(K2max,1);
%calculation of the couplings and fields contributions for maximal sets

tic
parpool('local',n_core);
%for cc=1:K2max
parfor cc=1:K2max
    
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
  
  
  
end
delete(gcp);
toc;


%%%%%%%%%%Iterations in the number of Couplings in the model%%%%%%%%%%%%

J2=[];
cc=0;
h=h0;
spK_train=zeros(K2max+1,2);
sigmaK_train=zeros(K2max+1,2);
r2K_train=zeros(K2max+1,2);
nunK_train=zeros(K2max+1,2);






for cl=1:K2max+1
    %the first iteration cl=1,cc=0 correspond to the independent model 
    if (cc>0)   
    ii=contacts(cc,1);
    jj=contacts(cc,2);
   
    h{ii}=h{ii}+xh{cc}(1:q_kept_gauge(ii))-h0{ii};  
    h{jj}=h{jj}+xh{cc}(q_kept_gauge(ii)+1:q_kept_gauge(ii)+q_kept_gauge(jj))-h0{jj};
    J2=J(1:cc);
    end
    
    % single site variance sigma1, variance due to pairs sigma2, number of
    % unseens due to pairs n_unK, number of pairs of aa unseen for each site
    
     contacts_sel=contacts(1:cc,:);
    [sigma1,sigma2,sigma_Kd,n_unKd,nun2]=sigmak(freqmsa,contacts_sel,wt,M_cut,q_kept,sigmah,sigmaj);
   
    
    %insert  the 0 corresponding to the fields and couplings of the diss in the
    %cut alignement
    [J_diss,h_diss]=read_j_cons_from_sce_cl2(J2,h,q_tot,disscut,N,contacts_sel);
    
    %put the parameters  in wt gauge
    
    [J_wt,h_wt]=map2wt4redt(wtcut,J_diss,h_diss,q_tot,N,contacts_sel);
    
    
    
    %Evaluate the energy of the mutants msa_test coupled and grouped model
    %Attenzione per validazione su train test
    is_mutated = zeros(N*q, 1);
    %train mantengo nomi vecchi
    enh=zeros(1,M_test);
    for kp=1:M_test
        %selection of all mutations
        k=kp;
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
                enh(kp)=h_wt{idd}(refi)+log(lambdah/frefi) ;
                % se e' nel regrouped symbol
            else
                enh(kp)=h_wt{idd}(refi)+log(freqmsa(align_test(k,idd),idd)/frefi);
            end
            
        else
            enh(kp)=h_wt{idd}(indi_mut);
        end
    end
    variance = mean(mean(sigma2(is_mutated>0)));
    n_unseens=n_un;
    if (cc>0)
    n_unseens=n_unKd+n_un;
    end
   
    %sigma_Kd
    
    fitnesscomp_train=zeros(M_test,2);
    fitnesscomp_train(:,2)=fitness;
    fitnesscomp_train(:,1)=-enh';
    %r2=corr2(fitnesscomp(:,1),fitnesscomp(:,2))^2;
    r2mod=corr(fitness,enh')^2;
    [rsp_train,pvalsp_train] = corr(fitness,enh','Type','Spearman');
    
    spK_train(cl,1)=cc;
    spK_train(cl,2)=rsp_train;
    
    r2K_train(cl,2)=r2mod;
    r2K_train(cl,1)=cc;
    
    sigmaK_train(cl,1)=cc;
    sigmaK_train(cl,2)=variance;
    
    nunK_train(cl,1)=cc;
    nunK_train(cl,2)=n_unseens;
    
    
    %%%%%%%%%%%test%%%%%%%%%%%%%%%%%%%%%%%%%
    cc=cc+1;
end

% %Indipendent model is the first point in the iteration 
% %Speraman
% 
% vx=[0,spvsk(:,1)'];
% vy=[rsind,spvsk(:,2)'];
% spvskf(:,1)=vx';
% spvskf(:,2)=vy';
% 
% %numero unseens
% nx=[0,n_unK(:,1)'];
% ny=[n_un,n_unK(:,2)'];
% n_unKf(:,1)=nx'
% n_unKf(:,2)=ny'
% %varianza
% sx=[0,sigma_K(:,1)'];
% sy=[sigma1,sigma_K(:,2)']
% sigma_Kf(:,1)=sx';
% sigma_Kf(:,2)=sy';



%%%%%% Output names changing  the model  and iterarions parameters %%%%%%%%%%%%%%
parname=[];
pariter=[];
%Recording Sperman train vs Number of Couplings in the model
%outfile_rs=[dirnameo,'rs',parname,pariter,'train'];
outfile_rs_train=[dirnameo,'rs',parname,pariter,'_',num2str(nbr-1)]
%save(outfile_rs,'spvskf','-ascii');
save(outfile_rs_train,'spK_train','-ascii');

outfile_r2_train=[dirnameo,'r2',parname,pariter,'_',num2str(nbr-1)]
%save(outfile_rs,'spvskf','-ascii');
save(outfile_r2_train,'r2K_train','-ascii');


% %plots spearman vs k
% fig3=figure(3)
% hold on
% plot(spvsk_train(:,1),spvsk_train(:,2),'b.','MarkerSize',20);
% plot(spvsk_test(:,1),spvsk_test(:,2),'r.','MarkerSize',20);
% 
% %txt=[name,'r2=',num2str(r2),'rs=',num2str(rsp)]
% %title(txt)
% xlabel('Number of couplings K')
% ylabel('Spearmann')
% %figname=[dirnameo,'spvsk-',protname,'_dk10_k1055_withalphah.fig']
% %saveas(fig3,figname,'fig');
% %savefig(figname)
% figname=[dirnameo,'spvsk_train',parname,pariter,'train.pdf'];
% saveas(fig3,figname,'pdf');
% hold off


% %Recording Unseens vs Number of Couplings in the model
% outfile_un=[dirnameo,'unseens_train',parname,pariter];
% save(outfile_un,'n_unK_train','-ascii')
%  
% %Recording Variance vs Number of Couplings in the model
 outfile_var=[dirnameo,'variance',parname,pariter,'_',num2str(nbr-1)];
  save(outfile_var,'sigmaK_train','-ascii');
%  
% 
% %Recording Unseens vs Number of Couplings in the model
% outfile_un=[dirnameo,'unseens_test',parname,pariter];
% save(outfile_un,'n_unK_test','-ascii')
%  
% %Recording Variance vs Number of Couplings in the model
% outfile_var=[dirnameo,'variance_test',parname,pariter];
%  save(outfile_var,'sigma_K_test','-ascii');
 


 %contribution to spearman of a coupling train
 diffsp_train=diff(spK_train(:,2));
 %contribution to r2 of a coupling train
 diffr2_train=diff(r2K_train(:,2));
 
% recording coupling indexes in order of spearman contribution 
outfile_indexJ_sort_train=[dirnameo,'J_sort',parname,pariter,num2str(nbr-1)];
fid = fopen(outfile_indexJ_sort_train,'w');
for i=1:K2max
      fprintf(fid,'%5d %5d %5d %9.5f %9.8f %9.5f',Idx_s(i),indl(Idx_s(i),1),indl(Idx_s(i),2),diffsp_train(i),Fapc(indl(Idx_s(i),1),indl(Idx_s(i),2)),diffr2_train(i));
      %fprintf(fid2,'%5d %5d %5d',IdxF2_sorted_train(i),indl(IdxF2_sorted_train(i),1),indl(IdxF2_sorted_train(i),2));
             fprintf(fid,'\n');
              %fprintf(fid2,'\n');
              %recording K2max
              
 
end 

 [diffsp_sort_train,Idx_diffsp_sort_train]=sort(diffsp_train,'descend');
 %K2max_train=length(find(diffsp_sort_train>0.0001));

%attention the independent value for the spearman is too small
 Idx_s=Idx_s(Idx_diffsp_sort_train);
 
if (noprune==false) 
 K2max=length(find(diffsp_sort_train>0))
end





% [diffr2_sort_train,Idx_diffr2_sort_train]=sort(diffr2_train,'descend');
 %K2max_train=length(find(diffsp_sort_train>0.0001));

% %attention the independent value for the spearman is too small
%  Idx_r2=Idx_s(Idx_diffr2_sort_train);
%  
% % recording coupling indexes in order of r2 contribution 
% outfile_indexJ_sort_train_r2=[dirnameo,'J_sorted',parname,pariter,'_r2_',num2str(nbr)];
% fid = fopen(outfile_indexJ_sort_train,'w');
% for i=1:K2max
%       fprintf(fid,'%5d %5d %5d %9.5f %9.8f',Idx_r2(i),indl(Idx_r2(i),1),indl(Idx_r2(i),2),diffr2_sort_train(i),Fapc(indl(Idx_r2(i),1),indl(Idx_r2(i),2)));
%       %fprintf(fid2,'%5d %5d %5d',IdxF2_sorted_train(i),indl(IdxF2_sorted_train(i),1),indl(IdxF2_sorted_train(i),2));
%              fprintf(fid,'\n');
%               %fprintf(fid2,'\n');
%               %recording K2max
%               
%  
% end 
% 

 [rstar,kstar]=max(spK_train(:,2));
      spmax05=rstar*(1-0.005);
      ind=find(spK_train(:,2)>spmax05);
      kstar2=ind(1);
      rstar2=spK_train(kstar2,2);
      
   [r2star,k2star]=max(r2K_train(:,2));
      r2max05=r2star*(1-0.005);
      indr2=find(r2K_train(:,2)>r2max05);
      k2star2=indr2(1);
      r2star2=r2K_train(k2star2,2) ;  
      
      
pariter_train_prev=['_',num2str(nbr-1)]
 k2maxsave=[dirnameo,'K2maxsp',parname,pariter_train_prev];
 fidk =fopen(k2maxsave,'w');
 fprintf(fidk,'%5d %9.8f', kstar2,rstar2,k2star2,spK_train(k2star2,2));
 fprintf(fidk,'\n');
 
 
 k2maxsaver2=[dirnameo,'K2maxr2',parname,pariter_train_prev];
 fidkr2 =fopen(k2maxsaver2,'w');
 fprintf(fidkr2,'%5d %9.8f',k2star2,r2star2,kstar2,r2K_train(kstar2,2));
 fprintf(fidkr2,'\n'); 
  



  
% % recording cut J
% outfile_indexcut=[dirnameo,'IndexJ_cut',parname,pariter];
% fid = fopen(outfile_indexcut,'w');
%    for i=1:size(indcuts,2)
%              fprintf(fid,'%5d %4.3f',indcuts(i),diffsp(i));
%              fprintf(fid,'\n');
%    end   
%    
% end
   

end
% outfile_indexJ_sort_train=[dirnameo,'J_sorted',parname,pariter,'_',num2str(nbr)];
% fid = fopen(outfile_indexJ_sort_train,'w');
% [diffsp_sort_train,Idx_diffsp_sort_train]=sort(diffsp_train,'descend');
% Idx_s=Idx_s(Idx_diffsp_sort_train);
% 
% for i=1:K2max
%       fprintf(fid,'%5d %5d %5d %9.5f %9.8f,%9.5f',Idx_s(i),indl(Idx_s(i),1),indl(Idx_s(i),2),diffsp_sort_train(i),Fapc(indl(Idx_s(i),1),indl(Idx_s(i),2)),diffr2_train(Idx_diffsp_sort_train));
%       %fprintf(fid2,'%5d %5d %5d',IdxF2_sorted_train(i),indl(IdxF2_sorted_train(i),1),indl(IdxF2_sorted_train(i),2));
%              fprintf(fid,'\n');
%               %fprintf(fid2,'\n');
%               %recording K2max
%               
 
end 
