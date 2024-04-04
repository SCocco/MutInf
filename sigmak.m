%corretta aprile 2020 non avevo commentato alternativa nella var a 1,
%controllare che da lo stesso
function[sigma_in,sigma2,sigma_K,n_unK,nun2]=sigmak(freqmsa,contacts,wt,M_cut,q_kept,sigmah,sigmaj)
N=size(freqmsa,2);
q=21;
% calculation of variance according to the formula (and I will correlate the 2)
sigmain=zeros(21,N);
sigma2=zeros(21,N);
sigma1=zeros(21,N);
%number of contacts
nbc=zeros(N,1);
nun2=zeros(N,1);

for i=1:N
  %frequence of the wild type in the msa
  freqwtmsa(i)=freqmsa(wt(i),i);
  
end

%repeted matrix of Nxq for the wt frequences
 freqwtmsam=repmat(freqwtmsa,[q,1]);
 %rep2m is useful for the variance of a pair
 reg2m=repmat(1/(sigmah*M_cut),[q,q]);
  %parametri modello indipendente totale
%total number of unseen (for a fully connected model)
%  n_un2t(d)=0
for i=1:N
  %frequence of the wild type in the msa
  
    sigma1(:,i)=(1./(freqmsa(:,i)+ 1/(sigmah*M_cut))+1./(freqwtmsam(:,i) +1/(sigmah*M_cut)))/M_cut;
    %alternativa:
 %sigma1(:,i)=((1-freqmsa(:,i))./(freqmsa(:,i)+ 1/(sigmah*M_cut))+(1-freqwtmsam(:,i))./(freqwtmsam(:,i) +1/(sigmah*M_cut)))/M_cut;
%    for j=i+1:N
%    n_un2full(d)= n_unK(d)+(21-q_kept(i))*(21-q_kept(j));  
%    end
end
%unseen on the pairs questo e' per un modello fully connected. A noi 
%interessano gli unseen nei contatti scelti
%sigma_full_un(d)= n_un2full(d)*sigmah/((N*q))

%  hmsa=log(freqmsa+lambdah)'-log(freqwtmsam)';

sigma_in = mean(mean((1./(freqmsa + 1/(sigmah*M_cut)) + 1./(freqwtmsam + 1/(sigmah*M_cut))),1))/M_cut;

%sigma_in_appx(d)=mean(mean(1./(freqmsa+ 1/(sigmah*M_cut)),1))/M_cut;

%calcolo varianza a due approssimata stesso discorso
n_unK=0;
%variance of a pair
%sigmapair=cell(size(contacts,1),1);
for nc=1:size(contacts,1)
    ii=contacts(nc,1);
    jj=contacts(nc,2);
    nbc(ii)=nbc(ii)+1;
    nbc(jj)=nbc(jj)+1;
    
    sigma2(:,ii)=sigma2(:,ii)+(1./(freqmsa(:,ii).*freqwtmsam(:,jj)+ 1/(sigmaj*M_cut)))/M_cut+(1./(freqwtmsam(:,ii).*freqwtmsam(:,jj)+ 1/(sigmaj*M_cut)))/M_cut;
    %alternativa lorenzo
    %sigma2(:,ii)=sigma2(:,ii)+((1-freqmsa(:,ii).*freqwtmsam(:,jj))./(freqmsa(:,ii).*freqwtmsam(:,jj)+ 1/(sigmaj*M_cut)))/M_cut+((1-freqwtmsam(:,ii).*freqwtmsam(:,jj))./(freqwtmsam(:,ii).*freqwtmsam(:,jj)+ 1/(sigmaj*M_cut)))/M_cut;
    
    sigma2(:,jj)=sigma2(:,jj)+(1./(freqmsa(:,jj).*freqwtmsam(:,ii)+ 1/(sigmaj*M_cut)))/M_cut+(1./(freqwtmsam(:,ii).*freqwtmsam(:,jj)+ 1/(sigmaj*M_cut)))/M_cut;
    %alternativa lorenzo
    %sigma2(:,jj)=sigma2(:,jj)+((1-(freqmsa(:,jj).*freqwtmsam(:,ii)))./(freqmsa(:,jj).*freqwtmsam(:,ii)+ 1/(sigmaj*M_cut)))/M_cut+((1-freqwtmsam(:,ii).*freqwtmsam(:,jj))./(freqwtmsam(:,ii).*freqwtmsam(:,jj)+ 1/(sigmaj*M_cut)))/M_cut;
    
    n_unK= n_unK+(21-q_kept(ii))*(21-q_kept(jj));  

    nun2(ii)=nun2(ii)+(21-q_kept(ii))*(21-q_kept(jj));
    nun2(jj)=nun2(jj)+(21-q_kept(ii))*(21-q_kept(jj));
    %variance of a pair
  %  sigmapair{nc}=(1./(freqmsa(:,ii)*freqmsa(:,jj)'+reg2m)/M_cut+1./(freqwtmsam(:,ii)*freqmsa(:,jj)'+ reg2m)/M_cut+1./(freqmsa(:,ii)*freqwtmsam(:,jj)'+ reg2m)/M_cut+1./(freqwtmsam(:,ii)*freqwtmsam(:,jj)'+ reg2m)/M_cut)
    
end


%substract sigma1 for the sites with at least a contact and add sigma1 for
%the ones with no contact

cc=nbc>0;
for ii=1:N
    sigma2(:,ii)=sigma2(:,ii)+cc(ii)*(nbc(ii)-1)*sigma1(:,ii)+(1-cc(ii))*sigma1(:,ii);

%sigma_K_un(d)= sigma_K_un(d)+((1-cc(ii)+cc(ii)*(nbc(ii)-1))*(21-q_kept(ii)))/(N*q*(2*M*lambdah))+nun2(ii)/(N*q*(2*M*lambdaJ)) 
%sigma_K_pair_un(d)=sigma_K_pair_un(d)+nun2(ii)*sigmaj/(N*q);
%sigma_K_un(d)=sigma_K_un(d)+cc(ii)*(nbc(ii)-1)*sigma_in_un(d)+(1-cc(ii))*sigma_in_un(d)+nun2(ii)*sigmaj/(N*q)
%sigma_K_un(d)=sigma_K_un(d)+cc(ii)*(nbc(ii)-1)*sigma1(:,ii)+(1-cc(ii))*sigma1(:,ii)+nun2(ii)/(N*q*(2*M*lambdaJ))
end

sigma_K = mean(mean(sigma2,1));

%unseen on the pairs questo e' per un modello fully connected. A noi 
%interessano gli unseen nei contatti scelti
%sigma_K_pair2_un(d)=n_unK(d)*sigmaj/(N*q)
%sigma_K_un(d)= sigma_K_pair2_un(d)+sigma_K_un(d)

%predizioni songole mutazioni
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


end

