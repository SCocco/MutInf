function msa_num = conversion_align(filename)
    
    % convert the output of the MC simulation into a matrix alignement with number 1..20
    % usage: msa_num=conversion_align('./../../msa_generation/848msa.dat')
    %bug corretto 12/07/2018 sulla dichiarazione di MSA num
    output=importdata(filename);
    
    [M,N]=size(output);
    
    
    
     X = fastaread(filename);
     N = size(X(1).Sequence,2);
     M = size(X,1);
     msa_num = zeros(M,N);
    
     for i=1:M
      for j=1:N
        msa_num(i,j)=letter2number(X(i).Sequence(j));
        
      end
     end

function x=letter2number(a)
switch(a)

    case '-'
         x=1;
    case 'A'    
        x=2;    
    case 'C'    
        x=3;
    case 'D'
        x=4;
    case 'E'  
        x=5;
    case 'F'
        x=6;
    case 'G'  
        x=7;
    case 'H'
        x=8;
    case 'I'  
        x=9;
    case 'K'
        x=10;
    case 'L'  
        x=11;
    case 'M'
        x=12;
    case 'N'  
        x=13;
    case 'P'
        x=14;
    case 'Q'
        x=15;
    case 'R'
        x=16;
    case 'S'  
        x=17;
    case 'T'
        x=18;
    case 'V'
        x=19;
    case 'W'
        x=20;
    case 'Y'
        x=21;
    otherwise
        x=1;
end
end
 
   

%     % frequency in alignment
%     q=21;
%     N=27;
%     [M,N]=size(msa_num);
%     % alignment msa21 on N*Q sites
%     msa01=zeros(M,N*q);
%     for i=1:N
%         for a=1:q
%             % a.a. are labeled from 1 to q
%             msa01(:,(i-1)*q+a)=(msa_num(:,i)==a);
%         end
%     end
%     p=sum(msa01)/M; % p = frequence of amino acid on site i
%     pmat=reshape(p,q,N);
% 
%     bbb=[aaa;[0 0 0];[1 1 0];[1 0 1];[0 1 1];[1 0 0];[0 1 0];[0 0 1]];
%     figure
%     hold on
%     plot(pmat(1,:),'d','color',bbb(1,:),'MarkerSize',15)
%     for i=2:14
%     plot(pmat(i,:),'o','color',bbb(i,:))
%     end
%     for i=15:19
%     plot(pmat(i,:),'*','color',bbb(i-15+1,:),'MarkerSize',10)
%     end
%     plot(pmat(20,:),'o','color',bbb(1,:))
%     legend('C','M','F','I','L','V','W','Y','A','G','T','S','N','Q','D','E','H','R','K','P');
%     xlabel('Site')
%     ylabel('frequency')
%     title('Struc. #625')

end