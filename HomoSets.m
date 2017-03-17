function [ HomoSetsALL ] = HomoSets( Adj,HGThresh,minCN,maxCN )
%HOMOSETS Take as an input a NxN directed adjecancy matrix and the pvalue
%threshold of the hypergeometric test (the smaller it is, the greater the
%homogeneity of the output sets). The input should contain the range of common neighbors to the output sets.
%The output is a cell array that contains in every row:
% [X neuron] [Y neuron] [Z neurons] [Triad type of Z neurons] [Set type]
% for example: Sets = HomoSets(Adj,0.05,5,99)
n=length(Adj);
Adj(Adj<1) = 0;
Adj=Adj-eye(n);
Adj(Adj<0)=0;
Adj=Adj~=0;
con=Adj;
con2=Adj';

conn=zeros(n);
conn2=zeros(n);
for i=1:n
    temp=find(Adj(:,i));
    temp2=find(Adj(i,:));
    if (numel(temp) ~=0)
        conn(i,1:length(temp)) = temp;
    end
    if (numel(temp2) ~=0)
        conn2(i,1:length(temp2)) = temp2;
    end
    i
end

CommoneighborMat = zeros(n);
CommoneighborMat(CommoneighborMat==0)=0;
tens=[];
% 
for i=1:n
    current_indexes=1:n;
    con_shifted = circshift(con,i,1);
    con2_shifted = circshift(con2,i,1);
    current_indexes=circshift(current_indexes,i,2);
    Q=con&con_shifted;
    W=con2&con2_shifted;
    E=con&con2_shifted;
    R=con2&con_shifted;
    SUM=Q|W|E|R;
    QWER=sum(SUM,2);
    if (find(QWER) ~= 0)
        for l = 1:n
            if(QWER(l)>=minCN && QWER(l) <= maxCN)
                CommoneighborMat(l,current_indexes(l)) = QWER(l);
            end
        end
    end
    i
end
save('CommoneighborMat50.mat','CommoneighborMat','-v7.3');
for i=1:n
    for j=1:n
        if (i<j && CommoneighborMat(i,j) >= minCN && CommoneighborMat(i,j) <= maxCN  && (Adj(i,j) || Adj(j,i)))
            tens(size(tens,1) + 1,:) = [i j];
        end
    end
    i
end
save('PairsWith50PlusCN_BlueBrain.mat','tens','-v7.3');

%% HYPERGEOMETRIC TEST
HomoSetsALL=[];
load('PairsWith50PlusCN_BlueBrain.mat')
Slices = [1:5000:length(tens) length(tens)];
tensALL=tens;
for i=1:length(Slices)-1
    tens = tensALL(Slices(i):Slices(i+1),:);
    List={};
    superow=1;
    supermat=zeros(n*size(tens,1),9);
    supermat2=zeros(n*size(tens,1),4);
    for i=1:size(tens,1)
        currentPair=[];
        Q=intersect(conn(tens(i,1),:),conn(tens(i,2),:));
        W=intersect(conn2(tens(i,1),:),conn2(tens(i,2),:));
        E=intersect(conn(tens(i,1),:),conn2(tens(i,2),:));
        R=intersect(conn2(tens(i,1),:),conn(tens(i,2),:));
        QWER=[Q W E R];
        QWER=unique(QWER);
        QWER(QWER==0) = [];
        for k=QWER
            if (tens(i,1)~=k && tens(i,2)~=k && ismember(k,QWER))
                supermat(superow,:)= [tens(i,1) tens(i,2) k 0 0 0 0 0 0];
                supermat2(superow,:)= [tens(i,1) tens(i,2) k 0];
                if (Adj(tens(i,1),tens(i,2)) && ~Adj(tens(i,2),tens(i,1)))
                    lala=1;
                elseif (~Adj(tens(i,1),tens(i,2)) && Adj(tens(i,2),tens(i,1)))
                    lala=2;
                elseif (Adj(tens(i,1),tens(i,2)) && Adj(tens(i,2),tens(i,1)))
                    lala=3;
                elseif (~Adj(tens(i,1),tens(i,2)) && ~Adj(tens(i,2),tens(i,1)))
                    lala=4;
                end
                za=ismember(k,intersect(conn2(tens(i,1),:),conn2(tens(i,2),:)));
                if (za)
                    supermat(superow,4) = 1;
                end
                zb=ismember(k,intersect(conn(tens(i,1),:),conn(tens(i,2),:)));
                if (zb)
                    supermat(superow,5) = 1;
                end
                zc=ismember(k,intersect(conn(tens(i,1),:),conn2(tens(i,2),:)));
                if (zc)
                    supermat(superow,6) = 1;
                end
                zd=ismember(k,intersect(conn2(tens(i,1),:),conn(tens(i,2),:)));
                if (zd)
                    supermat(superow,7) = 1;
                end
                ze=ismember(k,intersect(conn2(tens(i,1),:),conn(tens(i,1),:)));
                if (ze)
                    supermat(superow,8) = 1;
                end
                zf=ismember(k,intersect(conn2(tens(i,2),:),conn(tens(i,2),:)));
                if (zf)
                    supermat(superow,9) = 1;
                end
                if (supermat(superow,1) ~= 0)
                    if (za && ~zb && ~zc && ~zd && ~ze && ~zf && lala==1)
                        supermat2(superow,4)=1;
                    elseif (~za && ~zb && ~zc && zd && ~ze && ~zf && lala==1)
                        supermat2(superow,4)=2;
                    elseif (za && ~zb && ~zc && zd && ~ze && zf && lala==1)
                        supermat2(superow,4)=3;
                    elseif (~za && ~zb && zc && ~zd && ~ze && ~zf && lala==1)
                        supermat2(superow,4)=4;
                    elseif (~za && zb && ~zc && ~zd && ~ze && ~zf && lala==1)
                        supermat2(superow,4)=5;
                    elseif (~za && zb && zc && ~zd && ~ze && zf && lala==1)
                        supermat2(superow,4)=6;
                    elseif (~za && zb && ~zc && zd && ze && ~zf && lala==1)
                        supermat2(superow,4)=7;
                    elseif (za && ~zb && zc && ~zd && ze && ~zf && lala==1)
                        supermat2(superow,4)=8;
                    elseif (za && zb && zc && zd && ze && zf && lala==1)
                        supermat2(superow,4)=9;
                    elseif (za && ~zb && ~zc && ~zd && ~ze && ~zf && lala==2)
                        supermat2(superow,4)=1;
                        q=supermat2(superow,2);
                        supermat2(superow,2)=supermat2(superow,1);
                        supermat2(superow,1)=q;
                    elseif (~za && ~zb && ~zc && zd && ~ze && ~zf && lala==2)
                        supermat2(superow,4)=4;
                        q=supermat2(superow,2);
                        supermat2(superow,2)=supermat2(superow,1);
                        supermat2(superow,1)=q;
                    elseif (za && ~zb && ~zc && zd && ~ze && zf && lala==2)
                        supermat2(superow,4)=8;
                        q=supermat2(superow,2);
                        supermat2(superow,2)=supermat2(superow,1);
                        supermat2(superow,1)=q;
                    elseif (~za && ~zb && zc && ~zd && ~ze && ~zf && lala==2)
                        supermat2(superow,4)=2;
                        q=supermat2(superow,2);
                        supermat2(superow,2)=supermat2(superow,1);
                        supermat2(superow,1)=q;
                    elseif (~za && zb && ~zc && ~zd && ~ze && ~zf && lala==2)
                        supermat2(superow,4)=5;
                        q=supermat2(superow,2);
                        supermat2(superow,2)=supermat2(superow,1);
                        supermat2(superow,1)=q;
                    elseif (~za && zb && zc && ~zd && ~ze && zf && lala==2)
                        supermat2(superow,4)=7;
                        q=supermat2(superow,2);
                        supermat2(superow,2)=supermat2(superow,1);
                        supermat2(superow,1)=q;
                    elseif (~za && zb && ~zc && zd && ze && ~zf && lala==2)
                        supermat2(superow,4)=6;
                        q=supermat2(superow,2);
                        supermat2(superow,2)=supermat2(superow,1);
                        supermat2(superow,1)=q;
                    elseif (za && ~zb && zc && ~zd && ze && ~zf && lala==2)
                        supermat2(superow,4)=3;
                        q=supermat2(superow,2);
                        supermat2(superow,2)=supermat2(superow,1);
                        supermat2(superow,1)=q;
                    elseif (za && zb && zc && zd && ze && zf && lala==2)
                        supermat2(superow,4)=9;
                        q=supermat2(superow,2);
                        supermat2(superow,2)=supermat2(superow,1);
                        supermat2(superow,1)=q;
                    elseif (za && ~zb && ~zc && ~zd && ~ze && ~zf && lala==3)
                        supermat2(superow,4)=10;
                    elseif (~za && ~zb && ~zc && zd && ~ze && ~zf && lala==3)
                        supermat2(superow,4)=11;
                    elseif (za && ~zb && ~zc && zd && ~ze && zf && lala==3)
                        supermat2(superow,4)=12;
                    elseif (~za && ~zb && zc && ~zd && ~ze && ~zf && lala==3)
                        supermat2(superow,4)=11;
                        %                                         q=supermat2(superow,2);
                        %                     supermat2(superow,2)=supermat2(superow,1);
                        %                     supermat2(superow,1)=q;
                        %  supermat2(superow,4)=16;
                    elseif (~za && zb && ~zc && ~zd && ~ze && ~zf && lala==3)
                        supermat2(superow,4)=13;
                    elseif (~za && zb && zc && ~zd && ~ze && zf && lala==3)
                        supermat2(superow,4)=14;
                    elseif (~za && zb && ~zc && zd && ze && ~zf && lala==3)
                        supermat2(superow,4)=14;
                        %                                         q=supermat2(superow,2);
                        %                     supermat2(superow,2)=supermat2(superow,1);
                        %                     supermat2(superow,1)=q;
                        % supermat2(superow,4)=17;
                    elseif (za && ~zb && zc && ~zd && ze && ~zf && lala==3)
                        supermat2(superow,4)=12;
                        %                                         q=supermat2(superow,2);
                        %                     supermat2(superow,2)=supermat2(superow,1);
                        %                     supermat2(superow,1)=q;
                        % supermat2(superow,4)=18;
                    elseif (za && zb && zc && zd && ze && zf && lala==3)
                        supermat2(superow,4)=15;
                    elseif (za && ~zb && ~zc && ~zd && ~ze && ~zf && lala==4)
                        supermat2(superow,4)=51;
                    elseif (~za && ~zb && ~zc && zd && ~ze && ~zf && lala==4)
                        supermat2(superow,4)=52;
                    elseif (za && ~zb && ~zc && zd && ~ze && zf && lala==4)
                        supermat2(superow,4)=53;
                    elseif (~za && ~zb && zc && ~zd && ~ze && ~zf && lala==4)
                        supermat2(superow,4)=52;
                    elseif (~za && zb && ~zc && ~zd && ~ze && ~zf && lala==4)
                        supermat2(superow,4)=54;
                    elseif (~za && zb && zc && ~zd && ~ze && zf && lala==4)
                        supermat2(superow,4)=55;
                    elseif (~za && zb && ~zc && zd && ze && ~zf && lala==4)
                        supermat2(superow,4)=55;
                    elseif (za && ~zb && zc && ~zd && ze && ~zf && lala==4)
                        supermat2(superow,4)=53;
                    elseif (za && zb && zc && zd && ze && zf && lala==4)
                        supermat2(superow,4)=56;
                    else
                        disp('wtf');
                    end
                end
                currentPair= [currentPair supermat2(superow,4)];
                superow = superow + 1;
            end
        end
        if (currentPair(1) <16)
            List{length(List) + 1} = currentPair;
            Pairs(1:2,length(List)) = supermat2(superow-1,1:2);
        end
        i
    end
    NumOfEnP=zeros(1,1);
    EnTriadsP=zeros(1001,15);
    EnSets={};
    AllSetsALL=List;
    TriadsALL=cell2mat(AllSetsALL);
    TriadsALL=hist(TriadsALL,1:15);
    Triads=cell2mat(AllSetsALL);
    K=TriadsALL;
    Triads=hist(Triads,1:15);
    NUni=sum(TriadsALL(1,1:9));
    NBi=sum(TriadsALL(1,10:15));
    counte=1;
    AllSets=AllSetsALL;
    for i=1:length(AllSetsALL)
        if (AllSets{i}(1) < 10)
            for j=1:9
                k=sum(AllSets{i} == j);
                n=length(AllSets{i});
                pvalue=1-hygecdf(k,NUni,K(j),n);
                PercentFromSet=k/n;
                if (pvalue<HGThresh && k/n>0.5)
                    %              disp([num2str(AllSets{i}) ' enriced for ' num2str(j) ' with pvalue=' num2str(pvalue)])
                    NumOfEnP=NumOfEnP+1;
                    EnTriadsP(j)=EnTriadsP(j)+1;
                    EnSets{counte} = AllSets{i};
                    EnPairs(1:3,counte)=[Pairs(1:2,i); j] ;
                    counte=counte+1;
                end
            end
        else
            for j=10:15
                k=sum(AllSets{i} == j);
                n=length(AllSets{i});
                pvalue=1-hygecdf(k,NBi,K(j),n);
                PercentFromSet=k/n;
                if (pvalue<HGThresh && k/n>0.5)
                    %              disp([num2str(AllSets{i}) ' enriced for ' num2str(j) ' with pvalue=' num2str(pvalue)]);
                    NumOfEnP=NumOfEnP+1;
                    EnTriadsP(j)=EnTriadsP(j)+1;
                    EnSets{counte} = AllSets{i};
                    EnPairs(1:3,counte)=[Pairs(1:2,i); j] ;
                    counte=counte+1;
                end
            end
        end
    end
    
    supermat2( ~any(supermat2,2), : ) = [];
    supermat5En=[];
    counter=1;
    for i=1:length(EnPairs)
        for j=1:length(supermat2)
            if (EnPairs(1,i) == supermat2(j,1) && EnPairs(2,i) == supermat2(j,2))
                supermat5En(counter,:)=supermat2(j,:);
                counter=counter+1;
            end
        end
    end
    for i=1:NumOfEnP
        HomoSets{i,1} = EnPairs(1,i);
        HomoSets{i,2} = EnPairs(2,i);
        temp=[];
        temp2=[];
        for j=1:length(supermat5En)
            if (supermat5En(j,1) == EnPairs(1,i) && supermat5En(j,2) == EnPairs(2,i))
                temp=[temp supermat5En(j,3)];
                temp2=[temp2 supermat5En(j,4)];
            end
        end
        HomoSets{i,3} = temp;
        HomoSets{i,4} = temp2;
        HomoSets{i,5} = EnPairs(3,i);
        
    end
    HomoSetsALL=[HomoSetsALL; HomoSets];
    save('HomoSetsALL50.mat','HomoSetsALL');
end
end
