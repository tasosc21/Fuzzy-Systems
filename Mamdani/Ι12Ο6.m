%-----------------------------------------------------------------------
% Αναγνώριση συστήματος με ασαφές σύστημα τύπου Mamdani.
%-----------------------------------------------------------------------
clc;
close all;
clear all;
%-----------------------------------------------------------------------
x=0:1:360;
y=linspace(-0.4,1,361);
z=0.5*sin(2*pi*x/360)+0.5*(sin(2*pi*(x+30)/360)).^2;
count=1;
%-----------------------------------------------------------------------
for i=0.1:0.1:0.7


    %-----------------------------------------------------------------------
    % Για σωστή αρίμθηση των εικόνων
    if count == 1
        cc = 3;
    else
        cc = count + 2;
    end
    %-----------------------------------------------------------------------
    % Βασικές Παράμετροι (όχι όλες)
    ruleNumber = 12;
    inputSets = ruleNumber;
    outputSets = 6;
    %-----------------------------------------------------------------------
    % Ομοιόμορφη κατανομή του χώρου εισόδου
        [Alpha,Beta,Gamma]=tri_MF_partition(min(x),max(x),inputSets,i);
        for (k=1:inputSets)
            mf_x(k,:)=tri_MF(x,Alpha(k),Beta(k),Gamma(k));
        end
    if count==1
        figure(1);
        STR=sprintf('input space uniform partition');
        axis([min(x) max(x) 0 1.05]);   title(STR);   hold on
        for (k=1:inputSets)
            plot(x,mf_x(k,:))
        end
    end
    %-----------------------------------------------------------------------
    % Ομοιόμορφη κατανομή του χώρου εξόδου
        [Alpha,Beta,Gamma]=tri_MF_partition(min(y),max(y),outputSets,i);
        for (k=1:outputSets)
            mf_y(k,:)=tri_MF(y,Alpha(k),Beta(k),Gamma(k));
        end
    if count==1
        figure(2);
        STR=sprintf('output space uniform partition');
        axis([min(y) max(y) 0 1.05]);   title(STR);   hold on
        for (k=1:outputSets)
            plot(y,mf_y(k,:))
        end
    end
    %-----------------------------------------------------------------------
    % Βάση κανόνων, αντιστοίχηση των ασαφών συνόλων στα σύνολα των κανόνων
    In(1,:)=mf_x(1,:);    In(2,:)=mf_x(2,:);  In(3,:)=mf_x(3,:);  In(4,:)=mf_x(4,:);
    In(5,:)=mf_x(5,:);    In(6,:)=mf_x(6,:);  In(7,:)=mf_x(7,:); In(8,:)=mf_x(8,:); 
    In(9,:)=mf_x(9,:);    In(10,:)=mf_x(10,:); In(11,:)=mf_x(11,:); In(12,:)=mf_x(12,:); 
    Out(1,:)=mf_y(3,:);   Out(2,:)=mf_y(5,:); Out(3,:)=mf_y(6,:); Out(4,:)=mf_y(6,:);
    Out(5,:)=mf_y(4,:);   Out(6,:)=mf_y(3,:); Out(7,:)=mf_y(3,:); Out(8,:)=mf_y(3,:);
    Out(9,:)=mf_y(2,:);   Out(10,:)=mf_y(1,:); Out(11,:)=mf_y(1,:); Out(12,:)=mf_y(3,:);

    %-----------------------------------------------------------------------
    for j=1:length(x)
        for l=1:ruleNumber
            Bbar(l,:)=min(In(l,j),Out(l,:));
        end    
        total=Bbar(1,:);
        for l=2:ruleNumber
            total=max(total,Bbar(l,:));       
        end
        output(j)=sum(y.*total)/sum(total);
    end
    %-----------------------------------------------------------------------
    figure(cc);
    plot(x,z,'b',x,output,'r--');
    STR=sprintf('Overlapping coef. %.2f:   Blue solid line: Cos, Red dotted line: fuzzy system''s output',i);
    title(STR);
    diff=sum(abs(z-output));
    fprintf("g=%f->%f\n",i,diff)


   count=count+1;
end
%-----------------------------------------------------------------------
