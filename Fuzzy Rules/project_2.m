%%-------------------------------------------------------------------------
%               Ασαφή Συδτήματα & Εξελικτική Υπολογιστική
%                           Δεύτερη Εργασία
%               Ασαφείς κανόνες-βάσεις ασαφών κανόνων
%       εξαγωγή συμπεράσματος και έξοδος του ασαφούς συστήματος
%
% Ονοματεπώνυμο: Αναστάσιος Χριστοδουλίδης
% Α.Μ.: aivc22019
%%-------------------------------------------------------------------------

%-------------------------------------------------------------------------

close all; clear all; clc; format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------------------------------------------------
%% Μέρος 1
%%------------------------------------------------------------------------

% Ασαφής Κανόνας: IF x is A THEN y is B

x=linspace(0,20,2001); % Πεδίο Ορισμού x, βήμα 0.01
y=linspace(0,100,10001); % Πεδίο Ορισμού y, βήμα 0.01

% Δημιουργία Ασαφών Συνόλων Εισόδου και Εξόδου
A=tri_MF(x,4,8,13); % Ασαφές σύνολο εισόδου Α
B=tri_MF(y,45,55,75); % Ασαφές σύνολο εξόδου Β

% Α' = τραπεζοειδής με a=7, b=11, c=13 και d=17.
A_bar=trap_MF(x,7,11,13,17); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------------------------------------------------
%% Μέρος 2
%%------------------------------------------------------------------------

%Γραφική απεικόνιση συνόλων A & A'
figure('Name','Figure 1-Μέρος 2','NumberTitle','off') 
plot(x,A,'b',x,A_bar,'r'); 
title("Blue line->A,     Red line->A'");

%Γραφική απεικόνιση συνόλου Β & B'...
figure('Name','Figure 2-Μέρος 2/3','NumberTitle','off')
axis([min(y) max(y) 0 1.05]);
plot(y,B,'b');
title("Blue line->Β");
hold on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------------------------------------------------
%% Μέρος 3
%%------------------------------------------------------------------------

% Υπολογισμός w (Μέγιστο της τομής Α & Α')
w=max(min(A_bar,A)); 

% Υπολογισμός Ασαφούς συνόλου Β'
B_bar=min(w,B); % Ασαφές σύνολο Β'

%...Γραφική απεικόνιση συνόλου Β & B' (συνέχεια)
area(y,B_bar);
title("Blue line->Β,    Filled Area->B'");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------------------------------------------------
%% Μέρος 4.α
%%------------------------------------------------------------------------

%-------------------------------
% A'=Τριγωνική με a=4, b=8, c=13
%-------------------------------
A_bar=tri_MF(x,4,8,13); % Ασαφές σύνολο Α'
w=max(min(A_bar,A)); % Υπολογισμός w (Μέγιστο της τομής Α & Α')
B_bar=min(w,B); % Ασαφές σύνολο Β'

%Γραφική απεικόνιση συνόλων A & A'
figure('Name','Figure 3-Μέρος 4α','NumberTitle','off')
my_display_2(1,x,A,'A')
my_display_2(2,x,A_bar,"A'")

%Γραφική απεικόνιση συνόλου Β & B'
figure('Name','Figure 4-Μέρος 4α','NumberTitle','off')
plot(y,B,'b');
hold on;
area(y,B_bar);
title("Blue line->Β,    Filled Area->B'");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------------------------------------------------
%% Μέρος 4.β
%%------------------------------------------------------------------------

%---------------------------
% A'=Γκαουσιανή με m=3, σ=2
%---------------------------
A_bar=gauss_MF(x,3,2); % Ασαφές σύνολο Α'
w=max(min(A_bar,A)); % Υπολογισμός w (Μέγιστο της τομής Α & Α')
B_bar=min(w,B); % Ασαφές σύνολο Β'

%Γραφική απεικόνιση συνόλων A & A'
figure('Name','Figure 5-Μέρος 4β','NumberTitle','off')
plot(x,A,'b',x,A_bar,'r'); 
title("Blue line->A,     Red line->A'");

%Γραφική απεικόνιση συνόλου Β & B'
figure('Name','Figure 6-Μέρος 4β','NumberTitle','off')
plot(y,B,'b');
hold on;
area(y,B_bar);
title("Blue line->Β,    Filled Area->B'");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------------------------------------------------
%% Μέρος 5
%%------------------------------------------------------------------------

%-----------------------------------------------------------
% Ασαφείς κανόνες
% R1: IF x is Slow AND y is Very Small THEN z is Very Large
% R2: IF x is Moderate AND y is Large THEN z is Small
%-----------------------------------------------------------

% Βήμα 0.01 για τις μεταβλητές εισόδου
x=linspace(0,100,10001); % Π.Ο. Ταχύτητας
inputSets_x=4; % Ασαφή σύναλα Ταχύτητας (Slow, Moderate, Cruising, Fast)

y=linspace(0,30,3001); % Π.Ο. Επιτάχυνσης
inputSets_y=2; % Ασαφή σύναλα Επιτάχυνσης (Small, Large)

z=linspace(0,3,4); % Π.Ο. Δύναμη που εφαρμόζεται στο γκάζι, βήμα 1
z_optimal=linspace(0,3,301); % Π.Ο. Δύναμη που εφαρμόζεται στο γκάζι, βήμα 0.01
outputSets=2; % Ασαφή σύναλα Δύναμη που εφαρμόζεται στο γκάζι (Small, Large)

overlap=0.4; % Επικάλυψη

%-------------
% Ομοιόμορφη διαμέριση χώρων με τριγωνικές συναρτήσεις

% Διαμέριση και απεικόνιση x σε Slow, Moderate, Cruising, Fast 
[Alpha,Beta,Gamma]=tri_MF_partition(min(x),max(x),inputSets_x,overlap);

for (k=1:inputSets_x)
  mf_x(k,:)=tri_MF(x,Alpha(k),Beta(k),Gamma(k));  
end

% Απεικόνιση ασαφών συνόλων x (τριγωνικές)
figure('Name','Figure 7-Μέρος 5.2','NumberTitle','off')
STR=sprintf('Partition of x with Tri MF');
plot(x,mf_x); hold on;
q=overlap*ones(size(x));
plot(x,q,'k--');
legend('Slow','Moderate','Cruising','Fast', 'Overlap')
axis([min(x) max(x) 0 1.05]);   title(STR);

% Διαμέριση και απεικόνιση y σε Small, Large
[Alpha,Beta,Gamma]=tri_MF_partition(min(y),max(y),inputSets_y,overlap);

for (k=1:inputSets_y)
  mf_y(k,:)=tri_MF(y,Alpha(k),Beta(k),Gamma(k));
end

% Απεικόνιση ασαφών συνόλων y (τριγωνικές)
figure('Name','Figure 8-Μέρος 5.2','NumberTitle','off')
STR=sprintf('Partition of y with Tri MF');
plot(y,mf_y); hold on;
q=overlap*ones(size(y));
plot(y,q,'k--');
legend('Small','Large', 'Overlap')
axis([min(y) max(y) 0 1.05]);   title(STR);

% Διαμέριση z σε Small, Large (βήμα 1)
[Alpha,Beta,Gamma]=tri_MF_partition(min(z),max(z),outputSets,overlap);
for (k=1:outputSets)
  mf_z(k,:)=tri_MF(z,Alpha(k),Beta(k),Gamma(k));
end

% Διαμέριση z σε Small, Large (βήμα 0.01)
[Alpha,Beta,Gamma]=tri_MF_partition(min(z_optimal),max(z_optimal ...
    ),outputSets,overlap);
for (k=1:outputSets)
  mf_z_optimal(k,:)=tri_MF(z_optimal,Alpha(k),Beta(k),Gamma(k));
end

% Απεικόνιση ασαφών συνόλων z (τριγωνικές)
figure('Name','Figure 9-Μέρος 5.2','NumberTitle','off')
my_display_2(1,z,mf_z,'Partition of z with Tri MF(step 1)') % για βήμα 1
hold on;
q=overlap*ones(size(z));
plot(z,q,'k--');
legend('Small','Large', 'Overlap')
axis([min(z) max(z) 0 1.05]);

my_display_2(2,z_optimal,mf_z_optimal,'Partition of z with Tri MF(step 0.01)') % για βήμα 0.01
hold on;
q=overlap*ones(size(z));
plot(z,q,'k--');
legend('Small','Large', 'Overlap')
axis([min(z) max(z) 0 1.05]); 


% Υλοποίηση και Απεικόνιση μετασχηματισμού "Very Small" του y (τριγωνικές)
y_VerySmall=mf_y(1,:).^2;
figure('Name','Figure 10-Μέρος 5.3','NumberTitle','off')
STR=sprintf('y=Very Small');
plot(y,y_VerySmall);
axis([min(y) max(y) 0 1.05]);   title(STR);

% Υλοποίηση μετασχηματισμού "Very Large" του z(βήμα 1) (τριγωνικές)
z_VeryLarge=mf_z(2,:).^2;
% Υλοποίηση μετασχηματισμού "Very Large" του z(βήμα 0.01) (τριγωνικές)
z_VeryLarge_optimal=mf_z_optimal(2,:).^2;

% Απεικόνιση μετασχηματισμού "Very Large" z με βήμα 1 & 0.01 (τριγωνικές)
figure('Name','Figure 11-Μέρος 5.3','NumberTitle','off')
my_display_2(1,z,z_VeryLarge,'z=Very Large(step 1)')
my_display_2(2,z_optimal,z_VeryLarge_optimal,'z=Very Large(step 0.01)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------------------------------------------------
%% Μέρος 6
%%------------------------------------------------------------------------

%-------------
% Ομοιόμορφη διαμέριση χώρων με Γκαουσιανές συναρτήσεις συμμετοχής

% Διαμέριση x σε Slow, Moderate, Cruising, Fast 
A=min(x); B=max(x);
temp2=(B-A)/(inputSets_x-1);
meso=linspace(A,B,inputSets_x);
sigma=(meso(2)-meso(1))/(2*sqrt((-2)*log(overlap)));

for (k=1:inputSets_x)    
    mf_x_gauss(k,:)=exp(-0.5*(x-meso(k)).^2/(sigma^2));
end

% Απεικόνιση ασαφών συνόλων x (Γκαουσιανές)
figure('Name','Figure 12-Μέρος 6','NumberTitle','off')
STR=sprintf('Partition of x with Gaussian MF');
plot(x,mf_x_gauss); hold on;
q=overlap*ones(size(x));
plot(x,q,'k--');
legend('Slow','Moderate','Cruising','Fast', 'Overlap')
axis([min(x) max(x) 0 1.05]);   title(STR);

% Διαμέριση y σε Small, Large
A=min(y); B=max(y);
temp2=(B-A)/(inputSets_y-1);
meso=linspace(A,B,inputSets_y);
sigma=(meso(2)-meso(1))/(2*sqrt((-2)*log(overlap)));

for (k=1:inputSets_y)
  mf_y_gauss(k,:)=exp(-0.5*(y-meso(k)).^2/(sigma^2));
end

% Απεικόνιση ασαφών συνόλων y (Γκαουσιανές)
figure('Name','Figure 13-Μέρος 6','NumberTitle','off')
STR=sprintf('Partition of y with Gaussian MF');
plot(y,mf_y_gauss); hold on;
q=overlap*ones(size(y));
plot(y,q,'k--');
legend('Small','Large', 'Overlap')
axis([min(y) max(y) 0 1.05]);   title(STR);

% Διαμέριση z σε Small, Large (βήμα 1)
A=min(z); B=max(z);
temp2=(B-A)/(outputSets-1);
meso=linspace(A,B,outputSets);
sigma=(meso(2)-meso(1))/(2*sqrt((-2)*log(overlap)));

for (k=1:outputSets)    
    mf_z_gauss(k,:)=exp(-0.5*(z-meso(k)).^2/(sigma^2));  
end

% Διαμέριση z σε Small, Large (βήμα 0.01)
A=min(z_optimal); B=max(z_optimal);
temp2=(B-A)/(outputSets-1);
meso=linspace(A,B,outputSets);
sigma=(meso(2)-meso(1))/(2*sqrt((-2)*log(overlap)));

for (k=1:outputSets)
  mf_z_gauss_optimal(k,:)=exp(-0.5*(z_optimal-meso(k)).^2/(sigma^2));  
end

% Απεικόνιση ασαφών συνόλων z (Γκοαουσιανές)
figure('Name','Figure 14-Μέρος 6','NumberTitle','off')
my_display_2(1,z,mf_z_gauss,'Partition of z with Gauss MF(step 1)') % για βήμα 1
hold on;
q=overlap*ones(size(z));
plot(z,q,'k--');
legend('Small','Large', 'Overlap')
axis([min(z) max(z) 0 1.05]);

my_display_2(2,z_optimal,mf_z_gauss_optimal,'Partition of z with Gauss MF(step 0.01)') % για βήμα 0.01
hold on;
q=overlap*ones(size(z_optimal));
plot(z_optimal,q,'k--');
legend('Small','Large', 'Overlap')
axis([min(z) max(z) 0 1.05]); 

% Υλοποίηση μετασχηματισμού "Very Small" του y
y_VerySmall_gauss=mf_y_gauss(1,:).^2;

% Υλοποίηση μετασχηματισμού "Very Large" του z(βήμα 1)
z_VeryLarge_gauss=mf_z_gauss(2,:).^2;

% Υλοποίηση μετασχηματισμού "Very Large" του z(βήμα 0.1)
z_VeryLarge_gauss_optimal=mf_z_gauss_optimal(2,:).^2;

% Απεικόνιση μετασχηματισμού "Very Small" του y
figure('Name','Figure 15-Μέρος 6','NumberTitle','off')
STR=sprintf('y=Very Small with Gaussian MF');
plot(y,y_VerySmall_gauss);
axis([min(y) max(y) 0 1.05]);   title(STR);

% Απεικόνιση μετασχηματισμού "Very Large" του z βήμα 1 & 0.01
figure('Name','Figure 16-Μέρος 6','NumberTitle','off')
my_display_2(1,z,z_VeryLarge_gauss,'z=Very Large with Gaussian MF(step 1)') % για βήμα 1
my_display_2(2,z_optimal,z_VeryLarge_gauss_optimal,'z=Very Large with Gaussian MF(step 0.01)') % για βήμα 0.01


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------------------------------------------------
%% Μέρος 7
%%------------------------------------------------------------------------

fprintf("Μέρος 7\nΑσαφοποίηση ζεύγους (35,12) με Τριγωνικές ΣΣ\n\n")

%%------------------------------------------------------------------------
%% Μέρος 7.1
%%------------------------------------------------------------------------

% Χρήση Τριγωνικών ΣΣ και βήμα διακριτοποίησης 0.01 από Μέρος 5
fprintf('Μέρος 7.1\nΤριγωνικές ΣΣ(βήμα 0.01)\n\n')

% Ζεύγος εισόδων (35,12)
% Ασαφοποίηση με τριγωνικές ΣΣ
A_bar=tri_MF(x,15,35,55); % Ασαφές σύνολο Α'
B_bar=tri_MF(y,9,12,15); % Ασαφές σύνολο Β'

% Υπολογισμός w1 και w2
disp('Rule 1')
w1=w_2_inputs(mf_x(1,:),A_bar,y_VerySmall,B_bar);
disp('Rule 2')
w2=w_2_inputs(mf_x(2,:),A_bar,mf_y(2,:),B_bar);

% Υπολογισμός ασαφών συνόλων εξόδου C1' και C2' με βήμα 0.01
C1_bar= min(w1,z_VeryLarge_optimal);
C2_bar= min(w2,mf_z_optimal(1,:));

% Υπολογισμός Ασαφούς συνόλου εξόδου της βάσης κανόνων με βήμα 0.01
C=max(C1_bar,C2_bar);

% Απεικόνιση ασαφών συνόλών εξόδου C1' & C2', και
% συνολικού ασαφούς συνόλου εξόδου της βάσης κανόνων C με βήμα 0.01
figure('Name','Figure 17-Μέρος 7.1','NumberTitle','off')
my_display_2(1,z_optimal,C1_bar,"Fuzzy output set C1' & C2'(step 0.01)(Tri MF)")
hold on
plot(z_optimal,C2_bar);
legend("C1'","C2'")
my_display_2(2,z_optimal,C,"Fuzzy output set C'(step 0.01)(Tri MF)")

% Σαφής έξοδος του συστήματος με βήμα διακριτοποίησης 0.01
disp("Σαφής έξοδος του συστήματος με βήμα διακριτοποίησης 0.01")
output=sum(z_optimal.*C)/sum(C)

%--------------------

% Χρήση Γκαουσιανών ΣΣ και βήμα διακριτοποίησης 0.01 από Μέρος 6

fprintf('\nΜέρος 7.1\nΓκαουσιανές ΣΣ(βήμα 0.01)\n\n')

% Υπολογισμός w1 και w2
disp('Rule 1')
w1=w_2_inputs(mf_x_gauss(1,:),A_bar,y_VerySmall_gauss,B_bar);
disp('Rule 2')
w2=w_2_inputs(mf_x_gauss(2,:),A_bar,mf_y_gauss(2,:),B_bar);

% Υπολογισμός ασαφών συνόλων εξόδου C1' και C2' με βήμα 0.01
C1_bar= min(w1,z_VeryLarge_gauss_optimal);
C2_bar= min(w2,mf_z_gauss_optimal(1,:));

% Υπολογισμός Ασαφούς συνόλου εξόδου της βάσης κανόνων με βήμα 0.01
C=max(C1_bar,C2_bar);

% Απεικόνιση ασαφών συνόλών εξόδου C1' & C2', και
% συνολικού ασαφούς συνόλου εξόδου της βάσης κανόνων C με βήμα 0.01
figure('Name','Figure 18-Μέρος 7.1','NumberTitle','off')
my_display_2(1,z_optimal,C1_bar,"Fuzzy output set C1' & C2'(step 0.01)(Gauss MF)")
hold on
plot(z_optimal,C2_bar);
legend("C1'","C2'")
my_display_2(2,z_optimal,C,"Fuzzy output set C'(step 0.01)(Gauss MF)")

% Σαφής έξοδος του συστήματος με βήμα διακριτοποίησης 0.01
disp("Σαφής έξοδος του συστήματος με βήμα διακριτοποίησης 0.01")
output=sum(z_optimal.*C)/sum(C)


%%------------------------------------------------------------------------
%% Μέρος 7.2
%%------------------------------------------------------------------------

% Χρήση Τριγωνικών ΣΣ και βήμα διακριτοποίησης 1 από Μέρος 5

fprintf('------------------------------------------------\nΜέρος 7.2\nΤριγωνικές ΣΣ(βήμα 1)\n');

% Υπολογισμός w1 και w2
disp('Rule 1')
w1=w_2_inputs(mf_x(1,:),A_bar,y_VerySmall,B_bar);
disp('Rule 2')
w2=w_2_inputs(mf_x(2,:),A_bar,mf_y(2,:),B_bar);

% Υπολογισμός ασαφών συνόλων εξόδου C1' και C2' με βήμα 1
C1_bar= min(w1,z_VeryLarge);
C2_bar= min(w2,mf_z(1,:));

% Υπολογισμός Ασαφούς συνόλου εξόδου της βάσης κανόνων με βήμα 1
C=max(C1_bar,C2_bar);


% Απεικόνιση ασαφών συνόλών εξόδου C1' & C2', και
% συνολικού ασαφούς συνόλου εξόδου της βάσης κανόνων C με βήμα 1
figure('Name','Figure 19-Μέρος 7.2','NumberTitle','off')
my_display_2(1,z,C1_bar,"Fuzzy output set C1' & C2'(step 1)(Tri MF)")
hold on
plot(z,C2_bar);
legend("C1'","C2'")
my_display_2(2,z,C,"Fuzzy output set C'(step 1)(Tri MF)")

% Σαφής έξοδος του συστήματος με βήμα διακριτοποίησης 1
disp("Σαφής έξοδος του συστήματος με βήμα διακριτοποίησης 1")
output=sum(z.*C)/sum(C)

%---------------

% Χρήση Γκαουσιανών ΣΣ και βήμα διακριτοποίησης 1 από Μέρος 6

fprintf('\nΓκαουσιανές ΣΣ\n');

% Υπολογισμός w1 και w2
disp('Rule 1')
w1=w_2_inputs(mf_x_gauss(1,:),A_bar,y_VerySmall_gauss,B_bar);
disp('Rule 2')
w2=w_2_inputs(mf_x_gauss(2,:),A_bar,mf_y_gauss(2,:),B_bar);

% Υπολογισμός ασαφών συνόλων εξόδου C1' και C2' με βήμα 1
C1_bar= min(w1,z_VeryLarge_gauss);
C2_bar= min(w2,mf_z_gauss(1,:));

% Υπολογισμός Ασαφούς συνόλου εξόδου της βάσης κανόνων με βήμα 1
C=max(C1_bar,C2_bar);

% Απεικόνιση ασαφών συνόλών εξόδου C1' & C2', και
% συνολικού ασαφούς συνόλου εξόδου της βάσης κανόνων C με βήμα 0.01
figure('Name','Figure 20-Μέρος 7.2','NumberTitle','off')
my_display_2(1,z,C1_bar,"Fuzzy output set C1' & C2'(step 1)(Gauss MF)")
hold on
plot(z,C2_bar);
legend("C1'","C2'")
my_display_2(2,z,C,"Fuzzy output set C'(step 1)(Gauss MF)")

% Σαφής έξοδος του συστήματος με βήμα διακριτοποίησης 1
disp("Σαφής έξοδος του συστήματος με βήμα διακριτοποίησης 1")
output=sum(z.*C)/sum(C)

%%------------------------------------------------------------------------
%% Μέρος 7.3
%%------------------------------------------------------------------------

fprintf('------------------------------------------------\nΜέρος 7.3\nΤριγωνικές ΣΣ\n');

%Δημιουργία ασαφών singleton
x0=35;
y0=12;
A_bar=zeros(size(x));
A_bar(find(x==x0))=1;
B_bar=zeros(size(y));
B_bar(find(y==y0))=1;

% Υπολογισμός w1 και w2
disp('Rule 1')
w1=w_2_inputs(mf_x(1,:),A_bar,y_VerySmall,B_bar);
disp('Rule 2')
w2=w_2_inputs(mf_x(2,:),A_bar,mf_y(2,:),B_bar);

% Υπολογισμός ασαφών συνόλων εξόδου C1' και C2' με βήμα 0.01
C1_bar= min(w1,z_VeryLarge_optimal);
C2_bar= min(w2,mf_z_optimal(1,:));

% Υπολογισμός Ασαφούς συνόλου εξόδου της βάσης κανόνων με βήμα 0.01
C=max(C1_bar,C2_bar);

% Απεικόνιση ασαφών συνόλών εξόδου C1' & C2', και
% συνολικού ασαφούς συνόλου εξόδου της βάσης κανόνων C με βήμα 0.01
figure('Name','Figure 21-Μέρος 7.3','NumberTitle','off')
my_display_2(1,z_optimal,C1_bar,"Fuzzy output set C1' & C2'(step 0.01)(Tri MF)")
hold on
plot(z_optimal,C2_bar);
legend("C1'","C2'")
my_display_2(2,z_optimal,C,"Fuzzy output set C'(step 0.01)(Tri MF)")

% Σαφής έξοδος του συστήματος με βήμα διακριτοποίησης 0.01
disp("Σαφής έξοδος του συστήματος με βήμα διακριτοποίησης 0.01")
output=sum(z_optimal.*C)/sum(C)


%------------------------

fprintf('\nΓκαουσιανές ΣΣ\n');

% Υπολογισμός w1 και w2
disp('Rule 1')
w1=w_2_inputs(mf_x_gauss(1,:),A_bar,y_VerySmall_gauss,B_bar);
disp('Rule 2')
w2=w_2_inputs(mf_x_gauss(2,:),A_bar,mf_y_gauss(2,:),B_bar);

% Υπολογισμός ασαφών συνόλων εξόδου C1' και C2' με βήμα 0.01
C1_bar= min(w1,z_VeryLarge_gauss_optimal);
C2_bar= min(w2,mf_z_optimal(1,:));

% Υπολογισμός Ασαφούς συνόλου εξόδου της βάσης κανόνων με βήμα 0.01
C=max(C1_bar,C2_bar);

% Απεικόνιση ασαφών συνόλών εξόδου C1' & C2', και
% συνολικού ασαφούς συνόλου εξόδου της βάσης κανόνων C με βήμα 0.01
figure('Name','Figure 22-Μέρος 7.3','NumberTitle','off')
my_display_2(1,z_optimal,C1_bar,"Fuzzy output set C1' & C2'(step 0.01)(Gauss MF)")
hold on
plot(z_optimal,C2_bar);
legend("C1'","C2'")
my_display_2(2,z_optimal,C,"Fuzzy output set C'(step 0.01)(Gauss MF)")

% Σαφής έξοδος του συστήματος με βήμα διακριτοποίησης 0.01
disp("Σαφής έξοδος του συστήματος με βήμα διακριτοποίησης 0.01")
output=sum(z_optimal.*C)/sum(C)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------------------------------------------------
%% Μέρος 8
%%------------------------------------------------------------------------

fprintf("Μέρος 8\nΑσαφοποίηση ζεύγους (35,12) με Γκαουσιανές ΣΣ\n\n")

%Ζεύγος εισόδων (35,12)
%Ασαφοποίηση με Γκαουσιανές ΣΣ
A_bar=gauss_MF(x,35,8); % Ασαφές σύνολο Α'
B_bar=gauss_MF(y,12,5); % Ασαφές σύνολο Β'

%%------------------------------------------------------------------------
%% Μέρος 8.1
%%------------------------------------------------------------------------

% Χρήση Τριγωνικών ΣΣ και βήμα διακριτοποίησης 0.01 από Μέρος 5
fprintf('Μέρος 8.1\nΤριγωνικές ΣΣ(βήμα 0.01)\n\n')

% Υπολογισμός w1 και w2
disp('Rule 1')
w1=w_2_inputs(mf_x(1,:),A_bar,y_VerySmall,B_bar);
disp('Rule 2')
w2=w_2_inputs(mf_x(2,:),A_bar,mf_y(2,:),B_bar);

% Υπολογισμός ασαφών συνόλων εξόδου C1' και C2' με βήμα 0.01
C1_bar= min(w1,z_VeryLarge_optimal);
C2_bar= min(w2,mf_z_optimal(1,:));

% Υπολογισμός Ασαφούς συνόλου εξόδου της βάσης κανόνων με βήμα 0.01
C=max(C1_bar,C2_bar);

% Απεικόνιση ασαφών συνόλών εξόδου C1' & C2', και
% συνολικού ασαφούς συνόλου εξόδου της βάσης κανόνων C με βήμα 0.01
figure('Name','Figure 23-Μέρος 8.1','NumberTitle','off')
my_display_2(1,z_optimal,C1_bar,"Fuzzy output set C1' & C2'(step 0.01)(Tri MF)")
hold on
plot(z_optimal,C2_bar);
legend("C1'","C2'")
my_display_2(2,z_optimal,C,"Fuzzy output set C'(step 0.01)(Tri MF)")

% Σαφής έξοδος του συστήματος με βήμα διακριτοποίησης 0.01
disp("Σαφής έξοδος του συστήματος με βήμα διακριτοποίησης 0.01")
output=sum(z_optimal.*C)/sum(C)

%--------------------

% Χρήση Γκαουσιανών ΣΣ και βήμα διακριτοποίησης 0.01 από Μέρος 6
fprintf('\nΜέρος 8.1\nΓκαουσιανές ΣΣ(βήμα 0.01)\n\n')

% Υπολογισμός w1 και w2
disp('Rule 1')
w1=w_2_inputs(mf_x_gauss(1,:),A_bar,y_VerySmall_gauss,B_bar);
disp('Rule 2')
w2=w_2_inputs(mf_x_gauss(2,:),A_bar,mf_y_gauss(2,:),B_bar);

% Υπολογισμός ασαφών συνόλων εξόδου C1' και C2' με βήμα 0.01
C1_bar= min(w1,z_VeryLarge_gauss_optimal);
C2_bar= min(w2,mf_z_gauss_optimal(1,:));

% Υπολογισμός Ασαφούς συνόλου εξόδου της βάσης κανόνων με βήμα 0.01
C=max(C1_bar,C2_bar);

% Απεικόνιση ασαφών συνόλών εξόδου C1' & C2', και
% συνολικού ασαφούς συνόλου εξόδου της βάσης κανόνων C με βήμα 0.01
figure('Name','Figure 24-Μέρος 8.1','NumberTitle','off')
my_display_2(1,z_optimal,C1_bar,"Fuzzy output set C1' & C2'(step 0.01)(Gauss MF)")
hold on
plot(z_optimal,C2_bar);
legend("C1'","C2'")
my_display_2(2,z_optimal,C,"Fuzzy output set C'(step 0.01)(Gauss MF)")

% Σαφής έξοδος του συστήματος με βήμα διακριτοποίησης 0.01
disp("Σαφής έξοδος του συστήματος με βήμα διακριτοποίησης 0.01")
output=sum(z_optimal.*C)/sum(C)


%%------------------------------------------------------------------------
%% Μέρος 8.2
%%------------------------------------------------------------------------

% Χρήση Τριγωνικών ΣΣ και βήμα διακριτοποίησης 1 από Μέρος 5

fprintf('------------------------------------------------\nΜέρος 8.2\nΤριγωνικές ΣΣ(βήμα 1)\n');

% Υπολογισμός w1 και w2
disp('Rule 1')
w1=w_2_inputs(mf_x(1,:),A_bar,y_VerySmall,B_bar);
disp('Rule 2')
w2=w_2_inputs(mf_x(2,:),A_bar,mf_y(2,:),B_bar);

% Υπολογισμός ασαφών συνόλων εξόδου C1' και C2' με βήμα 1
C1_bar= min(w1,z_VeryLarge);
C2_bar= min(w2,mf_z(1,:));

% Υπολογισμός Ασαφούς συνόλου εξόδου της βάσης κανόνων με βήμα 1
C=max(C1_bar,C2_bar);


% Απεικόνιση ασαφών συνόλών εξόδου C1' & C2', και
% συνολικού ασαφούς συνόλου εξόδου της βάσης κανόνων C με βήμα 1
figure('Name','Figure 25-Μέρος 8.2','NumberTitle','off')
my_display_2(1,z,C1_bar,"Fuzzy output set C1' & C2'(step 1)(Tri MF)")
hold on
plot(z,C2_bar);
legend("C1'","C2'")
my_display_2(2,z,C,"Fuzzy output set C'(step 1)(Tri MF)")

% Σαφής έξοδος του συστήματος με βήμα διακριτοποίησης 1
disp("Σαφής έξοδος του συστήματος με βήμα διακριτοποίησης 1")
output=sum(z.*C)/sum(C)

%---------------

% Χρήση Γκαουσιανών ΣΣ και βήμα διακριτοποίησης 1 από Μέρος 6

fprintf('\nΓκαουσιανές ΣΣ\n');

% Υπολογισμός w1 και w2
disp('Rule 1')
w1=w_2_inputs(mf_x_gauss(1,:),A_bar,y_VerySmall_gauss,B_bar);
disp('Rule 2')
w2=w_2_inputs(mf_x_gauss(2,:),A_bar,mf_y_gauss(2,:),B_bar);

% Υπολογισμός ασαφών συνόλων εξόδου C1' και C2' με βήμα 1
C1_bar= min(w1,z_VeryLarge_gauss);
C2_bar= min(w2,mf_z_gauss(1,:));

% Υπολογισμός Ασαφούς συνόλου εξόδου της βάσης κανόνων με βήμα 1
C=max(C1_bar,C2_bar);

% Απεικόνιση ασαφών συνόλών εξόδου C1' & C2', και
% συνολικού ασαφούς συνόλου εξόδου της βάσης κανόνων C με βήμα 1
figure('Name','Figure 26-Μέρος 8.2','NumberTitle','off')
my_display_2(1,z,C1_bar,"Fuzzy output set C1' & C2'(step 1)(Gauss MF)")
hold on
plot(z,C2_bar);
legend("C1'","C2'")
my_display_2(2,z,C,"Fuzzy output set C'(step 1)(Gauss MF)")

% Σαφής έξοδος του συστήματος με βήμα διακριτοποίησης 1
disp("Σαφής έξοδος του συστήματος με βήμα διακριτοποίησης 1")
output=sum(z.*C)/sum(C)