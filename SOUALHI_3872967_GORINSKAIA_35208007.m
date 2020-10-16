clc 
clf
clear all 
close all 
%###########################################################################################
%
% SOUALHI TAKIEDDINE 3872967
% GORINSKAIA KATERINA 30 
%
% Ce Programme une fois executé, va calculer et afficher pour les
% differntes valeur de N : la solution analytique, la solution numerique et
% la l'erreur globale, et traçe aprés l'erreur globale on fonction de pas h
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



n=[5 10 20 40 80]; %differentes valeurs de N 
pas=[1/5 1/10 1/20 1/40 1/80]; %differentes valeurs de pas 
Err_globale=zeros(1,length(n)); %Initialisation de l'erreur globale 


for m=1:5;  

    
    
N=n(m); 
N_dim = N-1;
sprintf('la valeur de N =%d',N);

% contructuion de la matrice A de dim N ;

%contruction de la la diagonale 
A = 4*eye(N_dim);
A(2:N_dim,1:N_dim-1) = A(2:N_dim,1:N_dim-1) - eye(N_dim-1);
A(1:N_dim-1,2:N_dim) = A(1:N_dim-1,2:N_dim) - eye(N_dim-1);  
A=-A;


%mettre en forme la matrice A dans la matrice globale 
U = zeros((N-1)*N);

for  i=1:N
    M = U(i*N_dim-(N_dim-1): i*N_dim , i*N_dim-(N_dim-1): i*N_dim  ) + A;
    U(i*N_dim-(N_dim-1): i*N_dim , i*N_dim-(N_dim-1): i*N_dim  ) = M;
end

M = U(N*N_dim-(N_dim-1): N*N_dim , N*N_dim-(N_dim-1): N*N_dim  ) - (0.5*A);
U(N*N_dim-(N_dim-1): N*N_dim , N*N_dim-(N_dim-1): N*N_dim  ) = M;

%construction des sous-matrices identité 
I3 = eye(N_dim);

for  i=1:N_dim
    
   M = U(((i+1)*N_dim)-(N_dim-1): (i+1)*N_dim , (i)*N_dim-(N_dim-1): i*N_dim) + I3;
    U(((i+1)*N_dim)-(N_dim-1): (i+1)*N_dim , (i)*N_dim-(N_dim-1): i*N_dim) = M;
    
   M = U(((i)*N_dim)-(N_dim-1): (i)*N_dim , (i+1)*N_dim-(N_dim-1): (i+1)*N_dim) + I3;
    U(((i)*N_dim)-(N_dim-1): (i)*N_dim , (i+1)*N_dim-(N_dim-1): (i+1)*N_dim) = M;
    
end


% Resolutuion 
% matrices Tg Th Tb 

 Ub=0;
 Ug=0;
 Uh=1;
    
B=zeros(N*(N-1),1);
Tb=zeros(N*(N-1),1);
Th=zeros(N*(N-1),1);
Tg=zeros(N*(N-1),1);


for  i=1:N

    Th(i*N_dim) = Tg(i*N_dim) + Uh ; 
    
end 

% Construction du vecteur B
B = -Th -Tg -Tb;
B(end)=-0.5;

% matrice T Inconnue 
T = U\B ;

             
pas_h = 1/N;   % Pas d'espace en x et y
x_discretise = ([1:N+1]-1)*pas_h;   % Liste des valeurs de x, de 0 a 1 inclus
y_discretise = ([1:N+1]-1)*pas_h;   % Liste des valeurs de y, de 0 a 1 inclus

% Construction d'une grille de maillage 
% X_grid : contient les abcisses de tous les points de la grille
% Y_grid : contient les ordonnï¿½es de tous les points de la grille
[X_grid,Y_grid] = meshgrid(x_discretise,y_discretise);

%reformulation de la solution calculé numeriquement
T_grid = reshape(T,N-1,N);
V1=[0 ones(1,N)];
V2=zeros(N_dim,1);
V3=zeros(1,N+1);

%Solution Numerique Finale 
U_numerique=[V3; V2 T_grid ;V1];


%Calcul de la Solution Analytique en utilisant la formule fournie 
Dim=size(U_numerique);
U_analytique=zeros(Dim(1),Dim(2));
for i=1:Dim(1)
    for j=1:Dim(2) 
         for p=1:2:400
             U_analytique(i,j)=U_analytique(i,j)+(sin((pi/2)*p*X_grid(i,j))*sinh(p*(pi/2)*Y_grid(i,j)))/(p*sinh(p*(pi/2))); %implementation de la formule donnée 
         end 
    U_analytique(i,j)=U_analytique(i,j)*(4/pi); 
    end
end
U_analytique(Dim(1),:)=V1;

%Calcul de l'erreur locale
Err_locale=zeros(Dim(1),Dim(2)); %initialisation du vecteur a zero 
for i=1:Dim(1)
    for j=1:Dim(2) 
    Err_locale(i,j)= abs((U_numerique(i,j)-U_analytique(i,j)));
    end
end 

% calcul de l'erreur globale 
Err_globale(m)=sum(sum(power(Err_locale,2)));
Err_globale(m)=sqrt(Err_globale(m)/(N*(N-1)));


%Visualisation de la solution numerique, la solution analytique et l'erreur
%locale pour un N donné

figure(m)
hold on

% Visualisation solution numerique de u(x,y)pour un N donné
subplot(3,1,1)
contourf(X_grid,Y_grid,U_numerique)
xlabel('Abcisse x'), ylabel('Ordonnée y');
title(['Solution Numerique N=',num2str(n(m))]);

% Visualisation solution Analytique de u(x,y)pour un N donné
subplot(3,1,2)
contourf(X_grid,Y_grid,U_analytique)
xlabel('Abcisse x'), ylabel('Ordonnée y');
title(['Solution Analytique N=',num2str(n(m))]);

% Visualisation de l'erreur locale pour un N donné
subplot(3,1,3)
contourf(X_grid,Y_grid,Err_locale)
xlabel('Abcisse x'), ylabel('Ordonnée y ');
title(['Erreur Locale N=',num2str(n(m))]);
hold off 

end 


figure(6) 
% Visualisation de l'erreur globale en fonction du pas pour un les
% differentes valeur de N 
plot(pas,Err_globale)
xlabel('le Pas'), ylabel('Erreur Globale ')
title(['tracé de lerreur globale en fonction de pas ']);
figure(7) 
plot(n,Err_globale)
xlabel('N'), ylabel('Erreur Globale ')
title(['tracé de lerreur globale en fonction de nombre des noeuds ']);