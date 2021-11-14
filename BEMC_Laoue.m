% BE : Approche multi-échelle progressive de la rupture
% Auteur : Nathan Laoué
% v1.0 : 9/11/2020
%----------------------------------------------------------%
clear
clc
close all
%----------------------------------------------------------%

%%
%----------------------Données----------------------%
disp("BE Nathan Laoué v1.0 ")
disp("---------------------------------------------------------------------------------------------------------------------------------------------------")
fprintf("Consignes d'entrée des valeurs :\n-aucun espace\n-pas de virgules ou autres caractère spécial seulement des points pour les chiffres a virgule\n-pas de lettres, il faut seulement entrer les chiffres sous ce format Ex: 255000 ou 25.12\n"); 
disp("---------------------------------------------------------------------------------------------------------------------------------------------------")
disp('Propriétées du stratifié :')

%Toutes les données sont entrées, si on veut permettre l'entrée de données
%différentes on peut les modifier directement sur le code ou alors supprimmer les données prérentré dans le
%code et décommenter tous les input de la partie "données" afin de pouvoir entrer les
%données via la command window. Tout a été réaliser en MPa.

E11= 140000; disp(['E11[MPa] : ',num2str(E11)]) %input('E11[MPa] : '); %valeur de E11
E22= 10000; disp(['E22[MPa] : ',num2str(E22)])%input('E22 [MPa] : '); %valeur de E22
G12= 6000; disp(['G12[MPa] : ',num2str(G12)])%input('G12 [MPa] : '); %valeur de G12
nu12 = 0.3; disp(['nu12 : ',num2str(nu12)])%input('nu12 : '); % valeur de nu12
ep=0.256; disp(['Épaisseur [mm] : ',num2str(ep)])%input('Épaisseur [mm] : '); %valeur de l'épaisseur

fprintf('\n Contraintes à  rupture du pli : \n')
Xt=1990; disp(['Xt [Mpa]=',num2str(Xt)])%input('Xt [Mpa] : ') %valeur de Xt
Xc=-1500; disp(['Xc [Mpa]=',num2str(Xc)])%input('Xc [Mpa] : ') %valeur de Xc
Yt=38; disp(['Yt [Mpa]=',num2str(Yt)])%input('Yt [Mpa] : ') %valeur de Yt
Yc=-150; disp(['Yc [Mpa]=',num2str(Yc)])%input('Yc [Mpa] : ') %valeur de Yc
Sc=70; disp(['Sc [Mpa]=',num2str(Sc)])%input('Sc [Mpa] : ') %valeur de Sc

fprintf('\n Chargement appliqué au stratifié :\n')
N1=1598.796; disp(['Nx [N.mm^-1 : ',num2str(N1)])%input('N1 [N.mm^-1 : ') %valeur de N1
N2=0; disp(['Ny [N.mm^-1 : ',num2str(N2)])%input('N2 [N.mm^-1 : ') %valeur de N2
N3=0; disp(['Nz [N.mm^-1 : ',num2str(N3)])%input('N3 [N.mm^-1 : ') %valeur de N3
M1=0; disp(['Mx [N] : ',num2str(M1)])%input('M1 [N] : ') %valeur de M1
M2=0; disp(['My [N] : ',num2str(M2)])%input('M2 [N] : ') %valeur de M2
M3=0; disp(['Mz [N] : ',num2str(M3)])%input('M3 [N] : ') %valeur de M3
N=[N1 N2 N3]';
M=[M1 M2 M3]';

% Exemple :
% empilement=[0,45,135,90,90,135,45,0];
% taille=length(empilement);

Symetrie=input("La séquence d'empilement est symétrique ? (oui=1,non=2) \n"); %entrée de la séquence

if Symetrie == 1 %si la séquence est symétrique
    taille=input('Taille de la séquence : ');
    empilement=zeros(1,taille);
    for i=1:taille
        empilement(i)=input(['Entrer le ',num2str(i),' élément de la séquence (en Â°) : ']);
    end
    
    
    %augmentation de la taille du a la symétrie%
    n=taille;
    for i=(taille+1):taille*2
        empilement(i)=empilement(n);
        n=n-1;
    end
    taille=length(empilement);
    
else %si la séquence n'est pas symétrique
    taille=input('Taille de la séquence : ');
    empilement=zeros(1,taille);
    for i=1:taille
        empilement(i)=input(['Entrer le ',num2str(i),' élément de la séquence (en Â°) : ']);
    end
end

%%
%-----------------------calculs-----------------------%


%--calcul matrice de souplesse repère matériaux--%
S11 = 1/E11; S12 = -nu12/E11; S22=1/E22; S66=1/G12;
Smat = [S11 S12 0; S12 S22 0; 0 0 S66]; %en Mpa^-1

%--calcul matrice de rigidité réduite repère matériaux--%
Qmat = inv(Smat); %en Mpa

%--calcul des matrices de rigidité pour chaque pli--%
c=zeros(1,taille);
s=zeros(1,taille);
Tsig=zeros(3,3,taille);
Tepsi=zeros(3,3,taille);
Qglo=zeros(3,3,taille);

for i=1:taille
    c(i)=cosd(empilement(i));
    s(i)=sind(empilement(i));
    
    Tsig(:,:,i)=[c(i)^2 s(i)^2 -2*s(i)*c(i); s(i)^2 c(i)^2 2*s(i)*c(i); s(i)*c(i) -s(i)*c(i) c(i)^2-s(i)^2];
    Tepsi(:,:,i)=[c(i)^2 s(i)^2 -s(i)*c(i); s(i)^2 c(i)^2 s(i)*c(i); 2*s(i)*c(i) -2*s(i)*c(i) c(i)^2-s(i)^2];
    
    Qglo(:,:,i)=Tsig(:,:,i)*(Qmat/Tepsi(:,:,i));
end

%--calcul des matrices de souplesse pour chaque pli--%
Sglo=zeros(3,3,taille);

for i=1:taille
Sglo(:,:,i)=inv(Qglo(:,:,i));
end

%--calcul de A--%
A=zeros(3);
for i=1:taille
    A=A+Qglo(:,:,i); %en Mpa.mm
end
A=ep*A;

%--calcul de B--%
B=zeros(3);
htot=ep*taille;
hpli=-htot/2:ep:htot/2;
for i=1:(taille)
    coefh=((hpli(i+1))^2-(hpli(i))^2);
    B=B+coefh*Qglo(:,:,i); %en Mpa.mm^2
end
B=0.5*B;

%--calcul de D--%
D=zeros(3);
for i=1:taille
    coefh=((hpli(i+1))^3-(hpli(i))^3);
    D=D+coefh*Qglo(:,:,i); %en Mpa.mm^3
end
D=(1/3)*D;

%--calcul des déformations de membrane et courbures du stratifié--%
Matdefcont=[A B; B D]\[N;M]; %matrice déformation contraintes
Epsim=Matdefcont(1:3,:); %déformations de mebranes
k=Matdefcont(4:6,:); %courbures du stratifié

%--calcul des déformations et des contraintes dans chaque pli--%
Zk=zeros(1,taille);
for i=1:taille
    Zk(i)=1/2*(hpli(i)+hpli(i+1)); %calcul de Zk
end

epsipliglo=zeros(3,1,taille);
for i=1:taille
    epsipliglo(:,:,i)=Epsim+Zk(i)*k; %déformations de chaque plis dans le repère global
end

epsiplimat=zeros(3,1,taille);
for i=1:taille
    epsiplimat(:,:,i)=Tepsi(:,:,i)\epsipliglo(:,:,i); %déformations de chaque plis dans le repère local
end

sigmaplimat=zeros(3,1,taille);
for i=1:taille
    sigmaplimat(:,:,i)=Smat\epsiplimat(:,:,i); %contrainte dans chaque pli dans le repère local
end

sigmaglo=zeros(3,1,taille);
for i=1:taille
    sigmaglo(:,:,i)=Tsig(:,:,i)*sigmaplimat(:,:,i); %contrainte dans chaque pli dans le repère global
end

critXc=zeros(1,taille);
for i=1:taille
    critXc(i)=sigmaplimat(1,1,i)/Xc; %critère selon Xc
end

critXt=zeros(1,taille);
for i=1:taille
    critXt(i)=sigmaplimat(1,1,i)/Xt; %critère selon Xt
end

critYc=zeros(1,taille);
for i=1:taille
    critYc(i)=sigmaplimat(2,1,i)/Yc; %critère selon Yc
end

critYt=zeros(1,taille);
for i=1:taille
    critYt(i)=sigmaplimat(2,1,i)/Yt; %critère selon Yt
end

critSc=zeros(1,taille);
for i=1:taille
    critSc(i)=abs(sigmaplimat(3,1,i))/Sc; %critère selon Sc
end

test=critXt>critXc; test2=critXt<critXc;
longXt=test.*critXt; longXc=test2.*critXc;
long=longXt+longXc; %critère longitudianal

test=critYt>critYc; test2=critYt<critYc;
transYt=test.*critYt; transYc=test2.*critYc;
trans=transYt+transYc; %critère transverse

cisail=critSc; %critère cisaillement

Rtrans1=zeros(1,taille);
Rtrans2=zeros(1,taille);
Rlong1=zeros(1,taille);
Rlong2=zeros(1,taille);
Rcis1=zeros(1,taille);
Rcis2=zeros(1,taille);
for i=1:taille %calcul du reserve factor
    Rtrans1(i)=Xt/sigmaplimat(1,1,i);
    Rlong1(i)=Yt/sigmaplimat(2,1,i);
    Rcis1(i)=Sc/sigmaplimat(3,1,i);
    Rtrans2(i)=Xc/sigmaplimat(1,1,i);
    Rlong2(i)=Yc/sigmaplimat(2,1,i);
    Rcis2(i)=-Sc/sigmaplimat(3,1,i);
end

test1=Rtrans1>0; test2=Rtrans2>0;
Rtrans=test1.*Rtrans1+test2.*Rtrans2; %Reserve factor dans chaque pli en transversal

test1=Rlong1>0; test2Yc=Rlong2>0;
Rlong=test1.*Rlong1+test2Yc.*Rlong2; %Reserve factor dans chaque pli en longitudinal

test1=Rcis1>0; test2=Rcis2>0;
Rcis=test1.*Rcis1+test2.*Rcis2; %Reserve factor dans chaque pli en cisaillement


%%
%----------------------résultats----------------------%
rows=cell(1,taille);
for i=1:taille
    angle=empilement(i);
    rows(i)={['pli nÂ°',num2str(i),' à  ',num2str(angle),'Â°']}; %création du nom des lignes des tableaux
end

fprintf('\nTableau du critère de rupture en contrainte maximale : \n')
Tablecont=table(long',trans',cisail','RowName',rows); %tableau contrainte maximale
Tablecont.Properties.VariableNames = {'Longitudinal' 'Transverse' 'Cisaillement'};
disp(Tablecont)

MatallR=[Rlong',Rtrans',Rcis'];
Matallcont=[long',trans',cisail']; 
ordrecassure=zeros(1,taille);
typedecassure=zeros(1,taille);
n=0;
K=1;
while n~=1 %détermination de l'ordre de cassure des plis, le programme s'arrete des qu'il y a une rupture catastrophique
    [maxval,idx]=max(Matallcont(:));
    [rowmax,colmax]=ind2sub(size(Matallcont),idx);
    ordrecassure(K)=rowmax;
    typedecassure(K)=colmax;
% for j=1:taille
    if typedecassure(K)==1
        cassure='longitudinal, la cassure est catastrophique.';
        n=1;
    elseif typedecassure(K)==2 && test2Yc(rowmax)==1
        cassure='transverse, la cassure est en compression donc catastrophique.';
        n=1;
    elseif typedecassure(K)==2 && test2Yc(rowmax)==0
        cassure='transverse, la cassure est en traction donc n est pas catastrophique.';
    elseif typedecassure(K)==3
        cassure='cisaillement, la cassure n est pas catastrophique.';
    end
% end
Matallcont(rowmax,:)=0;

disp('----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------')
disp(['Le pli n°',num2str(ordrecassure(K)),' à  ',num2str(empilement(ordrecassure(K))),'° casse en ',num2str(K),' en ',cassure,' Il possède un Rfactor de ',num2str(1/maxval),])
disp('Cela implique que le pli casse a un chargement de : ')
disp(['Nx=',num2str(1/maxval*N1),' N/mm  Ny=',num2str(1/maxval*N2),' N/mm  Nz=',num2str(1/maxval*N3),' N/mm'])
disp(['Mx=',num2str(1/maxval*M1),' N  My=',num2str(1/maxval*M2),' N  Mz=',num2str(1/maxval*M3),' N'])
K=K+1;
end

%%
%--------------Approche détaillé--------------%
% clearvars -except E11 G12 E22 nu12 ep Xt Xc Yt Yc Sc taille N M empilement htot

%Voici une liste de fonction a disposition pour étudier plus en détail le
%programme et l'étude du stratifié. Les données du problème sont déja
%disponible dans le workspace, il suffit de les reprendre ou alors vous
%pouvez changer les valeurs vous màªme dans le code si dessous pour étudier
%d'autres propriétées.

%voici une utilisation possible de ces fonctions :
% S=Soup(E11,E22,nu12,G12);
% Q=Rig(E11,E22,nu12,G12);
% [A,B,D]=mat(htot,empilement,Q);
% [eps0,kap0]=deformationcourbure(A,B,D,N,M);
% epsm=deformation(eps0,kap0,htot,empilement);
% sigmam=contrainte(epsm,Q,empilement)
% contraintesMax(sigmam,Xc,Xt,Yc,Yt,Sc);
% TsaiHill(sigmam,Xc,Xt,Yc,Yt,Sc);



function M=Soup(E11,E22,mu12,G12) %Cette fonction permet de calculer la matrice de souplesse en repère matériaux
M=zeros(3,3);
M(1,1)=1/E11;
M(1,2)=-mu12/E11;
M(2,1)=M(1,2);
M(2,2)=1/E22;
M(3,3)=1/G12;
end

function M=Rig(E11,E22,mu12,G12) %Cette fonction permet de calculer la matrice de rigidité en repère matériaux
M=inv(Soup(E11,E22,mu12,G12));
end

function Qglob=rigiditeGlob(theta,QrM) %Cette fonction permet de calculer la matrice de rigidité globale
c=cosd(theta); s=sind(theta);
Tq=[c^4 s^4 2*c^2*s^2 4*c^2*s^2;s^4 c^4 2*c^2*s^2 4*c^2*s^2;c^2*s^2 c^2*s^2 c^4+s^4 -4*c^2*s^2;c^2*s^2 c^2*s^2 -2*c^2*s^2 (c^2-s^2)^2; c^3*s -c*s^3 -c*s*(c^2-s^2) -2*c*s*(c^2-s^2);c*s^3 -c^3*s c*s*(c^2-s^2) 2*c*s*(c^2-s^2)];
Qglob=zeros(3);
vectQ=Tq*[QrM(1,1);QrM(2,2);QrM(1,2);QrM(3,3)];
Qglob(1,1)=vectQ(1);
Qglob(2,2)=vectQ(2);
Qglob(1,2)=vectQ(3);
Qglob(2,1)=Qglob(1,2);
Qglob(3,3)=vectQ(4);
Qglob(1,3)=vectQ(5);
Qglob(2,3)=vectQ(6);
Qglob(3,1)=Qglob(1,3);
Qglob(3,2)=Qglob(2,3);
end

function [A,B,D]=mat(h,vectTheta,Qm) %Cette fonction permet de calculer les matrices A B et D
A=zeros(3);B=zeros(3);D=zeros(3);
n=length(vectTheta);
hk=linspace(-(h/2),(h/2),n+1);
for k=1:n
    Qk=rigiditeGlob(vectTheta(k),Qm);
    A=A+(hk(k+1)-hk(k))*Qk;
    B=B+(1/2)*(hk(k+1)^2-hk(k)^2)*Qk;
    D=D+(1/3)*(hk(k+1)^3-hk(k)^3)*Qk;
end
end

function [eps0,kap0]=deformationcourbure(A,B,D,N,M) %Cette fonction permet de calculer les déformations de membrane et de courbure
eps0=zeros(3,1);kap0=zeros(3,1);
res=[A B;B D]\[N;M];
eps0(1)=res(1);eps0(2)=res(2);eps0(3)=res(3);kap0(1)=res(4);kap0(2)=res(5);kap0(3)=res(6);
end

function epsm=deformation(eps0,kap0,h,vectTheta) %Cette fonction permet de calculer les déformations de chaque pli dans le repère matériaux
n=length(vectTheta);
eps=zeros(3,n);epsm=zeros(3,n);
hk=linspace(-(h/2),(h/2),n+1);
for k=1:n
    zk=0.5*(hk(k)+hk(k+1));
    eps(:,k)=eps0+zk*kap0;
    c=cosd(vectTheta(k));s=sind(vectTheta(k));
    Tk=[c^2 s^2 s*c;s^2 c^2 -s*c;-2*s*c 2*s*c c^2-s^2];
    epsm(:,k)=Tk*eps(:,k);
end
 
end

function sigmam=contrainte(epsm,Q,vectTheta) %Cette fonction permet de calculer les contraintes de chaque pli dans le repère matériaux
n=length(vectTheta);
sigmam=zeros(3,n);
for k=1:n
    sigmam(:,k)=Q*epsm(:,k);
end
end

function contraintesMax(sigmam,Xc,Xt,Yc,Yt,Sc) %Cette fonction permet d'afficher la surface de rupture pour le critère en contrainte maximale
figure();
[~,m]=size(sigmam);
subplot(1,2,1);
rectangle('Position',[Xc Yc -Xc+Xt -Yc+Yt],'EdgeColor','r')
hold on
for k=1:m
    plot(sigmam(1,k),sigmam(2,k),'b+')
end
hold off
subplot(1,2,2);
rectangle('Position',[Xc -Sc -Xc+Xt 2*Sc],'EdgeColor','r')
hold on
for k=1:m
    plot(sigmam(1,k),sigmam(3,k),'b+')
end
hold off
subplot(1,2,1);
hold on
plot([Xc Xt 0 0],[0 0 Yc Yt],'m*')
xlabel('\sigma_1','FontSize',20,'FontWeight','bold'); ylabel('\sigma_2','FontSize',20,'FontWeight','bold');
text(Xc+150,0,'X_c','color','magenta'); text(Xt+150,0,'X_t','color','magenta'); text(0,Yc+10,'Y_c','color','magenta'); text(0,Yt+10,'Y_t','color','magenta');
hold off
title('Critère de la contrainte maximale');
subplot(1,2,2);
hold on
plot([Xc Xt 0 0],[0 0 -Sc Sc],'m*')
xlabel('\sigma_1','FontSize',20,'FontWeight','bold'); ylabel('\sigma_6','FontSize',20,'FontWeight','bold');
text(Xc+150,0,'X_c','color','magenta'); text(Xt+150,0,'X_t','color','magenta'); text(0,-Sc+10,'-S_c','color','magenta'); text(0,Sc+10,'S_c','color','magenta');
hold off
end

function plotEllipse(x0,y0,a,b)
t=-pi:0.1:pi;
x=x0+a*cos(t);
y=y0+b*sin(t);
plot(x,y,'r','LineWidth',1)
end

function TsaiHill(sigmam,Xc,Xt,Yc,Yt,Sc) %Cette fonction permet d'afficher la surface de rupture pour le critère de Tsai-Hill
figure();
[~,m]=size(sigmam);
subplot(1,2,1);
plotEllipse(0,0,Xt,Yt);
hold on
for k=1:m
    plot(sigmam(1,k),sigmam(2,k),'b+')
end
hold off
subplot(1,2,2);
plotEllipse(0,0,Xt,Sc)
hold on
for k=1:m
    plot(sigmam(1,k),sigmam(3,k),'b+')
end
hold off
subplot(1,2,1);
hold on
plot([Xc Xt 0 0],[0 0 Yc Yt],'m*')
xlabel('\sigma_1','FontSize',20,'FontWeight','bold'); ylabel('\sigma_2','FontSize',20,'FontWeight','bold');
text(Xc+150,0,'X_c','color','magenta'); text(Xt+150,0,'X_t','color','magenta'); text(0,Yc+10,'Y_c','color','magenta'); text(0,Yt+10,'Y_t','color','magenta');
hold off
title('Critère de Tsai-Hill');
subplot(1,2,2);
hold on
plot([Xc Xt 0 0],[0 0 -Sc Sc],'m*')
xlabel('\sigma_1','FontSize',20,'FontWeight','bold'); ylabel('\sigma_6','FontSize',20,'FontWeight','bold');
text(Xc+150,0,'X_c','color','magenta'); text(Xt+150,0,'X_t','color','magenta'); text(0,-Sc+10,'-S_c','color','magenta'); text(0,Sc+10,'S_c','color','magenta');
hold off
end
