%
% Calcul des contraintes et de l'indice de fissuration
%  dans le cas du retrait endogène empéché
%  prise en compte du fluage propre
%   voir document dans le répertoire
%


% testes git amanda

clc; clear; close all; tic
%
% Données matériaux voir fichier Excel
%
% Laitier activé
s = 0.374  ; t0 =  0.114*24*3600.; % Parametres communs [] et temps prise [s]
E28 = 22.8E9 ; nE = 0.546 ; % [Pa]
Retrait_endogene28 = -204.E-6 ; nretrait_endogene = 4.44;
Retrait_dessiccation28 = 0. * (-457.E-6) ; kt = 0.0343 ; ks = 1.25 ;
ft28 = 3.01E6 ; nft = 0.808;

%%%%%%%%%%%%%%%%%%%%%%%%%% FLUAGE
% Chaine de Kelvin-Voigt
t28 = 28.*24 * 3600 ;
E_kv28 = 1.E20* 21.1141*1.E9 ; n_E_kv = 0.3 ;
Eta_kv28 = 1.E20*6.82*24*3600 * E_kv28 ; n_Eta_kv = 0.3 ;
% Amortisseur seul (valeur à 28 jours)
eta28 = 1.E20*2543.6*1E9*24 * 3600 ; % [Pa.s]
% Energie d'activation pour le fluage propre
Ea_R_creep = 16.4E3/8.314 ; % [K]
%%%%%%%%%%%%%%%%%%%%%%%%%% FIN FLUAGE

% Données thermo-chimiques
parameters_material = zeros(6,1); affinity = zeros(3,1) ; conditions = zeros(4,1) ;
parameters_material(1) = 2188 ;                % Masse volumique kg/m3
parameters_material(2) = 519.;                 % Capacite calorifique massique J/kg/K
parameters_material(3) = 3.1;           % Conductivite thermique W/m/K (k)
parameters_material(4) = 77.E3/8.314;  % Énergie d'activation/R en Kelvin
Ea_R_thermique = parameters_material(4);

% Paramètres du temps équivalent
param_teq(1) = 203.7 ; % paramètre C [J/g]
param_teq(2) = -99.7 ; % paramètre Q1 [J/g]
param_teq(3) = 13.5 ; % paramètre t1 [heures]
param_teq(4) = 1. ; % paramètre dt1 [heures]
param_teq(5) = -223 ; % paramètre Q2 [J/g]
param_teq(6) = -1.2 ; % paramètre t2 [heures]
param_teq(7) = 8.8 ; % paramètre dt2 [heures]
param_teq(8) = 298.15 ; % Température de référence Arrhenius
Tref = param_teq(8);
param_teq(9) = 300. ; % Quantité de liant [kg/m3]

% Géométrie et Conditions aux limites
epaisseur = 1. ; % Epaisseur en [metres]
conditions(1) = 20+273.15 ;            % Temperature initiale [K]
conditions(2) = 20+273.15 ;          % Temperature exterieure [K]
conditions(3) = 17.5;                    % Coefficient de convection/rayonnement en W/(m^2K)
% Paramètre spatial
conditions(4) = 300 ;               % Nombre de pas d'espace
% Coefficient de dilatation thermique
coef_dilat = 0.*10.E-6 ;

% Données géométriques et CL
Rayon_sechage = 2*70*70/(2*70);   % [mm] (h0 = 2 V / S)
humidite_relative = 0.5 ; tsechage = 24 * 3600. ; % temps où le séchage démarre
% Discrétization temporelle
tfin = 28 * 24 * 3600 ; % Final time [secondes]
npastemps = 200; dt = (tfin - t0) / npastemps;

% Initialisation des variables
% FLUAGE
sigma_fluage = zeros(npastemps + 1,1) ; sigma_elas = zeros(npastemps + 1,1) ; sigma_fluage_amor = zeros(npastemps + 1,1) ;
sigma_kv = zeros(npastemps + 1,1) ;
% DEFORMATION LIBRE
Retrait_endogene = zeros(npastemps + 1,1) ;
Retrait_dessiccation = zeros(npastemps + 1,1) ;
Deformation_thermique = zeros(npastemps + 1,1) ;
% PROPRIETES MECANIQUES
ft = zeros(npastemps + 1,1) ; E = zeros(npastemps + 1,1) ;
% Indice de Sensibilité à la Fissuration (Cracking Index)
ISF_fluage = zeros(npastemps + 1,1) ; ISF_fluage_amor = zeros(npastemps + 1,1) ;
ISF_elas = zeros(npastemps + 1,1) ;
% Autre vecteurs
temps_eq = zeros(npastemps + 1,1) ; temps_eq(1) = 0; % Temps équivalent
Tmoy = zeros(npastemps + 1,1) ; % Température moyenne
Real_Time = zeros(npastemps + 1,1) ; Real_Time(1) = 0; % Temps réel

% Calculs préliminaires (drying shrinkage)
kh = 1 - ((humidite_relative)^3) ; % HR < 0.98
taush = kt * ((ks * Rayon_sechage)^2) ; % Tau_ds

% Boucle pour le calcul des contraintes
% en conditions parfaitement restreintes
% Attention pour le calcul il est faux avant la prise

% Calcul thermique
T = resolution_initial(parameters_material,param_teq,tfin,epaisseur,conditions,npastemps) ;
Tmoy = mean(T);
Tmax = max(T);

for i=2:npastemps+1
    Real_Time(i) = Real_Time(i-1) + dt ;
    % Calcul des propriétés mécaniques et retrait
    Tmoy_moy = 273.15 + ((Tmoy(i) + Tmoy(i-1) )/2.) ;
    coef_activation_thermique = exp(Ea_R_thermique*( 1/Tref - 1/Tmoy_moy) ) ;
    temps_eq(i) = temps_eq(i-1) + (dt*coef_activation_thermique) ; % Temps équivalent
    t = (temps_eq(i) + temps_eq(i-1) ) / 2. ;
    temp1 = (0.5 * ( (abs (t - tsechage) ) + (t - tsechage) )) / (24*3600.) ;
    S = tanh( sqrt(temp1 / taush) ) ;
    
    
    if t>t0,
        beta = exp(s*(1 - sqrt( (28*24*3600 - t0) /(t - t0) ) ) );
        
        E(i) = E28 * (beta^nE) ;
        ft(i) = ft28 * (beta^nft) ;
        
        % Calcul des retraits et deformation thermique
        Retrait_endogene(i) = Retrait_endogene28 * (beta^nretrait_endogene) ;
        Retrait_dessiccation(i) = Retrait_dessiccation28 * kh * S;
        Deformation_thermique(i) = coef_dilat * ( Tmoy(i) - Tmoy(1) ) ;
        Retrait_total_1 = Retrait_endogene(i-1) + Retrait_dessiccation(i-1) +  Deformation_thermique(i-1) ;
        Retrait_total_2 = Retrait_endogene(i) + Retrait_dessiccation(i) +  Deformation_thermique(i);
        
        % Calcul des contraintes et de l'ISF élastique
        dsigma = -E(i) *  (Retrait_total_2 - Retrait_total_1);
        sigma_elas(i) = sigma_elas(i-1) + dsigma;
        ISF_elas(i) = sigma_elas(i)/ft(i) ;
        if ISF_elas(i)<0,ISF_elas(i)=0.;end;
        
        % Activation thermique du fluage
        coef_activation_fluage = exp(Ea_R_creep*( 1/Tref - 1/Tmoy_moy) ) ;
        
        % Fluage amortisseur (irreversible)
        eta = eta28 * (t / (t28) ) * coef_activation_fluage;
        
        % Calcul des contraintes et de l'ISF avec fluage (amortisseur seul)
        denom1 = (1/ E(i) + (0.5 * (dt /eta) ) ) ;
        numer1 = (sigma_fluage_amor(i-1) * dt/eta) + (Retrait_total_2 - Retrait_total_1) ;
        dsigma1 = - numer1 / denom1;
        sigma_fluage_amor(i) = sigma_fluage_amor(i-1) + dsigma1;
        ISF_fluage_amor(i) = sigma_fluage_amor(i)/ft(i) ;
        if ISF_fluage_amor(i)<0,ISF_fluage_amor(i)=0.;end;
        
        % Fluage Kelvin-Voigt (reversible)
        E_kv = E_kv28 * (beta^n_E_kv) ;
        Eta_kv = Eta_kv28 * (beta^n_Eta_kv) * coef_activation_fluage ;
        tau = Eta_kv / E_kv;
        beta = (dt * exp (-dt/(2*tau) ) ) / Eta_kv ;
        %    beta = ((exp(dt/(tau))) - 1) * (exp(-temps(i-1)/tau)) / E_kv ;
        beta_prime = beta * (sigma_fluage(i-1) - sigma_kv(i-1));
        
        denom2 = (1/ E(i) + (0.5 * (dt /eta) + (beta/2.) ) ) ;
        numer2 = (sigma_fluage(i-1) * dt/eta) + (Retrait_total_2 - Retrait_total_1) + (beta_prime) ;
        dsigma2 = - numer2 / denom2;
        sigma_fluage(i) = sigma_fluage(i-1) + dsigma2;
        ISF_fluage(i) = sigma_fluage(i)/ft(i) ;
        if ISF_fluage(i)<0,ISF_fluage(i)=0.;end;
        
        % Mise à jour de la contrainte dans la chaine de KV
        lambda = ( sigma_kv(i-1) - sigma_fluage(i-1) - (dsigma2 / 2.) ) * (exp(temps_eq(i-1)/tau)) ;
        sigma_kv(i) = sigma_fluage(i-1) + (dsigma2 / 2.) + (lambda * (exp(-temps_eq(i)/tau))) ;
        a = 2 ;
        
    else,
        E(i)=0;
        ft(i)=1E-10;
    end;
    
    
end;

toc

Tmoy = Tmoy' ;
Tmax = Tmax' ;

figure; hold on; plot(Real_Time/86400.,1.E6*Retrait_endogene,'r'); plot(Real_Time/86400.,1.E6*Retrait_dessiccation,'b'); plot(Real_Time/86400.,1.E6*Deformation_thermique,'g'); title('Retrait endogène (rouge), dessiccation (bleu) et thermique (vert) [µm/m]'); hold off;
figure; hold on; plot(Real_Time/86400.,E/1E9,'r'); plot(Real_Time/86400.,ft/1E6,'b');  title('Module de Young [GPa] = rouge - Résistance en traction [MPa] = bleu'); hold off;
figure; hold on; plot(Real_Time/86400.,sigma_elas/1.E6,'r'); plot(Real_Time/86400,sigma_fluage/1.E6,'b') ; plot(Real_Time/86400,sigma_fluage_amor/1.E6,'k') ; plot(Real_Time/86400.,ft/1E6,'g'); title('Contraintes [MPa]: sans fluage = rouge - avec fluage amor = noir - avec fluage total = bleu / Résistance en traction [MPa] = vert'); hold off;
figure; hold on; plot(Real_Time/86400,ISF_elas,'r'); plot(Real_Time/86400,ISF_fluage,'b') ; plot(Real_Time/86400,ISF_fluage_amor,'g') ;  title('Indice de fissuration : sans fluage = rouge - avec fluage amor = vert - avec fluage = bleu'); hold off;
figure; hold on; plot(Real_Time/86400,Tmoy,'r'); plot(Real_Time/86400,Tmax,'b'); title('average temperature = Rouge - maximum temperature = Bleu'); hold off;


