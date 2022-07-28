function T=resolution_initial(parameters_material,param_teq,tfin,epaisseur,conditions,npastemps) ;
% La fonction effectue la résolution de l'équation de la chaleur pour un jeu de paramètres donnés
% Température [K] T(i,j) : i = pas d'espace et j = pas de temps

%%%%%%%%%%%%%%%%% Définition des constantes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Géométrie et Conditions aux limites
epaisseur = 1.;           % Epaisseur [m]
T0 = conditions(1) ;        % Temperature initiale [K]
Text = conditions(2) ;      % Temperature exterieure [K]
h = conditions(3);                 % Coefficient d'échange par convection [W/(m².K)]
duree_etude = tfin;  % Temps final [s]
Tref = param_teq(8) ;           % Température de référence Arrhenius [K]

% Paramètres matériaux
rho = parameters_material(1) ;              % Masse volumique [kg/m3]
Cp = parameters_material(2);                % Capacite calorifique massique [J/kg/K]
lambda = parameters_material(3);             % Conductivite thermique [W/m/K]
EaR = parameters_material(4) ; % Énergie d'activation/R [K]

% Composition
Ciment = param_teq(9) ;           % Quantité de ciment [kg/m3 béton]

% Paramètre de dégagement de chaleur
% param_teq

% Paramètres numériques
Nombre_pas_espace = conditions(4); % Nombre de pas d'espace pour la demi-épaisseur
Nombre_pas_temps = npastemps + 1;
theta = 1. ; % Theta méthode (0 = explicite / 0.5 = Cranck-Nichololson / 1 = implicite) => il faut une valeur supérieure ou égale à 0.5 pour une convergence inconditionnement stable

pas_espace = epaisseur/(2*Nombre_pas_espace);     % en metres
pas_temps  = duree_etude/Nombre_pas_temps;        % en secondes

alpha = pas_temps*lambda/(pas_espace^2 *rho*Cp);    % pas d'unités
beta = lambda/(h*pas_espace);                       % pas d'unités

x1 = 0.:(epaisseur/(2.*Nombre_pas_espace)):(epaisseur/2.); % Définition de l'abscisse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Matrices bandes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = sparse(1:(Nombre_pas_espace + 1),1:(Nombre_pas_espace + 1),(1 + (2*alpha*theta) )*ones(1,(Nombre_pas_espace + 1)),(Nombre_pas_espace + 1),(Nombre_pas_espace + 1));
E = sparse(2:(Nombre_pas_espace + 1),1:(Nombre_pas_espace + 1)-1,(-alpha*theta)*ones(1,(Nombre_pas_espace + 1)-1),(Nombre_pas_espace + 1),(Nombre_pas_espace + 1));
S = E+D+transpose(E); % Forme générale
S(1,1) = 1; S(1,2) = -1; % Symétrie
S((Nombre_pas_espace + 1),(Nombre_pas_espace + 1)-1) = -beta; S((Nombre_pas_espace + 1),(Nombre_pas_espace + 1))= beta+1; % Conditions aux limites
A = full(S); % Représentation complète de la matrice

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialisation %%%%%%%%%%%%%%%%%%%%%%%
% Matrice de temperature : T(i,j) : i = pas d'espace et j = pas de temps
T = zeros((Nombre_pas_espace + 1),Nombre_pas_temps) ;
T(:,1) = T0 ;                                              % Température initiale au premier pas de temps
Tu = zeros((Nombre_pas_espace + 1),1) ;                    % Matrice de second membre - Initialisation

teq = zeros((Nombre_pas_espace + 1),Nombre_pas_temps) ;    % Temps équivalent initialisé à 0 [secondes]
Q = zeros((Nombre_pas_espace + 1),Nombre_pas_temps) ;      % Quantité de chaleur initialisée à 0 [J/m3]
temps = zeros(1,Nombre_pas_temps) ;                        % Temps (réel) initialisé à 0 [secondes]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Résolution du problème thermique %%%%%%%%%%%%%%%%%%%%%%%%

for n = 2:Nombre_pas_temps
    temps(n) = temps(n-1) + pas_temps; % Mise jour du temps réel
    teq(:,n) = teq(:,n-1) + ( (exp(-EaR* (1./T(:,n-1) -1./Tref) ) )*pas_temps ); % Calcul du temps équivalent
    
    for i=2:Nombre_pas_espace
        Q(i,n) = Quantite_chaleur(param_teq,Ciment,teq,i,n);
        Tu(i) = T(i,n-1)+ (Q(i,n)-Q(i,n-1))/(rho*Cp);                                     % Second membre ( théta = 1 )
        Tu(i) = Tu(i) + alpha*(1 - theta) * (T(i+1,n-1) + T(i-1,n-1) - (2.*T(i,n-1)) ) ;  % Second membre ( théta différent de 1 )
    end
    Q(1,n) = Quantite_chaleur(param_teq,Ciment,teq,1,n);
    Q((Nombre_pas_espace+1),n) = Quantite_chaleur(param_teq,Ciment,teq,(Nombre_pas_espace+1),n);
  
    Tu(1)= 0;                                 % Gestion de la symétrie
    Tu((Nombre_pas_espace + 1)) = Text;       % Gestion de la conditions aux limites convectives
    
    T(:,n) = A\Tu ;               % On calcule ensuite la température au nième pas de temps sur tous les pas d'espace
    
end

T = T - 273.15*ones((Nombre_pas_espace + 1),Nombre_pas_temps); % Température [K] => [°C]
