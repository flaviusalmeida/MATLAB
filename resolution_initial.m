function T=resolution_initial(parameters_material,param_teq,tfin,epaisseur,conditions,npastemps) ;
% La fonction effectue la r�solution de l'�quation de la chaleur pour un jeu de param�tres donn�s
% Temp�rature [K] T(i,j) : i = pas d'espace et j = pas de temps

%%%%%%%%%%%%%%%%% D�finition des constantes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% G�om�trie et Conditions aux limites
epaisseur = 1.;           % Epaisseur [m]
T0 = conditions(1) ;        % Temperature initiale [K]
Text = conditions(2) ;      % Temperature exterieure [K]
h = conditions(3);                 % Coefficient d'�change par convection [W/(m�.K)]
duree_etude = tfin;  % Temps final [s]
Tref = param_teq(8) ;           % Temp�rature de r�f�rence Arrhenius [K]

% Param�tres mat�riaux
rho = parameters_material(1) ;              % Masse volumique [kg/m3]
Cp = parameters_material(2);                % Capacite calorifique massique [J/kg/K]
lambda = parameters_material(3);             % Conductivite thermique [W/m/K]
EaR = parameters_material(4) ; % �nergie d'activation/R [K]

% Composition
Ciment = param_teq(9) ;           % Quantit� de ciment [kg/m3 b�ton]

% Param�tre de d�gagement de chaleur
% param_teq

% Param�tres num�riques
Nombre_pas_espace = conditions(4); % Nombre de pas d'espace pour la demi-�paisseur
Nombre_pas_temps = npastemps + 1;
theta = 1. ; % Theta m�thode (0 = explicite / 0.5 = Cranck-Nichololson / 1 = implicite) => il faut une valeur sup�rieure ou �gale � 0.5 pour une convergence inconditionnement stable

pas_espace = epaisseur/(2*Nombre_pas_espace);     % en metres
pas_temps  = duree_etude/Nombre_pas_temps;        % en secondes

alpha = pas_temps*lambda/(pas_espace^2 *rho*Cp);    % pas d'unit�s
beta = lambda/(h*pas_espace);                       % pas d'unit�s

x1 = 0.:(epaisseur/(2.*Nombre_pas_espace)):(epaisseur/2.); % D�finition de l'abscisse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Matrices bandes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = sparse(1:(Nombre_pas_espace + 1),1:(Nombre_pas_espace + 1),(1 + (2*alpha*theta) )*ones(1,(Nombre_pas_espace + 1)),(Nombre_pas_espace + 1),(Nombre_pas_espace + 1));
E = sparse(2:(Nombre_pas_espace + 1),1:(Nombre_pas_espace + 1)-1,(-alpha*theta)*ones(1,(Nombre_pas_espace + 1)-1),(Nombre_pas_espace + 1),(Nombre_pas_espace + 1));
S = E+D+transpose(E); % Forme g�n�rale
S(1,1) = 1; S(1,2) = -1; % Sym�trie
S((Nombre_pas_espace + 1),(Nombre_pas_espace + 1)-1) = -beta; S((Nombre_pas_espace + 1),(Nombre_pas_espace + 1))= beta+1; % Conditions aux limites
A = full(S); % Repr�sentation compl�te de la matrice

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialisation %%%%%%%%%%%%%%%%%%%%%%%
% Matrice de temperature : T(i,j) : i = pas d'espace et j = pas de temps
T = zeros((Nombre_pas_espace + 1),Nombre_pas_temps) ;
T(:,1) = T0 ;                                              % Temp�rature initiale au premier pas de temps
Tu = zeros((Nombre_pas_espace + 1),1) ;                    % Matrice de second membre - Initialisation

teq = zeros((Nombre_pas_espace + 1),Nombre_pas_temps) ;    % Temps �quivalent initialis� � 0 [secondes]
Q = zeros((Nombre_pas_espace + 1),Nombre_pas_temps) ;      % Quantit� de chaleur initialis�e � 0 [J/m3]
temps = zeros(1,Nombre_pas_temps) ;                        % Temps (r�el) initialis� � 0 [secondes]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% R�solution du probl�me thermique %%%%%%%%%%%%%%%%%%%%%%%%

for n = 2:Nombre_pas_temps
    temps(n) = temps(n-1) + pas_temps; % Mise jour du temps r�el
    teq(:,n) = teq(:,n-1) + ( (exp(-EaR* (1./T(:,n-1) -1./Tref) ) )*pas_temps ); % Calcul du temps �quivalent
    
    for i=2:Nombre_pas_espace
        Q(i,n) = Quantite_chaleur(param_teq,Ciment,teq,i,n);
        Tu(i) = T(i,n-1)+ (Q(i,n)-Q(i,n-1))/(rho*Cp);                                     % Second membre ( th�ta = 1 )
        Tu(i) = Tu(i) + alpha*(1 - theta) * (T(i+1,n-1) + T(i-1,n-1) - (2.*T(i,n-1)) ) ;  % Second membre ( th�ta diff�rent de 1 )
    end
    Q(1,n) = Quantite_chaleur(param_teq,Ciment,teq,1,n);
    Q((Nombre_pas_espace+1),n) = Quantite_chaleur(param_teq,Ciment,teq,(Nombre_pas_espace+1),n);
  
    Tu(1)= 0;                                 % Gestion de la sym�trie
    Tu((Nombre_pas_espace + 1)) = Text;       % Gestion de la conditions aux limites convectives
    
    T(:,n) = A\Tu ;               % On calcule ensuite la temp�rature au ni�me pas de temps sur tous les pas d'espace
    
end

T = T - 273.15*ones((Nombre_pas_espace + 1),Nombre_pas_temps); % Temp�rature [K] => [�C]
