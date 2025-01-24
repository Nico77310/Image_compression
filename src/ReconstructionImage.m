%%  Application de la SVD : compression d'images

clear all
close all

% Lecture de l'image
I = imread('BD_Asterix_1.png');
I = rgb2gray(I);
I = double(I);

[q, p] = size(I);

% Décomposition par SVD
fprintf('Décomposition en valeurs singulières\n')
tic
[U, S, V] = svd(I);
toc

l = min(p,q);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% On choisit de ne considérer que 200 vecteurs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vecteur pour stocker la différence entre l'image et l'image reconstuite
inter = 1:40:(200+40);
inter(end) = 200;
differenceSVD = zeros(size(inter,2), 1);

% images reconstruites en utilisant de 1 à 200 vecteurs (avec un pas de 40)
ti = 0;
td = 0;
for k = inter

    % Calcul de l'image de rang k
    Im_k = U(:, 1:k)*S(1:k, 1:k)*V(:, 1:k)';

    % Affichage de l'image reconstruite
    ti = ti+1;
    figure(ti)
    colormap('gray')
    imagesc(Im_k), axis equal
    
    % Calcul de la différence entre les 2 images
    td = td + 1;
    differenceSVD(td) = sqrt(sum(sum((I-Im_k).^2)));
    pause
end

% Figure des différences entre image réelle et image reconstruite
ti = ti+1;
figure(ti)
hold on 
plot(inter, differenceSVD, 'rx')
ylabel('RMSE')
xlabel('rank k')
pause


%% Plugger les différentes méthodes : eig, puissance itérée et les 4 versions de la "subspace iteration method" 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUELQUES VALEURS PAR DÉFAUT DE PARAMÈTRES, 
% VALEURS QUE VOUS POUVEZ/DEVEZ FAIRE ÉVOLUER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tolérance
eps = 1e-8;
% nombre d'itérations max pour atteindre la convergence
maxit = 10000;

% taille de l'espace de recherche (m)
search_space = 400 ;

% pourcentage que l'on se fixe
percentage = 0.998;

% p pour les versions 2 et 3 (attention p déjà utilisé comme taille)
puiss = 2;

%%%%%%%%%%%%%
% À COMPLÉTER
%%%%%%%%%%%%%

%
% calcul des couples propres
if q > p
    M = I'*I;
    [ vect_v, D, n_ev, ~, ~, ~] = subspace_iter_v3( M, search_space, percentage, puiss, eps, maxit );
else 
    M = I*I';
    [ vect_u, D, n_ev, ~, ~, ~] = subspace_iter_v3( M, search_space, percentage, puiss, eps, maxit );
end

%
% calcul des valeurs singulières
Singular_Values = sqrt(diag(D));

%
% calcul de l'autre ensemble de vecteurs
if q > p
    vect_u = repmat(1./Singular_Values',q,1) .* (I*vect_v);
else
    vect_v = repmat(1./Singular_Values',p,1) .* (I'*vect_u);
end

%
% calcul de I_k

Sigma = diag(Singular_Values);

k = 150; % définition de k 

I_k = vect_u(:,1:k)*Sigma(1:k,1:k)*vect_v(:,1:k)';

% Affichage de l'image reconstruite
ti = 0;
    ti = ti+1;
    figure(ti)
    colormap('gray')
    imagesc(I_k), axis equal

sqrt(sum(sum((I-I_k).^2)))
%%
% calcul des meilleures approximations de rang faible
%%
k_values = [];
low_rank_approximation = [];
for k = 1:10:n_ev
    k_values = [k_values, k];
    I_low_rank = vect_u(:,1:k)*Sigma(1:k,1:k)*vect_v(:,1:k)';
    low_rank_approximation = [low_rank_approximation, sqrt(sum(sum((I-I_low_rank).^2)))];
end
figure 
plot(k_values,low_rank_approximation)
ylabel('RMSE')
xlabel('rank k')
grid on
