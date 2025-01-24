n = 100; % Nombre de valeurs propres de A

figure;

for imat = 1:4
    [A, D, info] = matgen_csad(imat, n);
    
    eigenvalues = eig(A);
    eigenvalues = sort(eigenvalues, 'descend');
    
    cond_A = cond(A);
    
    subplot(2, 2, imat);
    scatter(1:n, eigenvalues);
    text(n/2, max(eigenvalues), ['Cond = ', num2str(cond_A)], 'HorizontalAlignment', 'center');
    title(['Type ', num2str(imat)]);
    xlabel('Indice');
    ylabel('Valeur propre');
end