function X = AXeqXB(listeA, listeB)

liste_uA = [];
liste_uB = [];

for ii=1:length(listeA)
    [~, uA] = r2thetau(listeA{ii});
    liste_uA = [liste_uA, uA];
end

for ii=1:length(listeB)
    [~, uB] = r2thetau(listeB{ii});
    liste_uB = [liste_uB, uB];
end

%liste_uA
%liste_uB

T = loc3D_opt(liste_uA, liste_uB)

Rx = T(1:3,1:3);

for ii=1:length(listeA)
    M(3*ii-2:3*ii,1:3) = listeA{ii}(1:3,1:3) - eye(3);
    N(3*ii-2:3*ii,1) = Rx*listeB{ii}(1:3,4) - listeA{ii}(1:3,4);
end

tx = pinv(M) * N;

X = [Rx tx; 0 0 0 1];

end

