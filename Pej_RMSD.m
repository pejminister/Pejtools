function E = Pej_RMSD(X,Y)
E = sqrt(mean((X-Y).^2));
end