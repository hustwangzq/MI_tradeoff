function H = SIM_diffraction_matrix(rxPos, txPos, lambda, At)
%DIFFRACTION_MATRIX Rayleigh-Sommerfeld near-field diffraction matrix.
% rxPos and txPos are located on parallel yz planes, and propagation is along x.

Nr = size(rxPos,1); Nt = size(txPos,1);
H = zeros(Nr,Nt);
for r = 1:Nr
    for t = 1:Nt
        vecRT = rxPos(r,:) - txPos(t,:);
        dist = norm(vecRT);
        coschi = abs(vecRT(1))/max(dist,1e-12);
        H(r,t) = (At*coschi/dist) * (1/(2*pi*dist) - 1j/lambda) * exp(1j*2*pi*dist/lambda);
    end
end
end
