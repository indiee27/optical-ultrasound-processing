function pASA_tot = propagate_ASA_multi(p0,dz,dx,dy,dt,c,oversample)
% size(p0) = [Nx,Ny,Nt];

if size(oversample) == 1
    oversample = [oversample,oversample,oversample];
end

p0 = cat(1,p0,zeros(size(p0).*[oversample(1)-1 1 1]));
p0 = cat(2,p0,zeros(size(p0).*[1 oversample(2)-1 1]));
p0 = cat(3,p0,zeros(size(p0).*[1 1 oversample(3)-1]));

p0 = fftn(p0);

[Nx,Ny,Nt] = size(p0);

dk  = 2*pi / (Nt*dt*c);
dkx = 2*pi / (Nx*dx);
dky = 2*pi / (Ny*dy);

k  = dk  * [0:ceil(Nt/2)-1,-floor(Nt/2):-1];
kx = dkx * [0:ceil(Nx/2)-1,-floor(Nx/2):-1];
ky = dky * [0:ceil(Ny/2)-1,-floor(Ny/2):-1];

kx = kx'*ones(1,Ny);
ky = ones(1,Nx)'*ky;
kx2ky2 = (kx.^2) + (ky.^2);
k2 = k.^2;
clear kx ky dk dkx dky dx dy;

Ndz = length(dz);
pASA = zeros(size(p0));
pASA_tot = zeros([size(p0),Ndz]);
w = waitbar(0,'Propagation progress:');
for dzcnt = 1:Ndz
    if dz == 0
        pASA_tot(:,:,:,dzcnt) = p0;
    else
        for fcnt = 1:Nt
            exponent = exp(-1i * sign(k(fcnt)) * dz(dzcnt) * conj(sqrt(k2(fcnt)-kx2ky2)) );
        %     exponent(kx2ky2 > k2(fcnt)) = 0;
            exponent(abs(exponent)>1) = 0;
            pASA(:,:,fcnt) = squeeze(p0(:,:,fcnt)) .* exponent;
        end
        pASA_tot(:,:,:,dzcnt) = real(ifftn(pASA));
    end
    waitbar(dzcnt/Ndz,w);
end
close(w);

Nx = Nx/oversample(1);
Ny = Ny/oversample(2);
Nt = Nt/oversample(3);
pASA_tot = pASA_tot(1:Nx,1:Ny,1:Nt,:);

clear P0 exponent fcnt k k2 kx2ky2 oversample;