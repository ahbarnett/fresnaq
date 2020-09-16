% figure out what makePupil does, since it has no docs. Barnett 9/15/20.

Nx = 100;
Ny = 100;
% see calling seq in makeStarshadeImage.m, l.78.
OR = 1;
IR = 0.417;  % somehow this gets set to this, can't figure out how
deltax = 0.3; deltay = -0.2;    % generic
pupil = makePupil(Nx, Ny, OR, IR, deltax, deltay);
figure; imagesc(pupil); axis equal;
