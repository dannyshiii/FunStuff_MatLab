function [vid, vc, vlow, vsparse] = DMD(U, S, V, vid, v2, tmp, time)
% This function contructs a Dynamic Mode Decomposition for 
% subtracting foreground and background of a video.
% First you need to convert your video file into a .mat file and 
% don't forget it needs to be in gray-scale.
% Then, you do a Singular Value Decomposition, and choose the corresponding
% mode / rank. Let U=U(:,1:mode),S=S(1:mode,1:mode),V=V(:,1:mode)
% Let the new U, S, V be the inputs here.
% tmp is the size of the frame. [r,c]=size(frame), tmp=r*c
% time is a specific frame you want to subtract at certain time.

Atilde = U.' * v2 * V / S;
[W,D] = eig(Atilde);
Phi = v2 * V/S * W;
t = linspace(0, 1, time);
dt = t(2)-t(1);
lambda = diag(D);
omega = log(lambda)/dt;
y = Phi \vid(:, 1);
%%
backg = find(abs(omega) == min(abs(omega)));
foreg = find(abs(omega) ~= min(abs(omega)));
omega1 = omega(backg);
omega2 = omega(foreg);
y1 = y(backg);
y2 = y(foreg);
%%
r1 = length(backg);
r2 = length(foreg);
uback = zeros(r1, length(t));
ufore = zeros(r2, length(t));
for k = 1 : length(t)
    uback(:, k) = (y1.*exp(omega1*(t(k))));
    ufore(:, k) = (y2.*exp(omega2*(t(k))));
end
%%
vlow = Phi(:, backg) * uback;
vsparse = Phi(:, foreg) * ufore;
%%
vc = vlow + vsparse;
vsparse = real(vc - abs(vlow));
R = zeros(tmp, time);
for k = 1:tmp
    idx = find(vsparse(k, :) < 0);
    R(k, idx) = vsparse(k, idx);
end
vlow = abs(vlow) + R;
vsparse = vsparse - R;
vc = real(vc); % since it's complex
end