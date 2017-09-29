function [F]=sumofThreeGauss(a, x)


F = a(1).*exp( -( (x-a(2)).^2 / (2.*a(3).^2) ) ) ...
    + a(4).*exp( -( (x-a(5)).^2 / (2.*a(6).^2) ) ) ...
    + a(7).*exp( -( (x-a(8)).^2 / (2.*a(9).^2) ) );

end