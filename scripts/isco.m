%0.68 is for an equal mass, nonspinning merger
af = 0.68;

Z1 = 1 + (1 - af^2)^(1/3)*((1 + af)^(1/3) + (1 - af)^(1/3));
Z2 = sqrt(3*af^2 + Z1^2);
risco = 3 + Z2 - sqrt((3 - Z1)*(3 + Z1 + 2*Z2))
Omisco = 1/(risco^(3/2) + af)