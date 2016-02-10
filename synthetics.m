TABLE = zeros(100,13);

STRIKE = 0 + (360).*rand(100,1);
DIP = 0 + 90.*rand(100,1);
RAKE = -180 + (180+180).*rand(100,1);

TABLE(:,3) = STRIKE;
TABLE(:,4) = DIP;
TABLE(:,5) = RAKE;

save('Synthetics','TABLE');