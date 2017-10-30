inf1 = 100000;
maxL = -inf1;
maxA = -inf1;
maxB = -inf1;
minL = inf1;
minA = inf1;
minB = inf1;

for r = 0:255
    for g = 255
        for b = 0:255
            l = rgb2lab([r,g,b]);
            maxL = max(maxL, l(1));
            maxA = max(maxA, l(2));
            maxB = max(maxB, l(3));
            minL = min(minL, l(1));
            minA = min(minA, l(2));
            minB = min(minB, l(3));
        end
    end
    r
end
