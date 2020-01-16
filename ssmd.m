function [SSMD] = ssmd(background, target)
    SSMD = (mean(target) - mean(background)) / sqrt(var(target) + var(background));
end