function [GSSMD] = gssmd(background, target)
    n = length(background);
    Range=linspace(min([background target]),max([background target]),round(1+log2(n)));
    % estimate probability density function
    pdf_b=hist(background, Range);
    pdf_t=hist(target, Range);    
    % estimate overlap
    OVL=sum(min([pdf_t./sum(pdf_t); pdf_b./sum(pdf_b)]));
    GSSMD = ((mean(target) - mean(background))/abs((mean(target) - mean(background))))*(1-OVL);   
    if (mean(target) == mean(background))
        GSSMD=0;
    end
end