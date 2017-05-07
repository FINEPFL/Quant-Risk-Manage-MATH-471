function ret = nll(pse, x, copulaname)
    ret = -sum(log(copulapdf(copulaname, pse, x)));
end