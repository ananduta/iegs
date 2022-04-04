function id = index_decision2(nu,nodeSet,h)

    sum_nu = 0;
    for ii = 1:length(nodeSet)
        i = nodeSet(ii);
        id_start(ii,1) = sum_nu + 1;
        id_end(ii,1) = sum_nu + nu(i)*h;
        sum_nu = sum_nu + nu(i)*h;
    end
    id =[id_start, id_end];
end