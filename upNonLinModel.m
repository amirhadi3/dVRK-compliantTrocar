function y = upNonLinModel(instObj,possig,Nn,ft,maxati,dofselector)
    instObj = Nn_to_wt(instObj,Nn,possig);
    wt = instObj.wt./maxati;
    y = (wt-ft).*dofselector;
end