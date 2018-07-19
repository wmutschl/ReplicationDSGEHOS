function ccc=varhac(datin,it1,it2,k1,k2,imax,ilag,imodel,imean)

% this is a MATLAB adaptation of Wouter den Haan's GAUSS procedure for 
% computing den Haan and Levin's VARHAC estimator

% The input is a data matrix "datin" where each row corresponds to a
% different observation

% it1:        first observation to be used
% it2:        last  observation to be used
% k1:         first column to be used
% k2:         last  column to be used
% imax:       maximum lag order considered  (if imax = 0, then no
%             correction for serial correlation will be made)
% ilag:       if equal to 1, then all elements enter with the same lag
%             if equal to 2, then the own lag can enter with a different lag
%             if equal to 3, then only the own lag enters
% imodel:     if equal to 1, then AIC is used
%             if equal to 2, then SCHWARTZ is used
%             if equal to 3, then a fixed lag order equal to imax is used
% imean:      if equal to 1, then the mean will be subtracted from each series

nt=it2-it1+1 ;
kdim=k2-k1+1 ;

dat=datin(it1:it2,k1:k2);

if imean==1
    ave=mean(dat);
    dat=dat-ones(nt,1)*ave ;
end
ddd=dat(imax+1:nt,1:kdim) ;
if imax~=0
    minpar=zeros(kdim,kdim*imax) ;
    minres=ddd' ;
    minorder=zeros(kdim,2) ;
    aic=((nt-imax-1)/(nt-imax))*std(ddd,1).^2 ;
    aic=log(aic) ;

    % DETERMINE HOW TO ESTIMATE THE PARAMETRIC MODEL 

    if imodel~=3
        if ilag==2

            % ILAG=2, THUS THE OWN ELEMENT CAN ENTER WITH A DIFFERENT LAG FOR EACH
            % L WE FIND THE MODEL THAT MINIMIZES THE MODEL SELECTION CRITERION

            for l=1:kdim

                % step 1: construct dependent and independent variables

                exo1=dat(1:nt,l)
                if l==1
                    exo2=dat(1:nt,2:kdim) ;
                end
                if l==kdim
                    exo2=dat(1:nt,1:kdim-1) ;
                end
                if l>1 
                    if l<kdim
                        temp1=dat(1:nt,1:l-1) ;
                        temp2=dat(1:nt,l+1:kdim) ;
                        exo2=[ temp1 temp2 ] ;
                    end
                end
                dep=dat(imax+1:nt,l:l) ;
                ss1=exo1(imax+1-1:nt-1,1) ;
                for iorder1=0:imax 
                    if iorder1>=2 
                        tt1=exo1(imax+1-iorder1:nt-iorder1,1) ;
                        ss1=[ ss1 tt1 ] ;
                    end
                    ss2=exo2(imax+1-1:nt-1,1:kdim-1) ;
                    for iorder2=0:imax
                        if (iorder1+iorder2)~=0
                            if iorder2>=2
                                tt2=exo2(imax+1-iorder2:nt-iorder2,1:kdim-1) ;
                                ss2=[ ss2 tt2 ];
                            end
                            if iorder2==0
                                ex1=ss1 ;
                            end
                            if iorder1==0
                                ex1=ss2 ;
                            end
                            if iorder1>0
                                if iorder2>0
                                    ex1=[ ss1 ss2 ] ;
                                end
                            end

                            % step 2: do the regression      

                            b=(ex1'*ex1)\(ex1'*dep) ;
                            resid=dep-ex1*b ;

                            % step 3: calculate the model selection criterion         

                            npar=cols(ex1) ;
                            aicnew=((nt-imax-1)/(nt-imax))*std(resid,1).^2 ;
                            if imodel==1 
                                aicnew=log(aicnew)+2*npar/(nt-imax) ;
                            else
                                aicnew=log(aicnew)+log(nt-imax)*npar/(nt-imax) ;
                            end

                            % step 4: if there is improvement store minpar, minres, and minorder 
     
                            if aicnew<aic(l) ;
                                aic(l)=aicnew ;
                                for t=1:kdim*imax 
                                    minpar(l,t)=0 ;
                                end
                                for iorder=1:iorder1 ;
                                    minpar(l,l+(iorder-1)*kdim)=b(iorder,1) ;
                                end
                                for iorder=1:iorder2
                                    if l==1 
                                        minpar(l,2+(iorder-1)*kdim:kdim+(iorder-1)*kdim)= ...
                                        b(iorder1+1+(iorder-1)*(kdim-1):iorder1+kdim-1+(iorder-1)*(kdim-1),1)' ;
                                    end
                                    if l==kdim
                                        minpar(l,1+(iorder-1)*kdim:kdim-1+(iorder-1)*kdim)= ...
                                        b(iorder1+1+(iorder-1)*(kdim-1):iorder1+kdim-1+(iorder-1)*(kdim-1),1)' ;
                                    end
                                    if l>1 
                                        if l<kdim
                                            minpar(l,1+(iorder-1)*kdim:(l-1)+(iorder-1)*kdim)= ...
                                             b(iorder1+1+(iorder-1)*(kdim-1):iorder1+(l-1)+(iorder-1)*(kdim-1),1)' ;
                                            minpar(l,l+1+(iorder-1)*kdim:kdim+(iorder-1)*kdim)= ...
                                             b(iorder1+l+(iorder-1)*(kdim-1):iorder1+kdim-1+(iorder-1)*(kdim-1),1)' ;
                                        end 
                                    end
                                end
                                minres(l,1:nt-imax)=resid' ; 
                                minorder(l,1)=iorder1 ;
                                minorder(l,2)=iorder2 ; 
                            end
                        end
                    end
                end
            end
        elseif ilag==3

            % ILAG=3, THUS ONLY THE OWN ELEMENT ENTERS FOR EACH L WE FIND
            % THE MODEL THAT MINIMIZES THE MODEL SELECTION CRITERION

            for l=1:kdim

                % step 1: construct dependent and independent variables 

                exo1=dat(1:nt,l) ;
                dep=dat(imax+1:nt,l:l) ;
                ss1=exo1(imax+1-1:nt-1,1) ;
                for iorder1=1:imax
                    if iorder1>=2
                        tt1=exo1(imax+1-iorder1:nt-iorder1,1) ;
                        ss1=[ ss1 tt1 ] ;
                    end
                    ex1=ss1 ;
        
                    % step 2: do the regression

                    b=(ex1'*ex1)\(ex1'*dep) ;
                    resid=dep-ex1*b ;

                    % step 3: calculate the model selection criterion

                    npar=cols(ex1) ;
                    aicnew=((nt-imax-1)/(nt-imax))*std(resid,1).^2 ;
                    if imodel==1 
                        aicnew=log(aicnew)+2*npar/(nt-imax) ;
                    else 
                        aicnew=log(aicnew)+log(nt-imax)*npar/(nt-imax) ;
                    end

                    % step 4: if there is improvement store minpar, minres, and minorder */

                    if aicnew < aic(l)
                        aic(l)=aicnew ;
                        for t=1:kdim*imax
                            minpar(l,t)=0 ;
                        end
                        for iorder=1:iorder1
                            minpar(l,l+(iorder-1)*kdim)=b(iorder,1) ;
                        end
                        minres(l,1:nt-imax)=resid' ;
                        minorder(l,1)=iorder1 ;
                    end
                end
            end
        else
            
            % ILAG = 1, THUS ALL VARIABLES ENTER WITH THE SAME LAG

            for l=1:kdim
                dep=dat(imax+1:nt,l:l) ;
                ex1=dat(imax+1-1:nt-1,1:kdim) ;
                for iorder=1:imax;
                    if iorder>=2 
                        ex2=dat(imax+1-iorder:nt-iorder,1:kdim) ;
                        ex1=[ ex1 ex2 ] ;
                    end
                    b=(ex1'*ex1)\(ex1'*dep) ;
                    resid=dep-ex1*b ;
                    aicnew=((nt-imax-1)/(nt-imax))*std(resid,1).^2 ;
                    if imodel==1
                        aicnew=log(aicnew)+2*kdim*iorder/(nt-imax) ;
                    else
                        aicnew=log(aicnew)+log(nt-imax)*kdim*iorder/(nt-imax) ;
                    end
                    if aicnew < aic(l)
                        aic(l)=aicnew ;
                        minpar(l,1:kdim*iorder)=b' ;
                        minres(l,1:nt-imax)=resid' ; 
                        minorder(l,1)=iorder ;
                    end
                end
            end
        end
    else
        
       % IMODEL=3, THUS A FIXED LAG FOR ALL ELEMENTS

       for l=1:kdim
            dep=dat(imax+1:nt,l:l) ;
            ex1=dat(imax+1-1:nt-1,1:kdim) ;
            for iorder=2:imax
                ex2=dat(imax+1-iorder:nt-iorder,1:kdim) ;
                ex1=[ ex1 ex2 ] ;
            end
            b=(ex1'*ex1)\(ex1'*dep) ;
            resid=dep-ex1*b ;
            minpar(l,1:kdim*imax)=b' ;
            minres(l,1:nt-imax)=resid' ; 
            minorder(l,1)=imax ;
            minorder(l,2)=imax ;
        end
    end
    
    % THE TIME-SERIES REPRESENTATION HAS NOW BEEN ESTIMATED, AND
    % WE CAN CALCULATE THE VARHAC ESTIMATOR

    covm=(minres*minres')/(nt-imax) ;
    bbb=eye(kdim) ;
    for iorder=1:imax
        bbb=bbb-minpar(1:kdim,(iorder-1)*kdim+1:iorder*kdim) ;
    end
    ccc=inv(bbb)*covm*inv(bbb') ;
else
    
    % IF NO CORRECTION FOR SERIAL CORRELATION IS REQUIRED, THEN
    %   THE VARIANCE IS EQUAL TO

    ccc=(ddd'*ddd)/nt ;
end