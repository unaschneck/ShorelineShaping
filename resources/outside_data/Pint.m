function [P,nbi,theta] = Pint(P1,P2,P3,eps)

% x=rand(1,3);
% y=rand(1,3);
% figure
% plot(x,y)
% x
% y
x0=P2(1);
y0=P2(2);
xp=P1(1);
yp=P1(2);
xn=P3(1);
yn=P3(2);
vp = [x0-xp,y0-yp]; % The vector from P1 to P2
vn = [xn-x0,yn-y0]; % The vector from P2 to P3
np = [vp(2),-vp(1)];
nn = [vn(2),-vn(1)];
np = np/sqrt(np(1)^2+np(2)^2); % vp rotated clockwise by 90 degrees and normalized
nn = nn/sqrt(nn(1)^2+nn(2)^2); % vn rotated clockwise by 90 degrees and normalized
% if any(isnan(nn))
%     disp(pause)% This happens when P3 = P2
% end
% pmid = [0.5*(xp+x0),0.5*(yp+y0)];
% nmid = [0.5*(xn+x0),0.5*(yn+y0)];
nbi=np+nn;
if ~any(nbi) % x0,y0 is on a promontory or a zero-width inlet, so nbi == [0,0] (true if all elements are zero)
    % ROSE: LAKES CAN HAVE PROMONTORIES BUT NO ZERO-WIDTH INLETS, RIGHT?
    % WOULD IT BE THE OTHER WAY AROUND FOR ISLANDS? HERE I HAVE ASSUMED
    % ONLY PROMONTORIES, SO WE CAN JUST EXTRAPOLATE ALONG THE PREVIOUS
    % SEGMENT TO GET THE OBSERVATION POINT
    nbi = vp;
    if vp == 0
        nbi = [x0-xn,y0-yn]; % The vector from P1 to P3 -- this occurs when P1 and P2 are the same point.
    end
%     % continue along the previous segment a distance eps
%     P = [x0,y0]+eps.*(vp./norm(vp));
% %     if norm(vp) == 0
% %         disp("we have a problem")
% %     end
end
nbi = nbi/sqrt(nbi(1)^2+nbi(2)^2);
P = [x0,y0]+eps*nbi;

end

% figure
% plot([xp,x0,xn],[yp,y0,yn])
% axis equal
% hold on
% 
% % plot([pmid(1),pmid(1)+np(1)],[pmid(2),pmid(2)+np(2)],'r')
% % plot([nmid(1),nmid(1)+nn(1)],[nmid(2),nmid(2)+nn(2)],'r')
% plot([x0,x0+nbi(1)],[y0,y0+nbi(2)],'r')
% plot(P(1),P(2),'ok')
