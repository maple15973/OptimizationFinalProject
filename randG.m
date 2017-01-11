function G = randG(LINK_NUM)
G_half = unifrnd(0,1,LINK_NUM,LINK_NUM);
G_other = (G_half+G_half')/100;
G_other(eye(LINK_NUM)~=0)=0;
G_diag = unifrnd(0.5,1,LINK_NUM,1);
G = G_other + diag(G_diag);