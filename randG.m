function G = randG(LINK_NUM)
G_half = unifrnd(0,1,LINK_NUM,LINK_NUM);
G = (G_half+G_half')/20 + diag(diag(G_half))*0.9;
