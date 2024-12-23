L=100e-3; H=20e-3; NumElemAlongSide=100;
[SysCoor,ElemNode,B1,B2,B3,B4] = get_mesh_filter(L,H,NumElemAlongSide);

SysCoor(:,2) = SysCoor(:,2) + H/2;

B0 = find(SysCoor(:,1) >= -1e-10 & SysCoor(:,2) < 1e-10);
B2 = find(SysCoor(:,2) > H-1e-10);
b0y = 2*B0; b0x = 2*B0-1; b2y = 2*B2; b2x = 2*B2-1;

BC=[b0y b0y*0];


