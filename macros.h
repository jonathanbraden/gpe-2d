#define SPECTRAL T

#define LAPLACIAN(x,y) fld(i+(x),j+(y),:)
#define GRAD2(x,y) (fld(i+(x),j+(y),:)-fld(i,j,:))**2
#define RANK0(O) (O(0,0))
#define RANK1(O) (O(-1,0) + O(1,0) + O(0,-1) + O(0,1))
#define RANK2(O) (O(-1,-1) + O(-1,1) + O(1,-1) + O(1,1))
#define STENCIL(C,O) ((C ## 0)*RANK0(O) + (C ## 1)*RANK1(O) + (C ## 2)*RANK2(O) )
