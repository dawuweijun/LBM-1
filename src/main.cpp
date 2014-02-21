#include "BGK2D9V.hpp"

#define N 100
#define M 100

int main(){
	BGK2D9V<N,M> bgk(1.3);
	
	char filename[64];
	for(int i=1;i<=100;i++){
		bgk.evolution(100);

		printf("%08d\n",i*100);
		sprintf(filename,"out/%08d.dat",i);
		FILE* fp=fopen(filename,"w");
		bgk.output(fp);
		fclose(fp);
	}

	return 0;
}
